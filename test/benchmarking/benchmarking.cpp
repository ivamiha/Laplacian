// File       : main.cpp
// Created    : Tue Mar 16 2021  9:46:31 am CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: Benchmarking suit for quantifying Laplacian kernel performance
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include "Cubism/Block/FieldLab.h"
#include "Cubism/Grid/Cartesian.h"
#include "Cubism/Mesh/StructuredUniform.h"
#include "Cubism/Block/Field.h"
#include "Cubism/IO/CartesianHDF.h"

#include "LaplacianSecondOrder.h"
#include "LaplacianFourthOrder.h"

#include <cstdio>
#include <vector>
#include <omp.h>
#include <iostream>
#include <fstream>

// enable fourth-order CDS, default is second-order CDS 
//#define USE_ACCUR 

using namespace Cubism;

/** 
 * @brief Get Roofline ceilings & performance based on benchmark parameters
 * @tparam T Datatype specifying arithmetic precision used in benchmark
 * @param time_zero Stored measured CPU time for each zero cache run
 * @param performance_zero Calculate and store zero cache performance GFlops/s
 * @param time_inft Stored measured CPU time for each inft cache run
 * @param performance_inft Calculate and store inft cache performance GFlops/s
 * @param roofCoords Calculate and store Roofline ceilings and ridge points 
 *
 * @rst Function template for calculating and storing Roofline-relevant 
 * data. Architecture-dependent parameters specified within function itself. 
 * @endrst
 * */
template <typename T>
void getRoofline(const std::vector<double> &time_zero,
                 std::vector<double> &performance_zero,
                 const std::vector<double> &time_inft,
                 std::vector<double> &performance_inft,
                 std::vector<double> &roofCoords)
{
    // specify architecture-dependent limits (Intel Xeon E5-2670 v3)
    const float maxFreq = 3.2;      // maximum CPU frequency [GHz] 
    const size_t maxCores = 12;     // maximum CPU cores available []
    const float maxBand = 68.0;     // maximum memory bandwidth [GB/s]
    const size_t FMA = 2;           // fused multiply-add effect []
    // specify vector-extensions (SIMD) - architecture has AVX2 as limit
    const size_t maxLanes = (256/8) / sizeof(T);

    // specify benchmarked kernel-dependent parameters
#ifdef USE_ACCUR
    const float opInt_zero = 0.180; // op. intensity no cache [Flops/Byte]
    const float opInt_inft = 1.438; // op. intensity infty cache [Flops/Byte]
    const size_t flopCell = 23;     // flops per processed cell [Flops/Cell]
#else
    const float opInt_zero = 0.175;
    const float opInt_inft = 0.875;
    const size_t flopCell = 14;
#endif /* USE_ACCUR */

    // compute measured performance for each benchmark run
    for (size_t i = 0; i < time_zero.size(); ++i) {
        performance_zero[i] = 1E-09 * std::pow(128,3) * std::pow(4,3)
                                                     * flopCell / time_zero[i];
        performance_inft[i] = 1E-09 * std::pow(32,3) * std::pow(4,3)
                                                     * flopCell / time_inft[i];
    }
    // sort performance vectors in ascending order
    std::sort(performance_zero.begin(), performance_zero.end());
    std::sort(performance_inft.begin(), performance_inft.end());

    // compute roofline ceilings
    const double ceil1 = maxFreq * FMA;
    const double ceil2 = maxFreq * FMA * maxLanes;
    const double ceil3 = maxFreq * FMA * maxLanes * maxCores;
    // compute average performance for both test scenarios
    double avg_zero = 0;
    double avg_inft = 0;
    for (size_t i = 0; i < performance_zero.size(); ++i) {
        avg_zero += performance_zero[i];
        avg_inft += performance_inft[i];
    }
    avg_zero /= performance_zero.size();
    avg_inft /= performance_inft.size();
    // store relevant roofline model coordinates
    roofCoords[0] = ceil1 / maxBand;
    roofCoords[1] = ceil1;
    roofCoords[2] = ceil2 / maxBand;
    roofCoords[3] = ceil2;
    roofCoords[4] = ceil3 / maxBand;
    roofCoords[5] = ceil3;
    roofCoords[6] = opInt_zero;
    roofCoords[7] = opInt_inft;
    roofCoords[8] = avg_zero;
    roofCoords[9] = avg_inft;
}

int main(int argc, char *argv[])
{
    // welcome & inform user what program configuration is running
    printf("--------------------------------------------------------------\n");
    printf("T E S T   L A P L A C I A N   C O M P U T E   K E R N E L\n");
    printf("==============================================================\n");
#ifdef USE_ACCUR
    printf("Benchmarking 4th-order CDS implementation.\n"); 
#else
    printf("Benchmarking 2nd-order CDS implementation.\n"); 
#endif /* USE_ACCUR */
    
    // identifers to be used for creating & managing meshes
    using IRange = Core::IndexRange<3>;
    using MIndex = typename IRange::MultiIndex;
    using Mesh = Mesh::StructuredUniform<double, IRange::Dim>;
    using PointType = typename Mesh::PointType;

    const MIndex nblocks(4);
    const MIndex block_cells_zero(128); // for zero cache benchmark
    const MIndex block_cells_inft(32);  // for infinite cache benchmark

    // identifiers for creating & managing scalar fields
    using SGrid = Grid::Cartesian<double, Mesh, EntityType::Cell, 0>;
    using DataType = typename SGrid::DataType;
    using FieldType = typename SGrid::BaseType;
    using FieldLab = Block::FieldLab<typename FieldType::FieldType>;
    using Stencil = typename FieldLab::StencilType;

    // define number of kernel calls to be executed in benchmarking
    size_t N = 100;
    // define number of threads to be utilized
    size_t nthreads = ((argc == 2) ? std::atoi(argv[1]) : 1);
    
    // map the blocks & cells on domain [0,1]^3 [m]
    SGrid sol_zero(nblocks, block_cells_zero);  // solution storage zero cache
    SGrid tmp_zero(nblocks, block_cells_zero);  // temporary storage zero cache
    SGrid sol_inft(nblocks, block_cells_inft);  // solution storage inft cache
    SGrid tmp_inft(nblocks, block_cells_inft);  // temporary storage inft cache

    // function for writing ICs utilizing input field (block) within grid
    auto fIC = [](FieldType &b) {
        const Mesh &bm = *b.getState().mesh;
        // loop over cells in the block's mesh using cell index (ci)
        for (auto &ci : bm[EntityType::Cell]) {
            const PointType x = bm.getCoordsCell(ci);
            b[ci] = std::sin(2.0 * M_PI * x[0]) *
                    std::cos(6.0 * M_PI * x[1]) *
                    std::sin(4.0 * M_PI * x[2]);
        }
    };

    // initialize solution storage grids by looping over blocks in the grid
    for (auto bf : sol_zero) {
        fIC(*bf);
    }
    for (auto bf : sol_inft) {
        fIC(*bf);
    }

    // dump the ICs to HDF5 files using single-precision
    IO::CartesianWriteHDF<float>("init_zero", "init", sol_zero, 0);
    IO::CartesianWriteHDF<float>("init_inft", "init", sol_inft, 0);

    // define numerical parameters based on chosen discretization
#ifdef USE_ACCUR
    const Stencil s(-2, 3, false);
#else
    const Stencil s(-1, 2, false);
#endif /* USE_ACCUR */

    // setup labs
    FieldLab flab_zero, flab_inft;
    flab_zero.allocate(s, sol_zero[0].getIndexRange());
    flab_inft.allocate(s, sol_inft[0].getIndexRange());

    // warm up - call Laplacian kernel once
    for (auto f : sol_zero) {
        const FieldType &bf = *f;
        const MIndex &bi = bf.getState().block_index;
        sol_zero.loadLab(bf, flab_zero);
        auto &tf = tmp_zero[bi];
        LaplacianSecondOrder(flab_zero, tf);
    }

    // setup required storage for timing & roofline
    std::vector<double> start(nthreads, 0);     // per-thread start time
    std::vector<double> time(nthreads, 0);      // per-thread measured time
    std::vector<double> time_zero(N, 0);        // time per zero cache run 
    std::vector<double> performance_zero(N);    // performance per zero run
    std::vector<double> time_inft(N, 0);        // time per inft cache run
    std::vector<double> performance_inft(N);    // performance per inft run
    std::vector<double> roofCoords(10);         // coords. for roofline model

    // run inft cache benchmark N times
    for (size_t i = 0; i < N; ++i)
    {
        // loop through blocks in the grid
#pragma omp parallel num_threads(nthreads) private(flab_inft)
        {
            const int tid = omp_get_thread_num();
            flab_inft.allocate(s, sol_inft[0].getIndexRange());
#pragma omp for nowait
            for (auto f : sol_inft)
            {
                const FieldType &bf = *f;
                const MIndex &bi = bf.getState().block_index;
                sol_inft.loadLab(bf, flab_inft);
                auto &tf = tmp_inft[bi];

                // benchmark selected Laplacian kernel
                start[tid] = omp_get_wtime();
#ifdef USE_ACCUR
                LaplacianFourthOrder(flab_inft, tf);
#else
                LaplacianSecondOrder(flab_inft, tf);
#endif /* USE_ACCUR */
                time[tid] += omp_get_wtime()- start[tid];
            }
        }
        // compute average execution time across all threads
        time_inft[i] = *max_element(std::begin(time), std::end(time));
        for (size_t t = 0; t < nthreads; ++t) time[t] = 0;
    }

    // run zero cache benchmark N times
    for (size_t i = 0; i < N; ++i) 
    {
        // loop through blocks in the grid
#pragma omp parallel num_threads(nthreads) private(flab_zero)
        {
            const int tid = omp_get_thread_num();
            flab_zero.allocate(s, sol_zero[0].getIndexRange());
#pragma omp for nowait
            for (auto f : sol_zero)
            {
                const FieldType &bf = *f;
                const MIndex &bi = bf.getState().block_index;
                sol_zero.loadLab(bf, flab_zero);
                auto &tf = tmp_zero[bi];

                // benchmark selected Laplacian kernel
                start[tid] = omp_get_wtime();
#ifdef USE_ACCUR
                LaplacianFourthOrder(flab_zero, tf);
#else
                LaplacianSecondOrder(flab_zero, tf);
#endif /* USE_ACCUR */
                time[tid] += omp_get_wtime() - start[tid];
            }
        }
        // compute average execution time across all threads
        time_zero[i] = *max_element(std::begin(time), std::end(time));
        for (size_t t = 0; t < nthreads; ++t) time[t] = 0;
    }

    // pass benchmarking measurements to getRoofline template function
    getRoofline<DataType>(time_zero, performance_zero,
                                     time_inft, performance_inft, roofCoords); 

    // define percentiles to be selected
    const int i1 = N * 0.9;
    const int i2 = N * 0.5;
    const int i3 = N * 0.1;

    // output 10th, 50th, and 90th percentiles
    printf("--------------------------------------------------------------\n"); 
    printf("Ranked peak performance for zero cache after %ld runs: \n", N);
    printf("10th percentile:\t%f GFlops/s\n", performance_zero[i3]);
    printf("50th percentile:\t%f GFlops/s\n", performance_zero[i2]);
    printf("90th percentile:\t%f GFlops/s\n", performance_zero[i1]);
    printf("Average performance:\t%f GFlops/s\n\n", roofCoords[8]);
    printf("Ranked peak performance for infinite cache after %ld runs: \n", N);
    printf("10th percentile:\t%f GFlops/s\n", performance_inft[i3]);
    printf("50th percentile:\t%f GFlops/s\n", performance_inft[i2]);
    printf("90th percentile:\t%f GFlops/s\n", performance_inft[i1]);
    printf("Average performance:\t%f GFlops/s\n", roofCoords[9]);
    printf("--------------------------------------------------------------\n"); 

    // write relevant data to file
    std::ofstream results;
    results.open("roof.txt");
    results << std::fixed << std::setprecision(2) << std::endl;
    for (size_t i = 0; i < roofCoords.size(); ++i) {
        results << roofCoords[i] << "\n";
    }
    results.close();

    // dump solutions to HDF5 files using single-precision
    IO::CartesianWriteHDF<float>("sol_zero", "sol", tmp_zero, 0);
    IO::CartesianWriteHDF<float>("sol_inft", "sol", tmp_inft, 0);

    return 0;
}
