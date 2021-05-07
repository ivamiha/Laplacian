// File       : main.cpp    
// Created    : Tue Mar 16 2021  9:46:31 am CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: Benchmarking suit for quantifying Laplacian kernel performance
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include "Cubism/Block/FieldLab.h"
#include "Cubism/Grid/Cartesian.h"
#include "Cubism/Mesh/StructuredUniform.h"
#include "Cubism/Block/Field.h"
#include "Cubism/Util/Timer.h"

#include "LaplacianSecondOrder.h"
#include "LaplacianFourthOrder.h"

#include <cstdio>
#include <vector>
#include <omp.h>
#include <iostream>
#include <fstream>

// enable fourth-order CDS, default is second-order CDS 
//#define USE_ACCUR 
// enable AVX vector extension based calculations, default is SSE
//#define USE_AVX
// enable OpenMP TLP, default is single core
//#define USE_TLP

using namespace Cubism;
using Util::Timer; 

/** 
 * @brief Get Roofline ceilings & performance based on benchmark parameters
 * @param T Datatype specifying arithmetic precision used in benchmark
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
    // specify vector extension-dependent limits (SIMD)
#ifdef USE_AVX 
    const size_t maxLanes = (256/8) / sizeof(T); 
#else
    const size_t maxLanes = (128/8) / sizeof(T); 
#endif /* USE_AVX */ 
    // specify execution-dependent values

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
        performance_zero[i] = 1E-09 * std::pow(64,3) * flopCell / time_zero[i];
        performance_inft[i] = 1E-09 * std::pow(32,3) * flopCell / time_inft[i];
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

    const MIndex nblocks(1); 
    const MIndex block_cells_zero(64);  // for zero cache benchmark
    const MIndex block_cells_inft(32);  // for infinite cache benchmark

    // identifiers for creating & managing scalar fields
    using SGrid = Grid::Cartesian<double, Mesh, EntityType::Cell, 0>;  
    using DataType = typename SGrid::DataType; 
    using FieldType = typename SGrid::BaseType; 
    using FieldLab = Block::FieldLab<typename FieldType::FieldType>; 
    using Stencil = typename FieldLab::StencilType;
   
    // define number of kernel calls to be executed in benchmarking
    size_t N = ((argc == 2) ? std::atoi(argv[1]) : 100);     
    
    // map the blocks & cells on domain [0,1]^3 [m] 
    SGrid sol_zero(nblocks, block_cells_zero);  // solution storage zero cache
    SGrid tmp_zero(nblocks, block_cells_zero);  // temporary storage zero cache
    SGrid sol_inft(nblocks, block_cells_inft);  // solution storage inft cache
    SGrid tmp_inft(nblocks, block_cells_inft);  // temporary storage inft cache

    // function for writing ICs utilizing input field (block) within grid     
    auto fIC = [](FieldType &b) {
        const DataType fac = 2 * M_PI;
        const Mesh &bm = *b.getState().mesh;
        // loop over cells in the block's mesh using cell index (ci) 
        for (auto &ci : bm[EntityType::Cell]) {
            const PointType x = bm.getCoordsCell(ci); 
            b[ci] = std::cos(fac * x[0]) * std::cos(fac * x[1]) *
                    std::sin(fac * x[2]);
        } 
    };

    // initialize solution storage grids by looping over blocks in the grid
    for (auto bf : sol_zero) {
        fIC(*bf); 
    }
    for (auto bf : sol_inft) {
        fIC(*bf); 
    }

    // define numerical parameters based on chosen discretization
#ifdef USE_ACCUR                      
    const Stencil s(-2, 3, true);                
#else                                          
    const Stencil s(-1, 2, false);              
#endif /* USE_ACCUR */
    
    // setup labs 
    FieldLab flab_zero, flab_inft; 
    flab_zero.allocate(s, sol_zero[0].getIndexRange()); 
    flab_inft.allocate(s, sol_inft[0].getIndexRange()); 

    // get block field index functors for periodic block accessing
    auto findex_zero_ = sol_zero.getIndexFunctor(0);
    auto findex_inft_ = sol_inft.getIndexFunctor(0); 

    // warm up - call Laplacian kernel once
    for (auto f : sol_zero) {
        const FieldType &bf = *f; 
        const MIndex &bi = bf.getState().block_index; 
        sol_zero.loadLab(bf, flab_zero); 
        auto &tf = tmp_zero[bi]; 
        LaplacianSecondOrder(flab_zero, tf); 
    }

    // setup timer & required storage 
    Timer t; 
    std::vector<double> time_zero(N);           // time per zero cache run 
    std::vector<double> performance_zero(N);    // performance per zero run
    std::vector<double> time_inft(N);           // time per inft cache run
    std::vector<double> performance_inft(N);    // performance per inft run
    std::vector<double> roofCoords(10);         // coords. for roofline model
    
    // run zero cache benchmark N times 
    for (size_t i = 0; i < N; ++i) 
    {   
#ifdef USE_TLP
        #pragma omp parallel for num_threads(4)
#endif /* USE TLP_ */
        // loop through blocks in the grid
        for (auto f : sol_zero) 
        {
            // reference fields & load data into flab object for current block
            const FieldType &bf = *f;  
            const MIndex &bi = bf.getState().block_index; 
            sol_zero.loadLab(bf, flab_zero); 
            auto &tf = tmp_zero[bi];
            
            t.start();
            // benchmark selected Laplacian kernel
#ifdef USE_ACCUR 
            LaplacianFourthOrder(flab_zero, tf); 
#else
            LaplacianSecondOrder(flab_zero, tf); 
#endif /* USE_ACCUR */
            time_zero[i] = t.stop(); 
        }
    }

    // run inft cache benchmark N times
    for (size_t i = 0; i < N; ++i) 
    {
#ifdef USE_TLP 
        #pragma omp parallel for num_threads(4)
#endif /* USE_TLP */
        // loop through blocks in the grid
        for (auto f : sol_inft)
        {
            const FieldType &bf = *f; 
            const MIndex &bi = bf.getState().block_index;
            sol_inft.loadLab(bf, flab_inft); 
            auto &tf = tmp_inft[bi];
           
            t.start(); 
            // benchmark selected Laplacian kernel 
#ifdef USE_ACCUR 
            LaplacianFourthOrder(flab_inft, tf); 
#else 
            LaplacianSecondOrder(flab_inft, tf); 
#endif /* USE_ACCUR */
            time_inft[i] = t.stop(); 
        }
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
    printf("Ranked peak performane for infinite cache after %ld runs: \n", N);
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

    return 0; 
}
