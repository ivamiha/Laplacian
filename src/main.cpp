// File       : main.cpp    
// Created    : Tue Mar 16 2021  9:46:31 am CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: File with main function for testing Laplacian function
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include "Cubism/Block/FieldLab.h"
#include "Cubism/Grid/Cartesian.h"
#include "Cubism/Mesh/StructuredUniform.h"
#include "Cubism/Block/Field.h"
#include "Cubism/Util/Timer.h"

#include "LaplacianSecondOrder.h"
#include "LaplacianFourthOrder.h"

#include <cstdio>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>

// enable fourth-order CDS, default is second-order CDS 
#define USE_ACCUR 
// enable AVX vector extension based calculations, default is SSE
//#define USE_AVX

using namespace Cubism;
using Util::Timer; 

/** 
 * @brief Get Roofline ceilings & performance based on benchmark parameters
 * @param T Datatype specifying arithmetic precision used in benchmark
 * @param time Stored measured CPU time for each benchmark run
 * @param performance Calculate and store benchmark performance GFlops/s
 * @param ceilings Calculate and store coordinates of Roofline ceilings 
 *
 * @rst Function template for calculating and storing Roofline-relevant 
 * data. Architecture-dependent parameters specified within function itself. 
 * @endrst
 * */
template <typename T>
void getRoofline(std::vector<double> &time, 
                 std::vector<double> &performance, 
                 std::vector<double> &roofCoords) 
{
    // specify architecture-dependent limits (Intel Xeon E5-2670 v3)
    const float maxFreq = 3.1;      // maximum CPU frequency [GHz] 
    const size_t maxCores = 12;     // maximum CPU cores available []
    const float maxBand = 68.0;     // maximum memory bandwidth [GB/s]
    const size_t FMA = 2;           // fused multiply-add effect []
    // specify vector extension-dependent limits (SIMD)
#ifdef USE_AVX 
    const size_t maxLanes = (256/8) / sizeof(T); 
#else
    const size_t maxLanes = (128/8) / sizeof(T); 
#endif /* USE_AVX */ 
    // specify actually used architecture-dependent parameters
    const size_t nCores = 1; 
    const size_t nLanes = 1; 

    // specify benchmarked kernel-dependent parameters
#ifdef USE_ACCUR
    const float opInt = 0.230;      // operational intensity [Flops/Byte]
    const size_t flopCell = 35;     // flops per processed cell [Flops/Cell]
#else 
    const float opInt = 0.192;
    const size_t flopCell = 20;  
#endif /* USE_ACCUR */

    // compute measured performance for each benchmark run
    for (size_t i = 0; i < time.size(); ++i) {
        performance[i] = 1.0E-09 * std::pow(32,3) * flopCell * FMA 
                       * nCores * nLanes / time[i];
    }
    // sort performance vector in ascending order 
    std::sort(performance.begin(), performance.end());  
    
    // compute roofline ceilings
    const double ceil1 = maxFreq * FMA;  
    const double ceil2 = maxFreq * FMA * maxLanes; 
    const double ceil3 = maxFreq * FMA * maxLanes * maxCores; 
    // store relevant roofline model coordinates
    roofCoords[0] = ceil1 / maxBand; 
    roofCoords[1] = ceil1; 
    roofCoords[2] = ceil2 / maxBand; 
    roofCoords[3] = ceil2; 
    roofCoords[4] = ceil3 / maxBand; 
    roofCoords[5] = ceil3; 
    roofCoords[6] = opInt; 
    roofCoords[7] = performance[performance.size()-1]; 
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
    const MIndex block_cells(32);

    // identifiers for creating & managing scalar fields
    using SGrid = Grid::Cartesian<double, Mesh, EntityType::Cell, 0>;  
    using DataType = typename SGrid::DataType; 
    using FieldType = typename SGrid::BaseType; 
    using FieldLab = Block::FieldLab<typename FieldType::FieldType>; 
    using Stencil = typename FieldLab::StencilType;
   
    // define number of kernel calls to be executed in benchmarking
    size_t N = ((argc == 2) ? std::atoi(argv[1]) : 100);     
    
    // map the blocks & cells on domain [0,1]^3 [m] 
    SGrid sol(nblocks, block_cells);        // solution storage
    SGrid tmp(nblocks, block_cells);        // temporary storage

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

    // initialize grid by looping over blocks in the grid
    for (auto bf : sol) {
        fIC(*bf); 
    }

    // define numerical parameters based on chosen discretization
#ifdef USE_ACCUR                      
    const Stencil s(-2, 3, true);                
#else                                          
    const Stencil s(-1, 2, false);              
#endif /* USE_ACCUR */
    
    // setup lab 
    FieldLab flab; 
    flab.allocate(s, sol[0].getIndexRange()); 

    // get block field index functor for periodic block accessing
    auto findex = sol.getIndexFunctor(0);

    // warm-up call (any Laplacian kernel will do) 
    for (auto f : sol) {
        const FieldType &bf = *f; 
        const MIndex &bi = bf.getState().block_index; 
        sol.loadLab(bf, flab); 
        auto &tf = tmp[bi]; 
        LaplacianSecondOrder(flab, tf); 
    }

    // setup timer & required storage 
    Timer timer;
    std::vector<double> time(N);            // time per benchmark run 
    std::vector<double> performance(N);     // performance per benchmark run
    std::vector<double> roofCoords(8);      // coords. for roofline model

    // run the benchmark N times 
    for (size_t i = 0; i < N; ++i) 
    {
        // loop through blocks in the grid  
        for (auto f : sol) 
        {
            // reference fields & load data into flab object for current block
            const FieldType &bf = *f;  
            const MIndex &bi = bf.getState().block_index; 
            sol.loadLab(bf, flab); 
            auto &tf = tmp[bi];

            // benchmark selected Laplacian kernel
            timer.start();
#ifdef USE_ACCUR 
            LaplacianFourthOrder(flab, tf); 
#else
            LaplacianSecondOrder(flab, tf); 
#endif /* USE_ACCUR */
            time[i] = timer.stop(); 
        }
    }
   
    // pass benchmarking measurements to getRoofline template function 
    getRoofline<DataType>(time, performance, roofCoords);     
    
    // define percentiles to be selected
    const int i1 = N * 0.9; 
    const int i2 = N * 0.5; 
    const int i3 = N * 0.1; 

    // output 10th, 50th, and 90th percentiles
    printf("--------------------------------------------------------------\n");     
    printf("Ranked peak performance after %ld runs: \n", N); 
    printf("10th percentile:\t%f GFLOP/s\n", performance[i3]); 
    printf("50th percentile:\t%f GFLOP/s\n", performance[i2]); 
    printf("90th percentile:\t%f GFLOP/s\n", performance[i1]); 
    printf("--------------------------------------------------------------\n");      

    // write relevant data to file 
    std::ofstream results; 
    results.open("roof.txt"); 
    for (size_t i = 0; i < roofCoords.size(); ++i) {
        results << roofCoords[i] << "\n"; 
    }
    results.close(); 
    // execute plotting commands
    system("wd=$PWD; cd ../../src/; python3 plot.py; cd $wd");  

    return 0; 
}
