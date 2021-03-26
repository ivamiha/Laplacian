// File       : main.cpp    
// Created    : Tue Mar 16 2021  9:46:31 am CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: File with main function for testing Laplacian function
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include "Cubism/Block/FieldLab.h"
#include "Cubism/Grid/Cartesian.h"
#include "Cubism/Mesh/StructuredUniform.h"
#include "Cubism/Block/Field.h"
#include "Cubism/Block/FieldOperator.h"
#include "Cubism/Util/Timer.h"

#include "Laplacian2f.h"
#include "Laplacian2s.h"
#include "Laplacian4f.h"
#include "Laplacian4s.h"

#include <cstdio>

// enable 4th-order CDS, default 2nd-order CDS 
#define USE_ACCUR 
// enable flat indexing, default spatial indexing
#define USE_FLAT

using namespace Cubism;
using Util::Timer; 

int main(int argc, char *argv[]) 
{    
    // welcome & inform user what program configuration is running
    printf("--------------------------------------------------------------\n");
    printf("T E S T   L A P L A C I A N   C O M P U T E   K E R N E L\n\n");
    printf("Optimized Laplacian compute kernel utilizing CubismNova\n");
    printf("==============================================================\n");
#ifdef USE_ACCUR
    printf("Benchmarking 4th-order CDS implementation with "); 
#else
    printf("Benchmarking 2nd-order CDS implementation with "); 
#endif /* USE_ACCUR */
#ifdef USE_FLAT 
    printf("flat indexing.\n"); 
#else
    printf("spatial indexing.\n"); 
#endif /* USE_FLAT */  

    // identifers to be used for creating & managing 3D scalar fields
    using Mesh = Mesh::StructuredUniform<double, 3>;
    using Point = typename Mesh::PointType; 
    using MIndex = typename Mesh::MultiIndex; 
    using SGrid = Grid::Cartesian<float, Mesh, EntityType::Cell, 0>;  
    using DataType = typename SGrid::DataType; 
    using FieldType = typename SGrid::BaseType; 
    using FieldLab = Block::FieldLab<typename FieldType::FieldType>; 
    using Stencil = typename FieldLab::StencilType;
   
    // define number of kernel calls to be executed in benchmarking
    size_t N = ((argc == 2) ? std::atoi(argv[1]) : 100);     
    // define number of blocks & cells per block, input in each dimension
    const MIndex nblocks(1);                        
    const MIndex block_cells(32);                   
    // map the blocks & cells on domain [0,1]^3 [m] 
    SGrid sol(nblocks, block_cells);                // solution storage
    SGrid tmp(nblocks, block_cells);                // temporary storage
    const Point h = sol.getMesh().getCellSize(0);   // grid spacing [m]

    // function for writing ICs utilizing input field (block) within grid     
    auto fIC = [h](FieldType &b) {
        const DataType fac = 2 * M_PI;
        const Mesh &bm = *b.getState().mesh;
        // loop over cells in the block's mesh using cell index (ci) 
        for (auto &ci : bm[EntityType::Cell]) {
            const Point x = bm.getCoordsCell(ci); 
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
    
    // warm-up Laplacian call (any will do)
    for (auto f : sol) 
    {
        const FieldType &bf = *f; 
        const MIndex &bi = bf.getState().block_index; 
        sol.loadLab(bf, flab); 
        auto &tf = tmp[bi]; 
        Laplacian2f(flab, tf);
    }
    
    // setup timer & time simulation duration 
    Timer timer; 
    double t0 = timer.stop(); 
#ifdef USE_ACCUR 
    // loop through time 
    for (size_t n = 0; n < N; ++n) 
    {
        timer.start(); 
        // loop through blocks in the grid  
        for (auto f : sol) 
        {
            // reference fields & load data into flab object for current block
            const FieldType &bf = *f;  
            const MIndex &bi = bf.getState().block_index; 
            sol.loadLab(bf, flab); 
            auto &tf = tmp[bi];

            // benchmark 4th-order Laplacian CDS
            timer.start(); 
            #ifdef USE_FLAT
                Laplacian4f(flab, tf);
            #else 
                Laplacian4s(flab, tf);
            #endif /* USE_FLAT */
            t0 += timer.stop(); 
        }
    }
#else 
    // loop through time 
    for (size_t n = 0; n < N; ++n) 
    { 
        // loop through blocks in the grid
        for (auto f : sol) 
        {
            // reference fields & load data into flab object for current block
            const FieldType &bf = *f;  
            const MIndex &bi = bf.getState().block_index; 
            sol.loadLab(bf, flab);
            auto &tf = tmp[bi];
                                                                             
            // benchmark 2nd-order Laplacian CDS
            timer.start(); 
            #ifdef USE_FLAT
                Laplacian2f(flab, tf);
            #else
                Laplacian2s(flab, tf);
            #endif /* USE_FLAT */
            t0 += timer.stop();  
        }
    }
#endif /* USE_ACCUR */
    
    // compute & communicate average kernel execution time
    t0 /= N; 
    printf("After %ld runs, average kernel execution time:\t%e [s].\n", N, t0); 

    return 0; 
}
