// File       : main.cpp    
// Created    : Tue Mar 16 2021  9:46:31 am CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: File with main function for testing Laplacian function
// Copyright 2021 ETH Zurich. All Rights Reserved.

// XXX: [fabianw@mavt.ethz.ch; 2021-03-25] In the current master branch
// (fa84c47) I have removed the `DataLab` type and renamed it to `FieldLab`
// (which is the higher level physical abstraction from the Data.h class where
// it inherits from.  Data.h is only concerned with data storage/layout.  Is is
// more reasonable to talk of a FieldLab (laboratory for all kind of numerical
// operations) instead of a DataLab.  I have replaced the occurrences below in
// case you want to merge.
#include "Cubism/Block/Field.h"
#include "Cubism/Block/FieldLab.h"
#include "Cubism/Common.h"
#include "Cubism/Grid/Cartesian.h"
#include "Cubism/Mesh/StructuredUniform.h"
#include "Cubism/Util/Timer.h"

#include "Laplacian2f.h"
#include "Laplacian2s.h"
#include "Laplacian4f.h"
#include "Laplacian4s.h"

#include <cstdio>

// enable 4th-order CDS, default 2nd-order CDS 
//#define USE_ACCUR 
// enable flat indexing, default spatial indexing
#define USE_FLAT  

using namespace Cubism;
using Util::Timer; 

int main(int argc, char *argv[]) 
{    
    // welcome & inform user what program configuration is running
    printf("\n------------------------------------------------------------\n");
    printf("T E S T   L A P L A C I A N   C O M P U T E   K E R N E L\n\n");
    printf("Optimized Laplacian compute kernel utilizing CubismNova\n");
    printf("============================================================\n");
#ifdef USE_ACCUR
    printf("Benchmarking 4th-order CDS implementation with "); 
#else
    printf("Benchmarking 2nd-order CDS implementation with "); 
#endif /* USE_ACCUR */
#ifdef USE_FLAT 
    printf("flat indexing.\n"); 
#else
    printf("spatial indexing\n"); 
#endif /* USE_FLAT */   
    
    // initialize simulation variables to be used throughout
    double time = 50.0;     // simulation duration [s]
    double D = 0.00002;     // diffusivity of fluid [m²/s]

    // identifers to be used for creating & managing 3D scalar fields
    using Mesh = Mesh::StructuredUniform<double, 3>;
    using Point = typename Mesh::PointType; 
    using MIndex = typename Mesh::MultiIndex; 
    using SGrid = Grid::Cartesian<float, Mesh, EntityType::Cell, 0>;  
    using DataType = typename SGrid::DataType; 
    using FieldType = typename SGrid::BaseType;
    using FieldLab = Block::FieldLab<typename FieldType::FieldType>;
    using Stencil = typename FieldLab::StencilType;

    // define number of blocks & cells per block, input in each dimension   
    const MIndex nblocks((argc == 2) ? std::atoi(argv[1]) : 8);                        
    const MIndex block_cells(25);                   
    // map the blocks & cells on domain [0,1]^3 [m] 
    SGrid grid(nblocks, block_cells);               // solution storage
    SGrid tmp(nblocks, block_cells);                // temporary storage
    const Point h = grid.getMesh().getCellSize(0);  // grid spacing [m]

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
    for (auto f : grid) {
        fIC(*f); 
    }

    // define numerical parameters based on chosen discretization
#ifdef USE_ACCUR                      
    const Stencil s(-2, 3, true);               // 4th-order CDS stencil 
    double dt = (h[0] * h[0]) / (2 * D);        // TODO: 4th-order stability
    double fac = (D * dt) / (12 * h[0] * h[0]); // 4th-order Laplacian factor
#else                                          
    const Stencil s(-1, 2, false);              // 2nd-order CDS stencil
    double dt = (h[0] * h[0]) / (2 * D);        // 2nd-order stability 
    double fac = (D * dt) / (h[0] * h[0]);      // 2nd-order Laplacian factor 
#endif /* USE_ACCUR */

    // setup lab 
    FieldLab flab; 
    flab.allocate(s, grid[0].getIndexRange()); 

    // get block field index functor for periodic block accessing
    auto findex = grid.getIndexFunctor(0); 

    // setup timer & time simulation duration 
    Timer timer; 
    double t0 = timer.stop(); 
    timer.start(); 
#ifdef USE_ACCUR 
    // loop through time 
    for (double t = 0.0; t < time; t += dt) 
    {
        // loop through blocks in the grid  
        for (auto f : grid) 
        {
            // reference fields & load data into flab object for current block
            const FieldType &bf = *f;  
            flab.loadData(bf.getState().block_index, findex);  
            auto &tf = tmp[bf.getState().block_index];

            // apply 4th-order central Laplacian discretization
        #ifdef USE_FLAT
            Laplacian4f(flab, tf);
        #else 
            Laplacian4s(flab, tf);
        #endif /* USE_FLAT */
        }

        // finalize computation via point-wise operations
        //tmp *= fac; 
        //flab += tmp;  
    }
#else 
    // loop through time 
    for (double t = 0.0; t < time; t += dt) 
    { 
        // loop through blocks in the grid
        for (auto f : grid) 
        {
            // reference fields & load data into flab object for current block
            const FieldType &bf = *f;  
            flab.loadData(bf.getState().block_index, findex);  
            auto &tf = tmp[bf.getState().block_index];
                                                                             
            // apply 2nd-order central Laplacian discretization
        #ifdef USE_FLAT
            Laplacian2f(flab, tf);
        #else
            Laplacian2s(flab, tf);
        #endif /* USE_FLAT */  
        }
                                                                             
        // finalize computation via point-wise operations  
        //tmp *= fac; 
        //flab += tmp; 
    }
#endif /* USE_ACCUR */  
    t0 += timer.stop();

    printf("Simulation execution time:\t%e [s].\n", t0);

    return 0; 
}
