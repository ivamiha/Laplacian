// File       : main.cpp    
// Created    : Tue Mar 16 2021  9:46:31 am CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: File with main function for testing Laplacian function
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include "Cubism/Block/DataLab.h"
#include "Cubism/Grid/Cartesian.h"
#include "Cubism/Mesh/StructuredUniform.h"
#include "Cubism/Block/Field.h"
#include "Cubism/Common.h"

#include "Laplacian2f.h"
#include "Laplacian2s.h"
#include "Laplacian4f.h"
#include "Laplacian4s.h"

#include <iostream>
#include <assert.h>

// enable 4th-order CDS, default 2nd-order CDS 
#define _ACCUR_ 
// enable benchmarking, default no benchmarking
// #define _BENCH_  

using namespace Cubism;

int main() 
{    
    // welcome & inform user what program configuration is running
    printf("\n------------------------------------------------------------\n");
    printf("T E S T   L A P L A C I A N   C O M P U T E   K E R N E L\n\n");
    printf("Optimized Laplacian compute kernel utilizing CubismNova\n");
    printf("============================================================\n");
#ifdef _BENCH_
    printf("Benchmarking "); 
#else
    printf("Running "); 
#endif /* BENCH */
#ifdef _ACCUR_ 
    printf("4th-order CDS implementation.\n"); 
#else
    printf("2nd-order CDS implementation.\n"); 
#endif /* ACCUR */   
    
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
    using DataLab = Block::DataLab<typename FieldType::FieldType>; 
    using Stencil = typename DataLab::StencilType;

    // define number of blocks & cells per block
    const MIndex nblocks(8);                        // 8^3 = 512 blocks
    const MIndex block_cells(16);                   // 16^3 = 4096 cells/block 
    // map the blocks & cells on domain [0,1]^3 [m] 
    SGrid grid(nblocks, block_cells);               // solution storage
    SGrid temp(nblocks, block_cells);               // temporary storage
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
#ifdef _ACCUR_                      
    const Stencil s(-2, 3, true);               // 4th-order CDS stencil 
    double dt = (h[0] * h[0]) / (2 * D);        // TODO: 4th-order stability
    double fac = (D * dt) / (12 * h[0] * h[0]); // 4th-order Laplacian factor
#else                                          
    const Stencil s(-1, 2, false);              // 2nd-order CDS stencil
    double dt = (h[0] * h[0]) / (2 * D);        // 2nd-order stability 
    double fac = (D * dt) / (h[0] * h[0]);      // 2nd-order Laplacian factor 
#endif /* _ACCUR_ */

    // setup lab 
    DataLab dlab; 
    dlab.allocate(s, grid[0].getIndexRange()); 

    // get block field index functor for periodic block accessing
    auto findex = grid.getIndexFunctor(0); 

#ifdef _ACCUR_ 
    // loop through time 
    for (double t = 0.0; t < time; t += dt) 
    {
        // loop through blocks in the grid  
        for (auto f : grid) 
        {
            // reference fields & load data into dlab object for current block
            const FieldType &bf = *f;  
            dlab.loadData(bf.getState().block_index, findex);  
            auto &tf = temp[bf.getState().block_index];

            // apply 4th-order central Laplacian discretization 
            Laplacian4f(dlab,tf);
        }

        // TODO: pass current solution into datalab
    }
#else 
    // loop through time 
    for (double t = 0.0; t < time; t += dt) 
    { 
        // loop through blocks in the grid
        for (auto f : grid) 
        {
            // reference fields & load data into dlab object for current block
            const FieldType &bf = *f;  
            dlab.loadData(bf.getState().block_index, findex);  
            auto &tf = temp[bf.getState().block_index];
                                                                             
            // apply 2nd-order central Laplacian discretization 
            Laplacian2f(dlab,tf);
        }
                                                                             
        // TODO: pass current solution into datalab
    }
#endif /* _ACCUR_ */

    return 0; 
}
