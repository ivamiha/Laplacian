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
#include "Cubism/Compiler.h"

#include "Laplacian2f.h"
#include "Laplacian2s.h"
#include "Laplacian4f.h"
#include "Laplacian4s.h"
#include <iostream>

using namespace Cubism;

int main() 
{    
    // welcome user & request input for simulation configurations
    int method;             // user defined sim configuration 
    printf("\n------------------------------------------------------------\n");
    printf("T E S T   L A P L A C I A N   C O M P U T E   K E R N E L\n\n");
    printf("Optimized Laplacian compute kernel utilizing CubismNova\n");
    printf("============================================================\n");
    printf("Selection menu:\n");
    printf("1) benchmark 2nd-order CDS implementation;\n");
    printf("2) benchmark 4th-order CDS implementation;\n");
    printf("3) run OVS for 2nd-order CDS implementation;\n");
    printf("4) run OVS for 4th-order CDS implementation.\n");
    printf("Please input desired simulation configuration:\n");
    std::cin >> method; 
    // ensure that selected method is indeed an integer in program menu
    method = round(method);
    while (method < 1 || method > 4) {
        printf("Please ensure input number is correct:\n");
        std::cin >> method;
        method = round(method);
    }
    
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
    // map the blocks & cells onto scalar grids on default domain [0,1]^3 [m] 
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

    // setup lab
    // TODO: define s4 & implement user "choice" in an intelligent way 
    DataLab dlab; 
    const Stencil s2(-1, 2, false);
    dlab.allocate(s2, grid[0].getIndexRange()); 

    // get block field index functor for periodic block accessing
    auto findex = grid.getIndexFunctor(0); 

    // compute simulation constants for 2nd-order CDS stability conditions
    // TODO: add 4th order stability conditions
    double dt = (h[0] * h[0]) / (2 * D);
    double fac = (D * dt) / (h[0] * h[0]); 

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

            // apply Laplacian discretization 
            // TODO: select intelligently which function is applied
            Laplacian2f(dlab,tf);
        }

        // TODO: pass current solution into datalab
    } 

    return 0; 
}
