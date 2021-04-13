// File       : ConvergenceRateTest.cpp
// Created    : Wed Mar 31 2021 10:08:03 am CEST (+0200)
// Author     : Ivan Mihajlovic Milin
// Description: Convergence rate test for Laplacian kernels
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include "Cubism/Block/Field.h"
#include "Cubism/Block/FieldLab.h"
#include "Cubism/Grid/Cartesian.h"
#include "Cubism/Mesh/StructuredUniform.h"
#include "gtest/gtest.h"

#include "Laplacian2f.h" 
#include "Laplacian2s.h"
#include "Laplacian4f.h"
#include "Laplacian2s.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

namespace 
{
using namespace Cubism;

void OVS() 
{
    // problem specifications for 3D concentration fields with cell entities
    using Mesh = Mesh::StructuredUniform<double, 3>; 
    using SGrid = Grid::Cartesian<double, Mesh, EntityType::Cell, 0>; 

    // useful identifiers
    using MIndex = Mesh::MultiIndex; 
    using PointType = typename Mesh::PointType; 
    using DataType = typename SGrid::DataType; 
    using FieldType = typename SGrid::BaseType; 
    using FieldLab = Block::FieldLab<typename FieldType::FieldType>; 
    using Stencil = typename FieldLab::StencilType; 

    // function for writing ICs utilizing input field (block) within grid
    auto fIC = [](FieldType &b) {
        const DataType fac = 2 * M_PI; 
        const Mesh &bm = *b.getState().mesh; 
        for (auto &ci : bm[EntityType::Cell]) {
            const PointType x = bm.getCoordsCell(ci); 
            b[ci] = std::cos(fac * x[0]) * std::cos(fac * x[1]) *
                    std::sin(fac * x[2]); 
        }
    };
    
    // function for writing exact solution utilizing input field within grid
    auto fEX = [](FieldType &b) {
        const DataType fac = 2 * M_PI; 
        const Mesh &bm = *b.getState().mesh; 
        for (auto &ci : bm[EntityType::Cell]) {
            const PointType x = bm.getCoordsCell(ci); 
            b[ci] = -12 * M_PI * M_PI * std::cos(fac * x[0]) *
                    std::cos(fac * x[1]) * std::sin(fac * x[2]);
        }
    };

    // define fundamental variables for conducting OVS
    const size_t n = 6;                 // number of grid refinements
    const float r = 0.5;                // grid refinement ratio 
    std::vector<DataType> h(n,0);       // grid spacing (uniform for OVS)  
    std::vector<DataType> L2(n,0);      // vector for storing computed L2 norms
    std::vector<DataType> OOA(n-1,0);   // vector for storing actual OOAs

    printf("----------------------------------------------------\n"); 

    // loop through the number of desired meshes in the OVS
    for (size_t i = 0; i < n; ++i) 
    {
        // define mesh parameters for current index
        const MIndex nblocks(pow(1/r,i)); 
        const MIndex block_cells(4);
        const size_t cells1D = pow(1/r,i) * 4; 
        const size_t cells = pow(cells1D, 3);
        h[i] = 1.0 / cells1D; 
        printf("Running study with h = %f [m] & %ld cells.\n", h[i], cells); 

        // map grid with range [0,1] [m] onto the defined mesh
        SGrid sol(nblocks, block_cells);    // solution storage
        SGrid tmp(nblocks, block_cells);    // temporary storage
        SGrid exa(nblocks, block_cells);    // exact solution storage

        // write ICs into the scalar grid for solution storage 
        for (auto bf : sol) {
            fIC(*bf); 
        }
        // write exact solution into scalar grid for later use
        for (auto bf : exa) {
            fEX(*bf); 
        }

        // define stencil 
        const Stencil s(-1, 2, false); 

        FieldLab flab; 
        flab.allocate(s, sol[0].getIndexRange()); 

        auto findex = sol.getIndexFunctor(0);

        // loop through blocks in the grid
        for (auto f : sol) 
        {
            // reference fields & load data into flab for current block
            const FieldType &bf = *f; 
            const MIndex &bi = bf.getState().block_index;
            sol.loadLab(bf, flab);  // solution block fieldLab 
            auto &tf = tmp[bi];     // temporary block field
            auto &ef = exa[bi];     // exact solution block field

            // apply Laplacian kernel discretization
            Laplacian2f(flab, tf);
                
            // compute L2-norm (diff. numerical & exact solutions)
            for (auto &ci : bf.getIndexRange()) {
                // compute local error 
                DataType e_loc = tf[ci] - ef[ci];
                // add local error contribution to global error (L2)
                L2[i] += e_loc * e_loc;           
            } 
        }
        // finalize global error computation through normalization (L2)
        L2[i] = std::sqrt(L2[i] / cells);    
        printf("Computed L2-norm: %f.\n", L2[i]); 
        printf("----------------------------------------------------\n"); 
    }

    // extract order of accuracy for all refinements & report them
    for (size_t i = 0; i < n-1; ++i) {
        OOA[i] = (std::log10(L2[i]) - std::log10(L2[i+1])) / 
                                    (std::log10(h[i]) - std::log10(h[i+1]));     
        printf("Computed actual order of accuracy: %f.\n", OOA[i]);
    }
    printf("----------------------------------------------------\n");  

    // write relevant data to file
    std::ofstream results; 
    results.open("ooa.txt"); 
    for (size_t i = 0; i < n-1; ++i) {
        results << std::log10(h[i+1]) << "\t" << OOA[i] << "\n"; 
    }
    results.close(); 

    // execute plotting commands
    system("wd=$PWD; cd ../../test/; python3 plot.py; cd $wd"); 
}
} // namespace

int main() 
{
    // TODO: use function pointer (& flexible stencil selection) to choose
    // which function (discretization) will be executed in OVS
    OVS(); 

    return 0; 
}
