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
#include <math.h> 

namespace 
{
using namespace Cubism;

int OVS() 
{
    // problem specifications for 3D concentration fields with cell entities
    using Mesh = Mesh::StructuredUniform<double, 3>; 
    using SGrid = Grid::Cartesian<float, Mesh, EntityType::Cell, 0>; 

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

    // number of mesh refinements to be included in the OVS
    size_t n = 5;
    // grid refinement ratio 
    size_t r = 2; 
    // vector for storing computed norms (global errors) 
    std::vector<DataType> L2(n);
    // vector for storing actual orders of accuracy between grid refinements
    std::vector<DataType> OOA(n-2); 

    // loop through the number of desired meshes in the OVS
    for (size_t i = 0; i < n; ++i) 
    {
        // define mesh parameters for current index
        const MIndex nblocks(pow(r,i)); 
        const MIndex block_cells(32); 

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

        size_t cells = 0; 

        // loop through blocks in the grid
        for (auto f : sol) 
        {
            // reference fields & load data into flab for current block
            const FieldType &bf = *f; 
            const MIndex &bi = bf.getState().block_index;
            sol.loadLab(bf, flab); 
            auto &tf = tmp[bi]; 
            auto &ef = exa[bi];

            // apply Laplacian kernel discretization
            Laplacian2f(flab, tf); 
            
            // compute part of global error corresponding to current block
            for (auto &ci : tf.getIndexRange()) {
                // compute local error
                tf[ci] = (tf[ci] - ef[ci]) / ef[ci]; 
                // add local error contribution to global error 
                L2[i] += sqrt(tf[ci]*tf[ci]) * sqrt(tf[ci]*tf[ci]);
                ++cells; 
            }
        }

        // finalize L2-norm computation
        L2[i] = sqrt(L2[i]/cells);  
        printf("For mesh with %f blocks/dimension, L2 norm is: %f.\n", 
                                                           powf(r,i), L2[i]); 
    }
    
    // extract order of accuracy & report it
    for (size_t i = 0; i < n-2; ++i) {
        OOA[i] = log((L2[i+2] - L2[i+1]) / (L2[i+1] - L2[i])) / log(r); 
        printf("Computed actual order of accuracy: %f.\n", OOA[i]); 
    }

    // average computed actual orders of accuracy & round average to 1DP
    
    // TODO: simplify & do old-school OVS 
    // - no TEST required
    // - feed some info in from main function to OVS (generalize) 
    // - generate plots along with slope of nominal OOA

    return 1; 
}

TEST(ConvergenceRate, Laplacian) 
{
    int result = OVS(); 
    EXPECT_EQ(result, 1);  
}

} // namespace

int main(int argc, char *argv[]) 
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS(); 
}
