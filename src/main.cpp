// File       : main.cpp
// Created    : Wed May 05 2021 03:48:25 PM (+0200)
// Author     : Ivan Mihajlovic Milin
// Description: Laplacian kernel application to 3D Gray-Scott system
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include "Cubism/Block/FieldLab.h"
#include "Cubism/Block/Field.h"
#include "Cubism/Grid/Cartesian.h"
#include "Cubism/Mesh/StructuredUniform.h" 
#include "Cubism/IO/CartesianHDF.h"
#include "Cubism/Util/Timer.h"

#include "LaplacianSecondOrder.h"
#include "LaplacianFourthOrder.h"

// enable fourth-order CDS, default is second-order CDS
// #define USE_ACCUR

#include <cstdio>
#include <omp.h>

using namespace Cubism;
using Util::Timer;

int main(int argc, char *argv[])
{
    // identifiers for creating & managing meshes
    using IRange = Core::IndexRange<3>; 
    using MIndex = typename IRange::MultiIndex; 
    using Mesh = Mesh::StructuredUniform<double, IRange::Dim>; 
    using PointType = typename Mesh::PointType;

    // identifiers for creating & managing scalar fields 
    using SGrid = Grid::Cartesian<double, Mesh, EntityType::Cell, 0>; 
    using DataType = typename SGrid::DataType; 
    using FieldType = typename SGrid::BaseType; 
    using FieldLab = Block::FieldLab<typename FieldType::FieldType>;
    using Stencil = typename FieldLab::StencilType;
  
    // define some critical simulation parameters
    const MIndex nblocks((argc == 2) ? std::atoi(argv[1]) : 8); 
    const MIndex block_cells(32); 
    const size_t nthreads = 4; 

    // map the blocks & cells onto default [0,1]^3 [m] domain
    SGrid u(nblocks, block_cells);      // field for current u solution
    SGrid v(nblocks, block_cells);      // field for current v solution 
    SGrid utmp(nblocks, block_cells);   // temporary field for new u solution
    SGrid vtmp(nblocks, block_cells);   // temporary field for new v solution

    // define simulation & Gray-Scott variables
    double time = 100.0;            // simulation duration [s]
    double F = 0.04;                // feed-rate (permeability to U) [m^2]
    double k = 0.06;                // feed-rate minus permeability to V [m^2]
    double Du = 0.00002;            // diffusivity of U species [m^2/s]
    double Dv = 0.00001;            // diffusivity of V species [m^2/s]
    double Fo = 0.8;                // mass Fourier number []

    Timer timer;

    // function for writing ICs into input block field within grid
    auto IC = [](FieldType &bf) {
        const DataType fac = 2 * M_PI; 
        const Mesh &bm = *bf.getState().mesh;
        // loop over cells in block mesh using cell index
        for (auto &ci : bm[EntityType::Cell]) {
            const PointType x = bm.getCoordsCell(ci); 
            bf[ci] = std::fabs( std::cos(fac * x[0]) * std::cos(fac * x[1]) * 
                     std::sin(fac * x[2]) );
        }
    };

    // initialize grid using IC function for u field
    for (auto bf : u) {
        IC(*bf); 
    }
    // initialize grid of v field to be 1 - u (u + v = 1 at all points)
    for (auto f : v) {
        const FieldType &bf = *f;
        const MIndex &bi = bf.getState().block_index;
        auto uf = u[bi]; 
        auto vf = v[bi]; 
        // loop over cells in block field using cell index
        for (auto &ci : bf.getIndexRange()) {
            vf[ci] = 1.0 - uf[ci];
        }
    }
    // dump the ICs into HDF5 files using single-precision
    timer.start(); 
    IO::CartesianWriteHDF<float>("init_u", "U", u, 0); 
    IO::CartesianWriteHDF<float>("init_v", "V", v, 0); 
    double tw = timer.stop(); 

    // define stencil and allocate FieldLabs
    FieldLab ulab, vlab;
#ifdef USE_ACCUR
    const Stencil s(-2, 3, true); 
#else
    const Stencil s(-1, 2, false);
#endif /* USE_ACCUR */
    ulab.allocate(s, u[0].getIndexRange());
    vlab.allocate(s, v[0].getIndexRange()); 
    auto ufindex = u.getIndexFunctor(); 
    auto vfindex = v.getIndexFunctor();

    // define timestep as half of minimum dt from stability analysis
    const PointType h = u.getMesh().getCellSize(0);
    const double dt = Fo * h[0] * h[0] / (4 * std::max(Du, Dv));
    
    timer.start(); 
    // loop through time
    for (double t = 0.0; t < time; t += dt) 
    {
        // inform user about time in current loop
        printf("Time:\t%f\n", t);
        
        // process block field in the u & v scalar fields
        #pragma omp parallel for num_threads(nthreads)
        for (auto f : u) 
        {
            // reference Fields & load FieldLabs for current block field
            const FieldType &bf = *f; 
            const MIndex &bi = bf.getState().block_index;
            ulab.loadData(bi, ufindex); // block field lab with u values 
            vlab.loadData(bi, vfindex); // block field lab with v values
            auto &utf = utmp[bi];       // block field for tmp u storage
            auto &vtf = vtmp[bi];       // block field for tmp v storage
            auto &uf  = u[bi];          // block field for u species
            auto &vf  = v[bi];          // block field for v species
            
            // apply Laplacian discretization & store it in temporary fields
#ifdef USE_ACCUR
            LaplacianFourthOrder(ulab, utf); 
            LaplacianFourthOrder(vlab, vtf); 
#else 
            LaplacianSecondOrder(ulab, utf); 
            LaplacianSecondOrder(vlab, vtf);
#endif /* USE_ACCUR */
            // loop over cells in block field for pointwise operations
            for (auto &ci : bf.getIndexRange()) {
                // add reaction & feed-rate terms to tmp fields
                utf[ci] = Du * utf[ci] 
                        + F * (1 - uf[ci]) - uf[ci] * vf[ci] * vf[ci];
                vtf[ci] = Dv * vtf[ci] 
                        + uf[ci] * vf[ci] * vf[ci] - vf[ci] * (F + k); 
                // advance current solutions
                uf[ci] += dt * utf[ci]; 
                vf[ci] += dt * vtf[ci]; 
            } 
        }
    } 
    double ts = timer.stop(); 

    // dump solutions into HDF5 files using single-precision 
    timer.start(); 
    IO::CartesianWriteHDF<float>("solu", "sol U", u, 0);
    IO::CartesianWriteHDF<float>("solv", "sol V", v, 0); 
    tw += timer.stop(); 

    printf("ncells:\t\t%f\nwrite_time:\t%f [s]\nsim_time:\t%f [s]\n",
            std::pow(static_cast<int>(nblocks.prod())/2, 3),
            tw,
            ts); 
    
    return 0; 
}
