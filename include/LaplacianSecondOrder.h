// File       : LaplacianSecondOrder.h
// Created    : Tue Mar 16 2021  9:50:39 am CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: Second-order Laplacian compute kernel
// Copyright 2021 ETH Zurich. All Rights Reserved.

#ifndef LAPLACIANSECONDORDER_H
#define LAPLACIANSECONDORDER_H

#include "ringBuff.h"
#include "LaplacianSecondOrderISPC.h"

// enable ISPC DLP, default is CubismNova's flat indexing
#define USE_DLP

/** 
 * @brief Second-order Laplacian compute kernel    
 * @param sol Current solution stored in data structure FieldLab 
 * @param tmp Temporary storage of new solution in data structure Field
 *
 * @rst Laplacian compute kernel utilizing second-order central discretization 
 * scheme (CDS). Data-level parallelism (DLP) implemented with ISPC (Intel® 
 * Implicit SPMD Program Compiler) and thread-level parallelism (TLP) 
 * implemented with OpenMP. 
 * @endrst 
 * */ 
template <typename FieldLab>
inline void LaplacianSecondOrder(FieldLab &sol, 
                                 typename FieldLab::FieldType &tmp)
{
    using DataType = typename FieldLab::FieldType::DataType; 
    using MIndex = typename FieldLab::MultiIndex;

    // extract mesh & relevant data
    const auto &bm = tmp.getState().mesh; 
    const auto h = bm->getCellSize(0);
    // compute inverse of grid spacing squared for given discretization
    const DataType ihx2 = 1.0 / (h[0] * h[0]);  
    const DataType ihy2 = 1.0 / (h[1] * h[1]); 
    const DataType ihz2 = 1.0 / (h[2] * h[2]);

#ifdef USE_DLP 
    // utilize ISPC-optimized second-order Laplacian compute kernel
    using Index = typename MIndex::DataType; 
    // extract extent of Field and FieldLab in each dimension
    const auto extent_f = tmp.getIndexRange().getExtent(); 
    const auto extent_l = sol.getIndexRange().getExtent();
    const size_t Nx = extent_f[0], Ny = extent_f[1], Nz = extent_f[2]; 
    const size_t Sx = extent_l[0], Sy = extent_l[1], Sz = extent_l[2];
    const size_t Nhalo = 1;
    // create ring buffer & load pointers to first two slices from FieldLab POV
    const size_t capacity = 3; 
    ringBuff_t* ring = ringCreate(capacity);
    const Index fx = (Sx - Nx) / 2;             // equal to x = 0 in Field 
    const Index fy = (Sy - Ny) / 2;             // equal to y = 0 in Field
    const Index fz = ((Sz - Nz) / 2) - Nhalo;   // equal to z = -Nhalo in Field 
    ringEnqueue(ring, &tmp(fx, fy, fz)); 
    ringEnqueue(ring, &tmp(fx, fy, fz + 1)); 
    
    // loop over slowest moving index to include all remaining slices
    // TODO: need to change limit once solve segmentation (Nz + Nhalo)
    for (Index iz = fz + 2; iz < Nz; ++iz) {
        // enqueue slice for latest z-coordinate
        // TODO: should index &sol rather than &tmp, but get segmentation
        ringEnqueue(ring, &tmp(fx, fy, iz)); 
        // call ISPC Laplacian kernel to process currently loaded slices
        ispc::LaplacianSecondOrderISPC(ring -> data, ring -> head, 
                                       ring -> capacity, Nx, Ny, Sx, Nhalo); 
    }
    
    // free dynamically allocated memory
    ringFree(ring); 
#else 
    // utilize naive second-order Laplacian compute kernel
    const MIndex ix{1, 0, 0}; 
    const MIndex iy{0, 1, 0};
    const MIndex iz{0, 0, 1};
    // loop over all elements in passed-in Field block  
    for (auto &i : tmp.getIndexRange()) {
        const DataType ddx = ihx2 * (sol[i + ix] - 2 * sol[i] + sol[i - ix]); 
        const DataType ddy = ihy2 * (sol[i + iy] - 2 * sol[i] + sol[i - iy]); 
        const DataType ddz = ihz2 * (sol[i + iz] - 2 * sol[i] + sol[i - iz]); 
        tmp[i] = ddx + ddy + ddz;         
    }
#endif /* USE_DLP */
}

#endif // LAPLACIANSECONDORDER_H
