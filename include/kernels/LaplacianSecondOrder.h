// File       : LaplacianSecondOrder.h
// Created    : Tue Mar 16 2021  9:50:39 am CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: Second-order Laplacian compute kernel
// Copyright 2021 ETH Zurich. All Rights Reserved.

#ifndef LAPLACIANSECONDORDER_H
#define LAPLACIANSECONDORDER_H

#include "sliceSecondOrder.h"

// enable ISPC DLP, default is CubismNova's flat indexing
//#define USE_DLP

/** 
 * @brief Second-order Laplacian compute kernel    
 * @param sol Current solution stored in data structure FieldLab 
 * @param tmp Temporary storage of new solution in data structure Field
 *
 * @rst Laplacian compute kernel utilizing second-order central discretization 
 * scheme (CDS). Data-level parallelism (DLP) implemented with ISPC (Intel® 
 * Implicit SPMD Program Compiler). 
 * @endrst 
 * */ 
template <typename FieldLab>
inline void LaplacianSecondOrder(FieldLab &sol, 
                                 typename FieldLab::FieldType &tmp)
{
    using DataType = typename FieldLab::FieldType::DataType; 
    using MIndex = typename FieldLab::MultiIndex;
    using Index = typename MIndex::DataType;

    // extract mesh & relevant data
    const auto &bm = tmp.getState().mesh; 
    const auto h = bm->getCellSize(0);
    // compute inverse of grid spacing squared and create pointer to array
    const DataType ih2[3] = {1.0 / (h[0] * h[0]), 
                             1.0 / (h[1] * h[1]), 
                             1.0 / (h[2] * h[2])}; 
    const DataType *pih2 = &ih2[0]; 
    // extract extents of Field and FieldLab which will be required
    const size_t Nx = tmp.getIndexRange().getExtent()[0]; 
    const size_t Ny = tmp.getIndexRange().getExtent()[1]; 
    const Index  Nz = tmp.getIndexRange().getExtent()[2]; 
    const size_t NxNy = Nx * Ny; 
    const size_t Sx = sol.getIndexRange().getExtent()[0]; 
    const size_t SxSy = Sx * sol.getIndexRange().getExtent()[1];
    // get pointers to first slices in both Field and FieldLab
    const MIndex zero{0, 0, 0}; 
    DataType *slice_tmp = &tmp[zero]; 
    DataType *slice_sol = sol.getInnerData(); 

#ifdef USE_DLP 
    // loop over slowest moving index to cover all slices in domain
    for (Index k = 0; k < Nz; ++k) { 
        // extract relevant slices & process them in ISPC kernel
        DataType *psol = &slice_sol[k * SxSy]; 
        DataType *ptmp = &slice_tmp[k * NxNy];  
        ispc::sliceSecondOrder(psol, ptmp, Nx, Ny, Sx, SxSy, pih2);
    }
#else 
    // loop over slowest moving index to cover all slices in domain 
    for (Index k = 0; k < Nz; ++k) {
        // extract relevant slices & proess them
        DataType *psol = &slice_sol[k * SxSy]; 
        DataType *ptmp = &slice_tmp[k * NxNy]; 
        // loop over slice indices
        for (size_t j = 0; j < Ny; ++j) {
            for (size_t i = 0; i < Nx; ++i) {
                // extract data for current indices
                const DataType c   = psol[i + j * Sx];          
                const DataType xp1 = psol[i + j * Sx + 1];      
                const DataType xm1 = psol[i + j * Sx - 1];      
                const DataType yp1 = psol[i + j * Sx + Sx];     
                const DataType ym1 = psol[i + j * Sx - Sx];     
                const DataType zp1 = psol[i + j * Sx + SxSy];   
                const DataType zm1 = psol[i + j * Sx - SxSy];   

                // apply second-order CDS
                const DataType ddx = pih2[0] * (xp1 - 2 * c + xm1); 
                const DataType ddy = pih2[1] * (yp1 - 2 * c + ym1); 
                const DataType ddz = pih2[2] * (zp1 - 2 * c + zm1);
                ptmp[i + j * Nx] = ddx + ddy + ddz; 
            }
        }
    }
#endif /* USE_DLP */
}

#endif // LAPLACIANSECONDORDER_H
