// File       : LaplacianFourthOrder.h
// Created    : Tue Mar 23 2021  5:45:28 pm CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: Fourth-order Laplacian compute kernel
// Copyright 2021 ETH Zurich. All Rights Reserved.

#ifndef LAPLACIANFOURTHORDER_H
#define LAPLACIANFOURTHORDER_H

#include "sliceFourthOrder.h"

// enable ISPC DLP, default is CubismNova's flat indexing
//#define USE_DLP 

/** 
 * @brief Fourth-order Laplacian compute kernel    
 * @param sol Current solution stored in data structure FieldLab
 * @param tmp Temporary storage of new solution in data structure Field
 *
 * @rst Laplacian compute kernel utilizing fourth-order central discretization 
 * scheme (CDS). Data-level parallelism (DLP) implemented with ISPC (Intel® 
 * Implicit SPMD Program Compiler). 
 * @endrst 
 * */ 
template <typename FieldLab>
inline void LaplacianFourthOrder(FieldLab &sol, 
                                 typename FieldLab::FieldType &tmp)
{
    using DataType = typename FieldLab::FieldType::DataType; 
    using MIndex = typename FieldLab::MultiIndex;           
    using Index = typename MIndex::DataType;
 
    // extract mesh & relevant data
    const auto &bm = tmp.getState().mesh; 
    const auto h = bm->getCellSize(0);
    // compute inverse of grid spacing squared and create pointer to array
    const DataType ih2[3] = {1.0 / (12 * h[0] * h[0]),
                             1.0 / (12 * h[1] * h[1]), 
                             1.0 / (12 * h[2] * h[2])};
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
        ispc::sliceFourthOrder(psol, ptmp, Nx, Ny, Sx, SxSy, pih2);
    }
#else 
    // loop over slowest moving index to cover all slices in domain
    for (Index k = 0; k < Nz; ++k) {
        // extract relevant slices & process them in ISPC kernel
        DataType *psol = &slice_sol[k * SxSy]; 
        DataType *ptmp = &slice_tmp[k * NxNy]; 
        // loop over slice indices 
        for (size_t j = 0; j < Ny; ++j) {
            for (size_t i = 0; i < Nx; ++i) {
                // extract data for relevant indices 
                const DataType c   = psol[i + j * Sx];            
                const DataType xp1 = psol[i + j * Sx + 1];        
                const DataType xm1 = psol[i + j * Sx - 1];        
                const DataType xp2 = psol[i + j * Sx + 2];        
                const DataType xm2 = psol[i + j * Sx - 2];        
                const DataType yp1 = psol[i + j * Sx + Sx];       
                const DataType ym1 = psol[i + j * Sx - Sx];       
                const DataType yp2 = psol[i + j * Sx + 2 * Sx];   
                const DataType ym2 = psol[i + j * Sx - 2 * Sx];   
                const DataType zp1 = psol[i + j * Sx + SxSy];     
                const DataType zm1 = psol[i + j * Sx - SxSy];     
                const DataType zp2 = psol[i + j * Sx + 2 * SxSy]; 
                const DataType zm2 = psol[i + j * Sx - 2 * SxSy]; 

                // apply fourth-order CDS
                const DataType ddx = pih2[0] * (16 * (xp1 + xm1) 
                                                        - xp2 - xm2 - 30 * c);
                const DataType ddy = pih2[1] * (16 * (yp1 + ym1) 
                                                        - yp2 - ym2 - 30 * c);
                const DataType ddz = pih2[2] * (16 * (zp1 + zm1)
                                                        - zp2 - zm2 - 30 * c);
                ptmp[i + j * Nx] = ddx + ddy + ddz; 
            }
        }
    }
#endif /* USE_DLP */
}

#endif // LAPLACIANFOURTHORDER_H
