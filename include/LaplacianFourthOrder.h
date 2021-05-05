// File       : LaplacianFourthOrder.h
// Created    : Tue Mar 23 2021  5:45:28 pm CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: Fourth-order Laplacian compute kernel
// Copyright 2021 ETH Zurich. All Rights Reserved.

#ifndef LAPLACIANFOURTHORDER_H
#define LAPLACIANFOURTHORDER_H

#include "sliceFourthOrder.h"

// enable ISPC DLP, default is CubismNova's flat indexing
#define USE_DLP 

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
 
    // extract mesh & relevant data
    const auto &bm = tmp.getState().mesh; 
    const auto h = bm->getCellSize(0);
    // compute inverse of grid spacing squared for given discretization
    const DataType ihx2 = 1.0 / (12 * h[0] * h[0]);  
    const DataType ihy2 = 1.0 / (12 * h[1] * h[1]); 
    const DataType ihz2 = 1.0 / (12 * h[2] * h[2]);

#ifdef USE_DLP
    // utilize ISPC-based implementation of fourth-order CDS
    using Index = typename MIndex::DataType; 
    // package inverse grid spacing squared values into array & create pointer
    const DataType ih2[3] = {ihx2, ihy2, ihz2}; 
    const DataType *pih2 = &ih2[0]; 
    // extract extents of Field and FieldLab which will be required 
    const size_t Nx = tmp.getIndexRange().getExtent()[0]; 
    const size_t Ny = tmp.getIndexRange().getExtent()[1]; 
    const Index  Nz = tmp.getIndexRange().getExtent()[2]; 
    const size_t NxNy = Nx * Ny; 
    const size_t Sx = sol.getIndexRange().getExtent()[0]; 
    const size_t SxSy = Sx * sol.getIndexRange().getExtent()[1]; 
    // get pointers to first slices in both FieldLab and Field 
    DataType *slice_sol = sol.getInnerData(); 
    const MIndex zero{0, 0, 0};
    DataType *slice_tmp = &tmp[zero];

    // loop over slowest moving index to cover all slices in domain
    for (Index iz = 0; iz < Nz; ++iz) {
        // extract relevant slices & process them in ISPC kernel
        DataType *psol = &slice_sol[iz * SxSy];
        DataType *ptmp = &slice_tmp[iz * NxNy]; 
        ispc::sliceFourthOrder(psol, ptmp, Nx, Ny, Sx, SxSy, pih2);
    }
#else 
    // utilize naive implementation of fourth-order CDS
    const MIndex ix{1, 0, 0}; 
    const MIndex iy{0, 1, 0};
    const MIndex iz{0, 0, 1};
    const MIndex iix{2, 0, 0};
    const MIndex iiy{0, 2, 0};
    const MIndex iiz{0, 0, 2};
    // loop over all elements in passed-in Field block  
    for (auto &i : tmp.getIndexRange()) {
        const DataType ddx = ihx2 * (16 * (sol[i + ix] + sol[i - ix]) 
                                - sol[i + iix] - sol[i - iix] - 30 * sol[i]);
        const DataType ddy = ihy2 * (16 * (sol[i + iy] + sol[i - iy])
                                - sol[i + iiy] - sol[i - iiy] - 30 * sol[i]);
        const DataType ddz = ihz2 * (16 * (sol[i + iz] + sol[i - iz]) 
                                - sol[i + iiz] - sol[i - iiz] - 30 * sol[i]); 
        tmp[i] = ddx + ddy + ddz;     
    }
#endif /* USE_DLP */
}

#endif // LAPLACIANFOURTHORDER_H
