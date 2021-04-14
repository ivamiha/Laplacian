// File       : Laplacian2s.h
// Created    : Tue Mar 23 2021 10:37:56 am CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: 2nd-order spatial-indexing Laplacian kernel specifications
// Copyright 2021 ETH Zurich. All Rights Reserved.

#ifndef LAPLACIAN2S_H
#define LAPLACIAN2S_H 

/** 
 * @brief Naive 2nd-order spatial-indexing Laplacian kernel    
 * @param sol Current solution stored in data structure DataLab 
 * @param tmp Temporary storage of new solution in data structure Field
 *
 * @rst Naive Laplacian compute kernel utilizing 2nd-order central 
 * discretization scheme (CDS). Implemented with CubismNova's spatial indexing.  
 * @endrst 
 * */ 
template <typename FieldLab>
inline void Laplacian2s(FieldLab &sol, 
                        typename FieldLab::FieldType &tmp) 
{
    using DataType = typename FieldLab::FieldType::DataType; 
    using MIndex = typename FieldLab::MultiIndex; 
    using Index = typename MIndex::DataType;
   
    // extract mesh & relevant data
    const auto &bm = tmp.getState().mesh; 
    const auto h = bm->getCellSize(0);
    // compute inverse of grid spacing squared  
    const DataType ihx2 = 1.0 / (h[0] * h[0]);  
    const DataType ihy2 = 1.0 / (h[1] * h[1]); 
    const DataType ihz2 = 1.0 / (h[2] * h[2]);

    // apply spatial-indexing & 2nd-order CDS 
    const auto extent = tmp.getIndexRange().getExtent(); 
    for (Index iz = 0; iz < extent[2]; ++iz) {
        for (Index iy = 0; iy < extent[1]; ++iy) {
            for (Index ix = 0; ix < extent[0]; ++ix) {
                const DataType ddx = 
                            ihx2 * (sol(ix + 1, iy, iz) - 2 * sol(ix, iy, iz) 
                                  + sol(ix - 1, iy, iz));
                const DataType ddy = 
                            ihy2 * (sol(ix, iy + 1, iz) - 2 * sol(ix, iy, iz) 
                                  + sol(ix, iy - 1, iz));
                const DataType ddz = 
                            ihz2 * (sol(ix, iy, iz + 1) - 2 * sol(ix, iy, iz) 
                                  + sol(ix, iy, iz - 1)); 
                tmp(ix, iy, iz) = ddx + ddy + ddz; 
            }
        }
    }
}

#endif // LAPLACIAN2S_H
