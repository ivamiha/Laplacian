// File       : Laplacian2f.h
// Created    : Tue Mar 16 2021  9:50:39 am CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: 2nd-order flat-indexing Laplacian kernel specifications
// Copyright 2021 ETH Zurich. All Rights Reserved.

#ifndef LAPLACIAN2F_H
#define LAPLACIAN2F_H 

/** 
 * @brief Naive 2nd-order flat-indexing Laplacian kernel    
 * @param sol Current solution stored in data structure DataLab 
 * @param tmp Temporary storage of new solution in data structure Field
 *
 * @rst Naive Laplacian compute kernel utilizing 2nd-order central
 * discretization scheme (CDS). Implemented with CubismNova's flat indexing.  
 * @endrst 
 * */ 
template <typename FieldLab>
void Laplacian2f(FieldLab &sol, 
                 typename FieldLab::FieldType &tmp)
{
    using DataType = typename FieldLab::FieldType::DataType; 
    using MIndex = typename FieldLab::MultiIndex;

    // extract mesh & relevant data
    const auto &bm = tmp.getState().mesh; 
    const auto h = bm->getCellSize(0);
    // compute inverse of grid spacing squared  
    const DataType ihx2 = 1.0 / (h[0] * h[0]);  
    const DataType ihy2 = 1.0 / (h[1] * h[1]); 
    const DataType ihz2 = 1.0 / (h[2] * h[2]);

    // apply flat-indexing & 2nd-order CDS   
    const MIndex ix{1, 0, 0}; 
    const MIndex iy{0, 1, 0};
    const MIndex iz{0, 0, 1};
    // loop over all elements in passed-in Field block  
    for (auto &i : tmp.getIndexRange()) {
        const DataType ddx = ihx2 * (sol[i + ix] - 2 * sol[i] + sol[i - ix]); 
        const DataType ddy = ihy2 * (sol[i + iy] - 2 * sol[i] + sol[i - iy]); 
        const DataType ddz = ihz2 * (sol[i + iz] - 2 * sol[i] + sol[i - iy]); 
        tmp[i] = ddx + ddy + ddz;         
    }
}

#endif // LAPLACIAN2F_H
