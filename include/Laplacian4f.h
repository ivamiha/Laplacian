// File       : Laplacian4f.h
// Created    : Tue Mar 23 2021  5:45:28 pm CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: 4th-order flat-indexing Laplacian kernel specifications
// Copyright 2021 ETH Zurich. All Rights Reserved.

#ifndef LAPLACIAN4F_H
#define LAPLACIAN4F_H 

/** 
 * @brief Naive 4th-order flat-indexing Laplacian kernel    
 * @param sol Current solution stored in data structure DataLab 
 * @param tmp Temporary storage of new solution in data structure Field
 *
 * @rst Naive Laplacian compute kernel utilizing 4th-order central
 * discretization scheme (CDS). Implemented with CubismNova's flat indexing.  
 * @endrst 
 * */ 
template <typename FieldLab>
void Laplacian4f(FieldLab &sol, 
                 typename FieldLab::FieldType &tmp)
{
    using DataType = typename FieldLab::FieldType::DataType; 
    using MIndex = typename FieldLab::MultiIndex;           
 
    // extract mesh & relevant data
    const auto &bm = tmp.getState().mesh; 
    const auto h = bm->getCellSize(0);
    // compute inverse of grid spacing squared  
    const DataType ihx2 = 1.0 / (12 * h[0] * h[0]);  
    const DataType ihy2 = 1.0 / (12 * h[1] * h[1]); 
    const DataType ihz2 = 1.0 / (12 * h[2] * h[2]);
    
    // apply flat-indexing & 4th-order CDS
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
}

#endif // LAPLACIAN4F_H
