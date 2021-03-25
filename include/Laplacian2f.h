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
 * @param fac Factor by which Laplacian discretization is multiplied, passed
 *        into function as a std::vector<DataType>
 *
 * @rst Naive Laplacian compute kernel utilizing 2nd-order central
 * discretization scheme (CDS). Implemented with CubismNova's flat indexing.  
 * @endrst 
 * */ 
template <typename FieldLab>
void Laplacian2f(FieldLab &sol, 
                 typename FieldLab::Data &tmp, 
                 std::vector<typename FieldLab::Data::DataType> &fac) 
{
    using DataType = typename FieldLab::Data::DataType; 
    using MIndex = typename FieldLab::MultiIndex;    
 
    // apply flat-indexing & 2nd-order CDS
    const MIndex ix{1, 0, 0}; 
    const MIndex iy{0, 1, 0};
    const MIndex iz{0, 0, 1};
    // loop over all elements in passed-in Field block  
    for (auto &i : tmp.getIndexRange()) {
        const DataType ddx = fac[0] * (sol[i + ix] - 2 * sol[i] + sol[i - ix]); 
        const DataType ddy = fac[1] * (sol[i + iy] - 2 * sol[i] + sol[i - iy]); 
        const DataType ddz = fac[2] * (sol[i + iz] - 2 * sol[i] + sol[i - iy]); 
        tmp[i] = ddx + ddy + ddz;         
    }
}

#endif // LAPLACIAN2F_H
