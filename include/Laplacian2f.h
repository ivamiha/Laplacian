// File       : Laplacian2f.h
// Created    : Tue Mar 16 2021  9:50:39 am CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: 2nd-order flat-indexing Laplacian kernel specifications
// Copyright 2021 ETH Zurich. All Rights Reserved.

#ifndef LAPLACIAN2f_H
#define LAPLACIAN2f_H 

using namespace Cubism;

/** 
 * @brief Naive 2nd-order flat-indexing Laplacian kernel    
 * @param sol Current solution stored in data structure DataLab 
 * @param tmp Temporary storage of new solution in data structure Field
 *
 * @rst Naive Laplacian compute kernel utilizing 2nd-order central
 * discretization scheme (CDS). Implemented with CubismNova's flat indexing.  
 * @endrst 
 * */ 
template <typename DataLab>
void Laplacian2f(DataLab &sol, 
                 typename DataLab::Data &tmp) 
{
    using DataType = typename DataLab::Data::DataType; 
    using Mesh = Mesh::StructuredUniform<DataType, 
                                         DataLab::Data::IndexRangeType::Dim>;
    using MIndex = typename Mesh::MultiIndex;           
 
    // apply flat-indexing & 2nd-order CDS
    const MIndex ix{1, 0, 0}; 
    const MIndex iy{0, 1, 0};
    const MIndex iz{0, 0, 1};
    // loop over all elements in passed-in Field block  
    for (auto &i : tmp.getIndexRange()) {
        tmp[i] =   sol[i + ix] + sol[i - ix] + sol[i + iy] + sol[i - iy];
                 + sol[i + iz] + sol[i - iz] - 6*sol[i];         
    }
}

#endif // LAPLACIAN2f_H
