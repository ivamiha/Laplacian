// File       : Laplacian2.h
// Created    : Tue Mar 16 2021  9:50:39 am CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: 2nd-order naive Laplacian compute kernel specifications
// Copyright 2021 ETH Zurich. All Rights Reserved.

#ifndef LAPLACIAN2_H
#define LAPLACIAN2_H 

#include <cstddef>
#include <cmath>

using namespace Cubism;

/** 
 * @brief Naive 2nd-order Laplacian compute kernel
 * @param sol Current solution stored in data structure DataLab 
 * @param tmp Temporary storage of new solution in data structure Field
 *
 * @rst Naive Laplacian compute kernel utilizing 2nd-order central
 * discretization scheme (CDS). Implemented with CubismNova.  
 * @endrst 
 * */ 
template <typename DataLab>
void Laplacian2(DataLab &sol, 
                typename DataLab::Data &tmp) 
{
    using IRange = Core::IndexRange<DataLab::Data::IndexRangeType::Dim>;
    using MIndex = typename IRange::MultiIndex; 
    using Index = typename MIndex::DataType;   
   
    for (Index &i : MIndex(0)) {
        for (Index &j : MIndex(1)) {
            for (Index &k : MIndex(2)) {
                tmp(i,j,k) = sol(i+1,j,k) + sol(i-1,j,k) + sol(i,j+1,k)
                                + sol(i,j-1,k) + sol(i,j,k+1) + sol(i,j,k-1)
                                - 6*sol(i,j,k);        
            }        
        }        
    }
}

#endif // LAPLACIAN2_H
