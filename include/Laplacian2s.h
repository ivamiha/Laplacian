// File       : Laplacian2s.h
// Created    : Tue Mar 23 2021 10:37:56 am CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: 2nd-order spatial-indexing Laplacian kernel specifications
// Copyright 2021 ETH Zurich. All Rights Reserved.

#ifndef LAPLACIAN2S_H
#define LAPLACIAN2S_H 

using namespace Cubism;

/** 
 * @brief Naive 2nd-order spatial-indexing Laplacian kernel    
 * @param sol Current solution stored in data structure DataLab 
 * @param tmp Temporary storage of new solution in data structure Field
 *
 * @rst Naive Laplacian compute kernel utilizing 2nd-order central 
 * discretization scheme (CDS). Implemented with CubismNova's spatial indexing.  
 * @endrst 
 * */ 
template <typename DataLab>
void Laplacian2s(DataLab &sol, 
                 typename DataLab::Data &tmp) 
{
    using DataType = typename DataLab::Data::DataType; 
    using Mesh = Mesh::StructuredUniform<DataType, 
                                         DataLab::Data::IndexRangeType::Dim>;
    using MIndex = typename Mesh::MultiIndex; 
    using Index = typename MIndex::DataType;
    
    // apply spatial-indexing & 2nd-order CDS 
    const auto extent = tmp.getIndexRange().getExtent(); 
    for (Index iz = 0; iz < extent[2]; ++iz) {
        for (Index iy = 0; iy < extent[1]; ++iy) {
            for (Index ix = 0; ix < extent[0]; ++ix) {
                tmp(ix, iy, iz) =   sol(ix + 1, iy, iz) + sol(ix - 1, iy, iz)
                                  + sol(ix, iy + 1, iz) + sol(ix, iy - 1, iz) 
                                  + sol(ix, iy, iz + 1) + sol(ix, iy, iz - 1)
                                  - 6*sol(ix, iy, iz);
            }
        }
    }
}

#endif // LAPLACIAN2S_H
