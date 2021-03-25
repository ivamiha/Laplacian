// File       : Laplacian4s.h
// Created    : Tue Mar 23 2021  6:12:27 pm CET (+0100) 
// Author     : Ivan Mihajlovic Milin
// Description: 4th-order spatial-indexing Laplacian kernel specifications
// Copyright 2021 ETH Zurich. All Rights Reserved.

#ifndef LAPLACIAN4S_H
#define LAPLACIAN4S_H

/** 
 * @brief Naive 4th-order spatial-indexing Laplacian kernel    
 * @param sol Current solution stored in data structure DataLab 
 * @param tmp Temporary storage of new solution in data structure Field
 *
 * @rst Naive Laplacian compute kernel utilizing 4th-order central
 * discretization scheme (CDS). Implemented with CubismNova's spatial indexing.  
 * @endrst 
 * */ 
template <typename FieldLab>
void Laplacian4s(FieldLab &sol, 
                 typename FieldLab::Data &tmp) 
{
    using DataType = typename FieldLab::Data::DataType; 
    using MIndex = typename FieldLab::MultiIndex;
    using Index = typename MIndex::DataType; 

    // apply spatial-indexing & 4th-order CDS 
    const auto extent = tmp.getIndexRange().getExtent(); 
    for (Index iz = 0; iz < extent[2]; ++iz) {
        for (Index iy = 0; iy < extent[1]; ++iy) {
            for (Index ix = 0; ix < extent[0]; ++ix) {
                tmp(ix, iy, iz) =   
                      16 * (sol(ix + 1, iy, iz) + sol(ix - 1, iy, iz)) 
                    + 16 * (sol(ix, iy + 1, iz) + sol(ix, iy - 1, iz)) 
                    + 16 * (sol(ix, iy, iz + 1) + sol(ix, iy, iz - 1))
                    - sol(ix + 2, iy, iz) - sol(ix - 2, iy, iz)
                    - sol(ix, iy + 2, iz) - sol(ix, iy -2, iz) 
                    - sol(ix, iy, iz + 2) - sol(ix, iy, iz - 2) 
                    - 90 * sol(ix, iy, iz);  
            }
        }
    }
}
#endif // LAPLACIAN4S_H
