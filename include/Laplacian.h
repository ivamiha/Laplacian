// File       : Laplacian.h
// Created    : Tue Mar 16 2021  9:50:39 am CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: Laplacian compute kernel specifications
// Copyright 2021 ETH Zurich. All Rights Reserved.

#ifndef LAPLACIAN_H
#define LAPLACIAN_H 

#include <cstddef>
#include <cmath>

/** 
 * @brief Optimized Laplacian compute kernel
 * @param dlab Data structure storing current solution 
 * @param field Data structure for storing interim solution
 *
 * @rst Some modifications may be required for non-CubismNova data types as the
 * current state of the template function is not enitrely agnostic. 
 * @endrst 
 * */ 
template <typename T, typename U> 
void Laplacian(T &dlab, U &field) {
    // get size of entire 3D data passed in 
    size_t cells = field.size();
    // assuming domain is a perfect cube, obtain cells per side
    size_t N = cbrt(cells); 

    // TODO: extend discretization possibilies to higher order CDSs
    // TODO: eliminate manual periodic BC treatment once DataLab updated 
    // apply 2nd-order CDS
    for (size_t i = 1; i < N-1; ++i) {
        for (size_t j = 1; j < N-1; ++j) {
            for (size_t k = 1; k < N-1; ++k) {
                field(i,j,k) =  dlab(i+1,j,k) + dlab(i-1,j,k) + dlab(i,j+1,k) 
                                + dlab(i,j-1,k) + dlab(i,j,k+1) + dlab(i,j,k-1)
                                - 6*dlab(i,j,k); 
                                 
            }
        }
    }
}

#endif // LAPLACIAN_H
