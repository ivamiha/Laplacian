// File       : Laplacian.h
// Created    : Tue Mar 16 2021  9:50:39 am CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: Laplacian compute kernel specifications
// Copyright 2021 ETH Zurich. All Rights Reserved.

#ifndef LAPLACIAN_HPP
#define LAPLACIAN_HPP 

#include "Cubism/Block/DataLab.h"
#include "Cubism/Block/Field.h"
#include "Cubism/Common.h"
#include "Cubism/Compiler.h"

#include <cstddef>

using namespace Cubism;

// definition: perform Laplacian discretization on input cell values
template <typename T, typename U> 
void Laplacian(T &dlab, U &field) {
    const ptrdiff_t N = 25;
    // apply 2nd-order CDS
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            for (size_t k = 0; k < N; ++k) {
                field(i,j,k) =  dlab(i+1,j,k) + dlab(i-1,j,k) + dlab(i,j+1,k)
                                + dlab(i,j-1,k) + dlab(i,j,k+1) + dlab(i,j,k-1)
                                - 6*dlab(i,j,k);
            }
        }
    }
}

#endif // LAPLACIAN_HPP
