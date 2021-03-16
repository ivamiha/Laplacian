// File       : Laplacian.cpp   
// Created    : Tue Mar 16 2021  9:48:33 am CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: Laplacian compute kernel
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include <cstddef>

template <typename T, typename U>
void Laplacian(T &dlab, U &field) {
    double N = 1000;
    
    // apply 2nd-order CDS
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            for (size_t k = 0; k < N; ++k) {
                 field(i,j,k) = dlab(i+1,j,k) + dlab(i-1,j,k) + dlab(i,j+1,k) 
                                + dlab(i,j-1,k) + dlab(i,j,k+1) + dlab(i,j,k-1) 
                                - 6*dlab(i,j,k);  
            }
        }
    } 
}
