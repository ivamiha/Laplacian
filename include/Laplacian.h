// File       : Laplacian.h
// Created    : Tue Mar 16 2021  9:50:39 am CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: Laplacian compute kernel specifications
// Copyright 2021 ETH Zurich. All Rights Reserved.

#ifndef LAPLACIAN_HPP
#define LAPLACIAN_HPP 

// definition: perform Laplacian discretization on input cell values
template <typename T, typename U> 
void Laplacian(T &dlab, U &field);

#endif // LAPLACIAN_HPP
