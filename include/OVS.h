// File       : OVS.h 
// Created    : Thu Mar 18 2021 10:02:24 am CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: Order verification study specifications
// Copyright 2021 ETH Zurich. All Rights Reserved.

#ifndef OVS_H
#define OVS_H 

#include "Laplacian.h"

#include <cstddef>
#include <cmath>

/** 
 * @brief Order verification study template function
 * @param dlab Data structure for storing current solution 
 * @param field Data structure for storing interim solution 
 *
 * @rst Order verification study carried out through the method of manufactured
 * solutions (MMS). Considering both temporal & spatial discretizations. 
 * @endrst
 * */
template<typename T, typename U>
void OVS(T &dlab, U &field) {
    // TODO: polish up the function footprint
    ;
}

#endif // OVS_H
