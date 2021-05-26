// File       : sliceSecondOrder.ispc
// Created    : Tue Apr 20 2021 05:33:51 PM (+0200)
// Author     : Ivan Mihajlovic Milin
// Description: ISPC implementation for processing slice with second-order CDS
// Copyright 2021 ETH Zurich. All Rights Reserved.

/**
 * @brief Sliced ISPC implementation of second-order Laplacian compute kernel
 * @param psol Pointer to currrent solution stored in FieldLab (use strides)
 * @param ptmp Pointer to temporary new solution storage in Field (use domain)
 * @param Nx Number of elements in computational domain in x-direction
 * @param Ny Number of elements in computational domain in y-direction
 * @param Sx Number of elements in FieldLab stride in x-direction
 * @param SxSy FieldLab x-direction * y-direction strides, use to change slices
 * @param pih2 Pointer to array with inverse grid spacing squared values 
 *
 * @rst Sliced ISPC implementation computes values of elements in a passed-in 
 * x-y plane (slice). Due to the problem being 3D, FD computations require 
 * both in-plane data and out-of-plane data (found in adjacent slices). SIMD
 * behavor introduced in innermost (fastest) looping index.
 * @endrst    
 * */
export inline void sliceSecondOrder(
                                const uniform double* uniform psol,
                                uniform double* uniform ptmp,
                                const uniform size_t Nx, 
                                const uniform size_t Ny, 
                                const uniform size_t Sx,
                                const uniform size_t SxSy, 
                                const uniform double* uniform pih2)
{
    for (uniform int j = 0; j < Ny; ++j) {
        foreach (i = 0 ... Nx) {
            // extract data for relevant indices in current SIMD lane
            double c   = psol[i + j * Sx];              // central cell
            double xp1 = psol[i + j * Sx + 1];          // x plus 1
            double xm1 = psol[i + j * Sx - 1];          // x minus 1
            double yp1 = psol[i + j * Sx + Sx];         // y plus 1
            double ym1 = psol[i + j * Sx - Sx];         // y minus 1
            double zp1 = psol[i + j * Sx + SxSy];       // z plus 1
            double zm1 = psol[i + j * Sx - SxSy];       // z minus 1 

            // apply second-order CDS
            const double ddx = pih2[0] * (xp1 - 2 * c + xm1); 
            const double ddy = pih2[1] * (yp1 - 2 * c + ym1); 
            const double ddz = pih2[2] * (zp1 - 2 * c + zm1);
            ptmp[i + j * Nx] = ddx + ddy + ddz; 
        }
    } 
} 
