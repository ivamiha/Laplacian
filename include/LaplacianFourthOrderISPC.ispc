// File       : LaplacianFourthOrderISPC.ispc
// Created    : Tue Apr 20 2021 05:42:48 PM (+0200)
// Author     : Ivan Mihajlovic Milin
// Description: Fourth-order Laplacian compute kernel ISPC implementation
// Copyright 2021 ETH Zurich. All Rights Reserved.
 
/**
 * @brief Sliced ISPC implementation of fourth-order Laplacian compute kernel
 * @param psol Pointer to current solution stored in FieldLab (use strides)
 * @param ptmp Pointer to temporary new solution storage in Field (use domain)
 * @param Nx Number of elements in computational domain in x-direction 
 * @param Ny Number of elements in computational domain in y-direction 
 * @param Sx Number of elements in FieldLab stride in x-direction 
 * @param SxSy FieldLab x-direction * y-direction strides, use to change slices
 * @param ih2 Inverse grid spacing squared to be used for CDS computation
 *
 * @rst Sliced ISPC implementation computes values of elements in a passed-in 
 * x-y plane (slice). Due to the problem being 3D, FD computations require 
 * both in-plane data and out-of-plane data (found in adjacent slices). SIMD 
 * behavior introduced in innermost (fastest) looping index. 
 * @endrst
 * */
export inline void LaplacianFourthOrderISPC(
                                        const uniform double* uniform psol,
                                        uniform double* uniform ptmp, 
                                        const uniform size_t Nx,
                                        const uniform size_t Ny,
                                        const uniform size_t Sx, 
                                        const uniform size_t SxSy,
                                        const uniform double ih2)
{
    for (uniform int j = 0; j < Ny; ++j) {
        foreach (i = 0 ... Nx) {
            // extract data for relevant indices in current SIMD lane
            double c   = psol[i + j * Sx];              // central cell
            double xp1 = psol[i + j * Sx + 1];          // x plus 1
            double xm1 = psol[i + j * Sx - 1];          // x minus 1
            double xp2 = psol[i + j * Sx + 2];          // x plus 2
            double xm2 = psol[i + j * Sx - 2];          // x minus 2
            double yp1 = psol[i + j * Sx + Sx];         // y plus 1
            double ym1 = psol[i + j * Sx - Sx];         // y minus 1
            double yp2 = psol[i + j * Sx + 2 * Sx];     // y plus 2
            double ym2 = psol[i + j * Sx - 2 * Sx];     // y minus 2
            double zp1 = psol[i + j * Sx + SxSy];       // z plus 1
            double zm1 = psol[i + j * Sx - SxSy];       // z minus 1
            double zp2 = psol[i + j * Sx + 2 * SxSy];   // z plus 2
            double zm2 = psol[i + j * Sx - 2 * SxSy];   // z minus 2

            // apply fourth-order CDS
            const double ddx = ih2 * (16 * (xp1 + xm1) - xp2 - xm2 - 30 * c);
            const double ddy = ih2 * (16 * (yp1 + ym1) - yp2 - ym2 - 30 * c);
            const double ddz = ih2 * (16 * (zp1 + zm1) - zp2 - zm2 - 30 * c);
            ptmp[i + j * Nx] = ddx + ddy + ddz;
        }
    }
}
