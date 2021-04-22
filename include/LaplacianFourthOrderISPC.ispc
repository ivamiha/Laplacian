// File       : LaplacianFourthOrderISPC.ispc
// Created    : Tue Apr 20 2021 05:42:48 PM (+0200)
// Author     : Ivan Mihajlovic Milin
// Description: Fourth-order Laplacian compute kernel ISPC implementation
// Copyright 2021 ETH Zurich. All Rights Reserved.
 
 /**
  * @brief Sliced ISPC implementation of fourth-order Laplacian compute kernel
  * @param ps Pointer to first slice required for current kernel call
  * @param Nx Number of elements in computational domain in x-direction
  * @param Ny Number of elements in computational domain in y-direction
  * @param Sx Number of elements in FieldLab stride in x-direction
  * @param Sy Number of elements in FieldLab stride in y-direction
  *
  * @rst Sliced ISPC implementation computes values of elements in a passed-in 
  * x-y plane (slice). Due to the problem being 3D, FD computations require 
  * both in-plane data and out-of-plane data (found in adjacent slices). SIMD 
  * behavior introduced in innermost (fastest) looping index. 
  * @endrst
  * */
export inline void LaplacianFourthOrderISPC(
                                        const uniform Real* uniform const ps,
                                        const uniform size_t Nx, 
                                        const uniform size_t Ny, 
                                        const uniform size_t Sx, 
                                        const uniform size_t Sy)
{
    const uniform size_t Nhalo = 2; 
    // loop over the y-axis with uniform 
    // loop over the x-axis with foreach
    ;
}
