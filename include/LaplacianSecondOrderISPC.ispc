// File       : LaplacianSecondOrderISPC.ispc
// Created    : Tue Apr 20 2021 05:33:51 PM (+0200)
// Author     : Ivan Mihajlovic Milin
// Description: Second-order Laplacian compute kernel ISPC implementation
// Copyright 2021 ETH Zurich. All Rights Reserved.

inline double* getSlice(uniform double** uniform, const uniform size_t, 
                                        const uniform size_t, uniform size_t);

/**
 * @brief Sliced ISPC implementation of second-order Laplacian compute kernel
 * @param data Pointer to array of pointers to double datatype (storing slices)
 * @param head Position of head currently in the ring buffer
 * @param capacity Maximum size of the created ring buffer
 * @param Nx Number of elements in computational domain in x-direction
 * @param Ny Number of elements in computational domain in y-direction
 * @param Sx Number of elements in FieldLab stride in x-direction
 * @param Nhalo Number of halo cells on each end of the Field for used stencil
 *
 * @rst Sliced ISPC implementation computes values of elements in a passed-in 
 * x-y plane (slice). Due to the problem being 3D, FD computations require 
 * both in-plane data and out-of-plane data (found in adjacent slices). SIMD
 * behavor introduced in innermost (fastest) looping index.
 * @endrst    
 * */
export void LaplacianSecondOrderISPC(
                                uniform double** uniform data,
                                const uniform size_t head, 
                                const uniform size_t capacity,
                                const uniform size_t Nx, 
                                const uniform size_t Ny, 
                                const uniform size_t Sx,
                                const uniform size_t Nhalo)
{
    // TODO: change how access values (from Field size to FieldLab stride) 
    // once have fixed the FieldLab segmentation
    
    // extract slices
    const uniform double *s1 = getSlice(data, head, capacity, -1);
    const uniform double *s2 = getSlice(data, head, capacity,  0);
    const uniform double *s3 = getSlice(data, head, capacity,  1); 
    
    // perform looping
    for (uniform int j = 0; j < Ny; ++j) {
        foreach (i = 0 ... Nx) {
            // extract pointers for indices in current SIMD lane
            double *xc  = &s2[i + j*Nx];
            double *xp1 = &xc[ 1];
            double *xm1 = &xc[-1];
            double *yp1 = &xc[ Nx];
            double *ym1 = &xc[-Nx];
            double *zp1 = &s3[i + j*Nx];
            double *zm1 = &s1[i + j*Nx];

            // apply second-order CDS 
            // TODO: store in some sort of temporary buffer (tmp once fix)
            *xc = *xp1 + *xm1 + *yp1 + *ym1 + *zp1 + *zm1 - 6*(*xc); 
        }
    } 
} 



/** 
 * @brief Conveniently index the array where slices have been stored 
 * @param data Pointer to array of pointers to double data types
 * @param head Position of head currently in the ring buffer
 * @param capacity Maximum size of the created ring buffer 
 * @param index Index specifying which slice is being accessed 
 *
 * @rst Convenient indexing is provided by the computation of the "realIndex"
 * which allows indexing according to the computational stencil being used.
 * For example, indeces -2, -1, 0, 1, 2 for fourth-order CDS. 
 * @endrst
 * */
inline double* getSlice(uniform double** uniform data, 
                 const uniform size_t head, 
                 const uniform size_t capacity,
                 uniform size_t index)
{
    uniform size_t offset = (capacity - 1) / 2; 
    uniform size_t realIndex = (head + offset + index) % capacity;

    return data[realIndex];
}
