// File       : ringBuff.h
// Created    : Wed Apr 21 2021 04:03:33 PM (+0200)
// Author     : Ivan Mihajlovic Milin
// Description: Bare-bone ring buffer implementation
// Copyright 2021 ETH Zurich. All Rights Reserved.

/** 
 * @brief Ring buffer (circular buffer) struct
 *
 * @rst Bare-bone implementation of a ring buffer. Several safety features 
 * which are commonly used (e.g. checking for NULL or full buffer) are 
 * ommitted. User never interacts with the struct, as only Laplacian kernels
 * access it -- therefore, this approach reduces struct memory footprint 
 * and speeds-up execution due to lack of if statements. 
 * @endrst 
 * */
// TODO: [imihajlovic@student.ethz.ch; Tue Apr 27 2021 04:30:02 PM (+0200)] 
// currently supports only double-precision, figure out how to extend to float
struct ringBuff {
    size_t head;        // current writing position    
    size_t capacity;    // maximum size of created ring buffer
    size_t offset;      // offset required to enable convenient indexing 
    double **data;      // pointer to array of pointers to doubles
}; 

typedef struct ringBuff ringBuff_t; 

/** 
 * @brief Allocate memory & initialize ring buffer of certain capacity
 * @param _capacity Capacity (maximum size) of ring buffer to be created 
 * */
ringBuff_t* ringCreate(size_t _capacity) {
    ringBuff_t *ring = new ringBuff_t; 

    ring -> head = 0; 
    ring -> capacity = _capacity; 
    ring -> offset = (_capacity - 1) / 2; 
    ring -> data = new double*[_capacity];

    return ring; 
}

/** 
 * @brief Enqueue/write item into buffer & advance head pointer 
 * @param ring Pointer to the created ringBuff_t struct
 * @param item Pointer to item which will be written to the buffer
 * */
void ringEnqueue(ringBuff_t *ring, double *item) {
    ring -> data[ring -> head] = item; 
    ring -> head = (ring -> head + 1) % ring -> capacity;
}

double* getSlice(ringBuff_t *ring, size_t index) {
    size_t realIndex = (ring -> head + ring -> offset + index) % 
                                                    ring -> capacity; 

    return (ring -> data[realIndex]); 
}

/** 
 * @brief Free memory that was dynamically allocated in ring buffer creation
 * @param ring Pointer to the created ringBuff_t struct
 * */
void ringFree(ringBuff_t* ring) {
    delete[] ring -> data; 
    delete ring;  
}
