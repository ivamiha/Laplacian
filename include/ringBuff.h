// File       : ringBuff.h
// Created    : Wed Apr 21 2021 04:03:33 PM (+0200)
// Author     : Ivan Mihajlovic Milin
// Description: Bare-bone ring buffer implementation
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include <cstdlib>
#include <cstddef>

/** 
 * @brief Ring buffer (circular buffer) structure template
 * @param head Current head position (write pointer)
 * @param tail Current tail position (read pointer)
 * @param capacity Maximum size of the ringBuff (array size)
 * @param data Pointer to an array where ringBuff data is stored
 *
 * @rst Bare-bone implementation of a ring buffer. Several "safety" features 
 * commonly included (e.g. checking for NULL full buffer) are ommitted, since
 * they require if statements which slow down execution. Only Laplacian kernels
 * utilize the ring buffer, such that the user is never exposed to them.  
 * @endrst 
 * */
template <typename DataType>
struct ringBuff {
    size_t head;        
    size_t tail;        
    size_t capacity;    
    DataType* data;     
}; 

// create global definition of ringBuff_t
typedef struct ringBuff<double> ringBuff_t; 

/** 
 * @brief Allocate memory & initialize ring buffer of certain capacity
 * @param _capacity Capacity (maximum size) of ring buffer to be created 
 * */
ringBuff_t* ringCreate(size_t _capacity) {
    // allocate memory of size of the ringBuff struct
    ringBuff_t* ring = (ringBuff_t*)malloc(sizeof(ringBuff_t)); 

    // initialize member varaibles of the ringBuff struct
    ring -> tail = 0; 
    ring -> head = 0; 
    ring -> capacity = _capacity; 
    ring -> data = (double*)malloc(_capacity * sizeof(double));

    return ring; 
}

/** 
 * @brief Enqueue/write item into buffer & advance head pointer 
 * @param ring Pointer to the created ringBuff_t structure
 * @param item Data to be written to the buffer 
 * */
void ringEnqueue(ringBuff_t* ring, double item) {
    // add item to buffer
    ring -> data[ring -> head] = item; 
    // advance head
    ring -> head = (ring -> head + 1) % ring -> capacity;
}

/** 
 * @brief Advance tail pointer to dequeue (forget) item first added to buffer
 * @param ring Pointer to the created ringBuff_t structure 
 * */
void ringDequeue(ringBuff_t* ring) {
    // advance tail
    ring -> tail = (ring -> tail + 1) % ring -> capacity; 
}

/** 
 * @brief Free memory that was dynamically allocated in ring buffer creation
 * @param ring Pointer to the created ringBuff_t structure
 * */
void ringFree(ringBuff_t* ring) {
    // free the dynamically allocated memory for the data array
    free(ring -> data); 
    // free the dynamically allocated memory for the ringBuff struct
    free(ring); 
}
