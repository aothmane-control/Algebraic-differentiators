#ifndef FUNCTIONS_ALGDIFF_H
#define FUNCTIONS_ALGDIFF_H


#include <stdio.h>
#include "g0.h"

// Structure to represent the ring buffer
typedef struct {
    double *buffer;
    int head;
    int length;
} RingBuffer;

/**
 * Initializes the given RingBuffer.
 *
 * @param rb The RingBuffer to be initialized.
 */
void init_buffer(RingBuffer *rb,int len);

/**
 * Inserts a value into the given RingBuffer.
 *
 * @param rb The RingBuffer to insert the value into.
 * @param value The value to be inserted.
 */
void insert_value(RingBuffer *rb, double value);

/**
 * Calculates the convolution of a RingBuffer with a given filter.
 *
 * @param rb The RingBuffer to convolve.
 * @param g The filter coefficients.
 * @param N The number of filter coefficients.
 * @return The result of the convolution.
 */
double convolution(RingBuffer *rb, const double *g, int N);

#endif /* FUNCTIONS_ALGDIFF_H */