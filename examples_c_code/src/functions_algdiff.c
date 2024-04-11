#include "functions_algdiff.h"
#include <stdlib.h>

// Function to initialize the ring buffer
void init_buffer(RingBuffer *rb,int len) {
    rb->buffer = malloc(len * sizeof(double));
    rb->head = 0;
    rb->length = len;
}

// Function to insert a new value into the buffer
void insert_value(RingBuffer *rb, double value) {
    rb->buffer[rb->head] = value;
    rb->head = (rb->head + 1) % rb->length;
}

// Function to perform discrete convolution with the given 'g' array
double convolution(RingBuffer *rb, const double *g, int N) {
    double result = 0.0;
    int i, j;

    // Perform convolution
    for (i = 0; i < N; i++) {
        j = (rb->head - i - 1 + rb->length) % rb->length;
        result += rb->buffer[j] * g[i];
    }

    return result;
}