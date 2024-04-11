#include <stdio.h>
#include <math.h>
#include "functions_algdiff.h"


int main() {
    // Initialize ring buffer
    RingBuffer buffer_measurements;
    init_buffer(&buffer_measurements,WINDOW_LENGTH_G0);

    // Insert new values and perform convolution
    double new_value;
    
    int N = 10000;
    double yApprox[N];
    double y[N];
    double t[N];
    t[0]= 0.0f;
    for (int i = 0; i < N; i++) {  // Insert N values for demonstration
        t[i] = i*TS_AlgDiff;
        new_value = sin(t[i]) /* Get new value */;
        insert_value(&buffer_measurements, new_value);
        y[i] = new_value;
        yApprox[i] = convolution(&buffer_measurements, AlgDiff_g0, WINDOW_LENGTH_G0);
    }

    // Files to store signals
    FILE *fileResult = fopen("result.txt", "w");
    FILE *fileSignal = fopen("signal.txt", "w");
    // Write the array elements to the file
    for (int i = 0; i < N; i++) {
        fprintf(fileResult, "%f,%f\n", t[i],yApprox[i]);
        fprintf(fileSignal, "%f,%f\n", t[i],y[i]);
    }
    
    // Close the file
    fclose(fileResult);
    fclose(fileSignal);
    
    printf("Array elements saved to file successfully.\n");
    return 0;
}