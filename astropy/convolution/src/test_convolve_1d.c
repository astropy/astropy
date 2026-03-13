/*Testing for the convolve.c convolution functions*/

#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "convolve.h"

void test_simple_array(){
    //Test simple 1d array
    // Input Signal 1 - A simple linear array of 5 elements
    size_t nx = 5;
    double f[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    
    // Output array
    double result[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

    // A moving sum kernel of 3 elements
    size_t nkx = 3;
    double g[] = {1.0, 1.0, 1.0};

    // Execute the convolution function
    convolve1d_c(result, f, nx, g, nkx, false, false, 1);

    // Verify the Output
    // For a kernel of size 3, the edges (index 0 and 4) are skipped by the inner loop.
    // Index 1 should be: (1.0 * 1.0) + (2.0 * 1.0) + (3.0 * 1.0) = 6.0
    // Index 2 should be: (2.0 * 1.0) + (3.0 * 1.0) + (4.0 * 1.0) = 9.0
    // Index 3 should be: (3.0 * 1.0) + (4.0 * 1.0) + (5.0 * 1.0) = 12.0
    
    printf("Computed Results:\n");
    for(int i = 0; i < nx; i++) {
        printf("result[%d] = %f\n", i, result[i]);
    }

    if (result[1] == 6.0 && result[2] == 9.0 && result[3] == 12.0) {
        printf("\n[PASS] The mathematical core is functioning correctly.\n");
        return;
    } else {
        printf("\n[FAIL] The output did not match the expected mathematical values.\n");
        return;
    }

}

int main() {
    //1D convolution tests
    printf("--- Starting 1D Convolution Tests ---\n");
    printf("Simple 1d array\n");
    void test_simple_array();

    return 0;
}