/*Testing for the convolve.c convolution 1d functions*/

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
    convolve1d_c(result, f, nx, g, nkx, false, true, 1);

    // Verify the Output
    // For kernel of size 3, the edges (index 0 and 4) are skipped.
    // Index 1 should be: (1.0 * 1.0) + (2.0 * 1.0) + (3.0 * 1.0) = 6.0
    // Index 2 should be: (2.0 * 1.0) + (3.0 * 1.0) + (4.0 * 1.0) = 9.0
    // Index 3 should be: (3.0 * 1.0) + (4.0 * 1.0) + (5.0 * 1.0) = 12.0
    
    printf("Computed Results:\n");
    for(int i = 0; i < nx; i++) {
        printf("result[%d] = %f\n", i, result[i]);
    }

    if (fabs(result[1] - 6.0) < 1e-6 && fabs(result[2] - 9.0) < 1e-6 
            && fabs(result[3] - 12.0) < 1e-6 && fabs(result[0] - 0.0) < 1e-6 
                && fabs(result[4] - 0.0) < 1e-6) {
        printf("\n[PASS] The mathematical core is functioning correctly.\n");
        return;
    } else {
        printf("\n[FAIL] The output did not match the expected mathematical"
             "values.\n");
        return;
    }
}

void test_simple_array_nan_interpolate() {
    //Test 1d array with NaN values
    // Input Signal 1 - A linear array of 5 elements including a NaN element
    size_t nx = 5;
    double f[] = {1.0, 2.0, NAN, 4.0, 5.0};
    
    // Output array
    double result[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

    // A moving sum kernel of 3 elements
    size_t nkx = 3;
    double g[] = {1.0, 1.0, 1.0};

    // Execute the convolution function
    convolve1d_c(result, f, nx, g, nkx, true, true, 1);

    /* Verify the Output
    * Index 1 should be: top/bot = 3.0/2.0 = 1.5, where  top = (1.0 * 1.0) + 
    * (2.0 * 1.0), bot = 1.0 + 1.0 
    * Index 2: top/bot = 6.0/2.0 = 3.0, where top = (2.0 * 1.0) + (4.0 * 1.0) = 
    * 6.0, bot = 2.0
    * Index 3 : top/bot = 9.0/2.0 = 4.5, where top = (4.0 * 1.0) + (5.0 * 1.0) = 
    * 9.0, bot = 2.0
    */

    printf("Computed Results:\n");
    for(int i = 0; i < nx; i++) {
        printf("result[%d] = %f\n", i, result[i]);
    }

    if (fabs(result[1] - 1.5) < 1e-9 && fabs(result[2] - 3.0) < 1e-9 
            && fabs(result[3] - 4.5) < 1e-9 && fabs(result[0] - 0.0) < 1e-9 
                && fabs(result[4] - 0.0) < 1e-9) {
        printf("\n[PASS] The mathematical core is functioning correctly.\n");
        return;
    } else {
        printf("\n[FAIL] The output did not match the expected mathematical values.\n");
        return;
    }

    return;
}

void test__bot0_nan_array() {
    //Test bot=0 fallback behavior
    // Input Signal 1 - A 1d array of 4 NaN elements
    size_t nx = 4;
    double f[] = {NAN, NAN, NAN, NAN};
    
    // Output array
    double result[5] = {0.0, 0.0, 0.0, 0.0};

    // A moving sum kernel of 3 elements
    size_t nkx = 3;
    double g[] = {1.0, 1.0, 1.0};

    // Execute the convolution function
    convolve1d_c(result, f, nx, g, nkx, true, true, 1);

    // Output should be same as input for all elements convoluted since bot = 0.
    //Index 0 and 3 are left 0.0 in output array

    printf("Computed Results for array containing only NaN:\n");
    for(int i = 0; i < nx; i++) {
        printf("result[%d] = %f\n", i, result[i]);
    }

    if (fabs(result[0] - 0.0) < 1e-9 && isnan(result[1]) && isnan(result[2]) 
        && fabs(result[3] - 0.0) < 1e-9) {
        printf("\n[PASS] The core is handling only NaN input correctly.\n");
        return;
    } else {
        printf("\n[FAIL] The output did not match the expected values.\n");
        return;
    }

    return;
}

int main() {
    //1D convolution tests
    printf("--- Starting 1D Convolution Tests ---\n");
    printf("Simple 1d array\n");
    test_simple_array();
    printf("Simple array containing NaN\n");
    test_simple_array_nan_interpolate();
    printf("Testing fallback when bot=0\n\n");
    test__bot0_nan_array();
    return 0;
}