/*Testing for the convolve.c convolution 1d functions*/

#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#include "convolve.h"

int test_simple_array()
{
    // Test simple 1d array
    //  Input Signal 1 - A simple linear array of 5 elements
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

    if (fabs(result[1] - 6.0) < 1e-6 && fabs(result[2] - 9.0) < 1e-6 &&
        fabs(result[3] - 12.0) < 1e-6 && fabs(result[0] - 0.0) < 1e-6 &&
        fabs(result[4] - 0.0) < 1e-6) {
        printf("\nPASS: The mathematical core is functioning correctly.\n");
        return 0;
    }
    else {
        fprintf(
            stderr,
            "\nFAIL: The output did not match the expected"
            "mathematical values.\n"
        );
        return 1;
    }
}

int test_simple_array_nan_interpolate()
{
    // Test 1d array with NaN values
    //  Input Signal 1 - A linear array of 5 elements including a NaN element
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

    if (fabs(result[1] - 1.5) < 1e-9 && fabs(result[2] - 3.0) < 1e-9 &&
        fabs(result[3] - 4.5) < 1e-9 && fabs(result[0] - 0.0) < 1e-9 &&
        fabs(result[4] - 0.0) < 1e-9) {
        printf("\nPASS: The mathematical core is functioning correctly.\n");
        return 0;
    }
    else {
        fprintf(
            stderr,
            "\nFAIL: The output did not match the expected"
            "mathematical values.\n"
        );
        return 1;
    }
}

int test_bot0_nan_array()
{
    // Test bot=0 fallback behavior
    //  Input Signal 1 - A 1d array of 4 NaN elements
    size_t nx = 4;
    double f[] = {NAN, NAN, NAN, NAN};

    // Output array
    double result[4] = {0.0, 0.0, 0.0, 0.0};

    // A moving sum kernel of 3 elements
    size_t nkx = 3;
    double g[] = {1.0, 1.0, 1.0};

    // Execute the convolution function
    convolve1d_c(result, f, nx, g, nkx, true, true, 1);

    // When bot=0 (all NaN input), the convolved elements are set to NaN.
    // Edge elements (index 0 and 3) aren't written and remain 0.0 from
    // initialization.

    if (fabs(result[0] - 0.0) < 1e-9 && isnan(result[1]) && isnan(result[2]) &&
        fabs(result[3] - 0.0) < 1e-9) {
        printf("\nPASS: The core is handling only NaN input correctly.\n");
        return 0;
    }
    else {
        fprintf(stderr, "\nFAIL: The output did not match the expected values.\n");
        return 1;
    }
}

int main()
{
    // 1D convolution tests
    int failures = 0;
    printf("--- Starting 1D Convolution Tests ---\n");
    printf("Simple 1d array\n");
    if (test_simple_array() != 0) {
        failures++;
    }
    printf("Simple array containing NaN\n");
    if (test_simple_array_nan_interpolate() != 0) {
        failures++;
    }
    printf("Testing fallback when bot=0\n\n");
    if (test_bot0_nan_array() != 0) {
        failures++;
    }

    if (failures > 0) {
        return 1;
    }
    else {
        return 0;
    }
}
