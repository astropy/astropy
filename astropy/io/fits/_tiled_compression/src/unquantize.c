// This file is copied/adapted from the imcompress.c file in CFITSIO, and
// includes the unquantize_* functions. These are included here because
// they are not exposed in the CFITSIO library so we need to compile against
// it even if we link against the system CFITSIO library.

# include <stdlib.h>
# include <stdio.h>

#define ZERO_VALUE -2147483646 /* value used to represent zero-valued pixels */
#define N_RANDOM 10000  /* DO NOT CHANGE THIS;  used when quantizing real numbers */
#define MEMORY_ALLOCATION 113  /* Could not allocate memory */
#define SUBTRACTIVE_DITHER_2 2

float *fits_rand_value = 0;

int fits_init_randoms(void) {

/* initialize an array of random numbers */

    int ii;
    double a = 16807.0;
    double m = 2147483647.0;
    double temp, seed;

    if (fits_rand_value) {
       return(0);  /* array is already initialized */
    }

    /* allocate array for the random number sequence */
    /* THIS MEMORY IS NEVER FREED */
    fits_rand_value = calloc(N_RANDOM, sizeof(float));

    if (!fits_rand_value) {
	return(MEMORY_ALLOCATION);
    }

    /*  We need a portable algorithm that anyone can use to generate this
        exact same sequence of random number.  The C 'rand' function is not
	suitable because it is not available to Fortran or Java programmers.
	Instead, use a well known simple algorithm published here:
	"Random number generators: good ones are hard to find", Communications of the ACM,
        Volume 31 ,  Issue 10  (October 1988) Pages: 1192 - 1201
    */

    /* initialize the random numbers */
    seed = 1;
    for (ii = 0; ii < N_RANDOM; ii++) {
        temp = a * seed;
	seed = temp -m * ((int) (temp / m) );
	fits_rand_value[ii] = (float) (seed / m);
    }

    /*
    IMPORTANT NOTE: the 10000th seed value must have the value 1043618065 if the
       algorithm has been implemented correctly */

    if ( (int) seed != 1043618065) {
        printf("fits_init_randoms generated incorrect random number sequence");
	return(1);
    } else {
        return(0);
    }
}

/*--------------------------------------------------------------------------*/
int unquantize_i1r4(long row, /* tile number = row number in table  */
            unsigned char *input, /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int dither_method,    /* I - dithering method to use             */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            unsigned char tnull,  /* I - value of FITS TNULLn keyword if any */
            float nullval,        /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            float *output,        /* O - array of converted pixels           */
            int *status)          /* IO - error status                       */
/*
    Unquantize byte values into the scaled floating point values
*/
{
    long ii;
    int nextrand, iseed;

    if (!fits_rand_value)
       if (fits_init_randoms()) return(MEMORY_ALLOCATION);

    /* initialize the index to the next random number in the list */
    iseed = (int) ((row - 1) % N_RANDOM);
    nextrand = (int) (fits_rand_value[iseed] * 500);

    if (nullcheck == 0)     /* no null checking required */
    {
             for (ii = 0; ii < ntodo; ii++)
            {
/*
		if (dither_method == SUBTRACTIVE_DITHER_2 && input[ii] == ZERO_VALUE)
		    output[ii] = 0.0;
		else
*/
                    output[ii] = (float) (((double) input[ii] - fits_rand_value[nextrand] + 0.5) * scale + zero);

	        nextrand++;
	        if (nextrand == N_RANDOM) {
	            iseed++;
		    if (iseed == N_RANDOM) iseed = 0;
		    nextrand = (int) (fits_rand_value[iseed] * 500);
	        }
            }
    }
    else        /* must check for null values */
    {
            for (ii = 0; ii < ntodo; ii++)
            {
                if (input[ii] == tnull)
                {
                    *anynull = 1;
                    if (nullcheck == 1)
                        output[ii] = nullval;
                    else
                        nullarray[ii] = 1;
                }
                else
                {
/*
		    if (dither_method == SUBTRACTIVE_DITHER_2 && input[ii] == ZERO_VALUE)
		        output[ii] = 0.0;
		    else
*/
                        output[ii] = (float) (((double) input[ii] - fits_rand_value[nextrand] + 0.5) * scale + zero);
                }

	        nextrand++;
	        if (nextrand == N_RANDOM) {
	            iseed++;
		    if (iseed == N_RANDOM) iseed = 0;
	            nextrand = (int) (fits_rand_value[iseed] * 500);
                }
            }
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int unquantize_i2r4(long row, /* seed for random values  */
            short *input,         /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int dither_method,    /* I - dithering method to use             */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            short tnull,          /* I - value of FITS TNULLn keyword if any */
            float nullval,        /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            float *output,        /* O - array of converted pixels           */
            int *status)          /* IO - error status                       */
/*
    Unquantize short integer values into the scaled floating point values
*/
{
    long ii;
    int nextrand, iseed;

    if (!fits_rand_value)
       if (fits_init_randoms()) return(MEMORY_ALLOCATION);

    /* initialize the index to the next random number in the list */
    iseed = (int) ((row - 1) % N_RANDOM);
    nextrand = (int) (fits_rand_value[iseed] * 500);

    if (nullcheck == 0)     /* no null checking required */
    {
           for (ii = 0; ii < ntodo; ii++)
            {
/*
		if (dither_method == SUBTRACTIVE_DITHER_2 && input[ii] == ZERO_VALUE)
		    output[ii] = 0.0;
		else
*/
                    output[ii] = (float) (((double) input[ii] - fits_rand_value[nextrand] + 0.5) * scale + zero);

	        nextrand++;
	        if (nextrand == N_RANDOM) {
	            iseed++;
		    if (iseed == N_RANDOM) iseed = 0;
		    nextrand = (int) (fits_rand_value[iseed] * 500);
	        }
            }
    }
    else        /* must check for null values */
    {
            for (ii = 0; ii < ntodo; ii++)
            {
                if (input[ii] == tnull)
                {
                    *anynull = 1;
                    if (nullcheck == 1)
                        output[ii] = nullval;
                    else
                        nullarray[ii] = 1;
                }
                else
                {
/*
                    if (dither_method == SUBTRACTIVE_DITHER_2 && input[ii] == ZERO_VALUE)
		        output[ii] = 0.0;
		    else
*/
                        output[ii] = (float) (((double) input[ii] - fits_rand_value[nextrand] + 0.5) * scale + zero);
                }

	        nextrand++;
	        if (nextrand == N_RANDOM) {
	            iseed++;
		    if (iseed == N_RANDOM) iseed = 0;
		    nextrand = (int) (fits_rand_value[iseed] * 500);
	        }
             }
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int unquantize_i4r4(long row, /* tile number = row number in table    */
            int *input,      /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int dither_method,    /* I - dithering method to use             */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            int tnull,       /* I - value of FITS TNULLn keyword if any */
            float nullval,        /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            float *output,        /* O - array of converted pixels           */
            int *status)          /* IO - error status                       */
/*
    Unquantize int integer values into the scaled floating point values
*/
{
    long ii;
    int nextrand, iseed;

    if (fits_rand_value == 0)
       if (fits_init_randoms()) return(MEMORY_ALLOCATION);

    /* initialize the index to the next random number in the list */
    iseed = (int) ((row - 1) % N_RANDOM);
    nextrand = (int) (fits_rand_value[iseed] * 500);

    if (nullcheck == 0)     /* no null checking required */
    {
            for (ii = 0; ii < ntodo; ii++)
            {
                if (dither_method == SUBTRACTIVE_DITHER_2 && input[ii] == ZERO_VALUE)
		    output[ii] = 0.0;
		else
                    output[ii] = (float) (((double) input[ii] - fits_rand_value[nextrand] + 0.5) * scale + zero);

	        nextrand++;
	        if (nextrand == N_RANDOM) {
	            iseed++;
		    if (iseed == N_RANDOM) iseed = 0;
		    nextrand = (int) (fits_rand_value[iseed] * 500);
	        }
            }
    }
    else        /* must check for null values */
    {
            for (ii = 0; ii < ntodo; ii++)
            {
                if (input[ii] == tnull)
                {
                    *anynull = 1;
                    if (nullcheck == 1)
                        output[ii] = nullval;
                    else
                        nullarray[ii] = 1;
                }
                else
                {
                    if (dither_method == SUBTRACTIVE_DITHER_2 && input[ii] == ZERO_VALUE)
		        output[ii] = 0.0;
		    else
                        output[ii] = (float) (((double) input[ii] - fits_rand_value[nextrand] + 0.5) * scale + zero);
                }

	        nextrand++;
	        if (nextrand == N_RANDOM) {
	            iseed++;
		    if (iseed == N_RANDOM) iseed = 0;
		    nextrand = (int) (fits_rand_value[iseed] * 500);
	        }
            }
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int unquantize_i1r8(long row, /* tile number = row number in table  */
            unsigned char *input, /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int dither_method,    /* I - dithering method to use             */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            unsigned char tnull,  /* I - value of FITS TNULLn keyword if any */
            double nullval,        /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            double *output,        /* O - array of converted pixels           */
            int *status)          /* IO - error status                       */
/*
    Unquantize byte values into the scaled floating point values
*/
{
    long ii;
    int nextrand, iseed;

    if (!fits_rand_value)
       if (fits_init_randoms()) return(MEMORY_ALLOCATION);

    /* initialize the index to the next random number in the list */
    iseed = (int) ((row - 1) % N_RANDOM);
    nextrand = (int) (fits_rand_value[iseed] * 500);

    if (nullcheck == 0)     /* no null checking required */
    {
            for (ii = 0; ii < ntodo; ii++)
            {
/*
                if (dither_method == SUBTRACTIVE_DITHER_2 && input[ii] == ZERO_VALUE)
		    output[ii] = 0.0;
		else
*/
                    output[ii] = (double) (((double) input[ii] - fits_rand_value[nextrand] + 0.5) * scale + zero);

	        nextrand++;
	        if (nextrand == N_RANDOM) {
	            iseed++;
		    if (iseed == N_RANDOM) iseed = 0;
		    nextrand = (int) (fits_rand_value[iseed] * 500);
	        }
            }
    }
    else        /* must check for null values */
    {
            for (ii = 0; ii < ntodo; ii++)
            {
                if (input[ii] == tnull)
                {
                    *anynull = 1;
                    if (nullcheck == 1)
                        output[ii] = nullval;
                    else
                        nullarray[ii] = 1;
                }
                else
                {
/*
                    if (dither_method == SUBTRACTIVE_DITHER_2 && input[ii] == ZERO_VALUE)
		        output[ii] = 0.0;
		    else
*/
                        output[ii] = (double) (((double) input[ii] - fits_rand_value[nextrand] + 0.5) * scale + zero);
                }

	        nextrand++;
	        if (nextrand == N_RANDOM) {
	            iseed++;
		    if (iseed == N_RANDOM) iseed = 0;
		    nextrand = (int) (fits_rand_value[iseed] * 500);
	        }
            }
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int unquantize_i2r8(long row, /* tile number = row number in table  */
            short *input,         /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int dither_method,    /* I - dithering method to use             */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            short tnull,          /* I - value of FITS TNULLn keyword if any */
            double nullval,        /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            double *output,        /* O - array of converted pixels           */
            int *status)          /* IO - error status                       */
/*
    Unquantize short integer values into the scaled floating point values
*/
{
    long ii;
    int nextrand, iseed;

    if (!fits_rand_value)
       if (fits_init_randoms()) return(MEMORY_ALLOCATION);

    /* initialize the index to the next random number in the list */
    iseed = (int) ((row - 1) % N_RANDOM);
    nextrand = (int) (fits_rand_value[iseed] * 500);

    if (nullcheck == 0)     /* no null checking required */
    {
           for (ii = 0; ii < ntodo; ii++)
            {
/*
                if (dither_method == SUBTRACTIVE_DITHER_2 && input[ii] == ZERO_VALUE)
		    output[ii] = 0.0;
		else
*/
                    output[ii] = (double) (((double) input[ii] - fits_rand_value[nextrand] + 0.5) * scale + zero);

	        nextrand++;
	        if (nextrand == N_RANDOM) {
	            iseed++;
		    if (iseed == N_RANDOM) iseed = 0;
		    nextrand = (int) (fits_rand_value[iseed] * 500);
	        }
            }
    }
    else        /* must check for null values */
    {
            for (ii = 0; ii < ntodo; ii++)
            {
                if (input[ii] == tnull)
                {
                    *anynull = 1;
                    if (nullcheck == 1)
                        output[ii] = nullval;
                    else
                        nullarray[ii] = 1;
                }
                else
                {
/*                    if (dither_method == SUBTRACTIVE_DITHER_2 && input[ii] == ZERO_VALUE)
		        output[ii] = 0.0;
		    else
*/
                        output[ii] = (double) (((double) input[ii] - fits_rand_value[nextrand] + 0.5) * scale + zero);
                }

	        nextrand++;
	        if (nextrand == N_RANDOM) {
	            iseed++;
		    if (iseed == N_RANDOM) iseed = 0;
		    nextrand = (int) (fits_rand_value[iseed] * 500);
	        }
            }
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int unquantize_i4r8(long row, /* tile number = row number in table    */
            int *input,      /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int dither_method,    /* I - dithering method to use             */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            int tnull,       /* I - value of FITS TNULLn keyword if any */
            double nullval,        /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            double *output,        /* O - array of converted pixels           */
            int *status)          /* IO - error status                       */
/*
    Unquantize int integer values into the scaled floating point values
*/
{
    long ii;
    int nextrand, iseed;

    if (fits_rand_value == 0)
       if (fits_init_randoms()) return(MEMORY_ALLOCATION);

    /* initialize the index to the next random number in the list */
    iseed = (int) ((row - 1) % N_RANDOM);
    nextrand = (int) (fits_rand_value[iseed] * 500);

    if (nullcheck == 0)     /* no null checking required */
    {
            for (ii = 0; ii < ntodo; ii++)
            {
                if (dither_method == SUBTRACTIVE_DITHER_2 && input[ii] == ZERO_VALUE)
		    output[ii] = 0.0;
		else
                    output[ii] = (double) (((double) input[ii] - fits_rand_value[nextrand] + 0.5) * scale + zero);

	        nextrand++;
	        if (nextrand == N_RANDOM) {
	            iseed++;
		    if (iseed == N_RANDOM) iseed = 0;
		    nextrand = (int) (fits_rand_value[iseed] * 500);
	        }
            }
    }
    else        /* must check for null values */
    {
            for (ii = 0; ii < ntodo; ii++)
            {
                if (input[ii] == tnull)
                {
                    *anynull = 1;
                    if (nullcheck == 1)
                        output[ii] = nullval;
                    else
                        nullarray[ii] = 1;
                }
                else
                {
                    if (dither_method == SUBTRACTIVE_DITHER_2 && input[ii] == ZERO_VALUE)
		        output[ii] = 0.0;
		    else
                        output[ii] = (double) (((double) input[ii] - fits_rand_value[nextrand] + 0.5) * scale + zero);
                }

	        nextrand++;
	        if (nextrand == N_RANDOM) {
	            iseed++;
		    if (iseed == N_RANDOM) iseed = 0;
		    nextrand = (int) (fits_rand_value[iseed] * 500);
	        }
            }
    }

    return(*status);
}
