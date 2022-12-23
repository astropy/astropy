int unquantize_i1r4(long row,
            unsigned char *input,         /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int dither_method,    /* I - which subtractive dither method to use */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            unsigned char tnull,          /* I - value of FITS TNULLn keyword if any */
            float nullval,        /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            float *output,        /* O - array of converted pixels           */
            int *status);          /* IO - error status                       */
int unquantize_i2r4(long row,
            short *input,         /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int dither_method,    /* I - which subtractive dither method to use */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            short tnull,          /* I - value of FITS TNULLn keyword if any */
            float nullval,        /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            float *output,        /* O - array of converted pixels           */
            int *status);          /* IO - error status                       */
int unquantize_i4r4(long row,
            int *input,      /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int dither_method,    /* I - which subtractive dither method to use */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            int tnull,       /* I - value of FITS TNULLn keyword if any */
            float nullval,        /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            float *output,        /* O - array of converted pixels           */
            int *status);          /* IO - error status                       */
int unquantize_i1r8(long row,
            unsigned char *input,         /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int dither_method,    /* I - which subtractive dither method to use */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            unsigned char tnull,          /* I - value of FITS TNULLn keyword if any */
            double nullval,        /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            double *output,        /* O - array of converted pixels           */
            int *status);          /* IO - error status                       */
int unquantize_i2r8(long row,
            short *input,         /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int dither_method,    /* I - which subtractive dither method to use */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            short tnull,          /* I - value of FITS TNULLn keyword if any */
            double nullval,        /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            double *output,        /* O - array of converted pixels           */
            int *status);          /* IO - error status                       */
int unquantize_i4r8(long row,
            int *input,      /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int dither_method,    /* I - which subtractive dither method to use */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            int tnull,       /* I - value of FITS TNULLn keyword if any */
            double nullval,        /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            double *output,        /* O - array of converted pixels           */
            int *status);          /* IO - error status                       */
