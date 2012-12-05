# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include <ctype.h>
# include <time.h>
# include "fitsio2.h"

#define NULL_VALUE -2147483647 /* value used to represent undefined pixels */

/* nearest integer function */
# define NINT(x)  ((x >= 0.) ? (int) (x + 0.5) : (int) (x - 0.5))

/* special quantize level value indicates that floating point image pixels */
/* should not be quantized and instead losslessly compressed (with GZIP) */
#define NO_QUANTIZE 9999


/* string array for storing the individual column compression stats */
char results[999][60];
float trans_ratio[999];

float *fits_rand_value = 0;

int imcomp_write_nocompress_tile(fitsfile *outfptr, long row, int datatype, 
    void *tiledata, long tilelen, int nullcheck, void *nullflagval, int *status);
int imcomp_convert_tile_tshort(fitsfile *outfptr, void *tiledata, long tilelen,
    int nullcheck, void *nullflagval, int nullval, int zbitpix, double scale,
    double zero, double actual_bzero, int *intlength, int *status);
int imcomp_convert_tile_tushort(fitsfile *outfptr, void *tiledata, long tilelen,
    int nullcheck, void *nullflagval, int nullval, int zbitpix, double scale,
    double zero, int *intlength, int *status);
int imcomp_convert_tile_tint(fitsfile *outfptr, void *tiledata, long tilelen,
    int nullcheck, void *nullflagval, int nullval, int zbitpix, double scale,
    double zero, int *intlength, int *status);
int imcomp_convert_tile_tuint(fitsfile *outfptr, void *tiledata, long tilelen,
    int nullcheck, void *nullflagval, int nullval, int zbitpix, double scale,
    double zero, int *intlength, int *status);
int imcomp_convert_tile_tbyte(fitsfile *outfptr, void *tiledata, long tilelen,
    int nullcheck, void *nullflagval, int nullval, int zbitpix, double scale,
    double zero, int *intlength, int *status);
int imcomp_convert_tile_tsbyte(fitsfile *outfptr, void *tiledata, long tilelen,
    int nullcheck, void *nullflagval, int nullval, int zbitpix, double scale,
    double zero, int *intlength, int *status);
int imcomp_convert_tile_tfloat(fitsfile *outfptr, long row, void *tiledata, long tilelen,
    long tilenx, long tileny, int nullcheck, void *nullflagval, int nullval, int zbitpix,
    double scale, double zero, int *intlength, int *flag, double *bscale, double *bzero,int *status);
int imcomp_convert_tile_tdouble(fitsfile *outfptr, long row, void *tiledata, long tilelen,
    long tilenx, long tileny, int nullcheck, void *nullflagval, int nullval, int zbitpix, 
    double scale, double zero, int *intlength, int *flag, double *bscale, double *bzero, int *status);

static int unquantize_i1r4(long row,
            unsigned char *input,         /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            unsigned char tnull,          /* I - value of FITS TNULLn keyword if any */
            float nullval,        /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            float *output,        /* O - array of converted pixels           */
            int *status);          /* IO - error status                       */
static int unquantize_i2r4(long row,
            short *input,         /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            short tnull,          /* I - value of FITS TNULLn keyword if any */
            float nullval,        /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            float *output,        /* O - array of converted pixels           */
            int *status);          /* IO - error status                       */
static int unquantize_i4r4(long row,
            INT32BIT *input,      /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            INT32BIT tnull,       /* I - value of FITS TNULLn keyword if any */
            float nullval,        /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            float *output,        /* O - array of converted pixels           */
            int *status);          /* IO - error status                       */
static int unquantize_i1r8(long row,
            unsigned char *input,         /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            unsigned char tnull,          /* I - value of FITS TNULLn keyword if any */
            double nullval,        /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            double *output,        /* O - array of converted pixels           */
            int *status);          /* IO - error status                       */
static int unquantize_i2r8(long row,
            short *input,         /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            short tnull,          /* I - value of FITS TNULLn keyword if any */
            double nullval,        /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            double *output,        /* O - array of converted pixels           */
            int *status);          /* IO - error status                       */
static int unquantize_i4r8(long row,
            INT32BIT *input,      /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            INT32BIT tnull,       /* I - value of FITS TNULLn keyword if any */
            double nullval,        /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            double *output,        /* O - array of converted pixels           */
            int *status);          /* IO - error status                       */
static int imcomp_float2nan(float *indata, long tilelen, int *outdata,
    float nullflagval,  int *status);
static int imcomp_double2nan(double *indata, long tilelen, LONGLONG *outdata,
    double nullflagval,  int *status);    
static int fits_read_write_compressed_img(fitsfile *fptr,   /* I - FITS file pointer */
            int  datatype,  /* I - datatype of the array to be returned      */
            LONGLONG  *infpixel, /* I - 'bottom left corner' of the subsection    */
            LONGLONG  *inlpixel, /* I - 'top right corner' of the subsection      */
            long  *ininc,    /* I - increment to be applied in each dimension */
            int  nullcheck,  /* I - 0 for no null checking                   */
                              /*     1: set undefined pixels = nullval       */
            void *nullval,    /* I - value for undefined pixels              */
            int  *anynul,     /* O - set to 1 if any values are null; else 0 */
            fitsfile *outfptr,   /* I - FITS file pointer                    */
            int  *status);

static int fits_shuffle_8bytes(char *heap, LONGLONG length, int *status);
static int fits_shuffle_4bytes(char *heap, LONGLONG length, int *status);
static int fits_shuffle_2bytes(char *heap, LONGLONG length, int *status);
static int fits_gzip_heap(fitsfile *infptr, fitsfile *outfptr, int *status);
static int fits_unshuffle_8bytes(char *heap, LONGLONG length, int *status);
static int fits_unshuffle_4bytes(char *heap, LONGLONG length, int *status);
static int fits_unshuffle_2bytes(char *heap, LONGLONG length, int *status);
static int fits_gunzip_heap(fitsfile *infptr, fitsfile *outfptr, int *status);

/*---------------------------------------------------------------------------*/
int fits_init_randoms(void) {

/* initialize an array of random numbers */

    int ii;
    double a = 16807.0;
    double m = 2147483647.0;
    double temp, seed;

    FFLOCK;
 
    if (fits_rand_value) {
       FFUNLOCK;
       return(0);  /* array is already initialized */
    }

    /* allocate array for the random number sequence */
    /* THIS MEMORY IS NEVER FREED */
    fits_rand_value = calloc(N_RANDOM, sizeof(float));

    if (!fits_rand_value) {
        FFUNLOCK;
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

    FFUNLOCK;

    /* 
    IMPORTANT NOTE: the 10000th seed value must have the value 1043618065 if the 
       algorithm has been implemented correctly */
    
    if ( (int) seed != 1043618065) {
        ffpmsg("fits_init_randoms generated incorrect random number sequence");
	return(1);
    } else {
        return(0);
    }
}
/*--------------------------------------------------------------------------*/
void bz_internal_error(int errcode)
{
    /* external function declared by the bzip2 code in bzlib_private.h */
    ffpmsg("bzip2 returned an internal error");
    ffpmsg("This should never happen");
    return;
}
/*--------------------------------------------------------------------------*/
int fits_set_compression_type(fitsfile *fptr,  /* I - FITS file pointer     */
       int ctype,    /* image compression type code;                        */
                     /* allowed values: RICE_1, GZIP_1, GZIP_2, PLIO_1,     */
                     /*  HCOMPRESS_1, BZIP2_1, and NOCOMPRESS               */
       int *status)  /* IO - error status                                   */
{
/*
   This routine specifies the image compression algorithm that should be
   used when writing a FITS image.  The image is divided into tiles, and
   each tile is compressed and stored in a row of at variable length binary
   table column.
*/
    (fptr->Fptr)->request_compress_type = ctype;
    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_set_tile_dim(fitsfile *fptr,  /* I - FITS file pointer             */
           int ndim,   /* number of dimensions in the compressed image      */
           long *dims, /* size of image compression tile in each dimension  */
                      /* default tile size = (NAXIS1, 1, 1, ...)            */
           int *status)         /* IO - error status                        */
{
/*
   This routine specifies the size (dimension) of the image
   compression  tiles that should be used when writing a FITS
   image.  The image is divided into tiles, and each tile is compressed
   and stored in a row of at variable length binary table column.
*/
    int ii;

    if (ndim < 0 || ndim > MAX_COMPRESS_DIM)
    {
        *status = BAD_DIMEN;
        return(*status);
    }

    for (ii = 0; ii < ndim; ii++)
    {
        (fptr->Fptr)->request_tilesize[ii] = dims[ii];
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_set_quantize_level(fitsfile *fptr,  /* I - FITS file pointer   */
           float qlevel,        /* floating point quantization level      */
           int *status)         /* IO - error status                */
{
/*
   This routine specifies the value of the quantization level, q,  that
   should be used when compressing floating point images.  The image is
   divided into tiles, and each tile is compressed and stored in a row
   of at variable length binary table column.
*/
    if (qlevel == 0.)
    {
        /* this means don't quantize the floating point values. Instead, */
	/* the floating point values will be losslessly compressed */
       (fptr->Fptr)->quantize_level = NO_QUANTIZE;
    } else {

        (fptr->Fptr)->quantize_level = qlevel;

    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_set_quantize_dither(fitsfile *fptr,  /* I - FITS file pointer   */
           int dither,        /* dither type      */
           int *status)         /* IO - error status                */
{
/*
   This routine specifies what type of dithering (randomization) should
   be performed when quantizing floating point images to integer prior to
   compression.   A value of -1 means do no dithering.  A value of 0 means
   used the default  SUBTRACTIVE_DITHER_1 (which is equivalent to dither = 1).
   A value of -1 means do not apply any dither.
*/

    if (dither == 0) dither = 1;
    (fptr->Fptr)->request_quantize_dither = dither;
    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_set_dither_offset(fitsfile *fptr,  /* I - FITS file pointer   */
           int offset,        /* random dithering offset value (1 to 10000) */
           int *status)         /* IO - error status                */
{
/*
   This routine specifies the value of the offset that should be applied when
   calculating the random dithering when quantizing floating point iamges.
   A random offset should be applied to each image to avoid quantization 
   effects when taking the difference of 2 images, or co-adding a set of
   images.  Without this random offset, the corresponding pixel in every image
   will have exactly the same dithering.
   
   offset = 0 means use the default random dithering based on system time
   offset = negative means randomly chose dithering based on 1st tile checksum
   offset = [1 - 10000] means use that particular dithering pattern

*/
    /* if positive, ensure that the value is in the range 1 to 10000 */
    if (offset > 0)
       (fptr->Fptr)->request_dither_offset = ((offset - 1) % 10000 ) + 1; 
    else
       (fptr->Fptr)->request_dither_offset = offset; 

    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_set_noise_bits(fitsfile *fptr,  /* I - FITS file pointer   */
           int noisebits,       /* noise_bits parameter value       */
                                /* (default = 4)                    */
           int *status)         /* IO - error status                */
{
/*
   ********************************************************************
   ********************************************************************
   THIS ROUTINE IS PROVIDED ONLY FOR BACKWARDS COMPATIBILITY;
   ALL NEW SOFTWARE SHOULD CALL fits_set_quantize_level INSTEAD
   ********************************************************************
   ********************************************************************

   This routine specifies the value of the noice_bits parameter that
   should be used when compressing floating point images.  The image is
   divided into tiles, and each tile is compressed and stored in a row
   of at variable length binary table column.

   Feb 2008:  the "noisebits" parameter has been replaced with the more
   general "quantize level" parameter.
*/
    float qlevel;

    if (noisebits < 1 || noisebits > 16)
    {
        *status = DATA_COMPRESSION_ERR;
        return(*status);
    }

    qlevel = (float) pow (2., (double)noisebits);
    fits_set_quantize_level(fptr, qlevel, status);
    
    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_set_hcomp_scale(fitsfile *fptr,  /* I - FITS file pointer   */
           float scale,       /* hcompress scale parameter value       */
                                /* (default = 0.)                    */
           int *status)         /* IO - error status                */
{
/*
   This routine specifies the value of the hcompress scale parameter that
  The image is
   divided into tiles, and each tile is compressed and stored in a row
   of at variable length binary table column.
*/
    (fptr->Fptr)->request_hcomp_scale = scale;
    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_set_hcomp_smooth(fitsfile *fptr,  /* I - FITS file pointer   */
           int smooth,       /* hcompress smooth parameter value       */
                                /* if scale > 1 and smooth != 0, then */
				/*  the image will be smoothed when it is */
				/* decompressed to remove some of the */
				/* 'blockiness' in the image produced */
				/* by the lossy compression    */
           int *status)         /* IO - error status                */
{
/*
   This routine specifies the value of the hcompress scale parameter that
  The image is
   divided into tiles, and each tile is compressed and stored in a row
   of at variable length binary table column.
*/

    (fptr->Fptr)->request_hcomp_smooth = smooth;
    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_set_lossy_int(fitsfile *fptr,  /* I - FITS file pointer   */
           int lossy_int,       /* I - True (!= 0) or False (0) */
           int *status)         /* IO - error status                */
{
/*
   This routine specifies whether images with integer pixel values should
   quantized and compressed the same way float images are compressed.
   The default is to not do this, and instead apply a lossless compression
   algorithm to integer images.
*/

    (fptr->Fptr)->request_lossy_int_compress = lossy_int;
    return(*status);
}/*--------------------------------------------------------------------------*/
int fits_get_compression_type(fitsfile *fptr,  /* I - FITS file pointer     */
       int *ctype,   /* image compression type code;                        */
                     /* allowed values:                                     */
		     /* RICE_1, GZIP_1, GZIP_2, PLIO_1, HCOMPRESS_1, BZIP2_1 */
       int *status)  /* IO - error status                                   */
{
/*
   This routine returns the image compression algorithm that should be
   used when writing a FITS image.  The image is divided into tiles, and
   each tile is compressed and stored in a row of at variable length binary
   table column.
*/
    *ctype = (fptr->Fptr)->request_compress_type;
    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_get_tile_dim(fitsfile *fptr,  /* I - FITS file pointer             */
           int ndim,   /* number of dimensions in the compressed image      */
           long *dims, /* size of image compression tile in each dimension  */
                       /* default tile size = (NAXIS1, 1, 1, ...)           */
           int *status)         /* IO - error status                        */
{
/*
   This routine returns the size (dimension) of the image
   compression  tiles that should be used when writing a FITS
   image.  The image is divided into tiles, and each tile is compressed
   and stored in a row of at variable length binary table column.
*/
    int ii;

    if (ndim < 0 || ndim > MAX_COMPRESS_DIM)
    {
        *status = BAD_DIMEN;
        return(*status);
    }

    for (ii = 0; ii < ndim; ii++)
    {
        dims[ii] = (fptr->Fptr)->request_tilesize[ii];
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_get_noise_bits(fitsfile *fptr,  /* I - FITS file pointer   */
           int *noisebits,       /* noise_bits parameter value       */
                                /* (default = 4)                    */
           int *status)         /* IO - error status                */
{
/*
   ********************************************************************
   ********************************************************************
   THIS ROUTINE IS PROVIDED ONLY FOR BACKWARDS COMPATIBILITY;
   ALL NEW SOFTWARE SHOULD CALL fits_set_quantize_level INSTEAD
   ********************************************************************
   ********************************************************************


   This routine returns the value of the noice_bits parameter that
   should be used when compressing floating point images.  The image is
   divided into tiles, and each tile is compressed and stored in a row
   of at variable length binary table column.

   Feb 2008: code changed to use the more general "quantize level" parameter
   rather than the "noise bits" parameter.  If quantize level is greater than
   zero, then the previous noisebits parameter is approximately given by
   
   noise bits = natural logarithm (quantize level) / natural log (2)
   
   This result is rounded to the nearest integer.
*/
    double qlevel;

    qlevel = (fptr->Fptr)->quantize_level;

    if (qlevel > 0. && qlevel < 65537. )
         *noisebits =  (int) ((log(qlevel) / log(2.0)) + 0.5);
    else 
        *noisebits = 0;

    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_get_quantize_level(fitsfile *fptr,  /* I - FITS file pointer   */
           float *qlevel,       /* quantize level parameter value       */
           int *status)         /* IO - error status                */
{
/*
   This routine returns the value of the noice_bits parameter that
   should be used when compressing floating point images.  The image is
   divided into tiles, and each tile is compressed and stored in a row
   of at variable length binary table column.
*/

    if ((fptr->Fptr)->quantize_level == NO_QUANTIZE) {
      *qlevel = 0;
    } else {
      *qlevel = (fptr->Fptr)->quantize_level;
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_get_dither_offset(fitsfile *fptr,  /* I - FITS file pointer   */
           int *offset,       /* dithering offset parameter value       */
           int *status)         /* IO - error status                */
{
/*
   This routine returns the value of the dithering offset parameter that
   is used when compressing floating point images.  The image is
   divided into tiles, and each tile is compressed and stored in a row
   of at variable length binary table column.
*/

    *offset = (fptr->Fptr)->request_dither_offset;
    return(*status);
}/*--------------------------------------------------------------------------*/
int fits_get_hcomp_scale(fitsfile *fptr,  /* I - FITS file pointer   */
           float *scale,          /* Hcompress scale parameter value       */
           int *status)         /* IO - error status                */

{
/*
   This routine returns the value of the noice_bits parameter that
   should be used when compressing floating point images.  The image is
   divided into tiles, and each tile is compressed and stored in a row
   of at variable length binary table column.
*/

    *scale = (fptr->Fptr)->request_hcomp_scale;
    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_get_hcomp_smooth(fitsfile *fptr,  /* I - FITS file pointer   */
           int *smooth,          /* Hcompress smooth parameter value       */
           int *status)         /* IO - error status                */

{
    *smooth = (fptr->Fptr)->request_hcomp_smooth;
    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_img_compress(fitsfile *infptr, /* pointer to image to be compressed */
                 fitsfile *outfptr, /* empty HDU for output compressed image */
                 int *status)       /* IO - error status               */

/*
   This routine initializes the output table, copies all the keywords,
   and  loops through the input image, compressing the data and
   writing the compressed tiles to the output table.
   
   This is a high level routine that is called by the fpack and funpack
   FITS compression utilities.
*/
{
    int bitpix, naxis;
    long naxes[MAX_COMPRESS_DIM];

    if (*status > 0)
        return(*status);

    /* get datatype and size of input image */
    if (fits_get_img_param(infptr, MAX_COMPRESS_DIM, &bitpix, 
                       &naxis, naxes, status) > 0)
        return(*status);

    if (naxis < 1 || naxis > MAX_COMPRESS_DIM)
    {
        ffpmsg("Image cannot be compressed: NAXIS out of range");
        return(*status = BAD_NAXIS);
    }

    /* if requested, treat integer images same as a float image. */
    /* Then the pixels will be quantized (lossy algorithm) to achieve */
    /* higher amounts of compression than with lossless algorithms */

    if ( (outfptr->Fptr)->request_lossy_int_compress != 0  && bitpix > 0) 
	bitpix = FLOAT_IMG;  /* compress integer images as if float */

    /* initialize output table */
    if (imcomp_init_table(outfptr, bitpix, naxis, naxes, 0, status) > 0)
        return (*status);    
    
    /* Copy the image header keywords to the table header. */
    if (imcomp_copy_img2comp(infptr, outfptr, status) > 0)
	    return (*status);

    /* turn off any intensity scaling (defined by BSCALE and BZERO */
    /* keywords) so that unscaled values will be read by CFITSIO */
    /* (except if quantizing an int image, same as a float image) */
    if ( (outfptr->Fptr)->request_lossy_int_compress == 0 && bitpix > 0) 
        ffpscl(infptr, 1.0, 0.0, status);

    /* force a rescan of the output file keywords, so that */
    /* the compression parameters will be copied to the internal */
    /* fitsfile structure used by CFITSIO */
    ffrdef(outfptr, status);

    /* turn off any intensity scaling (defined by BSCALE and BZERO */
    /* keywords) so that unscaled values will be written by CFITSIO */
    /* (except if quantizing an int image, same as a float image) */
    if ( (outfptr->Fptr)->request_lossy_int_compress == 0 && bitpix > 0) 
        ffpscl(outfptr, 1.0, 0.0, status);

    /* Read each image tile, compress, and write to a table row. */
    imcomp_compress_image (infptr, outfptr, status);

    /* force another rescan of the output file keywords, to */
    /* update PCOUNT and TFORMn = '1PB(iii)' keyword values. */
    ffrdef(outfptr, status);

    return (*status);
}
/*--------------------------------------------------------------------------*/
int fits_compress_img_OBSOLETE(fitsfile *infptr, /* pointer to image to be compressed */
                 fitsfile *outfptr, /* empty HDU for output compressed image */
                 int compress_type, /* compression type code               */
                                    /*  RICE_1, HCOMPRESS_1, etc.          */
                 long *intilesize,    /* size in each dimension of the tiles */
                                    /* NULL pointer means tile by rows */
		 int blocksize,     /* compression parameter: blocksize  */
                 int nbits,         /* compression parameter: nbits  */
                 int *status)       /* IO - error status               */

/*
  !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!!
   This routine is obsolete and should not be used.  The 
   ftools 'fimgzip' task used to call this routine (but that task has been deleted); 

  The name of the routine was changed 4/27/2011, to see if anyone complains.
  If not, then this routine should be deleted from the source code.
  !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!!
  
   This routine initializes the output table, copies all the keywords,
   and  loops through the input image, compressing the data and
   writing the compressed tiles to the output table.
*/
{
    int bitpix, naxis;
    long naxes[MAX_COMPRESS_DIM];

    if (*status > 0)
        return(*status);

    /* get datatype and size of input image */
    if (fits_get_img_param(infptr, MAX_COMPRESS_DIM, &bitpix, 
                       &naxis, naxes, status) > 0)
        return(*status);

    if (naxis < 1 || naxis > MAX_COMPRESS_DIM)
    {
        ffpmsg("Image cannot be compressed: NAXIS out of range");
        return(*status = BAD_NAXIS);
    }

    /* initialize output table */
    if (imcomp_init_table(outfptr, bitpix, naxis, naxes, 0, status) > 0)
        return (*status);

    /* Copy the image header keywords to the table header. */
    if (imcomp_copy_imheader(infptr, outfptr, status) > 0)
	    return (*status);

    /* turn off any intensity scaling (defined by BSCALE and BZERO */
    /* keywords) so that unscaled values will be read by CFITSIO */
    ffpscl(infptr, 1.0, 0.0, status);

    /* force a rescan of the output file keywords, so that */
    /* the compression parameters will be copied to the internal */
    /* fitsfile structure used by CFITSIO */
    ffrdef(outfptr, status);

    /* Read each image tile, compress, and write to a table row. */
    imcomp_compress_image (infptr, outfptr, status);

    /* force another rescan of the output file keywords, to */
    /* update PCOUNT and TFORMn = '1PB(iii)' keyword values. */
    ffrdef(outfptr, status);

    return (*status);
}
/*--------------------------------------------------------------------------*/
int imcomp_init_table(fitsfile *outfptr,
        int inbitpix,
        int naxis,
        long *naxes,
	int writebitpix,    /* write the ZBITPIX, ZNAXIS, and ZNAXES keyword? */
        int *status)
/* 
  create a BINTABLE extension for the output compressed image.
*/
{
    char keyname[FLEN_KEYWORD], zcmptype[12];
    int ii,  remain,  ncols, bitpix;
    long nrows;
    char *ttype[] = {"COMPRESSED_DATA", "ZSCALE", "ZZERO"};
    char *tform[3];
    char tf0[4], tf1[4], tf2[4];
    char *tunit[] = {"\0",            "\0",            "\0"  };
    char comm[FLEN_COMMENT];
    long actual_tilesize[MAX_COMPRESS_DIM]; /* Actual size to use for tiles */
    
    if (*status > 0)
        return(*status);

    /* check for special case of losslessly compressing floating point */
    /* images.  Only compression algorithm that supports this is GZIP */
    if ( (outfptr->Fptr)->quantize_level == NO_QUANTIZE) {
       if (((outfptr->Fptr)->request_compress_type != GZIP_1) &&
           ((outfptr->Fptr)->request_compress_type != GZIP_2)) {
         ffpmsg("Lossless compression of floating point images must use GZIP (imcomp_init_table)");
         return(*status = DATA_COMPRESSION_ERR);
       }
    }
    
    /* test for the 2 special cases that represent unsigned integers */
    if (inbitpix == USHORT_IMG)
        bitpix = SHORT_IMG;
    else if (inbitpix == ULONG_IMG)
        bitpix = LONG_IMG;
    else if (inbitpix == SBYTE_IMG)
        bitpix = BYTE_IMG;
    else 
        bitpix = inbitpix;

    /* reset default tile dimensions too if required */
    memcpy(actual_tilesize, outfptr->Fptr->request_tilesize, MAX_COMPRESS_DIM * sizeof(long));

    if ((outfptr->Fptr)->request_compress_type == HCOMPRESS_1) {

         if (naxis < 2 ) {
            ffpmsg("Hcompress cannot be used with 1-dimensional images (imcomp_init_table)");
            return(*status = DATA_COMPRESSION_ERR);

	 } else if  (naxes[0] < 4 || naxes[1] < 4) {
            ffpmsg("Hcompress minimum image dimension is 4 pixels (imcomp_init_table)");
            return(*status = DATA_COMPRESSION_ERR);
         }

         if ((actual_tilesize[0] == 0) &&
             (actual_tilesize[1] == 0) ){
	     
	    /* compress the whole image as a single tile */
             actual_tilesize[0] = naxes[0];
             actual_tilesize[1] = naxes[1];

              for (ii = 2; ii < naxis; ii++) {
	             /* set all higher tile dimensions = 1 */
                     actual_tilesize[ii] = 1;
	      }

         } else if ((actual_tilesize[0] == 0) &&
             (actual_tilesize[1] == 1) ){
	     
             /*
              The Hcompress algorithm is inherently 2D in nature, so the row by row
	      tiling that is used for other compression algorithms is not appropriate.
	      If the image has less than 30 rows, then the entire image will be compressed
	      as a single tile.  Otherwise the tiles will consist of 16 rows of the image. 
	      This keeps the tiles to a reasonable size, and it also includes enough rows
	      to allow good compression efficiency.  If the last tile of the image 
	      happens to contain less than 4 rows, then find another tile size with
	      between 14 and 30 rows (preferably even), so that the last tile has 
	      at least 4 rows
	     */ 
	      
             /* 1st tile dimension is the row length of the image */
             actual_tilesize[0] = naxes[0];

              if (naxes[1] <= 30) {  /* use whole image if it is small */
                   actual_tilesize[1] = naxes[1];
	      } else {
                /* look for another good tile dimension */
	          if        (naxes[1] % 16 == 0 || naxes[1] % 16 > 3) {
                      actual_tilesize[1] = 16;
		  } else if (naxes[1] % 24 == 0 || naxes[1] % 24 > 3) {
                      actual_tilesize[1] = 24;
		  } else if (naxes[1] % 20 == 0 || naxes[1] % 20 > 3) {
                      actual_tilesize[1] = 20;
		  } else if (naxes[1] % 30 == 0 || naxes[1] % 30 > 3) {
                      actual_tilesize[1] = 30;
		  } else if (naxes[1] % 28 == 0 || naxes[1] % 28 > 3) {
                      actual_tilesize[1] = 28;
		  } else if (naxes[1] % 26 == 0 || naxes[1] % 26 > 3) {
                      actual_tilesize[1] = 26;
		  } else if (naxes[1] % 22 == 0 || naxes[1] % 22 > 3) {
                      actual_tilesize[1] = 22;
		  } else if (naxes[1] % 18 == 0 || naxes[1] % 18 > 3) {
                      actual_tilesize[1] = 18;
		  } else if (naxes[1] % 14 == 0 || naxes[1] % 14 > 3) {
                      actual_tilesize[1] = 14;
		  } else  {
                      actual_tilesize[1] = 17;

		  }
	      }

        } else if (actual_tilesize[0] < 4 ||
                   actual_tilesize[1] < 4) {

            /* user-specified tile size is too small */
            ffpmsg("Hcompress minimum tile dimension is 4 pixels (imcomp_init_table)");
            return(*status = DATA_COMPRESSION_ERR);
	}
	
        /* check if requested tile size causes the last tile to to have less than 4 pixels */
        remain = naxes[0] % (actual_tilesize[0]);  /* 1st dimension */
        if (remain > 0 && remain < 4) {
            (actual_tilesize[0])++; /* try increasing tile size by 1 */
	   
            remain = naxes[0] % (actual_tilesize[0]);
            if (remain > 0 && remain < 4) {
                ffpmsg("Last tile along 1st dimension has less than 4 pixels (imcomp_init_table)");
                return(*status = DATA_COMPRESSION_ERR);	
            }        
        }

        remain = naxes[1] % (actual_tilesize[1]);  /* 2nd dimension */
        if (remain > 0 && remain < 4) {
            (actual_tilesize[1])++; /* try increasing tile size by 1 */
	   
            remain = naxes[1] % (actual_tilesize[1]);
            if (remain > 0 && remain < 4) {
                ffpmsg("Last tile along 2nd dimension has less than 4 pixels (imcomp_init_table)");
                return(*status = DATA_COMPRESSION_ERR);	
            }        
        }

    } /* end, if HCOMPRESS_1 */
    
    for (ii = 0; ii < naxis; ii++) {
        if (actual_tilesize[ii] <= 0) {
	    /* tile size of 0 means use the image size of that dimension */
            actual_tilesize[ii] = naxes[ii];
	}
    }

    /* ---- set up array of TFORM strings -------------------------------*/
    strcpy(tf0, "1PB");
    strcpy(tf1, "1D");
    strcpy(tf2, "1D");

    tform[0] = tf0;
    tform[1] = tf1;
    tform[2] = tf2;

    /* calculate number of rows in output table */
    nrows = 1;
    for (ii = 0; ii < naxis; ii++)
    {
        nrows = nrows * ((naxes[ii] - 1)/ (actual_tilesize[ii]) + 1);
    }

    /* determine the default  number of columns in the output table */
    if (bitpix < 0 && (outfptr->Fptr)->quantize_level != NO_QUANTIZE)  
        ncols = 3;  /* quantized and scaled floating point image */
    else
        ncols = 1; /* default table has just one 'COMPRESSED_DATA' column */

    if ((outfptr->Fptr)->request_compress_type == RICE_1)
    {
        strcpy(zcmptype, "RICE_1");
    }
    else if ((outfptr->Fptr)->request_compress_type == GZIP_1)
    {
        strcpy(zcmptype, "GZIP_1");
    }
    else if ((outfptr->Fptr)->request_compress_type == GZIP_2)
    {
        strcpy(zcmptype, "GZIP_2");
    }
    else if ((outfptr->Fptr)->request_compress_type == BZIP2_1)
    {
        strcpy(zcmptype, "BZIP2_1");
    }
    else if ((outfptr->Fptr)->request_compress_type == PLIO_1)
    {
        strcpy(zcmptype, "PLIO_1");
       /* the PLIO compression algorithm outputs short integers, not bytes */
       strcpy(tform[0], "1PI");
    }
    else if ((outfptr->Fptr)->request_compress_type == HCOMPRESS_1)
    {
        strcpy(zcmptype, "HCOMPRESS_1");
    }
    else if ((outfptr->Fptr)->request_compress_type == NOCOMPRESS)
    {
        strcpy(zcmptype, "NOCOMPRESS");
    }    
    else
    {
        ffpmsg("unknown compression type (imcomp_init_table)");
        return(*status = DATA_COMPRESSION_ERR);
    }

    /* create the bintable extension to contain the compressed image */
    ffcrtb(outfptr, BINARY_TBL, nrows, ncols, ttype, 
                tform, tunit, 0, status);

    /* Add standard header keywords. */
    ffpkyl (outfptr, "ZIMAGE", 1, 
           "extension contains compressed image", status);  

    if (writebitpix) {
        /*  write the keywords defining the datatype and dimensions of */
	/*  the uncompressed image.  If not, these keywords will be */
        /*  copied later from the input uncompressed image  */
	   
        ffpkyj (outfptr, "ZBITPIX", bitpix,
			"data type of original image", status);
        ffpkyj (outfptr, "ZNAXIS", naxis,
			"dimension of original image", status);

        for (ii = 0;  ii < naxis;  ii++)
        {
            sprintf (keyname, "ZNAXIS%d", ii+1);
            ffpkyj (outfptr, keyname, naxes[ii],
			"length of original image axis", status);
        }
    }
                      
    for (ii = 0;  ii < naxis;  ii++)
    {
        sprintf (keyname, "ZTILE%d", ii+1);
        ffpkyj (outfptr, keyname, actual_tilesize[ii],
			"size of tiles to be compressed", status);
    }

    if (bitpix < 0) {
       
	if ((outfptr->Fptr)->quantize_level == NO_QUANTIZE) {
	    ffpkys(outfptr, "ZQUANTIZ", "NONE", 
	      "Lossless compression without quantization", status);
	} else {
	    
	    /* Unless dithering has been specifically turned off by setting */
	    /* request_quantize_dither = -1, use dithering by default */
	    /* when quantizing floating point images. */
	
	    if ( (outfptr->Fptr)->request_quantize_dither == 0) 
              (outfptr->Fptr)->request_quantize_dither = SUBTRACTIVE_DITHER_1;
       
	    if ((outfptr->Fptr)->request_quantize_dither == SUBTRACTIVE_DITHER_1) {
	      ffpkys(outfptr, "ZQUANTIZ", "SUBTRACTIVE_DITHER_1", 
	        "Pixel Quantization Algorithm", status);

	      /* also write the associated ZDITHER0 keyword with a default value */
	      /* which may get updated later. */
              ffpky(outfptr, TINT, "ZDITHER0", &((outfptr->Fptr)->request_dither_offset), 
	       "dithering offset when quantizing floats", status);
            }
	}
    }

    ffpkys (outfptr, "ZCMPTYPE", zcmptype,
	          "compression algorithm", status);

    /* write any algorithm-specific keywords */
    if ((outfptr->Fptr)->request_compress_type == RICE_1)
    {
        ffpkys (outfptr, "ZNAME1", "BLOCKSIZE",
            "compression block size", status);

        /* for now at least, the block size is always 32 */
        ffpkyj (outfptr, "ZVAL1", 32,
			"pixels per block", status);

        ffpkys (outfptr, "ZNAME2", "BYTEPIX",
            "bytes per pixel (1, 2, 4, or 8)", status);

        if (bitpix == BYTE_IMG)
            ffpkyj (outfptr, "ZVAL2", 1,
			"bytes per pixel (1, 2, 4, or 8)", status);
        else if (bitpix == SHORT_IMG)
            ffpkyj (outfptr, "ZVAL2", 2,
			"bytes per pixel (1, 2, 4, or 8)", status);
        else 
            ffpkyj (outfptr, "ZVAL2", 4,
			"bytes per pixel (1, 2, 4, or 8)", status);

    }
    else if ((outfptr->Fptr)->request_compress_type == HCOMPRESS_1)
    {
        ffpkys (outfptr, "ZNAME1", "SCALE",
            "HCOMPRESS scale factor", status);
        ffpkye (outfptr, "ZVAL1", (outfptr->Fptr)->request_hcomp_scale,
		7, "HCOMPRESS scale factor", status);

        ffpkys (outfptr, "ZNAME2", "SMOOTH",
            "HCOMPRESS smooth option", status);
        ffpkyj (outfptr, "ZVAL2", (long) (outfptr->Fptr)->request_hcomp_smooth,
			"HCOMPRESS smooth option", status);
    }

    /* Write the BSCALE and BZERO keywords, if an unsigned integer image */
    if (inbitpix == USHORT_IMG)
    {
        strcpy(comm, "offset data range to that of unsigned short");
        ffpkyg(outfptr, "BZERO", 32768., 0, comm, status);
        strcpy(comm, "default scaling factor");
        ffpkyg(outfptr, "BSCALE", 1.0, 0, comm, status);
    }
    else if (inbitpix == SBYTE_IMG)
    {
        strcpy(comm, "offset data range to that of signed byte");
        ffpkyg(outfptr, "BZERO", -128., 0, comm, status);
        strcpy(comm, "default scaling factor");
        ffpkyg(outfptr, "BSCALE", 1.0, 0, comm, status);
    }
    else if (inbitpix == ULONG_IMG)
    {
        strcpy(comm, "offset data range to that of unsigned long");
        ffpkyg(outfptr, "BZERO", 2147483648., 0, comm, status);
        strcpy(comm, "default scaling factor");
        ffpkyg(outfptr, "BSCALE", 1.0, 0, comm, status);
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int imcomp_calc_max_elem (int comptype, int nx, int zbitpix, int blocksize)

/* This function returns the maximum number of bytes in a compressed
   image line.

    nx = maximum number of pixels in a tile
    blocksize is only relevant for RICE compression
*/
{    
    if (comptype == RICE_1)
    {
        if (zbitpix == 16)
            return (sizeof(short) * nx + nx / blocksize + 2 + 4);
	else
            return (sizeof(float) * nx + nx / blocksize + 2 + 4);
    }
    else if ((comptype == GZIP_1) || (comptype == GZIP_2))
    {
        /* gzip usually compressed by at least a factor of 2 for I*4 images */
        /* and somewhat less for I*2 images */
        /* If this size turns out to be too small, then the gzip */
        /* compression routine will allocate more space as required */
        /* to be on the safe size, allocate buffer same size as input */
	
        if (zbitpix == 16)
            return(nx * 2);
	else if (zbitpix == 8)
            return(nx);
	else
            return(nx * 4);
    }
    else if (comptype == BZIP2_1)
    {
        /* To guarantee that the compressed data will fit, allocate an output
	   buffer of size 1% larger than the uncompressed data, plus 600 bytes */

            return((int) (nx * 1.01 * zbitpix / 8. + 601.));
    }
     else if (comptype == HCOMPRESS_1)
    {
        /* Imperical evidence suggests in the worst case, 
	   the compressed stream could be up to 10% larger than the original
	   image.  Add 26 byte overhead, only significant for very small tiles
	   
         Possible improvement: may need to allow a larger size for 32-bit images */

        if (zbitpix == 16 || zbitpix == 8)
	
            return( (int) (nx * 2.2 + 26));   /* will be compressing 16-bit int array */
        else
            return( (int) (nx * 4.4 + 26));   /* will be compressing 32-bit int array */
    }
    else
        return(nx * sizeof(int));
}
/*--------------------------------------------------------------------------*/
int imcomp_compress_image (fitsfile *infptr, fitsfile *outfptr, int *status)

/* This routine does the following:
        - reads an image one tile at a time
        - if it is a float or double image, then it tries to quantize the pixels
          into scaled integers.
        - it then compressess the integer pixels, or if the it was not
	  possible to quantize the floating point pixels, then it losslessly
	  compresses them with gzip
	- writes the compressed byte stream to the output FITS file
*/
{
    double *tiledata;
    int anynul, gotnulls = 0, datatype;
    long ii, row;
    int naxis;
    double dummy = 0., dblnull = DOUBLENULLVALUE;
    float fltnull = FLOATNULLVALUE;
    long maxtilelen, tilelen, incre[] = {1, 1, 1, 1, 1, 1};
    long naxes[MAX_COMPRESS_DIM], fpixel[MAX_COMPRESS_DIM];
    long lpixel[MAX_COMPRESS_DIM], tile[MAX_COMPRESS_DIM];
    long tilesize[MAX_COMPRESS_DIM];
    long i0, i1, i2, i3, i4, i5;
    char card[FLEN_CARD];

    if (*status > 0)
        return(*status);

    maxtilelen = (outfptr->Fptr)->maxtilelen;

    /* 
     Allocate buffer to hold 1 tile of data; size depends on which compression 
     algorithm is used:

      Rice and GZIP will compress byte, short, or int arrays without conversion.
      PLIO requires 4-byte int values, so byte and short arrays must be converted to int.
      HCompress internally converts byte or short values to ints, and
         converts int values to 8-byte longlong integers.
    */
    
    if ((outfptr->Fptr)->zbitpix == FLOAT_IMG)
    {
        datatype = TFLOAT;

        if ( (outfptr->Fptr)->compress_type == HCOMPRESS_1) {
	    /* need twice as much scratch space (8 bytes per pixel) */
            tiledata = (double*) malloc (maxtilelen * 2 *sizeof (float));	
	} else {
            tiledata = (double*) malloc (maxtilelen * sizeof (float));
	}
    }
    else if ((outfptr->Fptr)->zbitpix == DOUBLE_IMG)
    {
        datatype = TDOUBLE;
        tiledata = (double*) malloc (maxtilelen * sizeof (double));
    }
    else if ((outfptr->Fptr)->zbitpix == SHORT_IMG)
    {
        datatype = TSHORT;
        if ( (outfptr->Fptr)->compress_type == RICE_1  ||
	     (outfptr->Fptr)->compress_type == GZIP_1  ||
	     (outfptr->Fptr)->compress_type == GZIP_2  ||
	     (outfptr->Fptr)->compress_type == BZIP2_1 ||
             (outfptr->Fptr)->compress_type == NOCOMPRESS) {
	    /* only need  buffer of I*2 pixels for gzip, bzip2, and Rice */

            tiledata = (double*) malloc (maxtilelen * sizeof (short));	
	} else {
 	    /*  need  buffer of I*4 pixels for Hcompress and PLIO */
            tiledata = (double*) malloc (maxtilelen * sizeof (int));
        }
    }
    else if ((outfptr->Fptr)->zbitpix == BYTE_IMG)
    {

        datatype = TBYTE;
        if ( (outfptr->Fptr)->compress_type == RICE_1  ||
	     (outfptr->Fptr)->compress_type == BZIP2_1 ||
	     (outfptr->Fptr)->compress_type == GZIP_1  ||
	     (outfptr->Fptr)->compress_type == GZIP_2) {
	    /* only need  buffer of I*1 pixels for gzip, bzip2, and Rice */

            tiledata = (double*) malloc (maxtilelen);	
	} else {
 	    /*  need  buffer of I*4 pixels for Hcompress and PLIO */
            tiledata = (double*) malloc (maxtilelen * sizeof (int));
        }
    }
    else if ((outfptr->Fptr)->zbitpix == LONG_IMG)
    {
        datatype = TINT;
        if ( (outfptr->Fptr)->compress_type == HCOMPRESS_1) {
	    /* need twice as much scratch space (8 bytes per pixel) */

            tiledata = (double*) malloc (maxtilelen * 2 * sizeof (int));	
	} else {
 	    /* only need  buffer of I*4 pixels for gzip, bzip2,  Rice, and PLIO */

            tiledata = (double*) malloc (maxtilelen * sizeof (int));
        }
    }
    else
    {
	ffpmsg("Bad image datatype. (imcomp_compress_image)");
	return (*status = MEMORY_ALLOCATION);
    }
    
    if (tiledata == NULL)
    {
	ffpmsg("Out of memory. (imcomp_compress_image)");
	return (*status = MEMORY_ALLOCATION);
    }

    /*  calculate size of tile in each dimension */
    naxis = (outfptr->Fptr)->zndim;
    for (ii = 0; ii < MAX_COMPRESS_DIM; ii++)
    {
        if (ii < naxis)
        {
             naxes[ii] = (outfptr->Fptr)->znaxis[ii];
             tilesize[ii] = (outfptr->Fptr)->tilesize[ii];
        }
        else
        {
            naxes[ii] = 1;
            tilesize[ii] = 1;
        }
    }
    row = 1;

    /* set up big loop over up to 6 dimensions */
    for (i5 = 1; i5 <= naxes[5]; i5 += tilesize[5])
    {
     fpixel[5] = i5;
     lpixel[5] = minvalue(i5 + tilesize[5] - 1, naxes[5]);
     tile[5] = lpixel[5] - fpixel[5] + 1;
     for (i4 = 1; i4 <= naxes[4]; i4 += tilesize[4])
     {
      fpixel[4] = i4;
      lpixel[4] = minvalue(i4 + tilesize[4] - 1, naxes[4]);
      tile[4] = lpixel[4] - fpixel[4] + 1;
      for (i3 = 1; i3 <= naxes[3]; i3 += tilesize[3])
      {
       fpixel[3] = i3;
       lpixel[3] = minvalue(i3 + tilesize[3] - 1, naxes[3]);
       tile[3] = lpixel[3] - fpixel[3] + 1;
       for (i2 = 1; i2 <= naxes[2]; i2 += tilesize[2])
       {
        fpixel[2] = i2;
        lpixel[2] = minvalue(i2 + tilesize[2] - 1, naxes[2]);
        tile[2] = lpixel[2] - fpixel[2] + 1;
        for (i1 = 1; i1 <= naxes[1]; i1 += tilesize[1])
        {
         fpixel[1] = i1;
         lpixel[1] = minvalue(i1 + tilesize[1] - 1, naxes[1]);
         tile[1] = lpixel[1] - fpixel[1] + 1;
         for (i0 = 1; i0 <= naxes[0]; i0 += tilesize[0])
         {
          fpixel[0] = i0;
          lpixel[0] = minvalue(i0 + tilesize[0] - 1, naxes[0]);
          tile[0] = lpixel[0] - fpixel[0] + 1;

          /* number of pixels in this tile */
          tilelen = tile[0];
          for (ii = 1; ii < naxis; ii++)
          {
             tilelen *= tile[ii];
          }

          /* read next tile of data from image */
	  anynul = 0;
          if (datatype == TFLOAT)
          {
              ffgsve(infptr, 1, naxis, naxes, fpixel, lpixel, incre, 
                  FLOATNULLVALUE, (float *) tiledata,  &anynul, status);
          }
          else if (datatype == TDOUBLE)
          {
              ffgsvd(infptr, 1, naxis, naxes, fpixel, lpixel, incre, 
                  DOUBLENULLVALUE, tiledata, &anynul, status);
          }
          else if (datatype == TINT)
          {
              ffgsvk(infptr, 1, naxis, naxes, fpixel, lpixel, incre, 
                  0, (int *) tiledata,  &anynul, status);
          }
          else if (datatype == TSHORT)
          {
              ffgsvi(infptr, 1, naxis, naxes, fpixel, lpixel, incre, 
                  0, (short *) tiledata,  &anynul, status);
          }
          else if (datatype == TBYTE)
          {
              ffgsvb(infptr, 1, naxis, naxes, fpixel, lpixel, incre, 
                  0, (unsigned char *) tiledata,  &anynul, status);
          }
          else 
          {
              ffpmsg("Error bad datatype of image tile to compress");
              free(tiledata);
              return (*status);
          }

          /* now compress the tile, and write to row of binary table */
          /*   NOTE: we don't have to worry about the presence of null values in the
	       array if it is an integer array:  the null value is simply encoded
	       in the compressed array just like any other pixel value.  
	       
	       If it is a floating point array, then we need to check for null
	       only if the anynul parameter returned a true value when reading the tile
	  */
          if (anynul && datatype == TFLOAT) {
              imcomp_compress_tile(outfptr, row, datatype, tiledata, tilelen,
                               tile[0], tile[1], 1, &fltnull, status);
          } else if (anynul && datatype == TDOUBLE) {
              imcomp_compress_tile(outfptr, row, datatype, tiledata, tilelen,
                               tile[0], tile[1], 1, &dblnull, status);
          } else {
              imcomp_compress_tile(outfptr, row, datatype, tiledata, tilelen,
                               tile[0], tile[1], 0, &dummy, status);
          }

          /* set flag if we found any null values */
          if (anynul)
              gotnulls = 1;

          /* check for any error in the previous operations */
          if (*status > 0)
          {
              ffpmsg("Error writing compressed image to table");
              free(tiledata);
              return (*status);
          }

	  row++;
         }
        }
       }
      }
     }
    }

    free (tiledata);  /* finished with this buffer */

    /* insert ZBLANK keyword if necessary; only for TFLOAT or TDOUBLE images */
    if (gotnulls)
    {
          ffgcrd(outfptr, "ZCMPTYPE", card, status);
          ffikyj(outfptr, "ZBLANK", COMPRESS_NULL_VALUE, 
             "null value in the compressed integer array", status);
    }

    return (*status);
}
/*--------------------------------------------------------------------------*/
int imcomp_compress_tile (fitsfile *outfptr,
    long row,  /* tile number = row in the binary table that holds the compressed data */
    int datatype, 
    void *tiledata, 
    long tilelen,
    long tilenx,
    long tileny,
    int nullcheck,
    void *nullflagval,
    int *status)

/*
   This is the main compression routine.

   This routine does the following to the input tile of pixels:
        - if it is a float or double image, then it quantizes the pixels
        - compresses the integer pixel values
        - writes the compressed byte stream to the FITS file.

   If the tile cannot be quantized than the raw float or double values
   are losslessly compressed with gzip and then written to the output table.
   
   This input array may be modified by this routine.  If the array is of type TINT
   or TFLOAT, and the compression type is HCOMPRESS, then it must have been 
   allocated to be twice as large (8 bytes per pixel) to provide scratch space.

  Note that this routine does not fully support the implicit datatype conversion that
  is supported when writing to normal FITS images.  The datatype of the input array
  must have the same datatype (either signed or unsigned) as the output (compressed)
  FITS image in some cases.
*/
{
    int *idata;		/* quantized integer data */
    int cn_zblank, zbitpix, nullval;
    int flag = 1;  /* true by default; only = 0 if float data couldn't be quantized */
    int intlength;      /* size of integers to be compressed */
    double scale, zero, actual_bzero;
    long ii;
    size_t clen;		/* size of cbuf */
    short *cbuf;	/* compressed data */
    int  nelem = 0;		/* number of bytes */
    size_t gzip_nelem = 0;
    unsigned int bzlen;
    int ihcompscale;
    float hcompscale;
    double noise2, noise3, noise5;
    double bscale[1] = {1.}, bzero[1] = {0.};	/* scaling parameters */
    long  hcomp_len;
    LONGLONG *lldata;

    if (*status > 0)
        return(*status);

    /* check for special case of losslessly compressing floating point */
    /* images.  Only compression algorithm that supports this is GZIP */
    if ( (outfptr->Fptr)->quantize_level == NO_QUANTIZE) {
       if (((outfptr->Fptr)->compress_type != GZIP_1) &&
           ((outfptr->Fptr)->compress_type != GZIP_2)) {
         ffpmsg("Lossless compression of floating point images must use GZIP");
         return(*status = DATA_COMPRESSION_ERR);
       }
    }

    /* free the previously saved tile if the input tile is for the same row */
    if ((outfptr->Fptr)->tilerow == row) {
        if ((outfptr->Fptr)->tiledata) {
            free((outfptr->Fptr)->tiledata);
        }
	  
        if ((outfptr->Fptr)->tilenullarray) {
            free((outfptr->Fptr)->tilenullarray);
        }

        (outfptr->Fptr)->tiledata = 0;
        (outfptr->Fptr)->tilenullarray = 0;
        (outfptr->Fptr)->tilerow = 0;
        (outfptr->Fptr)->tiledatasize = 0;
        (outfptr->Fptr)->tiletype = 0;
    }

    if ( (outfptr->Fptr)->compress_type == NOCOMPRESS) {
         /* Special case when using NOCOMPRESS for diagnostic purposes in fpack */
         if (imcomp_write_nocompress_tile(outfptr, row, datatype, tiledata, tilelen, 
	     nullcheck, nullflagval, status) > 0) {
             return(*status);
         }
         return(*status);
    }

    /* =========================================================================== */
    /* initialize various parameters */
    idata = (int *) tiledata;   /* may overwrite the input tiledata in place */

    /* zbitpix is the BITPIX keyword value in the uncompressed FITS image */
    zbitpix = (outfptr->Fptr)->zbitpix;

    /* if the tile/image has an integer datatype, see if a null value has */
    /* been defined (with the BLANK keyword in a normal FITS image).  */
    /* If so, and if the input tile array also contains null pixels, */
    /* (represented by pixels that have a value = nullflagval) then  */
    /* any pixels whose value = nullflagval, must be set to the value = nullval */
    /* before the pixel array is compressed.  These null pixel values must */
    /* not be inverse scaled by the BSCALE/BZERO values, if present. */

    cn_zblank = (outfptr->Fptr)->cn_zblank;
    nullval = (outfptr->Fptr)->zblank;

    if (zbitpix > 0 && cn_zblank != -1)  /* If the integer image has no defined null */
        nullcheck = 0;    /* value, then don't bother checking input array for nulls. */

    /* if the BSCALE and BZERO keywords exist, then the input values must */
    /* be inverse scaled by this factor, before the values are compressed. */
    /* (The program may have turned off scaling, which over rides the keywords) */
    
    scale = (outfptr->Fptr)->cn_bscale;
    zero  = (outfptr->Fptr)->cn_bzero;
    actual_bzero = (outfptr->Fptr)->cn_actual_bzero;

    /* =========================================================================== */
    /* prepare the tile of pixel values for compression */
    if (datatype == TSHORT) {
       imcomp_convert_tile_tshort(outfptr, tiledata, tilelen, nullcheck, nullflagval,
           nullval, zbitpix, scale, zero, actual_bzero, &intlength, status);
    } else if (datatype == TUSHORT) {
       imcomp_convert_tile_tushort(outfptr, tiledata, tilelen, nullcheck, nullflagval,
           nullval, zbitpix, scale, zero, &intlength, status);
    } else if (datatype == TBYTE) {
       imcomp_convert_tile_tbyte(outfptr, tiledata, tilelen, nullcheck, nullflagval,
           nullval, zbitpix, scale, zero,  &intlength, status);
    } else if (datatype == TSBYTE) {
       imcomp_convert_tile_tsbyte(outfptr, tiledata, tilelen, nullcheck, nullflagval,
           nullval, zbitpix, scale, zero,  &intlength, status);
    } else if (datatype == TINT) {
       imcomp_convert_tile_tint(outfptr, tiledata, tilelen, nullcheck, nullflagval,
           nullval, zbitpix, scale, zero, &intlength, status);
    } else if (datatype == TUINT) {
       imcomp_convert_tile_tuint(outfptr, tiledata, tilelen, nullcheck, nullflagval,
           nullval, zbitpix, scale, zero, &intlength, status);
    } else if (datatype == TLONG && sizeof(long) == 8) {
           ffpmsg("Integer*8 Long datatype is not supported when writing to compressed images");
           return(*status = BAD_DATATYPE);
    } else if (datatype == TULONG && sizeof(long) == 8) {
           ffpmsg("Unsigned integer*8 datatype is not supported when writing to compressed images");
           return(*status = BAD_DATATYPE);
    } else if (datatype == TFLOAT) {
        imcomp_convert_tile_tfloat(outfptr, row, tiledata, tilelen, tilenx, tileny, nullcheck,
        nullflagval, nullval, zbitpix, scale, zero, &intlength, &flag, bscale, bzero, status);
    } else if (datatype == TDOUBLE) {
       imcomp_convert_tile_tdouble(outfptr, row, tiledata, tilelen, tilenx, tileny, nullcheck,
       nullflagval, nullval, zbitpix, scale, zero, &intlength, &flag, bscale, bzero, status);
    } else {
          ffpmsg("unsupported image datatype (imcomp_compress_tile)");
          return(*status = BAD_DATATYPE);
    }

    if (*status > 0)
      return(*status);      /* return if error occurs */

    /* =========================================================================== */
    if (flag)   /* now compress the integer data array */
    {
        /* allocate buffer for the compressed tile bytes */
        clen = (outfptr->Fptr)->maxelem;
        cbuf = (short *) calloc (clen, sizeof (unsigned char));

        if (cbuf == NULL) {
            ffpmsg("Memory allocation failure. (imcomp_compress_tile)");
	    return (*status = MEMORY_ALLOCATION);
        }

        /* =========================================================================== */
        if ( (outfptr->Fptr)->compress_type == RICE_1)
        {
            if (intlength == 2) {
  	        nelem = fits_rcomp_short ((short *)idata, tilelen, (unsigned char *) cbuf,
                       clen, (outfptr->Fptr)->rice_blocksize);
            } else if (intlength == 1) {
  	        nelem = fits_rcomp_byte ((signed char *)idata, tilelen, (unsigned char *) cbuf,
                       clen, (outfptr->Fptr)->rice_blocksize);
            } else {
  	        nelem = fits_rcomp (idata, tilelen, (unsigned char *) cbuf,
                       clen, (outfptr->Fptr)->rice_blocksize);
            }

	    if (nelem < 0)  /* data compression error condition */
            {
	        free (cbuf);
                ffpmsg("error Rice compressing image tile (imcomp_compress_tile)");
                return (*status = DATA_COMPRESSION_ERR);
            }

	    /* Write the compressed byte stream. */
            ffpclb(outfptr, (outfptr->Fptr)->cn_compressed, row, 1,
                     nelem, (unsigned char *) cbuf, status);
        }

        /* =========================================================================== */
        else if ( (outfptr->Fptr)->compress_type == PLIO_1)
        {
              for (ii = 0; ii < tilelen; ii++)  {
                if (idata[ii] < 0 || idata[ii] > 16777215)
                {
                   /* plio algorithn only supports positive 24 bit ints */
                   ffpmsg("data out of range for PLIO compression (0 - 2**24)");
                   return(*status = DATA_COMPRESSION_ERR);
                }
              }

  	      nelem = pl_p2li (idata, 1, cbuf, tilelen);

	      if (nelem < 0)  /* data compression error condition */
              {
	        free (cbuf);
                ffpmsg("error PLIO compressing image tile (imcomp_compress_tile)");
                return (*status = DATA_COMPRESSION_ERR);
              }

	      /* Write the compressed byte stream. */
              ffpcli(outfptr, (outfptr->Fptr)->cn_compressed, row, 1,
                     nelem, cbuf, status);
        }

        /* =========================================================================== */
        else if ( ((outfptr->Fptr)->compress_type == GZIP_1) ||
                  ((outfptr->Fptr)->compress_type == GZIP_2) )   {

	    if ((outfptr->Fptr)->quantize_level == NO_QUANTIZE && datatype == TFLOAT) {
	      /* Special case of losslessly compressing floating point pixels with GZIP */
	      /* In this case we compress the input tile array directly */

#if BYTESWAPPED
               ffswap4((int*) tiledata, tilelen); 
#endif
               if ( (outfptr->Fptr)->compress_type == GZIP_2 )
		    fits_shuffle_4bytes((char *) tiledata, tilelen, status);

                compress2mem_from_mem((char *) tiledata, tilelen * sizeof(float),
                    (char **) &cbuf,  &clen, realloc, &gzip_nelem, status);

	    } else if ((outfptr->Fptr)->quantize_level == NO_QUANTIZE && datatype == TDOUBLE) {
	      /* Special case of losslessly compressing double pixels with GZIP */
	      /* In this case we compress the input tile array directly */

#if BYTESWAPPED
               ffswap8((double *) tiledata, tilelen); 
#endif
               if ( (outfptr->Fptr)->compress_type == GZIP_2 )
		    fits_shuffle_8bytes((char *) tiledata, tilelen, status);

                compress2mem_from_mem((char *) tiledata, tilelen * sizeof(double),
                    (char **) &cbuf,  &clen, realloc, &gzip_nelem, status);

	    } else {

	        /* compress the integer idata array */

#if BYTESWAPPED
	       if (intlength == 2)
                 ffswap2((short *) idata, tilelen); 
	       else if (intlength == 4)
                 ffswap4(idata, tilelen); 
#endif

               if (intlength == 2) {

                  if ( (outfptr->Fptr)->compress_type == GZIP_2 )
		    fits_shuffle_2bytes((char *) tiledata, tilelen, status);

                  compress2mem_from_mem((char *) idata, tilelen * sizeof(short),
                   (char **) &cbuf,  &clen, realloc, &gzip_nelem, status);

               } else if (intlength == 1) {

                  compress2mem_from_mem((char *) idata, tilelen * sizeof(unsigned char),
                   (char **) &cbuf,  &clen, realloc, &gzip_nelem, status);

               } else {

                  if ( (outfptr->Fptr)->compress_type == GZIP_2 )
		    fits_shuffle_4bytes((char *) tiledata, tilelen, status);

                  compress2mem_from_mem((char *) idata, tilelen * sizeof(int),
                   (char **) &cbuf,  &clen, realloc, &gzip_nelem, status);
               }
            }

	    /* Write the compressed byte stream. */
            ffpclb(outfptr, (outfptr->Fptr)->cn_compressed, row, 1,
                     gzip_nelem, (unsigned char *) cbuf, status);

        /* =========================================================================== */
        } else if ( (outfptr->Fptr)->compress_type == BZIP2_1) {

#if BYTESWAPPED
	   if (intlength == 2)
               ffswap2((short *) idata, tilelen); 
	   else if (intlength == 4)
               ffswap4(idata, tilelen); 
#endif

           bzlen = (unsigned int) clen;
	   
           /* call bzip2 with blocksize = 900K, verbosity = 0, and default workfactor */

/*  bzip2 is not supported in the public release.  This is only for test purposes.
           if (BZ2_bzBuffToBuffCompress( (char *) cbuf, &bzlen,
	         (char *) idata, (unsigned int) (tilelen * intlength), 9, 0, 0) ) 
*/
	   {
                   ffpmsg("bzip2 compression error");
                   return(*status = DATA_COMPRESSION_ERR);
           }

	    /* Write the compressed byte stream. */
            ffpclb(outfptr, (outfptr->Fptr)->cn_compressed, row, 1,
                     bzlen, (unsigned char *) cbuf, status);

        /* =========================================================================== */
        }  else if ( (outfptr->Fptr)->compress_type == HCOMPRESS_1)     {
	    /*
	      if hcompscale is positive, then we have to multiply
	      the value by the RMS background noise to get the 
	      absolute scale value.  If negative, then it gives the
	      absolute scale value directly.
	    */
            hcompscale = (outfptr->Fptr)->hcomp_scale;

	    if (hcompscale > 0.) {
	       fits_img_stats_int(idata, tilenx, tileny, nullcheck,
	                nullval, 0,0,0,0,0,0,&noise2,&noise3,&noise5,status);

		/* use the minimum of the 3 noise estimates */
		if (noise2 != 0. && noise2 < noise3) noise3 = noise2;
		if (noise5 != 0. && noise5 < noise3) noise3 = noise5;
		
		hcompscale = (float) (hcompscale * noise3);

	    } else if (hcompscale < 0.) {

		hcompscale = hcompscale * -1.0F;
	    }

	    ihcompscale = (int) (hcompscale + 0.5);

            hcomp_len = clen;  /* allocated size of the buffer */
	    
            if (zbitpix == BYTE_IMG || zbitpix == SHORT_IMG) {
                fits_hcompress(idata, tilenx, tileny, 
		  ihcompscale, (char *) cbuf, &hcomp_len, status);

            } else {
                 /* have to convert idata to an I*8 array, in place */
                 /* idata must have been allocated large enough to do this */
                lldata = (LONGLONG *) idata;
		
                for (ii = tilelen - 1; ii >= 0; ii--) {
		    lldata[ii] = idata[ii];
		}

                fits_hcompress64(lldata, tilenx, tileny, 
		  ihcompscale, (char *) cbuf, &hcomp_len, status);
            }

	    /* Write the compressed byte stream. */
            ffpclb(outfptr, (outfptr->Fptr)->cn_compressed, row, 1,
                     hcomp_len, (unsigned char *) cbuf, status);
        }

        /* =========================================================================== */
        if ((outfptr->Fptr)->cn_zscale > 0)
        {
              /* write the linear scaling parameters for this tile */
	      ffpcld (outfptr, (outfptr->Fptr)->cn_zscale, row, 1, 1,
                      bscale, status);
	      ffpcld (outfptr, (outfptr->Fptr)->cn_zzero,  row, 1, 1,
                      bzero,  status);
        }

        free(cbuf);  /* finished with this buffer */

    /* =========================================================================== */
    } else {    /* if flag == 0., floating point data couldn't be quantized */

	 /* losslessly compress the data with gzip. */

         /* if gzip2 compressed data column doesn't exist, create it */
         if ((outfptr->Fptr)->cn_gzip_data < 1) {
                 fits_insert_col(outfptr, 999, "GZIP_COMPRESSED_DATA", "1PB", status);

                 if (*status <= 0)  /* save the number of this column */
                       ffgcno(outfptr, CASEINSEN, "GZIP_COMPRESSED_DATA",
                                &(outfptr->Fptr)->cn_gzip_data, status);
         }
 
         if (datatype == TFLOAT)  {
               /* allocate buffer for the compressed tile bytes */
	       /* make it 10% larger than the original uncompressed data */
               clen = tilelen * sizeof(float) * 1.1;
               cbuf = (short *) calloc (clen, sizeof (unsigned char));

               if (cbuf == NULL)
               {
                   ffpmsg("Memory allocation error. (imcomp_compress_tile)");
	           return (*status = MEMORY_ALLOCATION);
               }

	       /* convert null values to NaNs in place, if necessary */
	       if (nullcheck == 1) {
	           imcomp_float2nan((float *) tiledata, tilelen, (int *) tiledata,
	               *(float *) (nullflagval), status);
	       }

#if BYTESWAPPED
               ffswap4((int*) tiledata, tilelen); 
#endif
               compress2mem_from_mem((char *) tiledata, tilelen * sizeof(float),
                    (char **) &cbuf,  &clen, realloc, &gzip_nelem, status);

         } else if (datatype == TDOUBLE) {

               /* allocate buffer for the compressed tile bytes */
	       /* make it 10% larger than the original uncompressed data */
               clen = tilelen * sizeof(double) * 1.1;
               cbuf = (short *) calloc (clen, sizeof (unsigned char));

               if (cbuf == NULL)
               {
                   ffpmsg("Memory allocation error. (imcomp_compress_tile)");
	           return (*status = MEMORY_ALLOCATION);
               }

	       /* convert null values to NaNs in place, if necessary */
	       if (nullcheck == 1) {
	           imcomp_double2nan((double *) tiledata, tilelen, (LONGLONG *) tiledata,
	               *(double *) (nullflagval), status);
	       }

#if BYTESWAPPED
               ffswap8((double*) tiledata, tilelen); 
#endif
               compress2mem_from_mem((char *) tiledata, tilelen * sizeof(double),
                    (char **) &cbuf,  &clen, realloc, &gzip_nelem, status);
        }

	/* Write the compressed byte stream. */
        ffpclb(outfptr, (outfptr->Fptr)->cn_gzip_data, row, 1,
             gzip_nelem, (unsigned char *) cbuf, status);

        free(cbuf);  /* finished with this buffer */
    }

    return(*status);
}

/*--------------------------------------------------------------------------*/
int imcomp_write_nocompress_tile(fitsfile *outfptr,
    long row,
    int datatype, 
    void *tiledata, 
    long tilelen,
    int nullcheck,
    void *nullflagval,
    int *status)
{
    char coltype[4];

    /* Write the uncompressed image tile pixels to the tile-compressed image file. */
    /* This is a special case when using NOCOMPRESS for diagnostic purposes in fpack. */ 
    /* Currently, this only supports a limited number of data types and */
    /* does not fully support null-valued pixels in the image. */

    if ((outfptr->Fptr)->cn_uncompressed < 1) {
        /* uncompressed data column doesn't exist, so append new column to table */
        if (datatype == TSHORT) {
	    strcpy(coltype, "1PI");
	} else if (datatype == TINT) {
	    strcpy(coltype, "1PJ");
	} else if (datatype == TFLOAT) {
	    strcpy(coltype, "1PE");
        } else {
	    ffpmsg("NOCOMPRESSION option only supported for int*2, int*4, and float*4 images");
            return(*status = DATA_COMPRESSION_ERR);
        }

        fits_insert_col(outfptr, 999, "UNCOMPRESSED_DATA", coltype, status); /* create column */
    }

    fits_get_colnum(outfptr, CASEINSEN, "UNCOMPRESSED_DATA",
                    &(outfptr->Fptr)->cn_uncompressed, status);  /* save col. num. */
    
    fits_write_col(outfptr, datatype, (outfptr->Fptr)->cn_uncompressed, row, 1,
                      tilelen, tiledata, status);  /* write the tile data */
    return (*status);
}
 /*--------------------------------------------------------------------------*/
int imcomp_convert_tile_tshort(
    fitsfile *outfptr,
    void *tiledata, 
    long tilelen,
    int nullcheck,
    void *nullflagval,
    int nullval,
    int zbitpix,
    double scale,
    double zero,
    double actual_bzero,
    int *intlength,
    int *status)
{
    /*  Prepare the input tile array of pixels for compression.
    /*  Convert input integer*2 tile array in place to 4 or 8-byte ints for compression, */
    /*  If needed, convert 4 or 8-byte ints and do null value substitution. */
    /*  Note that the calling routine must have allocated the input array big enough */
    /* to be able to do this.  */

    short *sbuff;
    int flagval, *idata;
    long ii;
    
       /* We only support writing this integer*2 tile data to a FITS image with 
          BITPIX = 16 and with BZERO = 0 and BSCALE = 1.  */
	  
       if (zbitpix != SHORT_IMG || scale != 1.0 || zero != 0.0) {
           ffpmsg("Datatype conversion/scaling is not supported when writing to compressed images");
           return(*status = DATA_COMPRESSION_ERR);
       } 

       sbuff = (short *) tiledata;
       idata = (int *) tiledata;
       
       if ( (outfptr->Fptr)->compress_type == RICE_1 || (outfptr->Fptr)->compress_type == GZIP_1
         || (outfptr->Fptr)->compress_type == GZIP_2 || (outfptr->Fptr)->compress_type == BZIP2_1 ) 
       {
           /* don't have to convert to int if using gzip, bzip2 or Rice compression */
           *intlength = 2;
             
           if (nullcheck == 1) {
               /* reset pixels equal to flagval to the FITS null value, prior to compression */
               flagval = *(short *) (nullflagval);
               if (flagval != nullval) {
                  for (ii = tilelen - 1; ii >= 0; ii--) {
	            if (sbuff[ii] == (short) flagval)
		       sbuff[ii] = (short) nullval;
                  }
               }
           }
       } else if ((outfptr->Fptr)->compress_type == HCOMPRESS_1) {
           /* have to convert to int if using HCOMPRESS */
           *intlength = 4;

           if (nullcheck == 1) {
               /* reset pixels equal to flagval to the FITS null value, prior to compression */
               flagval = *(short *) (nullflagval);
               for (ii = tilelen - 1; ii >= 0; ii--) {
	            if (sbuff[ii] == (short) flagval)
		       idata[ii] = nullval;
                    else
                       idata[ii] = (int) sbuff[ii];
               }
           } else {  /* just do the data type conversion to int */
               for (ii = tilelen - 1; ii >= 0; ii--) 
                   idata[ii] = (int) sbuff[ii];
           }
       } else {
           /* have to convert to int if using PLIO */
           *intlength = 4;
           if (zero == 0. && actual_bzero == 32768.) {
             /* Here we are compressing unsigned 16-bit integers that have */
	     /* been offset by -32768 using the standard FITS convention. */
	     /* Since PLIO cannot deal with negative values, we must apply */
	     /* the shift of 32786 to the values to make them all positive. */
	     /* The inverse negative shift will be applied in */
	     /* imcomp_decompress_tile when reading the compressed tile. */
             if (nullcheck == 1) {
               /* reset pixels equal to flagval to the FITS null value, prior to compression */
               flagval = *(short *) (nullflagval);
               for (ii = tilelen - 1; ii >= 0; ii--) {
	            if (sbuff[ii] == (short) flagval)
		       idata[ii] = nullval;
                    else
                       idata[ii] = (int) sbuff[ii] + 32768;
               }
             } else {  /* just do the data type conversion to int */
               for (ii = tilelen - 1; ii >= 0; ii--) 
                   idata[ii] = (int) sbuff[ii] + 32768;
             }
           } else {
	     /* This is not an unsigned 16-bit integer array, so process normally */
             if (nullcheck == 1) {
               /* reset pixels equal to flagval to the FITS null value, prior to compression */
               flagval = *(short *) (nullflagval);
               for (ii = tilelen - 1; ii >= 0; ii--) {
	            if (sbuff[ii] == (short) flagval)
		       idata[ii] = nullval;
                    else
                       idata[ii] = (int) sbuff[ii];
               }
             } else {  /* just do the data type conversion to int */
               for (ii = tilelen - 1; ii >= 0; ii--) 
                   idata[ii] = (int) sbuff[ii];
             }
           }
        }
        return(*status);
}
 /*--------------------------------------------------------------------------*/
int imcomp_convert_tile_tushort(
    fitsfile *outfptr,
    void *tiledata, 
    long tilelen,
    int nullcheck,
    void *nullflagval,
    int nullval,
    int zbitpix,
    double scale,
    double zero,
    int *intlength,
    int *status)
{
    /*  Prepare the input  tile array of pixels for compression.
    /*  Convert input integer*2 tile array in place to 4 or 8-byte ints for compression, */
    /*  If needed, convert 4 or 8-byte ints and do null value substitution. */
    /*  Note that the calling routine must have allocated the input array big enough */
    /* to be able to do this.  */

    unsigned short *usbuff;
    short *sbuff;
    int flagval, *idata;
    long ii;
    
       /* datatype of input array is unsigned short.  We only support writing this datatype
          to a FITS image with BITPIX = 16 and with BZERO = 0 and BSCALE = 32768.  */

       if (zbitpix != SHORT_IMG || scale != 1.0 || zero != 32768.) {
           ffpmsg("Implicit datatype conversion is not supported when writing to compressed images");
           return(*status = DATA_COMPRESSION_ERR);
       } 

       usbuff = (unsigned short *) tiledata;
       sbuff = (short *) tiledata;
       idata = (int *) tiledata;

       if ((outfptr->Fptr)->compress_type == RICE_1 || (outfptr->Fptr)->compress_type == GZIP_1
        || (outfptr->Fptr)->compress_type == GZIP_2 || (outfptr->Fptr)->compress_type == BZIP2_1) 
       {
           /* don't have to convert to int if using gzip, bzip2, or Rice compression */
           *intlength = 2;

          /* offset the unsigned value by -32768 to a signed short value. */
	  /* It is more efficient to do this by just flipping the most significant of the 16 bits */

           if (nullcheck == 1) {
               /* reset pixels equal to flagval to the FITS null value, prior to compression  */
               flagval = *(unsigned short *) (nullflagval);
               for (ii = tilelen - 1; ii >= 0; ii--) {
	            if (usbuff[ii] == (unsigned short) flagval)
		       sbuff[ii] = (short) nullval;
                    else
		       usbuff[ii] =  (usbuff[ii]) ^ 0x8000;
               }
           } else {
               /* just offset the pixel values by 32768 (by flipping the MSB */
               for (ii = tilelen - 1; ii >= 0; ii--)
		       usbuff[ii] =  (usbuff[ii]) ^ 0x8000;
           }
       } else {
           /* have to convert to int if using HCOMPRESS or PLIO */
           *intlength = 4;

           if (nullcheck == 1) {
               /* offset the pixel values by 32768, and */
               /* reset pixels equal to flagval to nullval */
               flagval = *(unsigned short *) (nullflagval);
               for (ii = tilelen - 1; ii >= 0; ii--) {
	            if (usbuff[ii] == (unsigned short) flagval)
		       idata[ii] = nullval;
                    else
		       idata[ii] = ((int) usbuff[ii]) - 32768;
               }
           } else {  /* just do the data type conversion to int */
               for (ii = tilelen - 1; ii >= 0; ii--)
		       idata[ii] = ((int) usbuff[ii]) - 32768;
           }
        }

        return(*status);
}
 /*--------------------------------------------------------------------------*/
int imcomp_convert_tile_tint(
    fitsfile *outfptr,
    void *tiledata, 
    long tilelen,
    int nullcheck,
    void *nullflagval,
    int nullval,
    int zbitpix,
    double scale,
    double zero,
    int *intlength,
    int *status)
{
    /*  Prepare the input tile array of pixels for compression.
    /*  Convert input integer*2 tile array in place to 4 or 8-byte ints for compression, */
    /*  If needed, convert 4 or 8-byte ints and do null value substitution. */
    /*  Note that the calling routine must have allocated the input array big enough */
    /* to be able to do this.  */

    int flagval, *idata;
    long ii;
    
 
        /* datatype of input array is int.  We only support writing this datatype
           to a FITS image with BITPIX = 32 and with BZERO = 0 and BSCALE = 1.  */

       if (zbitpix != LONG_IMG || scale != 1.0 || zero != 0.) {
           ffpmsg("Implicit datatype conversion is not supported when writing to compressed images");
           return(*status = DATA_COMPRESSION_ERR);
       } 

       idata = (int *) tiledata;
       *intlength = 4;

       if (nullcheck == 1) {
               /* no datatype conversion is required for any of the compression algorithms,
	         except possibly for HCOMPRESS (to I*8), which is handled later.
		 Just reset pixels equal to flagval to the FITS null value */
               flagval = *(int *) (nullflagval);
               if (flagval != nullval) {
                  for (ii = tilelen - 1; ii >= 0; ii--) {
	            if (idata[ii] == flagval)
		       idata[ii] = nullval;
                  }
               }
       }

       return(*status);
}
 /*--------------------------------------------------------------------------*/
int imcomp_convert_tile_tuint(
    fitsfile *outfptr,
    void *tiledata, 
    long tilelen,
    int nullcheck,
    void *nullflagval,
    int nullval,
    int zbitpix,
    double scale,
    double zero,
    int *intlength,
    int *status)
{
    /*  Prepare the input tile array of pixels for compression.
    /*  Convert input integer*2 tile array in place to 4 or 8-byte ints for compression, */
    /*  If needed, convert 4 or 8-byte ints and do null value substitution. */
    /*  Note that the calling routine must have allocated the input array big enough */
    /* to be able to do this.  */

    int flagval, *idata;
    unsigned int *uintbuff, uintflagval;
    long ii;
    
 
       /* datatype of input array is unsigned int.  We only support writing this datatype
          to a FITS image with BITPIX = 32 and with BZERO = 0 and BSCALE = 2147483648.  */

       if (zbitpix != LONG_IMG || scale != 1.0 || zero != 2147483648.) {
           ffpmsg("Implicit datatype conversion is not supported when writing to compressed images");
           return(*status = DATA_COMPRESSION_ERR);
       } 

       *intlength = 4;
       idata = (int *) tiledata;
       uintbuff = (unsigned int *) tiledata;

       /* offset the unsigned value by -2147483648 to a signed int value. */
       /* It is more efficient to do this by just flipping the most significant of the 32 bits */

       if (nullcheck == 1) {
               /* reset pixels equal to flagval to nullval and */
               /* offset the other pixel values (by flipping the MSB) */
               uintflagval = *(unsigned int *) (nullflagval);
               for (ii = tilelen - 1; ii >= 0; ii--) {
	            if (uintbuff[ii] == uintflagval)
		       idata[ii] = nullval;
                    else
		       uintbuff[ii] = (uintbuff[ii]) ^ 0x80000000;
               }
       } else {
               /* just offset the pixel values (by flipping the MSB) */
               for (ii = tilelen - 1; ii >= 0; ii--)
		       uintbuff[ii] = (uintbuff[ii]) ^ 0x80000000;
       }

       return(*status);
}
 /*--------------------------------------------------------------------------*/
int imcomp_convert_tile_tbyte(
    fitsfile *outfptr,
    void *tiledata, 
    long tilelen,
    int nullcheck,
    void *nullflagval,
    int nullval,
    int zbitpix,
    double scale,
    double zero,
    int *intlength,
    int *status)
{
    /*  Prepare the input tile array of pixels for compression.
    /*  Convert input integer*2 tile array in place to 4 or 8-byte ints for compression, */
    /*  If needed, convert 4 or 8-byte ints and do null value substitution. */
    /*  Note that the calling routine must have allocated the input array big enough */
    /* to be able to do this.  */

    int flagval, *idata;
    long ii;
    unsigned char *usbbuff;
        
       /* datatype of input array is unsigned byte.  We only support writing this datatype
          to a FITS image with BITPIX = 8 and with BZERO = 0 and BSCALE = 1.  */

       if (zbitpix != BYTE_IMG || scale != 1.0 || zero != 0.) {
           ffpmsg("Implicit datatype conversion is not supported when writing to compressed images");
           return(*status = DATA_COMPRESSION_ERR);
       } 

       idata = (int *) tiledata;
       usbbuff = (unsigned char *) tiledata;

       if ( (outfptr->Fptr)->compress_type == RICE_1 || (outfptr->Fptr)->compress_type == GZIP_1
         || (outfptr->Fptr)->compress_type == GZIP_2 || (outfptr->Fptr)->compress_type == BZIP2_1 ) 
       {
           /* don't have to convert to int if using gzip, bzip2, or Rice compression */
           *intlength = 1;
             
           if (nullcheck == 1) {
               /* reset pixels equal to flagval to the FITS null value, prior to compression */
               flagval = *(unsigned char *) (nullflagval);
               if (flagval != nullval) {
                  for (ii = tilelen - 1; ii >= 0; ii--) {
	            if (usbbuff[ii] == (unsigned char) flagval)
		       usbbuff[ii] = (unsigned char) nullval;
                    }
               }
           }
       } else {
           /* have to convert to int if using HCOMPRESS or PLIO */
           *intlength = 4;

           if (nullcheck == 1) {
               /* reset pixels equal to flagval to the FITS null value, prior to compression */
               flagval = *(unsigned char *) (nullflagval);
               for (ii = tilelen - 1; ii >= 0; ii--) {
	            if (usbbuff[ii] == (unsigned char) flagval)
		       idata[ii] = nullval;
                    else
                       idata[ii] = (int) usbbuff[ii];
               }
           } else {  /* just do the data type conversion to int */
               for (ii = tilelen - 1; ii >= 0; ii--) 
                   idata[ii] = (int) usbbuff[ii];
           }
       }

       return(*status);
}
 /*--------------------------------------------------------------------------*/
int imcomp_convert_tile_tsbyte(
    fitsfile *outfptr,
    void *tiledata, 
    long tilelen,
    int nullcheck,
    void *nullflagval,
    int nullval,
    int zbitpix,
    double scale,
    double zero,
    int *intlength,
    int *status)
{
    /*  Prepare the input tile array of pixels for compression.
    /*  Convert input integer*2 tile array in place to 4 or 8-byte ints for compression, */
    /*  If needed, convert 4 or 8-byte ints and do null value substitution. */
    /*  Note that the calling routine must have allocated the input array big enough */
    /* to be able to do this.  */

    int flagval, *idata;
    long ii;
    signed char *sbbuff;
 
       /* datatype of input array is signed byte.  We only support writing this datatype
          to a FITS image with BITPIX = 8 and with BZERO = 0 and BSCALE = -128.  */

       if (zbitpix != BYTE_IMG|| scale != 1.0 || zero != -128.) {
           ffpmsg("Implicit datatype conversion is not supported when writing to compressed images");
           return(*status = DATA_COMPRESSION_ERR);
       }

       idata = (int *) tiledata;
       sbbuff = (signed char *) tiledata;

       if ( (outfptr->Fptr)->compress_type == RICE_1 || (outfptr->Fptr)->compress_type == GZIP_1
         || (outfptr->Fptr)->compress_type == GZIP_2 || (outfptr->Fptr)->compress_type == BZIP2_1 ) 
       {
           /* don't have to convert to int if using gzip, bzip2 or Rice compression */
           *intlength = 1;
             
           if (nullcheck == 1) {
               /* reset pixels equal to flagval to the FITS null value, prior to compression */
               /* offset the other pixel values (by flipping the MSB) */

               flagval = *(signed char *) (nullflagval);
               for (ii = tilelen - 1; ii >= 0; ii--) {
	            if (sbbuff[ii] == (signed char) flagval)
		       sbbuff[ii] = (signed char) nullval;
                    else
		       sbbuff[ii] = (sbbuff[ii]) ^ 0x80;               }
           } else {  /* just offset the pixel values (by flipping the MSB) */
               for (ii = tilelen - 1; ii >= 0; ii--) 
		       sbbuff[ii] = (sbbuff[ii]) ^ 0x80;
           }

       } else {
           /* have to convert to int if using HCOMPRESS or PLIO */
           *intlength = 4;

           if (nullcheck == 1) {
               /* reset pixels equal to flagval to the FITS null value, prior to compression */
               flagval = *(signed char *) (nullflagval);
               for (ii = tilelen - 1; ii >= 0; ii--) {
	            if (sbbuff[ii] == (signed char) flagval)
		       idata[ii] = nullval;
                    else
                       idata[ii] = ((int) sbbuff[ii]) + 128;
               }
           } else {  /* just do the data type conversion to int */
               for (ii = tilelen - 1; ii >= 0; ii--) 
                   idata[ii] = ((int) sbbuff[ii]) + 128;
           }
       }
 
       return(*status);
}
 /*--------------------------------------------------------------------------*/
int imcomp_convert_tile_tfloat(
    fitsfile *outfptr,
    long row,
    void *tiledata, 
    long tilelen,
    long tilenx,
    long tileny,
    int nullcheck,
    void *nullflagval,
    int nullval,
    int zbitpix,
    double scale,
    double zero,
    int *intlength,
    int *flag,
    double *bscale,
    double *bzero,
    int *status)
{
    /*  Prepare the input tile array of pixels for compression.
    /*  Convert input integer*2 tile array in place to 4 or 8-byte ints for compression, */
    /*  If needed, convert 4 or 8-byte ints and do null value substitution. */
    /*  Note that the calling routine must have allocated the input array big enough */
    /* to be able to do this.  */

    int flagval, *idata;
    long irow, ii;
    float floatnull;
    unsigned char *usbbuff;
    unsigned long dithersum;
    int iminval = 0, imaxval = 0;  /* min and max quantized integers */

           *intlength = 4;
           idata = (int *) tiledata;

          /* if the tile-compressed table contains zscale and zzero columns */
          /* then scale and quantize the input floating point data.    */

          if ((outfptr->Fptr)->cn_zscale > 0) {
	    /* quantize the float values into integers */

            if (nullcheck == 1)
	      floatnull = *(float *) (nullflagval);
	    else
	      floatnull = FLOATNULLVALUE;  /* NaNs are represented by this, by default */

            if ((outfptr->Fptr)->quantize_dither == SUBTRACTIVE_DITHER_1) {
	      
	          /* see if the dithering offset value needs to be initialized */                  
	          if ((outfptr->Fptr)->request_dither_offset == 0 && (outfptr->Fptr)->dither_offset == 0) {

		     /* This means randomly choose the dithering offset based on the system time. */
		     /* The offset will have a value between 1 and 10000, inclusive. */
		     /* The time function returns an integer value that is incremented each second. */
		     /* The clock function returns the elapsed CPU time, in integer CLOCKS_PER_SEC units. */
		     /* The CPU time returned by clock is typically (on linux PC) only good to 0.01 sec */
		     /* Summing the 2 quantities may help avoid cases where 2 executions of the program */
		     /* (perhaps in a multithreaded environoment) end up with exactly the same dither_offset */
		     /* value.  The sum is incremented by the current HDU number in the file to provide */
		     /* further randomization.  This randomization is desireable if multiple compressed */
		     /* images will be summed (or differenced). In such cases, the benefits of dithering */
		     /* may be lost if all the images use exactly the same sequence of random numbers when */
		     /* calculating the dithering offsets. */	     
		     
		     (outfptr->Fptr)->dither_offset = 
		       (( (int)time(NULL) + ( (int) clock() / (int) (CLOCKS_PER_SEC / 100)) + (outfptr->Fptr)->curhdu) % 10000) + 1;
		     
                     /* update the header keyword with this new value */
		     fits_update_key(outfptr, TINT, "ZDITHER0", &((outfptr->Fptr)->dither_offset), 
	                        NULL, status);

	          } else if ((outfptr->Fptr)->request_dither_offset < 0 && (outfptr->Fptr)->dither_offset < 0) {

		     /* this means randomly choose the dithering offset based on some hash function */
		     /* of the first input tile of data to be quantized and compressed.  This ensures that */
                     /* the same offset value is used for a given image every time it is compressed. */

		     usbbuff = (unsigned char *) tiledata;
		     dithersum = 0;
		     for (ii = 0; ii < 4 * tilelen; ii++) {
		         dithersum += usbbuff[ii];  /* doesn't matter if there is an integer overflow */
	             }
		     (outfptr->Fptr)->dither_offset = ((int) (dithersum % 10000)) + 1;
		
                     /* update the header keyword with this new value */
		     fits_update_key(outfptr, TINT, "ZDITHER0", &((outfptr->Fptr)->dither_offset), 
	                        NULL, status);
		  }

                  /* subtract 1 to convert from 1-based to 0-based element number */
	          irow = row + (outfptr->Fptr)->dither_offset - 1; /* dither the quantized values */

	      } else {
	          irow = 0;  /* do not dither the quantized values */
              }
	      
              *flag = fits_quantize_float (irow, (float *) tiledata, tilenx, tileny,
                   nullcheck, floatnull, (outfptr->Fptr)->quantize_level, idata,
                   bscale, bzero, &iminval, &imaxval);

              if (*flag > 1)
		   return(*status = *flag);
          }
          else if ((outfptr->Fptr)->quantize_level != NO_QUANTIZE)
	  {
	    /* if floating point pixels are not being losslessly compressed, then */
	    /* input float data is implicitly converted (truncated) to integers */
            if ((scale != 1. || zero != 0.))  /* must scale the values */
	       imcomp_nullscalefloats((float *) tiledata, tilelen, idata, scale, zero,
	           nullcheck, *(float *) (nullflagval), nullval, status);
             else
	       imcomp_nullfloats((float *) tiledata, tilelen, idata,
	           nullcheck, *(float *) (nullflagval), nullval,  status);
          }
          else if ((outfptr->Fptr)->quantize_level == NO_QUANTIZE)
	  {
	      /* just convert null values to NaNs in place, if necessary, then do lossless gzip compression */
		if (nullcheck == 1) {
	            imcomp_float2nan((float *) tiledata, tilelen, (int *) tiledata,
	                *(float *) (nullflagval), status);
		}
          }

          return(*status);
}
 /*--------------------------------------------------------------------------*/
int imcomp_convert_tile_tdouble(
    fitsfile *outfptr,
    long row,
    void *tiledata, 
    long tilelen,
    long tilenx,
    long tileny,
    int nullcheck,
    void *nullflagval,
    int nullval,
    int zbitpix,
    double scale,
    double zero,
    int *intlength,
    int *flag,
    double *bscale,
    double *bzero,
    int *status)
{
    /*  Prepare the input tile array of pixels for compression.
    /*  Convert input integer*2 tile array in place to 4 or 8-byte ints for compression, */
    /*  If needed, convert 4 or 8-byte ints and do null value substitution. */
    /*  Note that the calling routine must have allocated the input array big enough */
    /* to be able to do this.  */

    int flagval, *idata;
    long irow, ii;
    double doublenull;
    unsigned char *usbbuff;
    unsigned long dithersum;
    int iminval = 0, imaxval = 0;  /* min and max quantized integers */

           *intlength = 4;
           idata = (int *) tiledata;

          /* if the tile-compressed table contains zscale and zzero columns */
          /* then scale and quantize the input floating point data.    */
          /* Otherwise, just truncate the floats to integers.          */

          if ((outfptr->Fptr)->cn_zscale > 0)
          {
            if (nullcheck == 1)
	      doublenull = *(double *) (nullflagval);
	    else
	      doublenull = DOUBLENULLVALUE;
	      
            /* quantize the double values into integers */
              if ((outfptr->Fptr)->quantize_dither == SUBTRACTIVE_DITHER_1) {

	          /* see if the dithering offset value needs to be initialized (see above) */                  
	          if ((outfptr->Fptr)->request_dither_offset == 0 && (outfptr->Fptr)->dither_offset == 0) {

		     (outfptr->Fptr)->dither_offset = 
		       (( (int)time(NULL) + ( (int) clock() / (int) (CLOCKS_PER_SEC / 100)) + (outfptr->Fptr)->curhdu) % 10000) + 1;
		     
                     /* update the header keyword with this new value */
		     fits_update_key(outfptr, TINT, "ZDITHER0", &((outfptr->Fptr)->dither_offset), 
	                        NULL, status);

	          } else if ((outfptr->Fptr)->request_dither_offset < 0 && (outfptr->Fptr)->dither_offset < 0) {

		     usbbuff = (unsigned char *) tiledata;
		     dithersum = 0;
		     for (ii = 0; ii < 8 * tilelen; ii++) {
		         dithersum += usbbuff[ii];
	             }
		     (outfptr->Fptr)->dither_offset = ((int) (dithersum % 10000)) + 1;
		
                     /* update the header keyword with this new value */
		     fits_update_key(outfptr, TINT, "ZDITHER0", &((outfptr->Fptr)->dither_offset), 
	                        NULL, status);
		  }

	          irow = row + (outfptr->Fptr)->dither_offset - 1; /* dither the quantized values */

	      } else {
	          irow = 0;  /* do not dither the quantized values */
	      }

            *flag = fits_quantize_double (irow, (double *) tiledata, tilenx, tileny,
               nullcheck, doublenull, (outfptr->Fptr)->quantize_level, idata,
               bscale, bzero, &iminval, &imaxval);

            if (*flag > 1)
		return(*status = *flag);
          }
          else if ((outfptr->Fptr)->quantize_level != NO_QUANTIZE)
	  {
	    /* if floating point pixels are not being losslessly compressed, then */
	    /* input float data is implicitly converted (truncated) to integers */
             if ((scale != 1. || zero != 0.))  /* must scale the values */
	       imcomp_nullscaledoubles((double *) tiledata, tilelen, idata, scale, zero,
	           nullcheck, *(double *) (nullflagval), nullval, status);
             else
	       imcomp_nulldoubles((double *) tiledata, tilelen, idata,
	           nullcheck, *(double *) (nullflagval), nullval,  status);
          }
          else if ((outfptr->Fptr)->quantize_level == NO_QUANTIZE)
	  {
	      /* just convert null values to NaNs in place, if necessary, then do lossless gzip compression */
		if (nullcheck == 1) {
	            imcomp_double2nan((double *) tiledata, tilelen, (LONGLONG *) tiledata,
	                *(double *) (nullflagval), status);
		}
          }
 
          return(*status);
}
/*---------------------------------------------------------------------------*/
int imcomp_nullscale(
     int *idata, 
     long tilelen,
     int nullflagval,
     int nullval,
     double scale,
     double zero,
     int *status)
/*
   do null value substitution AND scaling of the integer array.
   If array value = nullflagval, then set the value to nullval.
   Otherwise, inverse scale the integer value.
*/
{
    long ii;
    double dvalue;
    
    for (ii=0; ii < tilelen; ii++)
    {
        if (idata[ii] == nullflagval)
	    idata[ii] = nullval;
	else 
	{
            dvalue = (idata[ii] - zero) / scale;

            if (dvalue < DINT_MIN)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = INT32_MIN;
            }
            else if (dvalue > DINT_MAX)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = INT32_MAX;
            }
            else
            {
                if (dvalue >= 0)
                    idata[ii] = (int) (dvalue + .5);
                else
                    idata[ii] = (int) (dvalue - .5);
            }
        }
    }
    return(*status);
}
/*---------------------------------------------------------------------------*/
int imcomp_nullvalues(
     int *idata, 
     long tilelen,
     int nullflagval,
     int nullval,
     int *status)
/*
   do null value substitution.
   If array value = nullflagval, then set the value to nullval.
*/
{
    long ii;
    
    for (ii=0; ii < tilelen; ii++)
    {
        if (idata[ii] == nullflagval)
	    idata[ii] = nullval;
    }
    return(*status);
}
/*---------------------------------------------------------------------------*/
int imcomp_scalevalues(
     int *idata, 
     long tilelen,
     double scale,
     double zero,
     int *status)
/*
   do inverse scaling the integer values.
*/
{
    long ii;
    double dvalue;
    
    for (ii=0; ii < tilelen; ii++)
    {
            dvalue = (idata[ii] - zero) / scale;

            if (dvalue < DINT_MIN)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = INT32_MIN;
            }
            else if (dvalue > DINT_MAX)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = INT32_MAX;
            }
            else
            {
                if (dvalue >= 0)
                    idata[ii] = (int) (dvalue + .5);
                else
                    idata[ii] = (int) (dvalue - .5);
            }
    }
    return(*status);
}
/*---------------------------------------------------------------------------*/
int imcomp_nullscalei2(
     short *idata, 
     long tilelen,
     short nullflagval,
     short nullval,
     double scale,
     double zero,
     int *status)
/*
   do null value substitution AND scaling of the integer array.
   If array value = nullflagval, then set the value to nullval.
   Otherwise, inverse scale the integer value.
*/
{
    long ii;
    double dvalue;
    
    for (ii=0; ii < tilelen; ii++)
    {
        if (idata[ii] == nullflagval)
	    idata[ii] = nullval;
	else 
	{
            dvalue = (idata[ii] - zero) / scale;

            if (dvalue < DSHRT_MIN)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = SHRT_MIN;
            }
            else if (dvalue > DSHRT_MAX)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = SHRT_MAX;
            }
            else
            {
                if (dvalue >= 0)
                    idata[ii] = (int) (dvalue + .5);
                else
                    idata[ii] = (int) (dvalue - .5);
            }
        }
    }
    return(*status);
}
/*---------------------------------------------------------------------------*/
int imcomp_nullvaluesi2(
     short *idata, 
     long tilelen,
     short nullflagval,
     short nullval,
     int *status)
/*
   do null value substitution.
   If array value = nullflagval, then set the value to nullval.
*/
{
    long ii;
    
    for (ii=0; ii < tilelen; ii++)
    {
        if (idata[ii] == nullflagval)
	    idata[ii] = nullval;
    }
    return(*status);
}
/*---------------------------------------------------------------------------*/
int imcomp_scalevaluesi2(
     short *idata, 
     long tilelen,
     double scale,
     double zero,
     int *status)
/*
   do inverse scaling the integer values.
*/
{
    long ii;
    double dvalue;
    
    for (ii=0; ii < tilelen; ii++)
    {
            dvalue = (idata[ii] - zero) / scale;

            if (dvalue < DSHRT_MIN)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = SHRT_MIN;
            }
            else if (dvalue > DSHRT_MAX)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = SHRT_MAX;
            }
            else
            {
                if (dvalue >= 0)
                    idata[ii] = (int) (dvalue + .5);
                else
                    idata[ii] = (int) (dvalue - .5);
            }
    }
    return(*status);
}
/*---------------------------------------------------------------------------*/
int imcomp_nullfloats(
     float *fdata,
     long tilelen,
     int *idata, 
     int nullcheck,
     float nullflagval,
     int nullval,
     int *status)
/*
   do null value substitution  of the float array.
   If array value = nullflagval, then set the output value to FLOATNULLVALUE.
*/
{
    long ii;
    double dvalue;
    
    if (nullcheck == 1) /* must check for null values */
    {
      for (ii=0; ii < tilelen; ii++)
      {
        if (fdata[ii] == nullflagval)
	    idata[ii] = nullval;
	else 
	{
            dvalue = fdata[ii];

            if (dvalue < DINT_MIN)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = INT32_MIN;
            }
            else if (dvalue > DINT_MAX)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = INT32_MAX;
            }
            else
            {
                if (dvalue >= 0)
                    idata[ii] = (int) (dvalue + .5);
                else
                    idata[ii] = (int) (dvalue - .5);
            }
        }
      }
    }
    else  /* don't have to worry about null values */
    {
      for (ii=0; ii < tilelen; ii++)
      {
            dvalue = fdata[ii];

            if (dvalue < DINT_MIN)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = INT32_MIN;
            }
            else if (dvalue > DINT_MAX)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = INT32_MAX;
            }
            else
            {
                if (dvalue >= 0)
                    idata[ii] = (int) (dvalue + .5);
                else
                    idata[ii] = (int) (dvalue - .5);
            }
      }
    }
    return(*status);
}
/*---------------------------------------------------------------------------*/
int imcomp_nullscalefloats(
     float *fdata,
     long tilelen,
     int *idata, 
     double scale,
     double zero,
     int nullcheck,
     float nullflagval,
     int nullval,
     int *status)
/*
   do null value substitution  of the float array.
   If array value = nullflagval, then set the output value to FLOATNULLVALUE.
   Otherwise, inverse scale the integer value.
*/
{
    long ii;
    double dvalue;
    
    if (nullcheck == 1) /* must check for null values */
    {
      for (ii=0; ii < tilelen; ii++)
      {
        if (fdata[ii] == nullflagval)
	    idata[ii] = nullval;
	else 
	{
            dvalue = (fdata[ii] - zero) / scale;

            if (dvalue < DINT_MIN)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = INT32_MIN;
            }
            else if (dvalue > DINT_MAX)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = INT32_MAX;
            }
            else
            {
                if (dvalue >= 0.)
                    idata[ii] = (int) (dvalue + .5);
                else
                    idata[ii] = (int) (dvalue - .5);
            }
        }
      }
    }
    else  /* don't have to worry about null values */
    {
      for (ii=0; ii < tilelen; ii++)
      {
            dvalue = (fdata[ii] - zero) / scale;

            if (dvalue < DINT_MIN)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = INT32_MIN;
            }
            else if (dvalue > DINT_MAX)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = INT32_MAX;
            }
            else
            {
                if (dvalue >= 0.)
                    idata[ii] = (int) (dvalue + .5);
                else
                    idata[ii] = (int) (dvalue - .5);
            }
      }
    }
    return(*status);
}
/*---------------------------------------------------------------------------*/
int imcomp_nulldoubles(
     double *fdata,
     long tilelen,
     int *idata, 
     int nullcheck,
     double nullflagval,
     int nullval,
     int *status)
/*
   do null value substitution  of the float array.
   If array value = nullflagval, then set the output value to FLOATNULLVALUE.
   Otherwise, inverse scale the integer value.
*/
{
    long ii;
    double dvalue;
    
    if (nullcheck == 1) /* must check for null values */
    {
      for (ii=0; ii < tilelen; ii++)
      {
        if (fdata[ii] == nullflagval)
	    idata[ii] = nullval;
	else 
	{
            dvalue = fdata[ii];

            if (dvalue < DINT_MIN)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = INT32_MIN;
            }
            else if (dvalue > DINT_MAX)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = INT32_MAX;
            }
            else
            {
                if (dvalue >= 0.)
                    idata[ii] = (int) (dvalue + .5);
                else
                    idata[ii] = (int) (dvalue - .5);
            }
        }
      }
    }
    else  /* don't have to worry about null values */
    {
      for (ii=0; ii < tilelen; ii++)
      {
            dvalue = fdata[ii];

            if (dvalue < DINT_MIN)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = INT32_MIN;
            }
            else if (dvalue > DINT_MAX)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = INT32_MAX;
            }
            else
            {
                if (dvalue >= 0.)
                    idata[ii] = (int) (dvalue + .5);
                else
                    idata[ii] = (int) (dvalue - .5);
            }
      }
    }
    return(*status);
}
/*---------------------------------------------------------------------------*/
int imcomp_nullscaledoubles(
     double *fdata,
     long tilelen,
     int *idata, 
     double scale,
     double zero,
     int nullcheck,
     double nullflagval,
     int nullval,
     int *status)
/*
   do null value substitution  of the float array.
   If array value = nullflagval, then set the output value to FLOATNULLVALUE.
   Otherwise, inverse scale the integer value.
*/
{
    long ii;
    double dvalue;
    
    if (nullcheck == 1) /* must check for null values */
    {
      for (ii=0; ii < tilelen; ii++)
      {
        if (fdata[ii] == nullflagval)
	    idata[ii] = nullval;
	else 
	{
            dvalue = (fdata[ii] - zero) / scale;

            if (dvalue < DINT_MIN)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = INT32_MIN;
            }
            else if (dvalue > DINT_MAX)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = INT32_MAX;
            }
            else
            {
                if (dvalue >= 0.)
                    idata[ii] = (int) (dvalue + .5);
                else
                    idata[ii] = (int) (dvalue - .5);
            }
        }
      }
    }
    else  /* don't have to worry about null values */
    {
      for (ii=0; ii < tilelen; ii++)
      {
            dvalue = (fdata[ii] - zero) / scale;

            if (dvalue < DINT_MIN)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = INT32_MIN;
            }
            else if (dvalue > DINT_MAX)
            {
                *status = OVERFLOW_ERR;
                idata[ii] = INT32_MAX;
            }
            else
            {
                if (dvalue >= 0.)
                    idata[ii] = (int) (dvalue + .5);
                else
                    idata[ii] = (int) (dvalue - .5);
            }
      }
    }
    return(*status);
}
/*---------------------------------------------------------------------------*/
int fits_write_compressed_img(fitsfile *fptr,   /* I - FITS file pointer     */
            int  datatype,   /* I - datatype of the array to be written      */
            long  *infpixel, /* I - 'bottom left corner' of the subsection   */
            long  *inlpixel, /* I - 'top right corner' of the subsection     */
            int  nullcheck,  /* I - 0 for no null checking                   */
                             /*     1: pixels that are = nullval will be     */
                             /*     written with the FITS null pixel value   */
                             /*     (floating point arrays only)             */
            void *array,     /* I - array of values to be written            */
            void *nullval,   /* I - undefined pixel value                    */
            int  *status)    /* IO - error status                            */
/*
   Write a section of a compressed image.
*/
{
    int naxis[MAX_COMPRESS_DIM], tiledim[MAX_COMPRESS_DIM];
    long tilesize[MAX_COMPRESS_DIM], thistilesize[MAX_COMPRESS_DIM];
    long ftile[MAX_COMPRESS_DIM], ltile[MAX_COMPRESS_DIM];
    long tfpixel[MAX_COMPRESS_DIM], tlpixel[MAX_COMPRESS_DIM];
    long rowdim[MAX_COMPRESS_DIM], offset[MAX_COMPRESS_DIM],ntemp;
    long fpixel[MAX_COMPRESS_DIM], lpixel[MAX_COMPRESS_DIM];
    int ii, i5, i4, i3, i2, i1, i0, ndim, irow, pixlen, tilenul;
    int  tstatus, buffpixsiz;
    void *buffer;
    char *bnullarray = 0, card[FLEN_CARD];

    if (*status > 0) 
        return(*status);

    if (!fits_is_compressed_image(fptr, status) )
    {
        ffpmsg("CHDU is not a compressed image (fits_write_compressed_img)");
        return(*status = DATA_COMPRESSION_ERR);
    }

    /* reset position to the correct HDU if necessary */
    if (fptr->HDUposition != (fptr->Fptr)->curhdu)
        ffmahd(fptr, (fptr->HDUposition) + 1, NULL, status);

    /* rescan header if data structure is undefined */
    else if ((fptr->Fptr)->datastart == DATA_UNDEFINED)
        if ( ffrdef(fptr, status) > 0)               
            return(*status);


    /* ===================================================================== */


    if (datatype == TSHORT || datatype == TUSHORT)
    {
       pixlen = sizeof(short);
    }
    else if (datatype == TINT || datatype == TUINT)
    {
       pixlen = sizeof(int);
    }
    else if (datatype == TBYTE || datatype == TSBYTE)
    {
       pixlen = 1;
    }
    else if (datatype == TLONG || datatype == TULONG)
    {
       pixlen = sizeof(long);
    }
    else if (datatype == TFLOAT)
    {
       pixlen = sizeof(float);
    }
    else if (datatype == TDOUBLE)
    {
       pixlen = sizeof(double);
    }
    else
    {
        ffpmsg("unsupported datatype for compressing image");
        return(*status = BAD_DATATYPE);
    }

    /* ===================================================================== */

    /* allocate scratch space for processing one tile of the image */
    buffpixsiz = pixlen;  /* this is the minimum pixel size */
    
    if ( (fptr->Fptr)->compress_type == HCOMPRESS_1) { /* need 4 or 8 bytes per pixel */
        if ((fptr->Fptr)->zbitpix == BYTE_IMG ||
	    (fptr->Fptr)->zbitpix == SHORT_IMG )
                buffpixsiz = maxvalue(buffpixsiz, 4);
        else
	        buffpixsiz = 8;
    }
    else if ( (fptr->Fptr)->compress_type == PLIO_1) { /* need 4 bytes per pixel */
                buffpixsiz = maxvalue(buffpixsiz, 4);
    }
    else if ( (fptr->Fptr)->compress_type == RICE_1  ||
              (fptr->Fptr)->compress_type == GZIP_1 ||
              (fptr->Fptr)->compress_type == GZIP_2 ||
              (fptr->Fptr)->compress_type == BZIP2_1) {  /* need 1, 2, or 4 bytes per pixel */
        if ((fptr->Fptr)->zbitpix == BYTE_IMG)
            buffpixsiz = maxvalue(buffpixsiz, 1);
        else if ((fptr->Fptr)->zbitpix == SHORT_IMG)
            buffpixsiz = maxvalue(buffpixsiz, 2);
        else 
            buffpixsiz = maxvalue(buffpixsiz, 4);
    }
    else
    {
        ffpmsg("unsupported image compression algorithm");
        return(*status = BAD_DATATYPE);
    }
    
    /* cast to double to force alignment on 8-byte addresses */
    buffer = (double *) calloc ((fptr->Fptr)->maxtilelen, buffpixsiz);

    if (buffer == NULL)
    {
	    ffpmsg("Out of memory (fits_write_compress_img)");
	    return (*status = MEMORY_ALLOCATION);
    }

    /* ===================================================================== */

    /* initialize all the arrays */
    for (ii = 0; ii < MAX_COMPRESS_DIM; ii++)
    {
        naxis[ii] = 1;
        tiledim[ii] = 1;
        tilesize[ii] = 1;
        ftile[ii] = 1;
        ltile[ii] = 1;
        rowdim[ii] = 1;
    }

    ndim = (fptr->Fptr)->zndim;
    ntemp = 1;
    for (ii = 0; ii < ndim; ii++)
    {
        fpixel[ii] = infpixel[ii];
        lpixel[ii] = inlpixel[ii];

        /* calc number of tiles in each dimension, and tile containing */
        /* the first and last pixel we want to read in each dimension  */
        naxis[ii] = (fptr->Fptr)->znaxis[ii];
        if (fpixel[ii] < 1)
        {
            free(buffer);
            return(*status = BAD_PIX_NUM);
        }

        tilesize[ii] = (fptr->Fptr)->tilesize[ii];
        tiledim[ii] = (naxis[ii] - 1) / tilesize[ii] + 1;
        ftile[ii]   = (fpixel[ii] - 1)   / tilesize[ii] + 1;
        ltile[ii]   = minvalue((lpixel[ii] - 1) / tilesize[ii] + 1, 
                                tiledim[ii]);
        rowdim[ii]  = ntemp;  /* total tiles in each dimension */
        ntemp *= tiledim[ii];
    }

    /* support up to 6 dimensions for now */
    /* tfpixel and tlpixel are the first and last image pixels */
    /* along each dimension of the compression tile */
    for (i5 = ftile[5]; i5 <= ltile[5]; i5++)
    {
     tfpixel[5] = (i5 - 1) * tilesize[5] + 1;
     tlpixel[5] = minvalue(tfpixel[5] + tilesize[5] - 1, 
                            naxis[5]);
     thistilesize[5] = tlpixel[5] - tfpixel[5] + 1;
     offset[5] = (i5 - 1) * rowdim[5];
     for (i4 = ftile[4]; i4 <= ltile[4]; i4++)
     {
      tfpixel[4] = (i4 - 1) * tilesize[4] + 1;
      tlpixel[4] = minvalue(tfpixel[4] + tilesize[4] - 1, 
                            naxis[4]);
      thistilesize[4] = thistilesize[5] * (tlpixel[4] - tfpixel[4] + 1);
      offset[4] = (i4 - 1) * rowdim[4] + offset[5];
      for (i3 = ftile[3]; i3 <= ltile[3]; i3++)
      {
        tfpixel[3] = (i3 - 1) * tilesize[3] + 1;
        tlpixel[3] = minvalue(tfpixel[3] + tilesize[3] - 1, 
                              naxis[3]);
        thistilesize[3] = thistilesize[4] * (tlpixel[3] - tfpixel[3] + 1);
        offset[3] = (i3 - 1) * rowdim[3] + offset[4];
        for (i2 = ftile[2]; i2 <= ltile[2]; i2++)
        {
          tfpixel[2] = (i2 - 1) * tilesize[2] + 1;
          tlpixel[2] = minvalue(tfpixel[2] + tilesize[2] - 1, 
                                naxis[2]);
          thistilesize[2] = thistilesize[3] * (tlpixel[2] - tfpixel[2] + 1);
          offset[2] = (i2 - 1) * rowdim[2] + offset[3];
          for (i1 = ftile[1]; i1 <= ltile[1]; i1++)
          {
            tfpixel[1] = (i1 - 1) * tilesize[1] + 1;
            tlpixel[1] = minvalue(tfpixel[1] + tilesize[1] - 1, 
                                  naxis[1]);
            thistilesize[1] = thistilesize[2] * (tlpixel[1] - tfpixel[1] + 1);
            offset[1] = (i1 - 1) * rowdim[1] + offset[2];
            for (i0 = ftile[0]; i0 <= ltile[0]; i0++)
            {
              tfpixel[0] = (i0 - 1) * tilesize[0] + 1;
              tlpixel[0] = minvalue(tfpixel[0] + tilesize[0] - 1, 
                                    naxis[0]);
              thistilesize[0] = thistilesize[1] * (tlpixel[0] - tfpixel[0] + 1);
              /* calculate row of table containing this tile */
              irow = i0 + offset[1];

              /* read and uncompress this row (tile) of the table */
              /* also do type conversion and undefined pixel substitution */
              /* at this point */
              imcomp_decompress_tile(fptr, irow, thistilesize[0],
                    datatype, nullcheck, nullval, buffer, bnullarray, &tilenul,
                     status);

              if (*status == NO_COMPRESSED_TILE)
              {
                   /* tile doesn't exist, so initialize to zero */
                   memset(buffer, 0, pixlen * thistilesize[0]);
                   *status = 0;
              }

              /* copy the intersecting pixels to this tile from the input */
              imcomp_merge_overlap(buffer, pixlen, ndim, tfpixel, tlpixel, 
                     bnullarray, array, fpixel, lpixel, nullcheck, status);

              /* compress the tile again, and write it back to the FITS file */
              imcomp_compress_tile (fptr, irow, datatype, buffer, 
                                    thistilesize[0],
				    tlpixel[0] - tfpixel[0] + 1,
				    tlpixel[1] - tfpixel[1] + 1,
				    nullcheck, nullval, 
				    status);
            }
          }
        }
      }
     }
    }
    free(buffer);
    

    if ((fptr->Fptr)->zbitpix < 0 && nullcheck != 0) { 
/*
     This is a floating point FITS image with possible null values.
     It is too messy to test if any null values are actually written, so 
     just assume so.  We need to make sure that the
     ZBLANK keyword is present in the compressed image header.  If it is not
     there then we need to insert the keyword. 
*/   
        tstatus = 0;
        ffgcrd(fptr, "ZBLANK", card, &tstatus);

	if (tstatus) {   /* have to insert the ZBLANK keyword */
           ffgcrd(fptr, "ZCMPTYPE", card, status);
           ffikyj(fptr, "ZBLANK", COMPRESS_NULL_VALUE, 
                "null value in the compressed integer array", status);
	
           /* set this value into the internal structure; it is used if */
	   /* the program reads back the values from the array */
	 
          (fptr->Fptr)->zblank = COMPRESS_NULL_VALUE;
          (fptr->Fptr)->cn_zblank = -1;  /* flag for a constant ZBLANK */
        }  
    }  
    
    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_write_compressed_pixels(fitsfile *fptr, /* I - FITS file pointer   */
            int  datatype,  /* I - datatype of the array to be written      */
            LONGLONG   fpixel,  /* I - 'first pixel to write          */
            LONGLONG   npixel,  /* I - number of pixels to write      */
            int  nullcheck,  /* I - 0 for no null checking                   */
                             /*     1: pixels that are = nullval will be     */
                             /*     written with the FITS null pixel value   */
                             /*     (floating point arrays only)             */
            void *array,      /* I - array of values to write                */
            void *nullval,    /* I - value used to represent undefined pixels*/
            int  *status)     /* IO - error status                           */
/*
   Write a consecutive set of pixels to a compressed image.  This routine
   interpretes the n-dimensional image as a long one-dimensional array. 
   This is actually a rather inconvenient way to write compressed images in
   general, and could be rather inefficient if the requested pixels to be
   written are located in many different image compression tiles.    

   The general strategy used here is to write the requested pixels in blocks
   that correspond to rectangular image sections.  
*/
{
    int naxis, ii, bytesperpixel;
    long naxes[MAX_COMPRESS_DIM], nread;
    LONGLONG tfirst, tlast, last0, last1, dimsize[MAX_COMPRESS_DIM];
    long nplane, firstcoord[MAX_COMPRESS_DIM], lastcoord[MAX_COMPRESS_DIM];
    char *arrayptr;

    if (*status > 0)
        return(*status);

    arrayptr = (char *) array;

    /* get size of array pixels, in bytes */
    bytesperpixel = ffpxsz(datatype);

    for (ii = 0; ii < MAX_COMPRESS_DIM; ii++)
    {
        naxes[ii] = 1;
        firstcoord[ii] = 0;
        lastcoord[ii] = 0;
    }

    /*  determine the dimensions of the image to be written */
    ffgidm(fptr, &naxis, status);
    ffgisz(fptr, MAX_COMPRESS_DIM, naxes, status);

    /* calc the cumulative number of pixels in each successive dimension */
    dimsize[0] = 1;
    for (ii = 1; ii < MAX_COMPRESS_DIM; ii++)
         dimsize[ii] = dimsize[ii - 1] * naxes[ii - 1];

    /*  determine the coordinate of the first and last pixel in the image */
    /*  Use zero based indexes here */
    tfirst = fpixel - 1;
    tlast = tfirst + npixel - 1;
    for (ii = naxis - 1; ii >= 0; ii--)
    {
        firstcoord[ii] = (long) (tfirst / dimsize[ii]);
        lastcoord[ii]  = (long) (tlast / dimsize[ii]);
        tfirst = tfirst - firstcoord[ii] * dimsize[ii];
        tlast = tlast - lastcoord[ii] * dimsize[ii];
    }

    /* to simplify things, treat 1-D, 2-D, and 3-D images as separate cases */

    if (naxis == 1)
    {
        /* Simple: just write the requested range of pixels */

        firstcoord[0] = firstcoord[0] + 1;
        lastcoord[0] = lastcoord[0] + 1;
        fits_write_compressed_img(fptr, datatype, firstcoord, lastcoord,
            nullcheck, array, nullval, status);
        return(*status);
    }
    else if (naxis == 2)
    {
        nplane = 0;  /* write 1st (and only) plane of the image */
        fits_write_compressed_img_plane(fptr, datatype, bytesperpixel,
          nplane, firstcoord, lastcoord, naxes, nullcheck,
          array, nullval, &nread, status);
    }
    else if (naxis == 3)
    {
        /* test for special case: writing an integral number of planes */
        if (firstcoord[0] == 0 && firstcoord[1] == 0 &&
            lastcoord[0] == naxes[0] - 1 && lastcoord[1] == naxes[1] - 1)
        {
            for (ii = 0; ii < MAX_COMPRESS_DIM; ii++)
            {
                /* convert from zero base to 1 base */
                (firstcoord[ii])++;
                (lastcoord[ii])++;
            }

            /* we can write the contiguous block of pixels in one go */
            fits_write_compressed_img(fptr, datatype, firstcoord, lastcoord,
                nullcheck, array, nullval, status);
            return(*status);
        }

        /* save last coordinate in temporary variables */
        last0 = lastcoord[0];
        last1 = lastcoord[1];

        if (firstcoord[2] < lastcoord[2])
        {
            /* we will write up to the last pixel in all but the last plane */
            lastcoord[0] = naxes[0] - 1;
            lastcoord[1] = naxes[1] - 1;
        }

        /* write one plane of the cube at a time, for simplicity */
        for (nplane = firstcoord[2]; nplane <= lastcoord[2]; nplane++)
        {
            if (nplane == lastcoord[2])
            {
                lastcoord[0] = (long) last0;
                lastcoord[1] = (long) last1;
            }

            fits_write_compressed_img_plane(fptr, datatype, bytesperpixel,
              nplane, firstcoord, lastcoord, naxes, nullcheck,
              arrayptr, nullval, &nread, status);

            /* for all subsequent planes, we start with the first pixel */
            firstcoord[0] = 0;
            firstcoord[1] = 0;

            /* increment pointers to next elements to be written */
            arrayptr = arrayptr + nread * bytesperpixel;
        }
    }
    else
    {
        ffpmsg("only 1D, 2D, or 3D images are currently supported");
        return(*status = DATA_COMPRESSION_ERR);
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_write_compressed_img_plane(fitsfile *fptr, /* I - FITS file    */
            int  datatype,  /* I - datatype of the array to be written    */
            int  bytesperpixel, /* I - number of bytes per pixel in array */
            long   nplane,  /* I - which plane of the cube to write      */
            long *firstcoord, /* I coordinate of first pixel to write */
            long *lastcoord,  /* I coordinate of last pixel to write */
            long *naxes,     /* I size of each image dimension */
            int  nullcheck,  /* I - 0 for no null checking                   */
                             /*     1: pixels that are = nullval will be     */
                             /*     written with the FITS null pixel value   */
                             /*     (floating point arrays only)             */
            void *array,      /* I - array of values that are written        */
            void *nullval,    /* I - value for undefined pixels              */
            long *nread,      /* O - total number of pixels written          */
            int  *status)     /* IO - error status                           */

   /*
           in general we have to write the first partial row of the image,
           followed by the middle complete rows, followed by the last
           partial row of the image.  If the first or last rows are complete,
           then write them at the same time as all the middle rows.
    */
{
    /* bottom left coord. and top right coord. */
    long blc[MAX_COMPRESS_DIM], trc[MAX_COMPRESS_DIM]; 
    char *arrayptr;

    *nread = 0;

    arrayptr = (char *) array;

    blc[2] = nplane + 1;
    trc[2] = nplane + 1;

    if (firstcoord[0] != 0)
    { 
            /* have to read a partial first row */
            blc[0] = firstcoord[0] + 1;
            blc[1] = firstcoord[1] + 1;
            trc[1] = blc[1];  
            if (lastcoord[1] == firstcoord[1])
               trc[0] = lastcoord[0] + 1; /* 1st and last pixels in same row */
            else
               trc[0] = naxes[0];  /* read entire rest of the row */

            fits_write_compressed_img(fptr, datatype, blc, trc,
                nullcheck, arrayptr, nullval, status);

            *nread = *nread + trc[0] - blc[0] + 1;

            if (lastcoord[1] == firstcoord[1])
            {
               return(*status);  /* finished */
            }

            /* set starting coord to beginning of next line */
            firstcoord[0] = 0;
            firstcoord[1] += 1;
            arrayptr = arrayptr + (trc[0] - blc[0] + 1) * bytesperpixel;
    }

    /* write contiguous complete rows of the image, if any */
    blc[0] = 1;
    blc[1] = firstcoord[1] + 1;
    trc[0] = naxes[0];

    if (lastcoord[0] + 1 == naxes[0])
    {
            /* can write the last complete row, too */
            trc[1] = lastcoord[1] + 1;
    }
    else
    {
            /* last row is incomplete; have to read it separately */
            trc[1] = lastcoord[1];
    }

    if (trc[1] >= blc[1])  /* must have at least one whole line to read */
    {
        fits_write_compressed_img(fptr, datatype, blc, trc,
                nullcheck, arrayptr, nullval, status);

        *nread = *nread + (trc[1] - blc[1] + 1) * naxes[0];

        if (lastcoord[1] + 1 == trc[1])
               return(*status);  /* finished */

        /* increment pointers for the last partial row */
        arrayptr = arrayptr + (trc[1] - blc[1] + 1) * naxes[0] * bytesperpixel;

     }

    if (trc[1] == lastcoord[1] + 1)
        return(*status);           /* all done */

    /* set starting and ending coord to last line */

    trc[0] = lastcoord[0] + 1;
    trc[1] = lastcoord[1] + 1;
    blc[1] = trc[1];

    fits_write_compressed_img(fptr, datatype, blc, trc,
                nullcheck, arrayptr, nullval, status);

    *nread = *nread + trc[0] - blc[0] + 1;

    return(*status);
}

/* ######################################################################## */
/* ###                 Image Decompression Routines                     ### */
/* ######################################################################## */

/*--------------------------------------------------------------------------*/
int fits_img_decompress (fitsfile *infptr, /* image (bintable) to uncompress */
              fitsfile *outfptr,   /* empty HDU for output uncompressed image */
              int *status)         /* IO - error status               */

/* 
  This routine decompresses the whole image and writes it to the output file.
*/

{
    int ii, datatype = 0;
    int nullcheck, anynul;
    LONGLONG fpixel[MAX_COMPRESS_DIM], lpixel[MAX_COMPRESS_DIM];
    long inc[MAX_COMPRESS_DIM];
    long imgsize;
    float *nulladdr, fnulval;
    double dnulval;

    if (fits_img_decompress_header(infptr, outfptr, status) > 0)
    {
    	return (*status);
    }

    /* force a rescan of the output header keywords, then reset the scaling */
    /* in case the BSCALE and BZERO keywords are present, so that the       */
    /* decompressed values won't be scaled when written to the output image */
    ffrdef(outfptr, status);
    ffpscl(outfptr, 1.0, 0.0, status);
    ffpscl(infptr, 1.0, 0.0, status);

    /* initialize; no null checking is needed for integer images */
    nullcheck = 0;
    nulladdr =  &fnulval;

    /* determine datatype for image */
    if ((infptr->Fptr)->zbitpix == BYTE_IMG)
    {
        datatype = TBYTE;
    }
    else if ((infptr->Fptr)->zbitpix == SHORT_IMG)
    {
        datatype = TSHORT;
    }
    else if ((infptr->Fptr)->zbitpix == LONG_IMG)
    {
        datatype = TINT;
    }
    else if ((infptr->Fptr)->zbitpix == FLOAT_IMG)
    {
        /* In the case of float images we must check for NaNs  */
        nullcheck = 1;
        fnulval = FLOATNULLVALUE;
        nulladdr =  &fnulval;
        datatype = TFLOAT;
    }
    else if ((infptr->Fptr)->zbitpix == DOUBLE_IMG)
    {
        /* In the case of double images we must check for NaNs  */
        nullcheck = 1;
        dnulval = DOUBLENULLVALUE;
        nulladdr = (float *) &dnulval;
        datatype = TDOUBLE;
    }

    /* calculate size of the image (in pixels) */
    imgsize = 1;
    for (ii = 0; ii < (infptr->Fptr)->zndim; ii++)
    {
        imgsize *= (infptr->Fptr)->znaxis[ii];
        fpixel[ii] = 1;              /* Set first and last pixel to */
        lpixel[ii] = (infptr->Fptr)->znaxis[ii]; /* include the entire image. */
        inc[ii] = 1;
    }

    /* uncompress the input image and write to output image, one tile at a time */

    fits_read_write_compressed_img(infptr, datatype, fpixel, lpixel, inc,  
            nullcheck, nulladdr, &anynul, outfptr, status);

    return (*status);
}
/*--------------------------------------------------------------------------*/
int fits_decompress_img (fitsfile *infptr, /* image (bintable) to uncompress */
              fitsfile *outfptr,   /* empty HDU for output uncompressed image */
              int *status)         /* IO - error status               */

/* 
  THIS IS AN OBSOLETE ROUTINE.  USE fits_img_decompress instead!!!
  
  This routine decompresses the whole image and writes it to the output file.
*/

{
    double *data;
    int ii, datatype = 0, byte_per_pix = 0;
    int nullcheck, anynul;
    LONGLONG fpixel[MAX_COMPRESS_DIM], lpixel[MAX_COMPRESS_DIM];
    long inc[MAX_COMPRESS_DIM];
    long imgsize, memsize;
    float *nulladdr, fnulval;
    double dnulval;

    if (*status > 0)
        return(*status);

    if (!fits_is_compressed_image(infptr, status) )
    {
        ffpmsg("CHDU is not a compressed image (fits_decompress_img)");
        return(*status = DATA_DECOMPRESSION_ERR);
    }

    /* create an empty output image with the correct dimensions */
    if (ffcrim(outfptr, (infptr->Fptr)->zbitpix, (infptr->Fptr)->zndim, 
       (infptr->Fptr)->znaxis, status) > 0)
    {
        ffpmsg("error creating output decompressed image HDU");
    	return (*status);
    }
    /* Copy the table header to the image header. */
    if (imcomp_copy_imheader(infptr, outfptr, status) > 0)
    {
        ffpmsg("error copying header of compressed image");
    	return (*status);
    }

    /* force a rescan of the output header keywords, then reset the scaling */
    /* in case the BSCALE and BZERO keywords are present, so that the       */
    /* decompressed values won't be scaled when written to the output image */
    ffrdef(outfptr, status);
    ffpscl(outfptr, 1.0, 0.0, status);
    ffpscl(infptr, 1.0, 0.0, status);

    /* initialize; no null checking is needed for integer images */
    nullcheck = 0;
    nulladdr =  &fnulval;

    /* determine datatype for image */
    if ((infptr->Fptr)->zbitpix == BYTE_IMG)
    {
        datatype = TBYTE;
        byte_per_pix = 1;
    }
    else if ((infptr->Fptr)->zbitpix == SHORT_IMG)
    {
        datatype = TSHORT;
        byte_per_pix = sizeof(short);
    }
    else if ((infptr->Fptr)->zbitpix == LONG_IMG)
    {
        datatype = TINT;
        byte_per_pix = sizeof(int);
    }
    else if ((infptr->Fptr)->zbitpix == FLOAT_IMG)
    {
        /* In the case of float images we must check for NaNs  */
        nullcheck = 1;
        fnulval = FLOATNULLVALUE;
        nulladdr =  &fnulval;
        datatype = TFLOAT;
        byte_per_pix = sizeof(float);
    }
    else if ((infptr->Fptr)->zbitpix == DOUBLE_IMG)
    {
        /* In the case of double images we must check for NaNs  */
        nullcheck = 1;
        dnulval = DOUBLENULLVALUE;
        nulladdr = (float *) &dnulval;
        datatype = TDOUBLE;
        byte_per_pix = sizeof(double);
    }

    /* calculate size of the image (in pixels) */
    imgsize = 1;
    for (ii = 0; ii < (infptr->Fptr)->zndim; ii++)
    {
        imgsize *= (infptr->Fptr)->znaxis[ii];
        fpixel[ii] = 1;              /* Set first and last pixel to */
        lpixel[ii] = (infptr->Fptr)->znaxis[ii]; /* include the entire image. */
        inc[ii] = 1;
    }
    /* Calc equivalent number of double pixels same size as whole the image. */
    /* We use double datatype to force the memory to be aligned properly */
    memsize = ((imgsize * byte_per_pix) - 1) / sizeof(double) + 1;

    /* allocate memory for the image */
    data = (double*) calloc (memsize, sizeof(double));
    if (!data)
    { 
        ffpmsg("Couldn't allocate memory for the uncompressed image");
        return(*status = MEMORY_ALLOCATION);
    }

    /* uncompress the entire image into memory */
    /* This routine should be enhanced sometime to only need enough */
    /* memory to uncompress one tile at a time.  */
    fits_read_compressed_img(infptr, datatype, fpixel, lpixel, inc,  
            nullcheck, nulladdr, data, NULL, &anynul, status);

    /* write the image to the output file */
    if (anynul)
        fits_write_imgnull(outfptr, datatype, 1, imgsize, data, nulladdr, 
                          status);
    else
        fits_write_img(outfptr, datatype, 1, imgsize, data, status);

    free(data);
    return (*status);
}
/*--------------------------------------------------------------------------*/
int fits_img_decompress_header(fitsfile *infptr, /* image (bintable) to uncompress */
              fitsfile *outfptr,   /* empty HDU for output uncompressed image */
              int *status)         /* IO - error status               */

/* 
  This routine reads the header of the input tile compressed image and 
  converts it to that of a standard uncompress FITS image.
*/

{
    int writeprime = 0;
    int hdupos, inhdupos, numkeys;
    int nullprime = 0, copyprime = 0, norec = 0, tstatus;
    char card[FLEN_CARD];
    int ii, datatype = 0, naxis, bitpix;
    long naxes[MAX_COMPRESS_DIM];

    if (*status > 0)
        return(*status);
    else if (*status == -1) {
        *status = 0;
	writeprime = 1;
    }

    if (!fits_is_compressed_image(infptr, status) )
    {
        ffpmsg("CHDU is not a compressed image (fits_img_decompress)");
        return(*status = DATA_DECOMPRESSION_ERR);
    }

    /* get information about the state of the output file; does it already */
    /* contain any keywords and HDUs?  */
    fits_get_hdu_num(infptr, &inhdupos);  /* Get the current output HDU position */
    fits_get_hdu_num(outfptr, &hdupos);  /* Get the current output HDU position */
    fits_get_hdrspace(outfptr, &numkeys, 0, status);

    /* Was the input compressed HDU originally the primary array image? */
    tstatus = 0;
    if (!fits_read_card(infptr, "ZSIMPLE", card, &tstatus)) { 
      /* yes, input HDU was a primary array (not an IMAGE extension) */
      /* Now determine if we can uncompress it into the primary array of */
      /* the output file.  This is only possible if the output file */
      /* currently only contains a null primary array, with no addition */
      /* header keywords and with no following extension in the FITS file. */
      
      if (hdupos == 1) {  /* are we positioned at the primary array? */
            if (numkeys == 0) { /* primary HDU is completely empty */
	        nullprime = 1;
            } else {
                fits_get_img_param(outfptr, MAX_COMPRESS_DIM, &bitpix, &naxis, naxes, status);
	
	        if (naxis == 0) { /* is this a null image? */
                   nullprime = 1;

		   if (inhdupos == 2)  /* must be at the first extension */
		      copyprime = 1;
		}
           }
      }
    } 

    if (nullprime) {  
       /* We will delete the existing keywords in the null primary array
          and uncompress the input image into the primary array of the output.
	  Some of these keywords may be added back to the uncompressed image
	  header later.
       */

       for (ii = numkeys; ii > 0; ii--)
          fits_delete_record(outfptr, ii, status);

    } else  {

       /* if the ZTENSION keyword doesn't exist, then we have to 
          write the required keywords manually */
       tstatus = 0;
       if (fits_read_card(infptr, "ZTENSION", card, &tstatus)) {

          /* create an empty output image with the correct dimensions */
          if (ffcrim(outfptr, (infptr->Fptr)->zbitpix, (infptr->Fptr)->zndim, 
             (infptr->Fptr)->znaxis, status) > 0)
          {
             ffpmsg("error creating output decompressed image HDU");
    	     return (*status);
          }

	  norec = 1;  /* the required keywords have already been written */

       } else {  /* the input compressed image does have ZTENSION keyword */
       
          if (writeprime) {  /* convert the image extension to a primary array */
	      /* have to write the required keywords manually */

              /* create an empty output image with the correct dimensions */
              if (ffcrim(outfptr, (infptr->Fptr)->zbitpix, (infptr->Fptr)->zndim, 
                 (infptr->Fptr)->znaxis, status) > 0)
              {
                 ffpmsg("error creating output decompressed image HDU");
    	         return (*status);
              }

	      norec = 1;  /* the required keywords have already been written */

          } else {  /* write the input compressed image to an image extension */

              if (numkeys == 0) {  /* the output file is currently completely empty */
	  
	         /* In this case, the input is a compressed IMAGE extension. */
	         /* Since the uncompressed output file is currently completely empty, */
	         /* we need to write a null primary array before uncompressing the */
                 /* image extension */
	     
                 ffcrim(outfptr, 8, 0, naxes, status); /* naxes is not used */
	     
	         /* now create the empty extension to uncompress into */
                 if (fits_create_hdu(outfptr, status) > 0)
                 {
                      ffpmsg("error creating output decompressed image HDU");
    	              return (*status);
                 }
	  
	      } else {
                  /* just create a new empty extension, then copy all the required */
	          /* keywords into it.  */
                 fits_create_hdu(outfptr, status);
	      }
           }
       }

    }

    if (*status > 0)  {
        ffpmsg("error creating output decompressed image HDU");
    	return (*status);
    }

    /* Copy the table header to the image header. */

    if (imcomp_copy_comp2img(infptr, outfptr, norec, status) > 0)
    {
        ffpmsg("error copying header keywords from compressed image");
    }

    if (copyprime) {  
	/* append any unexpected keywords from the primary array.
	   This includes any keywords except SIMPLE, BITPIX, NAXIS,
	   EXTEND, COMMENT, HISTORY, CHECKSUM, and DATASUM.
	*/

        fits_movabs_hdu(infptr, 1, NULL, status);  /* move to primary array */
	
        /* do this so that any new keywords get written before any blank
	   keywords that may have been appended by imcomp_copy_comp2img  */
        fits_set_hdustruc(outfptr, status);

        if (imcomp_copy_prime2img(infptr, outfptr, status) > 0)
        {
            ffpmsg("error copying primary keywords from compressed file");
        }

        fits_movabs_hdu(infptr, 2, NULL, status); /* move back to where we were */
    }

    return (*status);
}
/*---------------------------------------------------------------------------*/
int fits_read_compressed_img(fitsfile *fptr,   /* I - FITS file pointer      */
            int  datatype,  /* I - datatype of the array to be returned      */
            LONGLONG  *infpixel, /* I - 'bottom left corner' of the subsection    */
            LONGLONG  *inlpixel, /* I - 'top right corner' of the subsection      */
            long  *ininc,    /* I - increment to be applied in each dimension */
            int  nullcheck,  /* I - 0 for no null checking                   */
                              /*     1: set undefined pixels = nullval       */
                              /*     2: set nullarray=1 for undefined pixels */
            void *nullval,    /* I - value for undefined pixels              */
            void *array,      /* O - array of values that are returned       */
            char *nullarray,  /* O - array of flags = 1 if nullcheck = 2     */
            int  *anynul,     /* O - set to 1 if any values are null; else 0 */
            int  *status)     /* IO - error status                           */
/*
   Read a section of a compressed image;  Note: lpixel may be larger than the 
   size of the uncompressed image.  Only the pixels within the image will be
   returned.
*/
{
    int naxis[MAX_COMPRESS_DIM], tiledim[MAX_COMPRESS_DIM];
    long tilesize[MAX_COMPRESS_DIM], thistilesize[MAX_COMPRESS_DIM];
    long ftile[MAX_COMPRESS_DIM], ltile[MAX_COMPRESS_DIM];
    long tfpixel[MAX_COMPRESS_DIM], tlpixel[MAX_COMPRESS_DIM];
    long rowdim[MAX_COMPRESS_DIM], offset[MAX_COMPRESS_DIM],ntemp;
    long fpixel[MAX_COMPRESS_DIM], lpixel[MAX_COMPRESS_DIM];
    long inc[MAX_COMPRESS_DIM];
    int ii, i5, i4, i3, i2, i1, i0, ndim, irow, pixlen, tilenul;
    void *buffer;
    char *bnullarray = 0;
    double testnullval = 0.;

    if (*status > 0) 
        return(*status);

    if (!fits_is_compressed_image(fptr, status) )
    {
        ffpmsg("CHDU is not a compressed image (fits_read_compressed_img)");
        return(*status = DATA_DECOMPRESSION_ERR);
    }

    /* get temporary space for uncompressing one image tile */
    if (datatype == TSHORT)
    {
       buffer =  malloc ((fptr->Fptr)->maxtilelen * sizeof (short)); 
       pixlen = sizeof(short);
       if (nullval)
           testnullval = *(short *) nullval;
    }
    else if (datatype == TINT)
    {
       buffer =  malloc ((fptr->Fptr)->maxtilelen * sizeof (int));
       pixlen = sizeof(int);
       if (nullval)
           testnullval = *(int *) nullval;
    }
    else if (datatype == TLONG)
    {
       buffer =  malloc ((fptr->Fptr)->maxtilelen * sizeof (long));
       pixlen = sizeof(long);
       if (nullval)
           testnullval = *(long *) nullval;
    }
    else if (datatype == TFLOAT)
    {
       buffer =  malloc ((fptr->Fptr)->maxtilelen * sizeof (float));
       pixlen = sizeof(float);
       if (nullval)
           testnullval = *(float *) nullval;
    }
    else if (datatype == TDOUBLE)
    {
       buffer =  malloc ((fptr->Fptr)->maxtilelen * sizeof (double));
       pixlen = sizeof(double);
       if (nullval)
           testnullval = *(double *) nullval;
    }
    else if (datatype == TUSHORT)
    {
       buffer =  malloc ((fptr->Fptr)->maxtilelen * sizeof (unsigned short));
       pixlen = sizeof(short);
       if (nullval)
           testnullval = *(unsigned short *) nullval;
    }
    else if (datatype == TUINT)
    {
       buffer =  malloc ((fptr->Fptr)->maxtilelen * sizeof (unsigned int));
       pixlen = sizeof(int);
       if (nullval)
           testnullval = *(unsigned int *) nullval;
    }
    else if (datatype == TULONG)
    {
       buffer =  malloc ((fptr->Fptr)->maxtilelen * sizeof (unsigned long));
       pixlen = sizeof(long);
       if (nullval)
           testnullval = *(unsigned long *) nullval;
    }
    else if (datatype == TBYTE || datatype == TSBYTE)
    {
       buffer =  malloc ((fptr->Fptr)->maxtilelen * sizeof (char));
       pixlen = 1;
       if (nullval)
           testnullval = *(unsigned char *) nullval;
    }
    else
    {
        ffpmsg("unsupported datatype for uncompressing image");
        return(*status = BAD_DATATYPE);
    }

    /* If nullcheck ==1 and nullval == 0, then this means that the */
    /* calling routine does not want to check for null pixels in the array */
    if (nullcheck == 1 && testnullval == 0.)
        nullcheck = 0;

    if (buffer == NULL)
    {
	    ffpmsg("Out of memory (fits_read_compress_img)");
	    return (*status = MEMORY_ALLOCATION);
    }
	
    /* allocate memory for a null flag array, if needed */
    if (nullcheck == 2)
    {
        bnullarray = calloc ((fptr->Fptr)->maxtilelen, sizeof (char));

        if (bnullarray == NULL)
        {
	    ffpmsg("Out of memory (fits_read_compress_img)");
            free(buffer);
	    return (*status = MEMORY_ALLOCATION);
        }
    }

    /* initialize all the arrays */
    for (ii = 0; ii < MAX_COMPRESS_DIM; ii++)
    {
        naxis[ii] = 1;
        tiledim[ii] = 1;
        tilesize[ii] = 1;
        ftile[ii] = 1;
        ltile[ii] = 1;
        rowdim[ii] = 1;
    }

    ndim = (fptr->Fptr)->zndim;
    ntemp = 1;
    for (ii = 0; ii < ndim; ii++)
    {
        /* support for mirror-reversed image sections */
        if (infpixel[ii] <= inlpixel[ii])
        {
           fpixel[ii] = (long) infpixel[ii];
           lpixel[ii] = (long) inlpixel[ii];
           inc[ii]    = ininc[ii];
        }
        else
        {
           fpixel[ii] = (long) inlpixel[ii];
           lpixel[ii] = (long) infpixel[ii];
           inc[ii]    = -ininc[ii];
        }

        /* calc number of tiles in each dimension, and tile containing */
        /* the first and last pixel we want to read in each dimension  */
        naxis[ii] = (fptr->Fptr)->znaxis[ii];
        if (fpixel[ii] < 1)
        {
            if (nullcheck == 2)
            {
                free(bnullarray);
            }
            free(buffer);
            return(*status = BAD_PIX_NUM);
        }

        tilesize[ii] = (fptr->Fptr)->tilesize[ii];
        tiledim[ii] = (naxis[ii] - 1) / tilesize[ii] + 1;
        ftile[ii]   = (fpixel[ii] - 1)   / tilesize[ii] + 1;
        ltile[ii]   = minvalue((lpixel[ii] - 1) / tilesize[ii] + 1, 
                                tiledim[ii]);
        rowdim[ii]  = ntemp;  /* total tiles in each dimension */
        ntemp *= tiledim[ii];
    }

    if (anynul)
       *anynul = 0;  /* initialize */

    /* support up to 6 dimensions for now */
    /* tfpixel and tlpixel are the first and last image pixels */
    /* along each dimension of the compression tile */
    for (i5 = ftile[5]; i5 <= ltile[5]; i5++)
    {
     tfpixel[5] = (i5 - 1) * tilesize[5] + 1;
     tlpixel[5] = minvalue(tfpixel[5] + tilesize[5] - 1, 
                            naxis[5]);
     thistilesize[5] = tlpixel[5] - tfpixel[5] + 1;
     offset[5] = (i5 - 1) * rowdim[5];
     for (i4 = ftile[4]; i4 <= ltile[4]; i4++)
     {
      tfpixel[4] = (i4 - 1) * tilesize[4] + 1;
      tlpixel[4] = minvalue(tfpixel[4] + tilesize[4] - 1, 
                            naxis[4]);
      thistilesize[4] = thistilesize[5] * (tlpixel[4] - tfpixel[4] + 1);
      offset[4] = (i4 - 1) * rowdim[4] + offset[5];
      for (i3 = ftile[3]; i3 <= ltile[3]; i3++)
      {
        tfpixel[3] = (i3 - 1) * tilesize[3] + 1;
        tlpixel[3] = minvalue(tfpixel[3] + tilesize[3] - 1, 
                              naxis[3]);
        thistilesize[3] = thistilesize[4] * (tlpixel[3] - tfpixel[3] + 1);
        offset[3] = (i3 - 1) * rowdim[3] + offset[4];
        for (i2 = ftile[2]; i2 <= ltile[2]; i2++)
        {
          tfpixel[2] = (i2 - 1) * tilesize[2] + 1;
          tlpixel[2] = minvalue(tfpixel[2] + tilesize[2] - 1, 
                                naxis[2]);
          thistilesize[2] = thistilesize[3] * (tlpixel[2] - tfpixel[2] + 1);
          offset[2] = (i2 - 1) * rowdim[2] + offset[3];
          for (i1 = ftile[1]; i1 <= ltile[1]; i1++)
          {
            tfpixel[1] = (i1 - 1) * tilesize[1] + 1;
            tlpixel[1] = minvalue(tfpixel[1] + tilesize[1] - 1, 
                                  naxis[1]);
            thistilesize[1] = thistilesize[2] * (tlpixel[1] - tfpixel[1] + 1);
            offset[1] = (i1 - 1) * rowdim[1] + offset[2];
            for (i0 = ftile[0]; i0 <= ltile[0]; i0++)
            {
              tfpixel[0] = (i0 - 1) * tilesize[0] + 1;
              tlpixel[0] = minvalue(tfpixel[0] + tilesize[0] - 1, 
                                    naxis[0]);
              thistilesize[0] = thistilesize[1] * (tlpixel[0] - tfpixel[0] + 1);
              /* calculate row of table containing this tile */
              irow = i0 + offset[1];

/*
printf("row %d, %d %d, %d %d, %d %d; %d\n",
              irow, tfpixel[0],tlpixel[0],tfpixel[1],tlpixel[1],tfpixel[2],tlpixel[2],
	      thistilesize[0]);
*/   
              /* read and uncompress this row (tile) of the table */
              /* also do type conversion and undefined pixel substitution */
              /* at this point */

              imcomp_decompress_tile(fptr, irow, thistilesize[0],
                    datatype, nullcheck, nullval, buffer, bnullarray, &tilenul,
                     status);

              if (tilenul && anynul)
                  *anynul = 1;  /* there are null pixels */
/*
printf(" pixlen=%d, ndim=%d, %d %d %d, %d %d %d, %d %d %d\n",
     pixlen, ndim, fpixel[0],lpixel[0],inc[0],fpixel[1],lpixel[1],inc[1],
     fpixel[2],lpixel[2],inc[2]);
*/
              /* copy the intersecting pixels from this tile to the output */
              imcomp_copy_overlap(buffer, pixlen, ndim, tfpixel, tlpixel, 
                     bnullarray, array, fpixel, lpixel, inc, nullcheck, 
                     nullarray, status);
            }
          }
        }
      }
     }
    }
    if (nullcheck == 2)
    {
        free(bnullarray);
    }
    free(buffer);

    return(*status);
}
/*---------------------------------------------------------------------------*/
int fits_read_write_compressed_img(fitsfile *fptr,   /* I - FITS file pointer      */
            int  datatype,  /* I - datatype of the array to be returned      */
            LONGLONG  *infpixel, /* I - 'bottom left corner' of the subsection    */
            LONGLONG  *inlpixel, /* I - 'top right corner' of the subsection      */
            long  *ininc,    /* I - increment to be applied in each dimension */
            int  nullcheck,  /* I - 0 for no null checking                   */
                              /*     1: set undefined pixels = nullval       */
            void *nullval,    /* I - value for undefined pixels              */
            int  *anynul,     /* O - set to 1 if any values are null; else 0 */
            fitsfile *outfptr,   /* I - FITS file pointer                    */
            int  *status)     /* IO - error status                           */
/*
   This is similar to fits_read_compressed_img, except that it writes
   the pixels to the output image, on a tile by tile basis instead of returning
   the array.
*/
{
    int naxis[MAX_COMPRESS_DIM], tiledim[MAX_COMPRESS_DIM];
    long tilesize[MAX_COMPRESS_DIM], thistilesize[MAX_COMPRESS_DIM];
    long ftile[MAX_COMPRESS_DIM], ltile[MAX_COMPRESS_DIM];
    long tfpixel[MAX_COMPRESS_DIM], tlpixel[MAX_COMPRESS_DIM];
    long rowdim[MAX_COMPRESS_DIM], offset[MAX_COMPRESS_DIM],ntemp;
    long fpixel[MAX_COMPRESS_DIM], lpixel[MAX_COMPRESS_DIM];
    long inc[MAX_COMPRESS_DIM];
    int ii, i5, i4, i3, i2, i1, i0, ndim, irow, pixlen, tilenul;
    void *buffer;
    char *bnullarray = 0;
    double testnullval = 0.;
    LONGLONG firstelem;

    if (*status > 0) 
        return(*status);

    if (!fits_is_compressed_image(fptr, status) )
    {
        ffpmsg("CHDU is not a compressed image (fits_read_compressed_img)");
        return(*status = DATA_DECOMPRESSION_ERR);
    }

    /* get temporary space for uncompressing one image tile */
    if (datatype == TSHORT)
    {
       buffer =  malloc ((fptr->Fptr)->maxtilelen * sizeof (short)); 
       pixlen = sizeof(short);
       if (nullval)
           testnullval = *(short *) nullval;
    }
    else if (datatype == TINT)
    {
       buffer =  malloc ((fptr->Fptr)->maxtilelen * sizeof (int));
       pixlen = sizeof(int);
       if (nullval)
           testnullval = *(int *) nullval;
    }
    else if (datatype == TLONG)
    {
       buffer =  malloc ((fptr->Fptr)->maxtilelen * sizeof (long));
       pixlen = sizeof(long);
       if (nullval)
           testnullval = *(long *) nullval;
    }
    else if (datatype == TFLOAT)
    {
       buffer =  malloc ((fptr->Fptr)->maxtilelen * sizeof (float));
       pixlen = sizeof(float);
       if (nullval)
           testnullval = *(float *) nullval;
    }
    else if (datatype == TDOUBLE)
    {
       buffer =  malloc ((fptr->Fptr)->maxtilelen * sizeof (double));
       pixlen = sizeof(double);
       if (nullval)
           testnullval = *(double *) nullval;
    }
    else if (datatype == TUSHORT)
    {
       buffer =  malloc ((fptr->Fptr)->maxtilelen * sizeof (unsigned short));
       pixlen = sizeof(short);
       if (nullval)
           testnullval = *(unsigned short *) nullval;
    }
    else if (datatype == TUINT)
    {
       buffer =  malloc ((fptr->Fptr)->maxtilelen * sizeof (unsigned int));
       pixlen = sizeof(int);
       if (nullval)
           testnullval = *(unsigned int *) nullval;
    }
    else if (datatype == TULONG)
    {
       buffer =  malloc ((fptr->Fptr)->maxtilelen * sizeof (unsigned long));
       pixlen = sizeof(long);
       if (nullval)
           testnullval = *(unsigned long *) nullval;
    }
    else if (datatype == TBYTE || datatype == TSBYTE)
    {
       buffer =  malloc ((fptr->Fptr)->maxtilelen * sizeof (char));
       pixlen = 1;
       if (nullval)
           testnullval = *(unsigned char *) nullval;
    }
    else
    {
        ffpmsg("unsupported datatype for uncompressing image");
        return(*status = BAD_DATATYPE);
    }

    /* If nullcheck ==1 and nullval == 0, then this means that the */
    /* calling routine does not want to check for null pixels in the array */
    if (nullcheck == 1 && testnullval == 0.)
        nullcheck = 0;

    if (buffer == NULL)
    {
	    ffpmsg("Out of memory (fits_read_compress_img)");
	    return (*status = MEMORY_ALLOCATION);
    }

    /* initialize all the arrays */
    for (ii = 0; ii < MAX_COMPRESS_DIM; ii++)
    {
        naxis[ii] = 1;
        tiledim[ii] = 1;
        tilesize[ii] = 1;
        ftile[ii] = 1;
        ltile[ii] = 1;
        rowdim[ii] = 1;
    }

    ndim = (fptr->Fptr)->zndim;
    ntemp = 1;
    for (ii = 0; ii < ndim; ii++)
    {
        /* support for mirror-reversed image sections */
        if (infpixel[ii] <= inlpixel[ii])
        {
           fpixel[ii] = (long) infpixel[ii];
           lpixel[ii] = (long) inlpixel[ii];
           inc[ii]    = ininc[ii];
        }
        else
        {
           fpixel[ii] = (long) inlpixel[ii];
           lpixel[ii] = (long) infpixel[ii];
           inc[ii]    = -ininc[ii];
        }

        /* calc number of tiles in each dimension, and tile containing */
        /* the first and last pixel we want to read in each dimension  */
        naxis[ii] = (fptr->Fptr)->znaxis[ii];
        if (fpixel[ii] < 1)
        {
            free(buffer);
            return(*status = BAD_PIX_NUM);
        }

        tilesize[ii] = (fptr->Fptr)->tilesize[ii];
        tiledim[ii] = (naxis[ii] - 1) / tilesize[ii] + 1;
        ftile[ii]   = (fpixel[ii] - 1)   / tilesize[ii] + 1;
        ltile[ii]   = minvalue((lpixel[ii] - 1) / tilesize[ii] + 1, 
                                tiledim[ii]);
        rowdim[ii]  = ntemp;  /* total tiles in each dimension */
        ntemp *= tiledim[ii];
    }

    if (anynul)
       *anynul = 0;  /* initialize */

    firstelem = 1;

    /* support up to 6 dimensions for now */
    /* tfpixel and tlpixel are the first and last image pixels */
    /* along each dimension of the compression tile */
    for (i5 = ftile[5]; i5 <= ltile[5]; i5++)
    {
     tfpixel[5] = (i5 - 1) * tilesize[5] + 1;
     tlpixel[5] = minvalue(tfpixel[5] + tilesize[5] - 1, 
                            naxis[5]);
     thistilesize[5] = tlpixel[5] - tfpixel[5] + 1;
     offset[5] = (i5 - 1) * rowdim[5];
     for (i4 = ftile[4]; i4 <= ltile[4]; i4++)
     {
      tfpixel[4] = (i4 - 1) * tilesize[4] + 1;
      tlpixel[4] = minvalue(tfpixel[4] + tilesize[4] - 1, 
                            naxis[4]);
      thistilesize[4] = thistilesize[5] * (tlpixel[4] - tfpixel[4] + 1);
      offset[4] = (i4 - 1) * rowdim[4] + offset[5];
      for (i3 = ftile[3]; i3 <= ltile[3]; i3++)
      {
        tfpixel[3] = (i3 - 1) * tilesize[3] + 1;
        tlpixel[3] = minvalue(tfpixel[3] + tilesize[3] - 1, 
                              naxis[3]);
        thistilesize[3] = thistilesize[4] * (tlpixel[3] - tfpixel[3] + 1);
        offset[3] = (i3 - 1) * rowdim[3] + offset[4];
        for (i2 = ftile[2]; i2 <= ltile[2]; i2++)
        {
          tfpixel[2] = (i2 - 1) * tilesize[2] + 1;
          tlpixel[2] = minvalue(tfpixel[2] + tilesize[2] - 1, 
                                naxis[2]);
          thistilesize[2] = thistilesize[3] * (tlpixel[2] - tfpixel[2] + 1);
          offset[2] = (i2 - 1) * rowdim[2] + offset[3];
          for (i1 = ftile[1]; i1 <= ltile[1]; i1++)
          {
            tfpixel[1] = (i1 - 1) * tilesize[1] + 1;
            tlpixel[1] = minvalue(tfpixel[1] + tilesize[1] - 1, 
                                  naxis[1]);
            thistilesize[1] = thistilesize[2] * (tlpixel[1] - tfpixel[1] + 1);
            offset[1] = (i1 - 1) * rowdim[1] + offset[2];
            for (i0 = ftile[0]; i0 <= ltile[0]; i0++)
            {
              tfpixel[0] = (i0 - 1) * tilesize[0] + 1;
              tlpixel[0] = minvalue(tfpixel[0] + tilesize[0] - 1, 
                                    naxis[0]);
              thistilesize[0] = thistilesize[1] * (tlpixel[0] - tfpixel[0] + 1);
              /* calculate row of table containing this tile */
              irow = i0 + offset[1];
 
              /* read and uncompress this row (tile) of the table */
              /* also do type conversion and undefined pixel substitution */
              /* at this point */

              imcomp_decompress_tile(fptr, irow, thistilesize[0],
                    datatype, nullcheck, nullval, buffer, bnullarray, &tilenul,
                     status);

               /* write the image to the output file */

              if (tilenul && anynul) {     
                   /* this assumes that the tiled pixels are in the same order
		      as in the uncompressed FITS image.  This is not necessarily
		      the case, but it almost alway is in practice.  
		      Note that null checking is not performed for integer images,
		      so this could only be a problem for tile compressed floating
		      point images that use an unconventional tiling pattern.
		   */
                   fits_write_imgnull(outfptr, datatype, firstelem, thistilesize[0],
		      buffer, nullval, status);
              } else {
                  fits_write_subset(outfptr, datatype, tfpixel, tlpixel, 
		      buffer, status);
              }

              firstelem += thistilesize[0];

            }
          }
        }
      }
     }
    }

    free(buffer);

    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_read_compressed_pixels(fitsfile *fptr, /* I - FITS file pointer    */
            int  datatype,  /* I - datatype of the array to be returned     */
            LONGLONG   fpixel, /* I - 'first pixel to read          */
            LONGLONG   npixel,  /* I - number of pixels to read      */
            int  nullcheck,  /* I - 0 for no null checking                   */
                              /*     1: set undefined pixels = nullval       */
                              /*     2: set nullarray=1 for undefined pixels */
            void *nullval,    /* I - value for undefined pixels              */
            void *array,      /* O - array of values that are returned       */
            char *nullarray,  /* O - array of flags = 1 if nullcheck = 2     */
            int  *anynul,     /* O - set to 1 if any values are null; else 0 */
            int  *status)     /* IO - error status                           */
/*
   Read a consecutive set of pixels from a compressed image.  This routine
   interpretes the n-dimensional image as a long one-dimensional array. 
   This is actually a rather inconvenient way to read compressed images in
   general, and could be rather inefficient if the requested pixels to be
   read are located in many different image compression tiles.    

   The general strategy used here is to read the requested pixels in blocks
   that correspond to rectangular image sections.  
*/
{
    int naxis, ii, bytesperpixel, planenul;
    long naxes[MAX_COMPRESS_DIM], nread;
    long nplane, inc[MAX_COMPRESS_DIM];
    LONGLONG tfirst, tlast, last0, last1, dimsize[MAX_COMPRESS_DIM];
    LONGLONG firstcoord[MAX_COMPRESS_DIM], lastcoord[MAX_COMPRESS_DIM];
    char *arrayptr, *nullarrayptr;

    if (*status > 0)
        return(*status);

    arrayptr = (char *) array;
    nullarrayptr = nullarray;

    /* get size of array pixels, in bytes */
    bytesperpixel = ffpxsz(datatype);

    for (ii = 0; ii < MAX_COMPRESS_DIM; ii++)
    {
        naxes[ii] = 1;
        firstcoord[ii] = 0;
        lastcoord[ii] = 0;
        inc[ii] = 1;
    }

    /*  determine the dimensions of the image to be read */
    ffgidm(fptr, &naxis, status);
    ffgisz(fptr, MAX_COMPRESS_DIM, naxes, status);

    /* calc the cumulative number of pixels in each successive dimension */
    dimsize[0] = 1;
    for (ii = 1; ii < MAX_COMPRESS_DIM; ii++)
         dimsize[ii] = dimsize[ii - 1] * naxes[ii - 1];

    /*  determine the coordinate of the first and last pixel in the image */
    /*  Use zero based indexes here */
    tfirst = fpixel - 1;
    tlast = tfirst + npixel - 1;
    for (ii = naxis - 1; ii >= 0; ii--)
    {
        firstcoord[ii] = tfirst / dimsize[ii];
        lastcoord[ii] =  tlast / dimsize[ii];
        tfirst = tfirst - firstcoord[ii] * dimsize[ii];
        tlast = tlast - lastcoord[ii] * dimsize[ii];
    }

    /* to simplify things, treat 1-D, 2-D, and 3-D images as separate cases */

    if (naxis == 1)
    {
        /* Simple: just read the requested range of pixels */

        firstcoord[0] = firstcoord[0] + 1;
        lastcoord[0] = lastcoord[0] + 1;
        fits_read_compressed_img(fptr, datatype, firstcoord, lastcoord, inc,
            nullcheck, nullval, array, nullarray, anynul, status);
        return(*status);
    }
    else if (naxis == 2)
    {
        nplane = 0;  /* read 1st (and only) plane of the image */

        fits_read_compressed_img_plane(fptr, datatype, bytesperpixel,
          nplane, firstcoord, lastcoord, inc, naxes, nullcheck, nullval,
          array, nullarray, anynul, &nread, status);
    }
    else if (naxis == 3)
    {
        /* test for special case: reading an integral number of planes */
        if (firstcoord[0] == 0 && firstcoord[1] == 0 &&
            lastcoord[0] == naxes[0] - 1 && lastcoord[1] == naxes[1] - 1)
        {
            for (ii = 0; ii < MAX_COMPRESS_DIM; ii++)
            {
                /* convert from zero base to 1 base */
                (firstcoord[ii])++;
                (lastcoord[ii])++;
            }

            /* we can read the contiguous block of pixels in one go */
            fits_read_compressed_img(fptr, datatype, firstcoord, lastcoord, inc,
                nullcheck, nullval, array, nullarray, anynul, status);

            return(*status);
        }

        if (anynul)
            *anynul = 0;  /* initialize */

        /* save last coordinate in temporary variables */
        last0 = lastcoord[0];
        last1 = lastcoord[1];

        if (firstcoord[2] < lastcoord[2])
        {
            /* we will read up to the last pixel in all but the last plane */
            lastcoord[0] = naxes[0] - 1;
            lastcoord[1] = naxes[1] - 1;
        }

        /* read one plane of the cube at a time, for simplicity */
        for (nplane = (long) firstcoord[2]; nplane <= lastcoord[2]; nplane++)
        {
            if (nplane == lastcoord[2])
            {
                lastcoord[0] = last0;
                lastcoord[1] = last1;
            }

            fits_read_compressed_img_plane(fptr, datatype, bytesperpixel,
              nplane, firstcoord, lastcoord, inc, naxes, nullcheck, nullval,
              arrayptr, nullarrayptr, &planenul, &nread, status);

            if (planenul && anynul)
               *anynul = 1;  /* there are null pixels */

            /* for all subsequent planes, we start with the first pixel */
            firstcoord[0] = 0;
            firstcoord[1] = 0;

            /* increment pointers to next elements to be read */
            arrayptr = arrayptr + nread * bytesperpixel;
            if (nullarrayptr && (nullcheck == 2) )
                nullarrayptr = nullarrayptr + nread;
        }
    }
    else
    {
        ffpmsg("only 1D, 2D, or 3D images are currently supported");
        return(*status = DATA_DECOMPRESSION_ERR);
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_read_compressed_img_plane(fitsfile *fptr, /* I - FITS file   */
            int  datatype,  /* I - datatype of the array to be returned      */
            int  bytesperpixel, /* I - number of bytes per pixel in array */
            long   nplane,  /* I - which plane of the cube to read      */
            LONGLONG *firstcoord,  /* coordinate of first pixel to read */
            LONGLONG *lastcoord,   /* coordinate of last pixel to read */
            long *inc,         /* increment of pixels to read */
            long *naxes,      /* size of each image dimension */
            int  nullcheck,  /* I - 0 for no null checking                   */
                              /*     1: set undefined pixels = nullval       */
                              /*     2: set nullarray=1 for undefined pixels */
            void *nullval,    /* I - value for undefined pixels              */
            void *array,      /* O - array of values that are returned       */
            char *nullarray,  /* O - array of flags = 1 if nullcheck = 2     */
            int  *anynul,     /* O - set to 1 if any values are null; else 0 */
            long *nread,      /* O - total number of pixels read and returned*/
            int  *status)     /* IO - error status                           */

   /*
           in general we have to read the first partial row of the image,
           followed by the middle complete rows, followed by the last
           partial row of the image.  If the first or last rows are complete,
           then read them at the same time as all the middle rows.
    */
{
     /* bottom left coord. and top right coord. */
    LONGLONG blc[MAX_COMPRESS_DIM], trc[MAX_COMPRESS_DIM]; 
    char *arrayptr, *nullarrayptr;
    int tnull;

    if (anynul)
        *anynul = 0;

    *nread = 0;

    arrayptr = (char *) array;
    nullarrayptr = nullarray;

    blc[2] = nplane + 1;
    trc[2] = nplane + 1;

    if (firstcoord[0] != 0)
    { 
            /* have to read a partial first row */
            blc[0] = firstcoord[0] + 1;
            blc[1] = firstcoord[1] + 1;
            trc[1] = blc[1];  
            if (lastcoord[1] == firstcoord[1])
               trc[0] = lastcoord[0] + 1; /* 1st and last pixels in same row */
            else
               trc[0] = naxes[0];  /* read entire rest of the row */

            fits_read_compressed_img(fptr, datatype, blc, trc, inc,
                nullcheck, nullval, arrayptr, nullarrayptr, &tnull, status);

            *nread = *nread + (long) (trc[0] - blc[0] + 1);

            if (tnull && anynul)
               *anynul = 1;  /* there are null pixels */

            if (lastcoord[1] == firstcoord[1])
            {
               return(*status);  /* finished */
            }

            /* set starting coord to beginning of next line */
            firstcoord[0] = 0;
            firstcoord[1] += 1;
            arrayptr = arrayptr + (trc[0] - blc[0] + 1) * bytesperpixel;
            if (nullarrayptr && (nullcheck == 2) )
                nullarrayptr = nullarrayptr + (trc[0] - blc[0] + 1);

    }

    /* read contiguous complete rows of the image, if any */
    blc[0] = 1;
    blc[1] = firstcoord[1] + 1;
    trc[0] = naxes[0];

    if (lastcoord[0] + 1 == naxes[0])
    {
            /* can read the last complete row, too */
            trc[1] = lastcoord[1] + 1;
    }
    else
    {
            /* last row is incomplete; have to read it separately */
            trc[1] = lastcoord[1];
    }

    if (trc[1] >= blc[1])  /* must have at least one whole line to read */
    {
        fits_read_compressed_img(fptr, datatype, blc, trc, inc,
                nullcheck, nullval, arrayptr, nullarrayptr, &tnull, status);

        *nread = *nread + (long) ((trc[1] - blc[1] + 1) * naxes[0]);

        if (tnull && anynul)
           *anynul = 1;

        if (lastcoord[1] + 1 == trc[1])
               return(*status);  /* finished */

        /* increment pointers for the last partial row */
        arrayptr = arrayptr + (trc[1] - blc[1] + 1) * naxes[0] * bytesperpixel;
        if (nullarrayptr && (nullcheck == 2) )
                nullarrayptr = nullarrayptr + (trc[1] - blc[1] + 1) * naxes[0];
     }

    if (trc[1] == lastcoord[1] + 1)
        return(*status);           /* all done */

    /* set starting and ending coord to last line */

    trc[0] = lastcoord[0] + 1;
    trc[1] = lastcoord[1] + 1;
    blc[1] = trc[1];

    fits_read_compressed_img(fptr, datatype, blc, trc, inc,
                nullcheck, nullval, arrayptr, nullarrayptr, &tnull, status);

    if (tnull && anynul)
       *anynul = 1;

    *nread = *nread + (long) (trc[0] - blc[0] + 1);

    return(*status);
}
/*--------------------------------------------------------------------------*/
int imcomp_get_compressed_image_par(fitsfile *infptr, int *status)
 
/* 
    This routine reads keywords from a BINTABLE extension containing a
    compressed image.
*/
{
    char keyword[FLEN_KEYWORD];
    char value[FLEN_VALUE];
    int ii, tstatus, doffset;
    long expect_nrows, maxtilelen;

    if (*status > 0)
        return(*status);

    /* Copy relevant header keyword values to structure */
    if (ffgky (infptr, TSTRING, "ZCMPTYPE", value, NULL, status) > 0)
    {
        ffpmsg("required ZCMPTYPE compression keyword not found in");
        ffpmsg(" imcomp_get_compressed_image_par");
        return(*status);
    }

    (infptr->Fptr)->zcmptype[0] = '\0';
    strncat((infptr->Fptr)->zcmptype, value, 11);

    if (!FSTRCMP(value, "RICE_1") )
        (infptr->Fptr)->compress_type = RICE_1;
    else if (!FSTRCMP(value, "HCOMPRESS_1") )
        (infptr->Fptr)->compress_type = HCOMPRESS_1;
    else if (!FSTRCMP(value, "GZIP_1") )
        (infptr->Fptr)->compress_type = GZIP_1;
    else if (!FSTRCMP(value, "GZIP_2") )
        (infptr->Fptr)->compress_type = GZIP_2;
    else if (!FSTRCMP(value, "BZIP2_1") )
        (infptr->Fptr)->compress_type = BZIP2_1;
    else if (!FSTRCMP(value, "PLIO_1") )
        (infptr->Fptr)->compress_type = PLIO_1;
    else if (!FSTRCMP(value, "NOCOMPRESS") )
        (infptr->Fptr)->compress_type = NOCOMPRESS;
    else
    {
        ffpmsg("Unknown image compression type:");
        ffpmsg(value);
	return (*status = DATA_DECOMPRESSION_ERR);
    }

    /* get the floating point to integer quantization type, if present. */
    /* FITS files produced before 2009 will not have this keyword */
    tstatus = 0;
    if (ffgky(infptr, TSTRING, "ZQUANTIZ", value, NULL, &tstatus) > 0)
    {
        (infptr->Fptr)->quantize_dither = 0;
    } else {
        if (!FSTRCMP(value, "NONE") )
            (infptr->Fptr)->quantize_level = NO_QUANTIZE;
        else if (!FSTRCMP(value, "SUBTRACTIVE_DITHER_1") )
            (infptr->Fptr)->quantize_dither = SUBTRACTIVE_DITHER_1;
        else
            (infptr->Fptr)->quantize_dither = 0;
    }

    /* get the floating point quantization dithering offset, if present. */
    /* FITS files produced before October 2009 will not have this keyword */
    tstatus = 0;
    if (ffgky(infptr, TINT, "ZDITHER0", &doffset, NULL, &tstatus) > 0)
    {
	/* by default start with 1st element of random sequence */
        (infptr->Fptr)->dither_offset = 1;  
    } else {
        (infptr->Fptr)->dither_offset = doffset;
    }

    if (ffgky (infptr, TINT,  "ZBITPIX",  &(infptr->Fptr)->zbitpix,  
               NULL, status) > 0)
    {
        ffpmsg("required ZBITPIX compression keyword not found");
        return(*status);
    }

    if (ffgky (infptr,TINT, "ZNAXIS", &(infptr->Fptr)->zndim, NULL, status) > 0)
    {
        ffpmsg("required ZNAXIS compression keyword not found");
        return(*status);
    }

    if ((infptr->Fptr)->zndim < 1)
    {
        ffpmsg("Compressed image has no data (ZNAXIS < 1)");
	return (*status = BAD_NAXIS);
    }

    if ((infptr->Fptr)->zndim > MAX_COMPRESS_DIM)
    {
        ffpmsg("Compressed image has too many dimensions");
        return(*status = BAD_NAXIS);
    }

    expect_nrows = 1;
    maxtilelen = 1;
    for (ii = 0;  ii < (infptr->Fptr)->zndim;  ii++)
    {
        /* get image size */
        sprintf (keyword, "ZNAXIS%d", ii+1);
	ffgky (infptr, TLONG,keyword, &(infptr->Fptr)->znaxis[ii],NULL,status);

        if (*status > 0)
        {
            ffpmsg("required ZNAXISn compression keyword not found");
            return(*status);
        }

        /* get compression tile size */
	sprintf (keyword, "ZTILE%d", ii+1);

        /* set default tile size in case keywords are not present */
        if (ii == 0)
            (infptr->Fptr)->tilesize[0] = (infptr->Fptr)->znaxis[0];
        else
            (infptr->Fptr)->tilesize[ii] = 1;

        tstatus = 0;
	ffgky (infptr, TLONG, keyword, &(infptr->Fptr)->tilesize[ii], NULL, 
               &tstatus);

        expect_nrows *= (((infptr->Fptr)->znaxis[ii] - 1) / 
                  (infptr->Fptr)->tilesize[ii]+ 1);
        maxtilelen *= (infptr->Fptr)->tilesize[ii];
    }

    /* check number of rows */
    if (expect_nrows != (infptr->Fptr)->numrows)
    {
        ffpmsg(
        "number of table rows != the number of tiles in compressed image");
        return (*status = DATA_DECOMPRESSION_ERR);
    }

    /* read any algorithm specific parameters */
    if ((infptr->Fptr)->compress_type == RICE_1 )
    {
        if (ffgky(infptr, TINT,"ZVAL1", &(infptr->Fptr)->rice_blocksize,
                  NULL, status) > 0)
        {
            ffpmsg("required ZVAL1 compression keyword not found");
            return(*status);
        }

        tstatus = 0;
        if (ffgky(infptr, TINT,"ZVAL2", &(infptr->Fptr)->rice_bytepix,
                  NULL, &tstatus) > 0)
        {
            (infptr->Fptr)->rice_bytepix = 4;  /* default value */
        }

        if ((infptr->Fptr)->rice_blocksize < 16 &&
	    (infptr->Fptr)->rice_bytepix > 8) {
	     /* values are reversed */
	     tstatus = (infptr->Fptr)->rice_bytepix;
	     (infptr->Fptr)->rice_bytepix = (infptr->Fptr)->rice_blocksize;
	     (infptr->Fptr)->rice_blocksize = tstatus;
        }
    } else if ((infptr->Fptr)->compress_type == HCOMPRESS_1 ) {

        if (ffgky(infptr, TFLOAT,"ZVAL1", &(infptr->Fptr)->hcomp_scale,
                  NULL, status) > 0)
        {
            ffpmsg("required ZVAL1 compression keyword not found");
            return(*status);
        }

        tstatus = 0;
        ffgky(infptr, TINT,"ZVAL2", &(infptr->Fptr)->hcomp_smooth,
                  NULL, &tstatus);
    }    

    /* store number of pixels in each compression tile, */
    /* and max size of the compressed tile buffer */
    (infptr->Fptr)->maxtilelen = maxtilelen;

    (infptr->Fptr)->maxelem = 
           imcomp_calc_max_elem ((infptr->Fptr)->compress_type, maxtilelen, 
               (infptr->Fptr)->zbitpix, (infptr->Fptr)->rice_blocksize);

    /* Get Column numbers. */
    if (ffgcno(infptr, CASEINSEN, "COMPRESSED_DATA",
         &(infptr->Fptr)->cn_compressed, status) > 0)
    {
        ffpmsg("couldn't find COMPRESSED_DATA column (fits_get_compressed_img_par)");
        return(*status = DATA_DECOMPRESSION_ERR);
    }

    ffpmrk(); /* put mark on message stack; erase any messages after this */

    tstatus = 0;
    ffgcno(infptr,CASEINSEN, "UNCOMPRESSED_DATA",
          &(infptr->Fptr)->cn_uncompressed, &tstatus);

    tstatus = 0;
    ffgcno(infptr,CASEINSEN, "GZIP_COMPRESSED_DATA",
          &(infptr->Fptr)->cn_gzip_data, &tstatus);

    tstatus = 0;
    if (ffgcno(infptr, CASEINSEN, "ZSCALE", &(infptr->Fptr)->cn_zscale,
              &tstatus) > 0)
    {
        /* CMPSCALE column doesn't exist; see if there is a keyword */
        tstatus = 0;
        if (ffgky(infptr, TDOUBLE, "ZSCALE", &(infptr->Fptr)->zscale, NULL, 
                 &tstatus) <= 0)
            (infptr->Fptr)->cn_zscale = -1;  /* flag for a constant ZSCALE */
    }

    tstatus = 0;
    if (ffgcno(infptr, CASEINSEN, "ZZERO", &(infptr->Fptr)->cn_zzero,
               &tstatus) > 0)
    {
        /* CMPZERO column doesn't exist; see if there is a keyword */
        tstatus = 0;
        if (ffgky(infptr, TDOUBLE, "ZZERO", &(infptr->Fptr)->zzero, NULL, 
                  &tstatus) <= 0)
            (infptr->Fptr)->cn_zzero = -1;  /* flag for a constant ZZERO */
    }

    tstatus = 0;
    if (ffgcno(infptr, CASEINSEN, "ZBLANK", &(infptr->Fptr)->cn_zblank,
               &tstatus) > 0)
    {
        /* ZBLANK column doesn't exist; see if there is a keyword */
        tstatus = 0;
        if (ffgky(infptr, TINT, "ZBLANK", &(infptr->Fptr)->zblank, NULL,
                  &tstatus) <= 0)  {
            (infptr->Fptr)->cn_zblank = -1;  /* flag for a constant ZBLANK */

        } else {
           /* ZBLANK keyword doesn't exist; see if there is a BLANK keyword */
           tstatus = 0;
           if (ffgky(infptr, TINT, "BLANK", &(infptr->Fptr)->zblank, NULL,
                  &tstatus) <= 0)  
              (infptr->Fptr)->cn_zblank = -1;  /* flag for a constant ZBLANK */
        }
    }

    /* read the conventional BSCALE and BZERO scaling keywords, if present */
    tstatus = 0;
    if (ffgky (infptr, TDOUBLE, "BSCALE", &(infptr->Fptr)->cn_bscale, 
        NULL, &tstatus) > 0)
    {
        (infptr->Fptr)->cn_bscale = 1.0;
    }

    tstatus = 0;
    if (ffgky (infptr, TDOUBLE, "BZERO", &(infptr->Fptr)->cn_bzero, 
        NULL, &tstatus) > 0)
    {
        (infptr->Fptr)->cn_bzero = 0.0;
        (infptr->Fptr)->cn_actual_bzero = 0.0;
    } else {
        (infptr->Fptr)->cn_actual_bzero = (infptr->Fptr)->cn_bzero;
    }

    ffcmrk();  /* clear any spurious error messages, back to the mark */
    return (*status);
}
/*--------------------------------------------------------------------------*/
int imcomp_copy_imheader(fitsfile *infptr, fitsfile *outfptr, int *status)
/*
    This routine reads the header keywords from the input image and
    copies them to the output image;  the manditory structural keywords
    and the checksum keywords are not copied. If the DATE keyword is copied,
    then it is updated with the current date and time.
*/
{
    int nkeys, ii, keyclass;
    char card[FLEN_CARD];	/* a header record */

    if (*status > 0)
        return(*status);

    ffghsp(infptr, &nkeys, NULL, status); /* get number of keywords in image */

    for (ii = 5; ii <= nkeys; ii++)  /* skip the first 4 keywords */
    {
        ffgrec(infptr, ii, card, status);

	keyclass = ffgkcl(card);  /* Get the type/class of keyword */

        /* don't copy structural keywords or checksum keywords */
        if ((keyclass <= TYP_CMPRS_KEY) || (keyclass == TYP_CKSUM_KEY))
	    continue;

        if (FSTRNCMP(card, "DATE ", 5) == 0) /* write current date */
        {
            ffpdat(outfptr, status);
        }
        else if (FSTRNCMP(card, "EXTNAME ", 8) == 0) 
        {
            /* don't copy default EXTNAME keyword from a compressed image */
            if (FSTRNCMP(card, "EXTNAME = 'COMPRESSED_IMAGE'", 28))
            {
                /* if EXTNAME keyword already exists, overwrite it */
                /* otherwise append a new EXTNAME keyword */
                ffucrd(outfptr, "EXTNAME", card, status);
            }
        }
        else
        {
            /* just copy the keyword to the output header */
	    ffprec (outfptr, card, status);
        }

        if (*status > 0)
           return (*status);
    }
    return (*status);
}
/*--------------------------------------------------------------------------*/
int imcomp_copy_img2comp(fitsfile *infptr, fitsfile *outfptr, int *status)
/*
    This routine copies the header keywords from the uncompressed input image 
    and to the compressed image (in a binary table) 
*/
{
    char card[FLEN_CARD], card2[FLEN_CARD];	/* a header record */
    int nkeys, nmore, ii, jj, tstatus, bitpix;

    /* tile compressed image keyword translation table  */
    /*                        INPUT      OUTPUT  */
    /*                       01234567   01234567 */
    char *patterns[][2] = {{"SIMPLE",  "ZSIMPLE" },  
			   {"XTENSION", "ZTENSION" },
			   {"BITPIX",  "ZBITPIX" },
			   {"NAXIS",   "ZNAXIS"  },
			   {"NAXISm",  "ZNAXISm" },
			   {"EXTEND",  "ZEXTEND" },
			   {"BLOCKED", "ZBLOCKED"},
			   {"PCOUNT",  "ZPCOUNT" },  
			   {"GCOUNT",  "ZGCOUNT" },

			   {"CHECKSUM","ZHECKSUM"},  /* save original checksums */
			   {"DATASUM", "ZDATASUM"},
			   
			   {"*",       "+"       }}; /* copy all other keywords */
    int npat;

    if (*status > 0)
        return(*status);

    /* write a default EXTNAME keyword if it doesn't exist in input file*/
    fits_read_card(infptr, "EXTNAME", card, status);
    
    if (*status) {
       *status = 0;
       strcpy(card, "EXTNAME = 'COMPRESSED_IMAGE'");
       fits_write_record(outfptr, card, status);
    }

    /* copy all the keywords from the input file to the output */
    npat = sizeof(patterns)/sizeof(patterns[0][0])/2;
    fits_translate_keywords(infptr, outfptr, 1, patterns, npat,
			    0, 0, 0, status);


    if ( (outfptr->Fptr)->request_lossy_int_compress != 0) { 

	/* request was made to compress integer images as if they had float pixels. */
	/* If input image has positive bitpix value, then reset the output ZBITPIX */
	/* value to -32. */

	fits_read_key(infptr, TINT, "BITPIX", &bitpix, NULL, status);

	if (*status <= 0 && bitpix > 0) {
	    fits_modify_key_lng(outfptr, "ZBITPIX", -32, NULL, status);

	    /* also delete the BSCALE, BZERO, and BLANK keywords */
	    tstatus = 0;
	    fits_delete_key(outfptr, "BSCALE", &tstatus);
	    tstatus = 0;
	    fits_delete_key(outfptr, "BZERO", &tstatus);
	    tstatus = 0;
	    fits_delete_key(outfptr, "BLANK", &tstatus);
	}
    }

   /*
     For compatibility with software that uses an older version of CFITSIO,
     we must make certain that the new ZQUANTIZ keyword, if it exists, must
     occur after the other peudo-required keywords (e.g., ZSIMPLE, ZBITPIX,
     etc.).  Do this by trying to delete the keyword.  If that succeeds (and
     thus the keyword did exist) then rewrite the keyword at the end of header.
     In principle this should not be necessary once all software has upgraded
     to a newer version of CFITSIO (version number greater than 3.181, newer
     than August 2009).
     
     Do the same for the new ZDITHER0 keyword.
   */

   tstatus = 0;
   if (fits_read_card(outfptr, "ZQUANTIZ", card, &tstatus) == 0)
   {
        fits_delete_key(outfptr, "ZQUANTIZ", status);

        /* rewrite the deleted keyword at the end of the header */
        fits_write_record(outfptr, card, status);

	fits_write_history(outfptr, 
	    "Image was compressed by CFITSIO using scaled integer quantization:", status);
	sprintf(card2, "  q = %f / quantized level scaling parameter", 
	    (outfptr->Fptr)->quantize_level);
	fits_write_history(outfptr, card2, status); 
	fits_write_history(outfptr, card+10, status); 
   }

   tstatus = 0;
   if (fits_read_card(outfptr, "ZDITHER0", card, &tstatus) == 0)
   {
        fits_delete_key(outfptr, "ZDITHER0", status);

        /* rewrite the deleted keyword at the end of the header */
        fits_write_record(outfptr, card, status);
   }


    ffghsp(infptr, &nkeys, &nmore, status); /* get number of keywords in image */

    nmore = nmore / 36;  /* how many completely empty header blocks are there? */
     
     /* preserve the same number of spare header blocks in the output header */
     
    for (jj = 0; jj < nmore; jj++)
       for (ii = 0; ii < 36; ii++)
          fits_write_record(outfptr, "    ", status);

    return (*status);
}
/*--------------------------------------------------------------------------*/
int imcomp_copy_comp2img(fitsfile *infptr, fitsfile *outfptr, 
                          int norec, int *status)
/*
    This routine copies the header keywords from the compressed input image 
    and to the uncompressed image (in a binary table) 
*/
{
    char card[FLEN_CARD];	/* a header record */
    char *patterns[40][2];
    char negative[] = "-";
    int ii,jj, npat, nreq, nsp, tstatus = 0;
    int nkeys, nmore;
    
    /* tile compressed image keyword translation table  */
    /*                        INPUT      OUTPUT  */
    /*                       01234567   01234567 */

    /*  only translate these if required keywords not already written */
    char *reqkeys[][2] = {  
			   {"ZSIMPLE",   "SIMPLE" },  
			   {"ZTENSION", "XTENSION"},
			   {"ZBITPIX",   "BITPIX" },
			   {"ZNAXIS",    "NAXIS"  },
			   {"ZNAXISm",   "NAXISm" },
			   {"ZEXTEND",   "EXTEND" },
			   {"ZBLOCKED",  "BLOCKED"},
			   {"ZPCOUNT",   "PCOUNT" },  
			   {"ZGCOUNT",   "GCOUNT" },
			   {"ZHECKSUM",  "CHECKSUM"},  /* restore original checksums */
			   {"ZDATASUM",  "DATASUM"}}; 

    /* other special keywords */
    char *spkeys[][2] = {
			   {"XTENSION", "-"      },
			   {"BITPIX",  "-"       },
			   {"NAXIS",   "-"       },
			   {"NAXISm",  "-"       },
			   {"PCOUNT",  "-"       },
			   {"GCOUNT",  "-"       },
			   {"TFIELDS", "-"       },
			   {"TTYPEm",  "-"       },
			   {"TFORMm",  "-"       },
			   {"ZIMAGE",  "-"       },
			   {"ZQUANTIZ", "-"      },
			   {"ZDITHER0", "-"      },
			   {"ZTILEm",  "-"       },
			   {"ZCMPTYPE", "-"      },
			   {"ZBLANK",  "-"       },
			   {"ZNAMEm",  "-"       },
			   {"ZVALm",   "-"       },

			   {"CHECKSUM","-"       },  /* delete checksums */
			   {"DATASUM", "-"       },
			   {"EXTNAME", "+"       },  /* we may change this, below */
			   {"*",       "+"      }};  


    if (*status > 0)
        return(*status);
	
    nreq = sizeof(reqkeys)/sizeof(reqkeys[0][0])/2;
    nsp = sizeof(spkeys)/sizeof(spkeys[0][0])/2;

    /* construct translation patterns */

    for (ii = 0; ii < nreq; ii++) {
        patterns[ii][0] = reqkeys[ii][0];
	
        if (norec) 
            patterns[ii][1] = negative;
        else
            patterns[ii][1] = reqkeys[ii][1];
    }
    
    for (ii = 0; ii < nsp; ii++) {
        patterns[ii+nreq][0] = spkeys[ii][0];
        patterns[ii+nreq][1] = spkeys[ii][1];
    }

    npat = nreq + nsp;
    
    /* see if the EXTNAME keyword should be copied or not */
    fits_read_card(infptr, "EXTNAME", card, &tstatus);

    if (tstatus == 0) {
      if (!strncmp(card, "EXTNAME = 'COMPRESSED_IMAGE'", 28)) 
        patterns[npat-2][1] = negative;
    }
    
    /* translate and copy the keywords from the input file to the output */
    fits_translate_keywords(infptr, outfptr, 1, patterns, npat,
			    0, 0, 0, status);

    ffghsp(infptr, &nkeys, &nmore, status); /* get number of keywords in image */

    nmore = nmore / 36;  /* how many completely empty header blocks are there? */
     
    /* preserve the same number of spare header blocks in the output header */
     
    for (jj = 0; jj < nmore; jj++)
       for (ii = 0; ii < 36; ii++)
          fits_write_record(outfptr, "    ", status);


    return (*status);
}
/*--------------------------------------------------------------------------*/
int imcomp_copy_prime2img(fitsfile *infptr, fitsfile *outfptr, int *status)
/*
    This routine copies any unexpected keywords from the primary array
    of the compressed input image into the header of the uncompressed image
    (which is the primary array of the output file). 
*/
{
    int  nsp;

    /* keywords that will not be copied */
    char *spkeys[][2] = {
			   {"SIMPLE", "-"      },
			   {"BITPIX",  "-"       },
			   {"NAXIS",   "-"       },
			   {"NAXISm",  "-"       },
			   {"PCOUNT",  "-"       },
			   {"EXTEND",  "-"       },
			   {"GCOUNT",  "-"       },
			   {"CHECKSUM","-"       }, 
			   {"DATASUM", "-"       },
			   {"EXTNAME", "-"       },
			   {"HISTORY", "-"       },
			   {"COMMENT", "-"       },
			   {"*",       "+"      }};  

    if (*status > 0)
        return(*status);
	
    nsp = sizeof(spkeys)/sizeof(spkeys[0][0])/2;

    /* translate and copy the keywords from the input file to the output */
    fits_translate_keywords(infptr, outfptr, 1, spkeys, nsp,
			    0, 0, 0, status);

    return (*status);
}
/*--------------------------------------------------------------------------*/
int imcomp_decompress_tile (fitsfile *infptr,
          int nrow,            /* I - row of table to read and uncompress */
          int tilelen,         /* I - number of pixels in the tile        */
          int datatype,        /* I - datatype to be returned in 'buffer' */
          int nullcheck,       /* I - 0 for no null checking */
          void *nulval,        /* I - value to be used for undefined pixels */
          void *buffer,        /* O - buffer for returned decompressed values */
          char *bnullarray,    /* O - buffer for returned null flags */
          int *anynul,         /* O - any null values returned?  */
          int *status)

/* This routine decompresses one tile of the image */
{
    int *idata = 0;
    int tiledatatype, pixlen;          /* uncompressed integer data */
    size_t idatalen, tilebytesize;
    int ii, tnull;        /* value in the data which represents nulls */
    unsigned char *cbuf; /* compressed data */
    unsigned char charnull = 0;
    short snull = 0;
    int blocksize;
    float fnulval=0;
    float *tempfloat = 0;
    double dnulval=0;
    double bscale, bzero, actual_bzero, dummy = 0;    /* scaling parameters */
    long nelem = 0, offset = 0, tilesize;      /* number of bytes */
    int smooth, nx, ny, scale;  /* hcompress parameters */

    if (*status > 0)
       return(*status);

    /* **************************************************************** */
    /* check if this tile was cached; if so, just copy it out */
    if (nrow == (infptr->Fptr)->tilerow && datatype == (infptr->Fptr)->tiletype ) {

         memcpy(buffer, (infptr->Fptr)->tiledata, (infptr->Fptr)->tiledatasize);
	 
	 if (nullcheck == 2)
             memcpy(bnullarray, (infptr->Fptr)->tilenullarray, tilelen);

         *anynul = (infptr->Fptr)->tileanynull;
         return(*status);
    }

    /* **************************************************************** */
    /* get length of the compressed byte stream */
    ffgdes (infptr, (infptr->Fptr)->cn_compressed, nrow, &nelem, &offset, 
            status);

    /* EOF error here indicates that this tile has not yet been written */
    if (*status == END_OF_FILE)
           return(*status = NO_COMPRESSED_TILE);
      
    /* **************************************************************** */
    if (nelem == 0)  /* special case: tile was not compressed normally */
    {
        if ((infptr->Fptr)->cn_uncompressed >= 1 ) {

	    /* This option of writing the uncompressed floating point data */
	    /* to the tile compressed file was used until about May 2011. */
	    /* This was replaced by the more efficient option of gzipping the */
	    /* floating point data before writing it to the tile-compressed file */
	    
            /* no compressed data, so simply read the uncompressed data */
            /* directly from the UNCOMPRESSED_DATA column */   
            ffgdes (infptr, (infptr->Fptr)->cn_uncompressed, nrow, &nelem,
               &offset, status);

            if (nelem == 0 && offset == 0)  /* this should never happen */
	        return (*status = NO_COMPRESSED_TILE);

            if (nullcheck <= 1) { /* set any null values in the array = nulval */
                fits_read_col(infptr, datatype, (infptr->Fptr)->cn_uncompressed,
                  nrow, 1, nelem, nulval, buffer, anynul, status);
            } else  { /* set the bnullarray = 1 for any null values in the array */
                fits_read_colnull(infptr, datatype, (infptr->Fptr)->cn_uncompressed,
                  nrow, 1, nelem, buffer, bnullarray, anynul, status);
            }
        } else if ((infptr->Fptr)->cn_gzip_data >= 1) {

            /* This is the newer option, that was introduced in May 2011 */
            /* floating point data was not quantized,  so read the losslessly */
	    /* compressed data from the GZIP_COMPRESSED_DATA column */   

            ffgdes (infptr, (infptr->Fptr)->cn_gzip_data, nrow, &nelem,
               &offset, status);

            if (nelem == 0 && offset == 0) /* this should never happen */
	        return (*status = NO_COMPRESSED_TILE);

	    /* allocate memory for the compressed tile of data */
            cbuf = (unsigned char *) malloc (nelem);  
            if (cbuf == NULL) {
	        ffpmsg("error allocating memory for gzipped tile (imcomp_decompress_tile)");
	        return (*status = MEMORY_ALLOCATION);
            }

            /* read array of compressed bytes */
            if (fits_read_col(infptr, TBYTE, (infptr->Fptr)->cn_gzip_data, nrow,
                 1, nelem, &charnull, cbuf, NULL, status) > 0) {
                ffpmsg("error reading compressed byte stream from binary table");
	        free (cbuf);
                return (*status);
            }

            /* size of the returned (uncompressed) data buffer, in bytes */
            if ((infptr->Fptr)->zbitpix == FLOAT_IMG) {
	         idatalen = tilelen * sizeof(float);
            } else if ((infptr->Fptr)->zbitpix == DOUBLE_IMG) {
	         idatalen = tilelen * sizeof(double);
            } else {
                /* this should never happen! */
                ffpmsg("incompatible data type in gzipped floating-point tile-compressed image");
                free (cbuf);
                return (*status = DATA_DECOMPRESSION_ERR);
            }

            if (datatype == TDOUBLE && (infptr->Fptr)->zbitpix == FLOAT_IMG) {  
                /*  have to allocat a temporary buffer for the uncompressed data in the */
                /*  case where a gzipped "float" tile is returned as a "double" array   */
                tempfloat = (float*) malloc (idatalen); 

                if (tempfloat == NULL) {
	            ffpmsg("Memory allocation failure for tempfloat. (imcomp_decompress_tile)");
                    free (cbuf);
	            return (*status = MEMORY_ALLOCATION);
                }

                /* uncompress the data into temp buffer */
                if (uncompress2mem_from_mem ((char *)cbuf, nelem,
                     (char **) &tempfloat, &idatalen, NULL, &tilebytesize, status)) {
                    ffpmsg("failed to gunzip the image tile");
                    free (tempfloat);
                    free (cbuf);
                    return (*status);
                }
            } else {

                /* uncompress the data directly into the output buffer in all other cases */
                if (uncompress2mem_from_mem ((char *)cbuf, nelem,
                  (char **) &buffer, &idatalen, NULL, &tilebytesize, status)) {
                    ffpmsg("failed to gunzip the image tile");
                    free (cbuf);
                    return (*status);
                }
            }

            free(cbuf);

            /* do byte swapping and null value substitution for the tile of pixels */
            if (tilebytesize == 4 * tilelen) {  /* float pixels */

#if BYTESWAPPED
                if (tempfloat)
                    ffswap4((int *) tempfloat, tilelen);
                else
                    ffswap4((int *) buffer, tilelen);
#endif
               if (datatype == TFLOAT) {
                  if (nulval) {
		    fnulval = *(float *) nulval;
  		  }

                  fffr4r4((float *) buffer, (long) tilelen, 1., 0., nullcheck,   
                        fnulval, bnullarray, anynul,
                        (float *) buffer, status);
                } else if (datatype == TDOUBLE) {
                  if (nulval) {
		    dnulval = *(double *) nulval;
		  }

                  /* note that the R*4 data are in the tempfloat array in this case */
                  fffr4r8((float *) tempfloat, (long) tilelen, 1., 0., nullcheck,   
                   dnulval, bnullarray, anynul,
                    (double *) buffer, status);            
                  free(tempfloat);

                } else {
                  ffpmsg("implicit data type conversion is not supported for gzipped image tiles");
                  return (*status = DATA_DECOMPRESSION_ERR);
                }
            } else if (tilebytesize == 8 * tilelen) { /* double pixels */

#if BYTESWAPPED
                ffswap8((double *) buffer, tilelen);
#endif
                if (datatype == TFLOAT) {
                  if (nulval) {
		    fnulval = *(float *) nulval;
  		  }

                  fffr8r4((double *) buffer, (long) tilelen, 1., 0., nullcheck,   
                        fnulval, bnullarray, anynul,
                        (float *) buffer, status);
                } else if (datatype == TDOUBLE) {
                  if (nulval) {
		    dnulval = *(double *) nulval;
		  }

                  fffr8r8((double *) buffer, (long) tilelen, 1., 0., nullcheck,   
                   dnulval, bnullarray, anynul,
                    (double *) buffer, status);            
                } else {
                  ffpmsg("implicit data type conversion is not supported in tile-compressed images");
                  return (*status = DATA_DECOMPRESSION_ERR);
                }
	    } else {
                ffpmsg("error: uncompressed tile has wrong size");
                return (*status = DATA_DECOMPRESSION_ERR);
            }

          /* end of special case of losslessly gzipping a floating-point image tile */
        } else {  /* this should never happen */
	   *status = NO_COMPRESSED_TILE;
        }

        return(*status);
    }
   
    /* **************************************************************** */
    /* deal with the normal case of a compressed tile of pixels */
    if (nullcheck == 2)  {
        for (ii = 0; ii < tilelen; ii++)  /* initialize the null flage array */
            bnullarray[ii] = 0;
    }

    if (anynul)
       *anynul = 0;

    /* get linear scaling and offset values, if they exist */
    actual_bzero = (infptr->Fptr)->cn_actual_bzero;
    if ((infptr->Fptr)->cn_zscale == 0) {
         /* set default scaling, if scaling is not defined */
         bscale = 1.;
         bzero = 0.;
    } else if ((infptr->Fptr)->cn_zscale == -1) {
        bscale = (infptr->Fptr)->zscale;
        bzero  = (infptr->Fptr)->zzero;
    } else {
        /* read the linear scale and offset values for this row */
	ffgcvd (infptr, (infptr->Fptr)->cn_zscale, nrow, 1, 1, 0.,
				&bscale, NULL, status);
	ffgcvd (infptr, (infptr->Fptr)->cn_zzero, nrow, 1, 1, 0.,
				&bzero, NULL, status);
        if (*status > 0)
        {
          ffpmsg("error reading scaling factor and offset for compressed tile");
          return (*status);
        }

        /* test if floating-point FITS image also has non-default BSCALE and  */
	/* BZERO keywords.  If so, we have to combine the 2 linear scaling factors. */
	
	if ( ((infptr->Fptr)->zbitpix == FLOAT_IMG || 
	      (infptr->Fptr)->zbitpix == DOUBLE_IMG )
	    &&  
	      ((infptr->Fptr)->cn_bscale != 1.0 ||
	       (infptr->Fptr)->cn_bzero  != 0.0 )    ) 
	    {
	       bscale = bscale * (infptr->Fptr)->cn_bscale;
	       bzero  = bzero  * (infptr->Fptr)->cn_bscale + (infptr->Fptr)->cn_bzero;
	    }
    }

    if (bscale == 1.0 && bzero == 0.0 ) {
      /* if no other scaling has been specified, try using the values
         given by the BSCALE and BZERO keywords, if any */

        bscale = (infptr->Fptr)->cn_bscale;
        bzero  = (infptr->Fptr)->cn_bzero;
    }

    /* ************************************************************* */
    /* get the value used to represent nulls in the int array */
    if ((infptr->Fptr)->cn_zblank == 0) {
        nullcheck = 0;  /* no null value; don't check for nulls */
    } else if ((infptr->Fptr)->cn_zblank == -1) {
        tnull = (infptr->Fptr)->zblank;  /* use the the ZBLANK keyword */
    } else {
        /* read the null value for this row */
	ffgcvk (infptr, (infptr->Fptr)->cn_zblank, nrow, 1, 1, 0,
				&tnull, NULL, status);
        if (*status > 0) {
            ffpmsg("error reading null value for compressed tile");
            return (*status);
        }
    }

    /* ************************************************************* */
    /* allocate memory for the uncompressed array of tile integers */
    /* The size depends on the datatype and the compression type. */
    
    if ((infptr->Fptr)->compress_type == HCOMPRESS_1 &&
          ((infptr->Fptr)->zbitpix != BYTE_IMG &&
	   (infptr->Fptr)->zbitpix != SHORT_IMG) ) {

           idatalen = tilelen * sizeof(LONGLONG);  /* 8 bytes per pixel */

    } else if ( (infptr->Fptr)->compress_type == RICE_1 &&
               (infptr->Fptr)->zbitpix == BYTE_IMG && 
	       (infptr->Fptr)->rice_bytepix == 1) {

           idatalen = tilelen * sizeof(char); /* 1 byte per pixel */
    } else if ( ( (infptr->Fptr)->compress_type == GZIP_1  ||
                  (infptr->Fptr)->compress_type == GZIP_2  ||
                  (infptr->Fptr)->compress_type == BZIP2_1 ) &&
               (infptr->Fptr)->zbitpix == BYTE_IMG ) {

           idatalen = tilelen * sizeof(char); /* 1 byte per pixel */
    } else if ( (infptr->Fptr)->compress_type == RICE_1 &&
               (infptr->Fptr)->zbitpix == SHORT_IMG && 
	       (infptr->Fptr)->rice_bytepix == 2) {

           idatalen = tilelen * sizeof(short); /* 2 bytes per pixel */
    } else if ( ( (infptr->Fptr)->compress_type == GZIP_1  ||
                  (infptr->Fptr)->compress_type == GZIP_2  ||
                  (infptr->Fptr)->compress_type == BZIP2_1 )  &&
               (infptr->Fptr)->zbitpix == SHORT_IMG ) {

           idatalen = tilelen * sizeof(short); /* 2 bytes per pixel */
    } else if ( ( (infptr->Fptr)->compress_type == GZIP_1  ||
                  (infptr->Fptr)->compress_type == GZIP_2  ||
                  (infptr->Fptr)->compress_type == BZIP2_1 ) &&
               (infptr->Fptr)->zbitpix == DOUBLE_IMG ) {

           idatalen = tilelen * sizeof(double); /* 8 bytes per pixel  */
    } else {
           idatalen = tilelen * sizeof(int);  /* all other cases have int pixels */
    }

    idata = (int*) malloc (idatalen); 
    if (idata == NULL) {
	    ffpmsg("Memory allocation failure for idata. (imcomp_decompress_tile)");
	    return (*status = MEMORY_ALLOCATION);
    }

    /* ************************************************************* */
    /* allocate memory for the compressed bytes */

    if ((infptr->Fptr)->compress_type == PLIO_1) {
        cbuf = (unsigned char *) malloc (nelem * sizeof (short));
    } else {
        cbuf = (unsigned char *) malloc (nelem);
    }
    if (cbuf == NULL) {
	ffpmsg("Out of memory for cbuf. (imcomp_decompress_tile)");
        free(idata);
	return (*status = MEMORY_ALLOCATION);
    }
    
    /* ************************************************************* */
    /* read the compressed bytes from the FITS file */

    if ((infptr->Fptr)->compress_type == PLIO_1) {
        fits_read_col(infptr, TSHORT, (infptr->Fptr)->cn_compressed, nrow,
             1, nelem, &snull, (short *) cbuf, NULL, status);
    } else {
       fits_read_col(infptr, TBYTE, (infptr->Fptr)->cn_compressed, nrow,
             1, nelem, &charnull, cbuf, NULL, status);
    }

    if (*status > 0) {
        ffpmsg("error reading compressed byte stream from binary table");
	free (cbuf);
        free(idata);
        return (*status);
    }

    /* ************************************************************* */
    /*  call the algorithm-specific code to uncompress the tile */

    if ((infptr->Fptr)->compress_type == RICE_1) {

        blocksize = (infptr->Fptr)->rice_blocksize;

        if ((infptr->Fptr)->rice_bytepix == 1 ) {
            *status = fits_rdecomp_byte (cbuf, nelem, (unsigned char *)idata,
                        tilelen, blocksize);
            tiledatatype = TBYTE;
        } else if ((infptr->Fptr)->rice_bytepix == 2 ) {
            *status = fits_rdecomp_short (cbuf, nelem, (unsigned short *)idata,
                        tilelen, blocksize);
            tiledatatype = TSHORT;
        } else {
            *status = fits_rdecomp (cbuf, nelem, (unsigned int *)idata,
                         tilelen, blocksize);
            tiledatatype = TINT;
        }

    /* ************************************************************* */
    } else if ((infptr->Fptr)->compress_type == HCOMPRESS_1)  {

        smooth = (infptr->Fptr)->hcomp_smooth;

        if ( ((infptr->Fptr)->zbitpix == BYTE_IMG || (infptr->Fptr)->zbitpix == SHORT_IMG)) {
            *status = fits_hdecompress(cbuf, smooth, idata, &nx, &ny,
	        &scale, status);
        } else {  /* zbitpix = LONG_IMG (32) */
            /* idata must have been allocated twice as large for this to work */
            *status = fits_hdecompress64(cbuf, smooth, (LONGLONG *) idata, &nx, &ny,
	        &scale, status);
        }       

        tiledatatype = TINT;

    /* ************************************************************* */
    } else if ((infptr->Fptr)->compress_type == PLIO_1) {

        pl_l2pi ((short *) cbuf, 1, idata, tilelen);  /* uncompress the data */
        tiledatatype = TINT;

    /* ************************************************************* */
    } else if ( ((infptr->Fptr)->compress_type == GZIP_1) ||
                ((infptr->Fptr)->compress_type == GZIP_2) ) {

        uncompress2mem_from_mem ((char *)cbuf, nelem,
             (char **) &idata, &idatalen, realloc, &tilebytesize, status);

        /* determine the data type of the uncompressed array, and */
	/*  do byte unshuffling and unswapping if needed */
	if (tilebytesize == (size_t) (tilelen * 2)) {
	    /* this is a short I*2 array */
            tiledatatype = TSHORT;

            if ( (infptr->Fptr)->compress_type == GZIP_2 )
		    fits_unshuffle_2bytes((char *) idata, tilelen, status);

#if BYTESWAPPED
            ffswap2((short *) idata, tilelen);
#endif

	} else if (tilebytesize == (size_t) (tilelen * 4)) {
	    /* this is a int I*4 array (or maybe R*4) */
            tiledatatype = TINT;

            if ( (infptr->Fptr)->compress_type == GZIP_2 )
		    fits_unshuffle_4bytes((char *) idata, tilelen, status);

#if BYTESWAPPED
            ffswap4(idata, tilelen);
#endif

	} else if (tilebytesize == (size_t) (tilelen * 8)) {
	    /* this is a R*8 double array */
            tiledatatype = TDOUBLE;

            if ( (infptr->Fptr)->compress_type == GZIP_2 )
		    fits_unshuffle_8bytes((char *) idata, tilelen, status);
#if BYTESWAPPED
            ffswap8((double *) idata, tilelen);
#endif

        } else if (tilebytesize == (size_t) tilelen) {
	    
	    /* this is an unsigned char I*1 array */
            tiledatatype = TBYTE;

        } else {
            ffpmsg("error: uncompressed tile has wrong size");
            free(idata);
            return (*status = DATA_DECOMPRESSION_ERR);
        }

    /* ************************************************************* */
    } else if ((infptr->Fptr)->compress_type == BZIP2_1) {

/*  BZIP2 is not supported in the public release; this is only for test purposes 

        if (BZ2_bzBuffToBuffDecompress ((char *) idata, &idatalen, 
		(char *)cbuf, (unsigned int) nelem, 0, 0) )
*/
        {
            ffpmsg("bzip2 decompression error");
            free(idata);
            free (cbuf);
            return (*status = DATA_DECOMPRESSION_ERR);
        }

        if ((infptr->Fptr)->zbitpix == BYTE_IMG) {
	     tiledatatype = TBYTE;
        } else if ((infptr->Fptr)->zbitpix == SHORT_IMG) {
  	     tiledatatype = TSHORT;
#if BYTESWAPPED
            ffswap2((short *) idata, tilelen);
#endif
	} else {
  	     tiledatatype = TINT;
#if BYTESWAPPED
            ffswap4(idata, tilelen);
#endif
	}

    /* ************************************************************* */
    } else {
        ffpmsg("unknown compression algorithm");
        free(idata);
        return (*status = DATA_DECOMPRESSION_ERR);
    }

    free(cbuf);
    if (*status)  {  /* error uncompressing the tile */
            free(idata);
            return (*status);
    }

    /* ************************************************************* */
    /* copy the uncompressed tile data to the output buffer, doing */
    /* null checking, datatype conversion and linear scaling, if necessary */

    if (nulval == 0)
         nulval = &dummy;  /* set address to dummy value */

    if (datatype == TSHORT)
    {
        pixlen = sizeof(short);

	if ((infptr->Fptr)->quantize_level == NO_QUANTIZE) {
	 /* the floating point pixels were losselessly compressed with GZIP */
	 /* Just have to copy the values to the output array */
	 
          if (tiledatatype == TINT) {
              fffr4i2((float *) idata, tilelen, bscale, bzero, nullcheck,   
                *(short *) nulval, bnullarray, anynul,
                (short *) buffer, status);
          } else {
              fffr8i2((double *) idata, tilelen, bscale, bzero, nullcheck,   
                *(short *) nulval, bnullarray, anynul,
                (short *) buffer, status);
          }
        } else if (tiledatatype == TINT)
          if ((infptr->Fptr)->compress_type == PLIO_1 &&
	    bzero == 0. && actual_bzero == 32768.) {
	    /* special case where unsigned 16-bit integers have been */
	    /* offset by +32768 when using PLIO */
            fffi4i2(idata, tilelen, bscale, -32768., nullcheck, tnull,
             *(short *) nulval, bnullarray, anynul,
            (short *) buffer, status);
          } else {
            fffi4i2(idata, tilelen, bscale, bzero, nullcheck, tnull,
             *(short *) nulval, bnullarray, anynul,
            (short *) buffer, status);
          }
        else if (tiledatatype == TSHORT)
          fffi2i2((short *)idata, tilelen, bscale, bzero, nullcheck, (short) tnull,
           *(short *) nulval, bnullarray, anynul,
          (short *) buffer, status);
        else if (tiledatatype == TBYTE)
          fffi1i2((unsigned char *)idata, tilelen, bscale, bzero, nullcheck, (unsigned char) tnull,
           *(short *) nulval, bnullarray, anynul,
          (short *) buffer, status);
    }
    else if (datatype == TINT)
    {
        pixlen = sizeof(int);

	if ((infptr->Fptr)->quantize_level == NO_QUANTIZE) {
	 /* the floating point pixels were losselessly compressed with GZIP */
	 /* Just have to copy the values to the output array */
	 
          if (tiledatatype == TINT) {
              fffr4int((float *) idata, tilelen, bscale, bzero, nullcheck,   
                *(int *) nulval, bnullarray, anynul,
                (int *) buffer, status);
          } else {
              fffr8int((double *) idata, tilelen, bscale, bzero, nullcheck,   
                *(int *) nulval, bnullarray, anynul,
                (int *) buffer, status);
          }
        } else if (tiledatatype == TINT)
          fffi4int(idata, (long) tilelen, bscale, bzero, nullcheck, tnull,
           *(int *) nulval, bnullarray, anynul,
           (int *) buffer, status);
        else if (tiledatatype == TSHORT)
          fffi2int((short *)idata, tilelen, bscale, bzero, nullcheck, (short) tnull,
           *(int *) nulval, bnullarray, anynul,
           (int *) buffer, status);
        else if (tiledatatype == TBYTE)
          fffi1int((unsigned char *)idata, tilelen, bscale, bzero, nullcheck, (unsigned char) tnull,
           *(int *) nulval, bnullarray, anynul,
           (int *) buffer, status);
    }
    else if (datatype == TLONG)
    {
        pixlen = sizeof(long);

	if ((infptr->Fptr)->quantize_level == NO_QUANTIZE) {
	 /* the floating point pixels were losselessly compressed with GZIP */
	 /* Just have to copy the values to the output array */
	 
          if (tiledatatype == TINT) {
              fffr4i4((float *) idata, tilelen, bscale, bzero, nullcheck,   
                *(long *) nulval, bnullarray, anynul,
                (long *) buffer, status);
          } else {
              fffr8i4((double *) idata, tilelen, bscale, bzero, nullcheck,   
                *(long *) nulval, bnullarray, anynul,
                (long *) buffer, status);
          }
        } else if (tiledatatype == TINT)
          fffi4i4(idata, tilelen, bscale, bzero, nullcheck, tnull,
           *(long *) nulval, bnullarray, anynul,
            (long *) buffer, status);
        else if (tiledatatype == TSHORT)
          fffi2i4((short *)idata, tilelen, bscale, bzero, nullcheck, (short) tnull,
           *(long *) nulval, bnullarray, anynul,
            (long *) buffer, status);
        else if (tiledatatype == TBYTE)
          fffi1i4((unsigned char *)idata, tilelen, bscale, bzero, nullcheck, (unsigned char) tnull,
           *(long *) nulval, bnullarray, anynul,
            (long *) buffer, status);
    }
    else if (datatype == TFLOAT)
    {
        pixlen = sizeof(float);
        if (nulval) {
	      fnulval = *(float *) nulval;
	}
 
	if ((infptr->Fptr)->quantize_level == NO_QUANTIZE) {
	 /* the floating point pixels were losselessly compressed with GZIP */
	 /* Just have to copy the values to the output array */
	 
          if (tiledatatype == TINT) {
              fffr4r4((float *) idata, tilelen, bscale, bzero, nullcheck,   
                fnulval, bnullarray, anynul,
                (float *) buffer, status);
          } else {
              fffr8r4((double *) idata, tilelen, bscale, bzero, nullcheck,   
                fnulval, bnullarray, anynul,
                (float *) buffer, status);
          }
	
        } else if ((infptr->Fptr)->quantize_dither == SUBTRACTIVE_DITHER_1) {

         /* use the new dithering algorithm (introduced in July 2009) */

         if (tiledatatype == TINT)
          unquantize_i4r4(nrow + (infptr->Fptr)->dither_offset - 1, idata, 
	   tilelen, bscale, bzero, nullcheck, tnull,
           fnulval, bnullarray, anynul,
            (float *) buffer, status);
         else if (tiledatatype == TSHORT)
          unquantize_i2r4(nrow + (infptr->Fptr)->dither_offset - 1, (short *)idata, 
	   tilelen, bscale, bzero, nullcheck, (short) tnull,
           fnulval, bnullarray, anynul,
            (float *) buffer, status);
         else if (tiledatatype == TBYTE)
          unquantize_i1r4(nrow + (infptr->Fptr)->dither_offset - 1, (unsigned char *)idata, 
	   tilelen, bscale, bzero, nullcheck, (unsigned char) tnull,
           fnulval, bnullarray, anynul,
            (float *) buffer, status);

        } else {  /* use the old "round to nearest level" quantization algorithm */

         if (tiledatatype == TINT)
          fffi4r4(idata, tilelen, bscale, bzero, nullcheck, tnull,  
           fnulval, bnullarray, anynul,
            (float *) buffer, status);
         else if (tiledatatype == TSHORT)
          fffi2r4((short *)idata, tilelen, bscale, bzero, nullcheck, (short) tnull,  
           fnulval, bnullarray, anynul,
            (float *) buffer, status);
         else if (tiledatatype == TBYTE)
          fffi1r4((unsigned char *)idata, tilelen, bscale, bzero, nullcheck, (unsigned char) tnull,
           fnulval, bnullarray, anynul,
            (float *) buffer, status);
	}
    }
    else if (datatype == TDOUBLE)
    {
        pixlen = sizeof(double);
        if (nulval) {
	     dnulval = *(double *) nulval;
	}

	if ((infptr->Fptr)->quantize_level == NO_QUANTIZE) {
	 /* the floating point pixels were losselessly compressed with GZIP */
	 /* Just have to copy the values to the output array */

          if (tiledatatype == TINT) {
              fffr4r8((float *) idata, tilelen, bscale, bzero, nullcheck,   
                dnulval, bnullarray, anynul,
                (double *) buffer, status);
          } else {
              fffr8r8((double *) idata, tilelen, bscale, bzero, nullcheck,   
                dnulval, bnullarray, anynul,
                (double *) buffer, status);
          }
	
	} else if ((infptr->Fptr)->quantize_dither == SUBTRACTIVE_DITHER_1) {

         /* use the new dithering algorithm (introduced in July 2009) */
         if (tiledatatype == TINT)
          unquantize_i4r8(nrow + (infptr->Fptr)->dither_offset - 1, idata,
	   tilelen, bscale, bzero, nullcheck, tnull,
           dnulval, bnullarray, anynul,
            (double *) buffer, status);
         else if (tiledatatype == TSHORT)
          unquantize_i2r8(nrow + (infptr->Fptr)->dither_offset - 1, (short *)idata,
	   tilelen, bscale, bzero, nullcheck, (short) tnull,
           dnulval, bnullarray, anynul,
            (double *) buffer, status);
         else if (tiledatatype == TBYTE)
          unquantize_i1r8(nrow + (infptr->Fptr)->dither_offset - 1, (unsigned char *)idata,
	   tilelen, bscale, bzero, nullcheck, (unsigned char) tnull,
           dnulval, bnullarray, anynul,
            (double *) buffer, status);

        } else {  /* use the old "round to nearest level" quantization algorithm */

         if (tiledatatype == TINT)
          fffi4r8(idata, tilelen, bscale, bzero, nullcheck, tnull,
           dnulval, bnullarray, anynul,
            (double *) buffer, status);
         else if (tiledatatype == TSHORT)
          fffi2r8((short *)idata, tilelen, bscale, bzero, nullcheck, (short) tnull,
           dnulval, bnullarray, anynul,
            (double *) buffer, status);
         else if (tiledatatype == TBYTE)
          fffi1r8((unsigned char *)idata, tilelen, bscale, bzero, nullcheck, (unsigned char) tnull,
           dnulval, bnullarray, anynul,
            (double *) buffer, status);
	}
    }
    else if (datatype == TBYTE)
    {
        pixlen = sizeof(char);
        if (tiledatatype == TINT)
          fffi4i1(idata, tilelen, bscale, bzero, nullcheck, tnull,
           *(unsigned char *) nulval, bnullarray, anynul,
            (unsigned char *) buffer, status);
        else if (tiledatatype == TSHORT)
          fffi2i1((short *)idata, tilelen, bscale, bzero, nullcheck, (short) tnull,
           *(unsigned char *) nulval, bnullarray, anynul,
            (unsigned char *) buffer, status);
        else if (tiledatatype == TBYTE)
          fffi1i1((unsigned char *)idata, tilelen, bscale, bzero, nullcheck, (unsigned char) tnull,
           *(unsigned char *) nulval, bnullarray, anynul,
            (unsigned char *) buffer, status);
    }
    else if (datatype == TSBYTE)
    {
        pixlen = sizeof(char);
        if (tiledatatype == TINT)
          fffi4s1(idata, tilelen, bscale, bzero, nullcheck, tnull,
           *(signed char *) nulval, bnullarray, anynul,
            (signed char *) buffer, status);
        else if (tiledatatype == TSHORT)
          fffi2s1((short *)idata, tilelen, bscale, bzero, nullcheck, (short) tnull,
           *(signed char *) nulval, bnullarray, anynul,
            (signed char *) buffer, status);
        else if (tiledatatype == TBYTE)
          fffi1s1((unsigned char *)idata, tilelen, bscale, bzero, nullcheck, (unsigned char) tnull,
           *(signed char *) nulval, bnullarray, anynul,
            (signed char *) buffer, status);
    }
    else if (datatype == TUSHORT)
    {
        pixlen = sizeof(short);

	if ((infptr->Fptr)->quantize_level == NO_QUANTIZE) {
	 /* the floating point pixels were losselessly compressed with GZIP */
	 /* Just have to copy the values to the output array */
	 
          if (tiledatatype == TINT) {
              fffr4u2((float *) idata, tilelen, bscale, bzero, nullcheck,   
                *(unsigned short *) nulval, bnullarray, anynul,
                (unsigned short *) buffer, status);
          } else {
              fffr8u2((double *) idata, tilelen, bscale, bzero, nullcheck,   
                *(unsigned short *) nulval, bnullarray, anynul,
                (unsigned short *) buffer, status);
          }
        } else if (tiledatatype == TINT)
          fffi4u2(idata, tilelen, bscale, bzero, nullcheck, tnull,
           *(unsigned short *) nulval, bnullarray, anynul,
            (unsigned short *) buffer, status);
        else if (tiledatatype == TSHORT)
          fffi2u2((short *)idata, tilelen, bscale, bzero, nullcheck, (short) tnull,
           *(unsigned short *) nulval, bnullarray, anynul,
            (unsigned short *) buffer, status);
        else if (tiledatatype == TBYTE)
          fffi1u2((unsigned char *)idata, tilelen, bscale, bzero, nullcheck, (unsigned char) tnull,
           *(unsigned short *) nulval, bnullarray, anynul,
            (unsigned short *) buffer, status);
    }
    else if (datatype == TUINT)
    {
        pixlen = sizeof(int);

	if ((infptr->Fptr)->quantize_level == NO_QUANTIZE) {
	 /* the floating point pixels were losselessly compressed with GZIP */
	 /* Just have to copy the values to the output array */
	 
          if (tiledatatype == TINT) {
              fffr4uint((float *) idata, tilelen, bscale, bzero, nullcheck,   
                *(unsigned int *) nulval, bnullarray, anynul,
                (unsigned int *) buffer, status);
          } else {
              fffr8uint((double *) idata, tilelen, bscale, bzero, nullcheck,   
                *(unsigned int *) nulval, bnullarray, anynul,
                (unsigned int *) buffer, status);
          }
        } else         if (tiledatatype == TINT)
          fffi4uint(idata, tilelen, bscale, bzero, nullcheck, tnull,
           *(unsigned int *) nulval, bnullarray, anynul,
            (unsigned int *) buffer, status);
        else if (tiledatatype == TSHORT)
          fffi2uint((short *)idata, tilelen, bscale, bzero, nullcheck, (short) tnull,
           *(unsigned int *) nulval, bnullarray, anynul,
            (unsigned int *) buffer, status);
        else if (tiledatatype == TBYTE)
          fffi1uint((unsigned char *)idata, tilelen, bscale, bzero, nullcheck, (unsigned char) tnull,
           *(unsigned int *) nulval, bnullarray, anynul,
            (unsigned int *) buffer, status);
    }
    else if (datatype == TULONG)
    {
        pixlen = sizeof(long);

	if ((infptr->Fptr)->quantize_level == NO_QUANTIZE) {
	 /* the floating point pixels were losselessly compressed with GZIP */
	 /* Just have to copy the values to the output array */
	 
          if (tiledatatype == TINT) {
              fffr4u4((float *) idata, tilelen, bscale, bzero, nullcheck,   
                *(unsigned long *) nulval, bnullarray, anynul,
                (unsigned long *) buffer, status);
          } else {
              fffr8u4((double *) idata, tilelen, bscale, bzero, nullcheck,   
                *(unsigned long *) nulval, bnullarray, anynul,
                (unsigned long *) buffer, status);
          }
        } else if (tiledatatype == TINT)
          fffi4u4(idata, tilelen, bscale, bzero, nullcheck, tnull,
           *(unsigned long *) nulval, bnullarray, anynul, 
            (unsigned long *) buffer, status);
        else if (tiledatatype == TSHORT)
          fffi2u4((short *)idata, tilelen, bscale, bzero, nullcheck, (short) tnull,
           *(unsigned long *) nulval, bnullarray, anynul, 
            (unsigned long *) buffer, status);
        else if (tiledatatype == TBYTE)
          fffi1u4((unsigned char *)idata, tilelen, bscale, bzero, nullcheck, (unsigned char) tnull,
           *(unsigned long *) nulval, bnullarray, anynul, 
            (unsigned long *) buffer, status);
    }
    else
         *status = BAD_DATATYPE;

    free(idata);  /* don't need the uncompressed tile any more */

    /* **************************************************************** */
    /* cache the tile, in case the application wants it again  */

    /*   Don't cache the tile if tile is a single row of the image; 
         it is less likely that the cache will be used in this cases,
	 so it is not worth the time and the memory overheads.
    */
    if ((infptr->Fptr)->znaxis[0]   != (infptr->Fptr)->tilesize[0] ||
        (infptr->Fptr)->tilesize[1] != 1 )
    {
      tilesize = pixlen * tilelen;

      /* check that tile size/type has not changed */
      if (tilesize != (infptr->Fptr)->tiledatasize ||
        datatype != (infptr->Fptr)->tiletype )  {

        if ((infptr->Fptr)->tiledata) {
            free((infptr->Fptr)->tiledata);	    
        }
	
        (infptr->Fptr)->tiledata = 0;

        if ((infptr->Fptr)->tilenullarray) {
            free((infptr->Fptr)->tilenullarray);
        }
	
        (infptr->Fptr)->tilenullarray = 0;
        (infptr->Fptr)->tilerow = 0;
        (infptr->Fptr)->tiledatasize = 0;
        (infptr->Fptr)->tiletype = 0;

        /* allocate new array(s) */
	(infptr->Fptr)->tiledata = malloc(tilesize);
	if ((infptr->Fptr)->tiledata == 0)
	   return (*status);

        if (nullcheck == 2) {  /* also need array of null pixel flags */
	    (infptr->Fptr)->tilenullarray = malloc(tilelen);
	    if ((infptr->Fptr)->tilenullarray == 0)
	        return (*status);
        }

        (infptr->Fptr)->tiledatasize = tilesize;
        (infptr->Fptr)->tiletype = datatype;
      }

      /* copy the tile array(s) into cache buffer */
      memcpy((infptr->Fptr)->tiledata, buffer, tilesize);

      if (nullcheck == 2) {
	    if ((infptr->Fptr)->tilenullarray == 0)  {
       	      (infptr->Fptr)->tilenullarray = malloc(tilelen);
            }
            memcpy((infptr->Fptr)->tilenullarray, bnullarray, tilelen);
      }

      (infptr->Fptr)->tilerow = nrow;
      (infptr->Fptr)->tileanynull = *anynul;
    }

    return (*status);
}
/*--------------------------------------------------------------------------*/
int imcomp_copy_overlap (
    char *tile,         /* I - multi dimensional array of tile pixels */
    int pixlen,         /* I - number of bytes in each tile or image pixel */
    int ndim,           /* I - number of dimension in the tile and image */
    long *tfpixel,      /* I - first pixel number in each dim. of the tile */
    long *tlpixel,      /* I - last pixel number in each dim. of the tile */
    char *bnullarray,   /* I - array of null flags; used if nullcheck = 2 */
    char *image,        /* O - multi dimensional output image */
    long *fpixel,       /* I - first pixel number in each dim. of the image */
    long *lpixel,       /* I - last pixel number in each dim. of the image */
    long *ininc,        /* I - increment to be applied in each image dimen. */
    int nullcheck,      /* I - 0, 1: do nothing; 2: set nullarray for nulls */
    char *nullarray, 
    int *status)

/* 
  copy the intersecting pixels from a decompressed tile to the output image. 
  Both the tile and the image must have the same number of dimensions. 
*/
{
    long imgdim[MAX_COMPRESS_DIM]; /* product of preceding dimensions in the */
                                   /* output image, allowing for inc factor */
    long tiledim[MAX_COMPRESS_DIM]; /* product of preceding dimensions in the */
                                 /* tile, array;  inc factor is not relevant */
    long imgfpix[MAX_COMPRESS_DIM]; /* 1st img pix overlapping tile: 0 base, */
                                    /*  allowing for inc factor */
    long imglpix[MAX_COMPRESS_DIM]; /* last img pix overlapping tile 0 base, */
                                    /*  allowing for inc factor */
    long tilefpix[MAX_COMPRESS_DIM]; /* 1st tile pix overlapping img 0 base, */
                                    /*  allowing for inc factor */
    long inc[MAX_COMPRESS_DIM]; /* local copy of input ininc */
    long i1, i2, i3, i4;   /* offset along each axis of the image */
    long it1, it2, it3, it4;
    long im1, im2, im3, im4;  /* offset to image pixel, allowing for inc */
    long ipos, tf, tl;
    long t2, t3, t4;   /* offset along each axis of the tile */
    long tilepix, imgpix, tilepixbyte, imgpixbyte;
    int ii, overlap_bytes, overlap_flags;

    if (*status > 0)
        return(*status);

    for (ii = 0; ii < MAX_COMPRESS_DIM; ii++)
    {
        /* set default values for higher dimensions */
        inc[ii] = 1;
        imgdim[ii] = 1;
        tiledim[ii] = 1;
        imgfpix[ii] = 0;
        imglpix[ii] = 0;
        tilefpix[ii] = 0;
    }

    /* ------------------------------------------------------------ */
    /* calc amount of overlap in each dimension; if there is zero   */
    /* overlap in any dimension then just return  */
    /* ------------------------------------------------------------ */
    
    for (ii = 0; ii < ndim; ii++)
    {
        if (tlpixel[ii] < fpixel[ii] || tfpixel[ii] > lpixel[ii])
            return(*status);  /* there are no overlapping pixels */

        inc[ii] = ininc[ii];

        /* calc dimensions of the output image section */
        imgdim[ii] = (lpixel[ii] - fpixel[ii]) / labs(inc[ii]) + 1;
        if (imgdim[ii] < 1)
            return(*status = NEG_AXIS);

        /* calc dimensions of the tile */
        tiledim[ii] = tlpixel[ii] - tfpixel[ii] + 1;
        if (tiledim[ii] < 1)
            return(*status = NEG_AXIS);

        if (ii > 0)
           tiledim[ii] *= tiledim[ii - 1];  /* product of dimensions */

        /* first and last pixels in image that overlap with the tile, 0 base */
        tf = tfpixel[ii] - 1;
        tl = tlpixel[ii] - 1;

        /* skip this plane if it falls in the cracks of the subsampled image */
        while ((tf-(fpixel[ii] - 1)) % labs(inc[ii]))
        {
           tf++;
           if (tf > tl)
             return(*status);  /* no overlapping pixels */
        }

        while ((tl-(fpixel[ii] - 1)) % labs(inc[ii]))
        {
           tl--;
           if (tf > tl)
             return(*status);  /* no overlapping pixels */
        }
        imgfpix[ii] = maxvalue((tf - fpixel[ii] +1) / labs(inc[ii]) , 0);
        imglpix[ii] = minvalue((tl - fpixel[ii] +1) / labs(inc[ii]) ,
                               imgdim[ii] - 1);

        /* first pixel in the tile that overlaps with the image (0 base) */
        tilefpix[ii] = maxvalue(fpixel[ii] - tfpixel[ii], 0);

        while ((tfpixel[ii] + tilefpix[ii] - fpixel[ii]) % labs(inc[ii]))
        {
           (tilefpix[ii])++;
           if (tilefpix[ii] >= tiledim[ii])
              return(*status);  /* no overlapping pixels */
        }
/*
printf("ii tfpixel, tlpixel %d %d %d \n",ii, tfpixel[ii], tlpixel[ii]);
printf("ii, tf, tl, imgfpix,imglpix, tilefpix %d %d %d %d %d %d\n",ii,
 tf,tl,imgfpix[ii], imglpix[ii],tilefpix[ii]);
*/
        if (ii > 0)
           imgdim[ii] *= imgdim[ii - 1];  /* product of dimensions */
    }

    /* ---------------------------------------------------------------- */
    /* calc number of pixels in each row (first dimension) that overlap */
    /* multiply by pixlen to get number of bytes to copy in each loop   */
    /* ---------------------------------------------------------------- */

    if (inc[0] != 1)
       overlap_flags = 1;  /* can only copy 1 pixel at a time */
    else
       overlap_flags = imglpix[0] - imgfpix[0] + 1;  /* can copy whole row */

    overlap_bytes = overlap_flags * pixlen;

    /* support up to 5 dimensions for now */
    for (i4 = 0, it4=0; i4 <= imglpix[4] - imgfpix[4]; i4++, it4++)
    {
     /* increment plane if it falls in the cracks of the subsampled image */
     while (ndim > 4 &&  (tfpixel[4] + tilefpix[4] - fpixel[4] + it4)
                          % labs(inc[4]) != 0)
        it4++;

       /* offset to start of hypercube */
       if (inc[4] > 0)
          im4 = (i4 + imgfpix[4]) * imgdim[3];
       else
          im4 = imgdim[4] - (i4 + 1 + imgfpix[4]) * imgdim[3];

      t4 = (tilefpix[4] + it4) * tiledim[3];
      for (i3 = 0, it3=0; i3 <= imglpix[3] - imgfpix[3]; i3++, it3++)
      {
       /* increment plane if it falls in the cracks of the subsampled image */
       while (ndim > 3 &&  (tfpixel[3] + tilefpix[3] - fpixel[3] + it3)
                            % labs(inc[3]) != 0)
          it3++;

       /* offset to start of cube */
       if (inc[3] > 0)
          im3 = (i3 + imgfpix[3]) * imgdim[2] + im4;
       else
          im3 = imgdim[3] - (i3 + 1 + imgfpix[3]) * imgdim[2] + im4;

       t3 = (tilefpix[3] + it3) * tiledim[2] + t4;

       /* loop through planes of the image */
       for (i2 = 0, it2=0; i2 <= imglpix[2] - imgfpix[2]; i2++, it2++)
       {
          /* incre plane if it falls in the cracks of the subsampled image */
          while (ndim > 2 &&  (tfpixel[2] + tilefpix[2] - fpixel[2] + it2)
                               % labs(inc[2]) != 0)
             it2++;

          /* offset to start of plane */
          if (inc[2] > 0)
             im2 = (i2 + imgfpix[2]) * imgdim[1] + im3;
          else
             im2 = imgdim[2] - (i2 + 1 + imgfpix[2]) * imgdim[1] + im3;

          t2 = (tilefpix[2] + it2) * tiledim[1] + t3;

          /* loop through rows of the image */
          for (i1 = 0, it1=0; i1 <= imglpix[1] - imgfpix[1]; i1++, it1++)
          {
             /* incre row if it falls in the cracks of the subsampled image */
             while (ndim > 1 &&  (tfpixel[1] + tilefpix[1] - fpixel[1] + it1)
                                  % labs(inc[1]) != 0)
                it1++;

             /* calc position of first pixel in tile to be copied */
             tilepix = tilefpix[0] + (tilefpix[1] + it1) * tiledim[0] + t2;

             /* offset to start of row */
             if (inc[1] > 0)
                im1 = (i1 + imgfpix[1]) * imgdim[0] + im2;
             else
                im1 = imgdim[1] - (i1 + 1 + imgfpix[1]) * imgdim[0] + im2;
/*
printf("inc = %d %d %d %d\n",inc[0],inc[1],inc[2],inc[3]);
printf("im1,im2,im3,im4 = %d %d %d %d\n",im1,im2,im3,im4);
*/
             /* offset to byte within the row */
             if (inc[0] > 0)
                imgpix = imgfpix[0] + im1;
             else
                imgpix = imgdim[0] - 1 - imgfpix[0] + im1;
/*
printf("tilefpix0,1, imgfpix1, it1, inc1, t2= %d %d %d %d %d %d\n",
       tilefpix[0],tilefpix[1],imgfpix[1],it1,inc[1], t2);
printf("i1, it1, tilepix, imgpix %d %d %d %d \n", i1, it1, tilepix, imgpix);
*/
             /* loop over pixels along one row of the image */
             for (ipos = imgfpix[0]; ipos <= imglpix[0]; ipos += overlap_flags)
             {
               if (nullcheck == 2)
               {
                   /* copy overlapping null flags from tile to image */
                   memcpy(nullarray + imgpix, bnullarray + tilepix,
                          overlap_flags);
               }

               /* convert from image pixel to byte offset */
               tilepixbyte = tilepix * pixlen;
               imgpixbyte  = imgpix  * pixlen;
/*
printf("  tilepix, tilepixbyte, imgpix, imgpixbyte= %d %d %d %d\n",
          tilepix, tilepixbyte, imgpix, imgpixbyte);
*/
               /* copy overlapping row of pixels from tile to image */
               memcpy(image + imgpixbyte, tile + tilepixbyte, overlap_bytes);

               tilepix += (overlap_flags * labs(inc[0]));
               if (inc[0] > 0)
                 imgpix += overlap_flags;
               else
                 imgpix -= overlap_flags;
            }
          }
        }
      }
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
int imcomp_merge_overlap (
    char *tile,         /* O - multi dimensional array of tile pixels */
    int pixlen,         /* I - number of bytes in each tile or image pixel */
    int ndim,           /* I - number of dimension in the tile and image */
    long *tfpixel,      /* I - first pixel number in each dim. of the tile */
    long *tlpixel,      /* I - last pixel number in each dim. of the tile */
    char *bnullarray,   /* I - array of null flags; used if nullcheck = 2 */
    char *image,        /* I - multi dimensional output image */
    long *fpixel,       /* I - first pixel number in each dim. of the image */
    long *lpixel,       /* I - last pixel number in each dim. of the image */
    int nullcheck,      /* I - 0, 1: do nothing; 2: set nullarray for nulls */
    int *status)

/* 
  Similar to imcomp_copy_overlap, except it copies the overlapping pixels from
  the 'image' to the 'tile'.
*/
{
    long imgdim[MAX_COMPRESS_DIM]; /* product of preceding dimensions in the */
                                   /* output image, allowing for inc factor */
    long tiledim[MAX_COMPRESS_DIM]; /* product of preceding dimensions in the */
                                 /* tile, array;  inc factor is not relevant */
    long imgfpix[MAX_COMPRESS_DIM]; /* 1st img pix overlapping tile: 0 base, */
                                    /*  allowing for inc factor */
    long imglpix[MAX_COMPRESS_DIM]; /* last img pix overlapping tile 0 base, */
                                    /*  allowing for inc factor */
    long tilefpix[MAX_COMPRESS_DIM]; /* 1st tile pix overlapping img 0 base, */
                                    /*  allowing for inc factor */
    long inc[MAX_COMPRESS_DIM]; /* local copy of input ininc */
    long i1, i2, i3, i4;   /* offset along each axis of the image */
    long it1, it2, it3, it4;
    long im1, im2, im3, im4;  /* offset to image pixel, allowing for inc */
    long ipos, tf, tl;
    long t2, t3, t4;   /* offset along each axis of the tile */
    long tilepix, imgpix, tilepixbyte, imgpixbyte;
    int ii, overlap_bytes, overlap_flags;

    if (*status > 0)
        return(*status);

    for (ii = 0; ii < MAX_COMPRESS_DIM; ii++)
    {
        /* set default values for higher dimensions */
        inc[ii] = 1;
        imgdim[ii] = 1;
        tiledim[ii] = 1;
        imgfpix[ii] = 0;
        imglpix[ii] = 0;
        tilefpix[ii] = 0;
    }

    /* ------------------------------------------------------------ */
    /* calc amount of overlap in each dimension; if there is zero   */
    /* overlap in any dimension then just return  */
    /* ------------------------------------------------------------ */
    
    for (ii = 0; ii < ndim; ii++)
    {
        if (tlpixel[ii] < fpixel[ii] || tfpixel[ii] > lpixel[ii])
            return(*status);  /* there are no overlapping pixels */

        /* calc dimensions of the output image section */
        imgdim[ii] = (lpixel[ii] - fpixel[ii]) / labs(inc[ii]) + 1;
        if (imgdim[ii] < 1)
            return(*status = NEG_AXIS);

        /* calc dimensions of the tile */
        tiledim[ii] = tlpixel[ii] - tfpixel[ii] + 1;
        if (tiledim[ii] < 1)
            return(*status = NEG_AXIS);

        if (ii > 0)
           tiledim[ii] *= tiledim[ii - 1];  /* product of dimensions */

        /* first and last pixels in image that overlap with the tile, 0 base */
        tf = tfpixel[ii] - 1;
        tl = tlpixel[ii] - 1;

        /* skip this plane if it falls in the cracks of the subsampled image */
        while ((tf-(fpixel[ii] - 1)) % labs(inc[ii]))
        {
           tf++;
           if (tf > tl)
             return(*status);  /* no overlapping pixels */
        }

        while ((tl-(fpixel[ii] - 1)) % labs(inc[ii]))
        {
           tl--;
           if (tf > tl)
             return(*status);  /* no overlapping pixels */
        }
        imgfpix[ii] = maxvalue((tf - fpixel[ii] +1) / labs(inc[ii]) , 0);
        imglpix[ii] = minvalue((tl - fpixel[ii] +1) / labs(inc[ii]) ,
                               imgdim[ii] - 1);

        /* first pixel in the tile that overlaps with the image (0 base) */
        tilefpix[ii] = maxvalue(fpixel[ii] - tfpixel[ii], 0);

        while ((tfpixel[ii] + tilefpix[ii] - fpixel[ii]) % labs(inc[ii]))
        {
           (tilefpix[ii])++;
           if (tilefpix[ii] >= tiledim[ii])
              return(*status);  /* no overlapping pixels */
        }
/*
printf("ii tfpixel, tlpixel %d %d %d \n",ii, tfpixel[ii], tlpixel[ii]);
printf("ii, tf, tl, imgfpix,imglpix, tilefpix %d %d %d %d %d %d\n",ii,
 tf,tl,imgfpix[ii], imglpix[ii],tilefpix[ii]);
*/
        if (ii > 0)
           imgdim[ii] *= imgdim[ii - 1];  /* product of dimensions */
    }

    /* ---------------------------------------------------------------- */
    /* calc number of pixels in each row (first dimension) that overlap */
    /* multiply by pixlen to get number of bytes to copy in each loop   */
    /* ---------------------------------------------------------------- */

    if (inc[0] != 1)
       overlap_flags = 1;  /* can only copy 1 pixel at a time */
    else
       overlap_flags = imglpix[0] - imgfpix[0] + 1;  /* can copy whole row */

    overlap_bytes = overlap_flags * pixlen;

    /* support up to 5 dimensions for now */
    for (i4 = 0, it4=0; i4 <= imglpix[4] - imgfpix[4]; i4++, it4++)
    {
     /* increment plane if it falls in the cracks of the subsampled image */
     while (ndim > 4 &&  (tfpixel[4] + tilefpix[4] - fpixel[4] + it4)
                          % labs(inc[4]) != 0)
        it4++;

       /* offset to start of hypercube */
       if (inc[4] > 0)
          im4 = (i4 + imgfpix[4]) * imgdim[3];
       else
          im4 = imgdim[4] - (i4 + 1 + imgfpix[4]) * imgdim[3];

      t4 = (tilefpix[4] + it4) * tiledim[3];
      for (i3 = 0, it3=0; i3 <= imglpix[3] - imgfpix[3]; i3++, it3++)
      {
       /* increment plane if it falls in the cracks of the subsampled image */
       while (ndim > 3 &&  (tfpixel[3] + tilefpix[3] - fpixel[3] + it3)
                            % labs(inc[3]) != 0)
          it3++;

       /* offset to start of cube */
       if (inc[3] > 0)
          im3 = (i3 + imgfpix[3]) * imgdim[2] + im4;
       else
          im3 = imgdim[3] - (i3 + 1 + imgfpix[3]) * imgdim[2] + im4;

       t3 = (tilefpix[3] + it3) * tiledim[2] + t4;

       /* loop through planes of the image */
       for (i2 = 0, it2=0; i2 <= imglpix[2] - imgfpix[2]; i2++, it2++)
       {
          /* incre plane if it falls in the cracks of the subsampled image */
          while (ndim > 2 &&  (tfpixel[2] + tilefpix[2] - fpixel[2] + it2)
                               % labs(inc[2]) != 0)
             it2++;

          /* offset to start of plane */
          if (inc[2] > 0)
             im2 = (i2 + imgfpix[2]) * imgdim[1] + im3;
          else
             im2 = imgdim[2] - (i2 + 1 + imgfpix[2]) * imgdim[1] + im3;

          t2 = (tilefpix[2] + it2) * tiledim[1] + t3;

          /* loop through rows of the image */
          for (i1 = 0, it1=0; i1 <= imglpix[1] - imgfpix[1]; i1++, it1++)
          {
             /* incre row if it falls in the cracks of the subsampled image */
             while (ndim > 1 &&  (tfpixel[1] + tilefpix[1] - fpixel[1] + it1)
                                  % labs(inc[1]) != 0)
                it1++;

             /* calc position of first pixel in tile to be copied */
             tilepix = tilefpix[0] + (tilefpix[1] + it1) * tiledim[0] + t2;

             /* offset to start of row */
             if (inc[1] > 0)
                im1 = (i1 + imgfpix[1]) * imgdim[0] + im2;
             else
                im1 = imgdim[1] - (i1 + 1 + imgfpix[1]) * imgdim[0] + im2;
/*
printf("inc = %d %d %d %d\n",inc[0],inc[1],inc[2],inc[3]);
printf("im1,im2,im3,im4 = %d %d %d %d\n",im1,im2,im3,im4);
*/
             /* offset to byte within the row */
             if (inc[0] > 0)
                imgpix = imgfpix[0] + im1;
             else
                imgpix = imgdim[0] - 1 - imgfpix[0] + im1;
/*
printf("tilefpix0,1, imgfpix1, it1, inc1, t2= %d %d %d %d %d %d\n",
       tilefpix[0],tilefpix[1],imgfpix[1],it1,inc[1], t2);
printf("i1, it1, tilepix, imgpix %d %d %d %d \n", i1, it1, tilepix, imgpix);
*/
             /* loop over pixels along one row of the image */
             for (ipos = imgfpix[0]; ipos <= imglpix[0]; ipos += overlap_flags)
             {
               /* convert from image pixel to byte offset */
               tilepixbyte = tilepix * pixlen;
               imgpixbyte  = imgpix  * pixlen;
/*
printf("  tilepix, tilepixbyte, imgpix, imgpixbyte= %d %d %d %d\n",
          tilepix, tilepixbyte, imgpix, imgpixbyte);
*/
               /* copy overlapping row of pixels from image to tile */
               memcpy(tile + tilepixbyte, image + imgpixbyte,  overlap_bytes);

               tilepix += (overlap_flags * labs(inc[0]));
               if (inc[0] > 0)
                 imgpix += overlap_flags;
               else
                 imgpix -= overlap_flags;
            }
          }
        }
      }
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
static int unquantize_i1r4(long row, /* tile number = row number in table  */
            unsigned char *input, /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
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
static int unquantize_i2r4(long row, /* seed for random values  */
            short *input,         /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
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
static int unquantize_i4r4(long row, /* tile number = row number in table    */
            INT32BIT *input,      /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            INT32BIT tnull,       /* I - value of FITS TNULLn keyword if any */
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
static int unquantize_i1r8(long row, /* tile number = row number in table  */
            unsigned char *input, /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
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
static int unquantize_i2r8(long row, /* tile number = row number in table  */
            short *input,         /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
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
static int unquantize_i4r8(long row, /* tile number = row number in table    */
            INT32BIT *input,      /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            INT32BIT tnull,       /* I - value of FITS TNULLn keyword if any */
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
static int imcomp_float2nan(float *indata, 
    long tilelen,
    int *outdata,
    float nullflagval, 
    int *status)
/*
  convert pixels that are equal to nullflag to NaNs.
  Note that indata and outdata point to the same location.
*/
{
    int ii;
    
    for (ii = 0; ii < tilelen; ii++) {

      if (indata[ii] == nullflagval)
        outdata[ii] = -1;  /* integer -1 has the same bit pattern as a real*4 NaN */
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
static int imcomp_double2nan(double *indata, 
    long tilelen,
    LONGLONG *outdata,
    double nullflagval, 
    int *status)
/*
  convert pixels that are equal to nullflag to NaNs.
  Note that indata and outdata point to the same location.
*/
{
    int ii;
    
    for (ii = 0; ii < tilelen; ii++) {

      if (indata[ii] == nullflagval)
        outdata[ii] = -1;  /* integer -1 has the same bit pattern as a real*8 NaN */
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_transpose_table(fitsfile *infptr, fitsfile *outfptr, int *status)

/*
  Transpose the elements in the input table columns from row-major order into
  column-major order, and write to the output table (which may be the same as
  the input table).   For example, a table with 10000 rows and 2 '1I' columns
  will be transformed into a 1 row table with 2 '10000I' columns.
  
  In addition, compress each column of data and write as a 1-row variable length
  array column.
*/
{ 
    LONGLONG nrows, incolwidth[999], inrepeat[999], outcolstart[1000], outbytespan[999];
    LONGLONG headstart, datastart, dataend, startbyte, jj, kk, naxis1;
    long repeat, width, pcount;
    int ii, ncols, coltype, hdutype, ltrue = 1;
    char *buffer, *cptr, keyname[9], tform[40], colcode[999], colname[999][50];
    char comm[FLEN_COMMENT], *compressed_data;
    size_t dlen, datasize;
    float cratio[999];
    
    if (*status > 0)
        return(*status);
    
    fits_get_hdu_type(infptr, &hdutype, status);
    if (hdutype != BINARY_TBL) {
        *status = NOT_BTABLE;
        return(*status);
    }
        
    fits_get_num_rowsll(infptr, &nrows, status);
    fits_get_num_cols(infptr, &ncols, status);
    fits_read_key(infptr, TLONGLONG, "NAXIS1", &naxis1, NULL, status);
    if (*status > 0)
        return(*status);

    if (nrows < 1  || ncols < 1) {
	/* just copy the HDU if the table has 0 columns or rows */
	if (infptr != outfptr) {  /* copy input header to the output */
		fits_copy_hdu (infptr, outfptr, 0, status);
	}
	return(*status);
    }
 
    /* allocate space for the transposed table */
    buffer = calloc((size_t) naxis1, (size_t) nrows);
    if (!buffer) {
        ffpmsg("Could not allocate buffer for transformed table");
        *status = MEMORY_ALLOCATION;
        return(*status);
    }

    if (infptr != outfptr) {  /* copy input header to the output */
	fits_copy_header(infptr, outfptr, status);
    }

    outcolstart[0] = 0;
 
    /* do initial setup for each column */
    for (ii = 0; ii < ncols; ii++) {

	/* get the column name */
	fits_make_keyn("TTYPE", ii+1, keyname, status);
	fits_read_key(outfptr, TSTRING, keyname, colname[ii], comm, status);

	/* get the column type, repeat count, and unit width */
	fits_make_keyn("TFORM", ii+1, keyname, status);
	fits_read_key(outfptr, TSTRING, keyname, tform, comm, status);

	/* preserve the original TFORM value and comment string */
	keyname[0] = 'Z';
	fits_write_key(outfptr, TSTRING, keyname, tform, comm, status);
	keyname[0] = 'T';
 
        fits_binary_tform(tform, &coltype, &repeat, &width, status);

   	/* BIT columns are a difficult case */
	 /* round up to a multiple of 8 bits */
/*
	if (coltype == TBIT) {  
	    repeat = (repeat + 7) / 8 * 8; 
	}
*/

	cptr = tform;
	while(isdigit(*cptr)) cptr++;
	colcode[ii] = *cptr; /* save the column type code */

        /* all columns are now VLAs */
	fits_modify_key_str(outfptr, keyname, "1PB", "&", status);

	if (coltype == TBIT) {
	    repeat = (repeat + 7) / 8;  /* convert from bits to bytes */
	} else if (coltype == TSTRING) {
	    width = 1;  /* ignore the optional 'w' in 'rAw' format */
	} else if (coltype < 0) {  /* pointer to variable length array */
	    width = 8;
	    if (colcode[ii] == 'Q') width = 16;  /* this is a 'Q' not a 'P' column */
	    repeat = 1;
	}

	inrepeat[ii] = repeat;
	
	/* width (in bytes) of each element and field in the INPUT row-major table */
	incolwidth[ii] = repeat * width;
	
	/* starting offset of each field in the OUTPUT column-major table */
	outcolstart[ii + 1] = outcolstart[ii] + incolwidth[ii] * nrows;

	/* length of each sequence of bytes, after sorting them in signicant order */
	outbytespan[ii] = (incolwidth[ii] * nrows) / width;
    }

    /* the transformed table has only 1 row */
    /* output table width 8 bytes per column */
    fits_modify_key_lng(outfptr, "NAXIS2", 1, "&", status);
    fits_modify_key_lng(outfptr, "NAXIS1", ncols * 8, "&", status);

    /* move to the start of the input table */
    fits_get_hduaddrll(infptr, &headstart, &datastart, &dataend, status);
    ffmbyt(infptr, datastart, 0, status);

    /* now transpose the table into an array in memory */
    for (jj = 0; jj < nrows; jj++)   {    /* loop over rows */
      for (ii = 0; ii < ncols; ii++) {  /* loop over columns */
      
	kk = 0;	

	    cptr = buffer + (outcolstart[ii] + (jj * incolwidth[ii]));   /* addr to copy to */

	    startbyte = (infptr->Fptr)->bytepos;  /* save the starting byte location */

	    ffgbyt(infptr, incolwidth[ii], cptr, status);  /* copy all the bytes */

	    if (incolwidth[ii] >= MINDIRECT) { /* have to explicitly move to next byte */
		ffmbyt(infptr, startbyte + incolwidth[ii], 0, status);
	    }
      }
    }

    fits_set_hdustruc(outfptr, status);
    
    /* now compress each column with GZIP and write out to output table */
    for (ii = 0; ii < ncols; ii++) {  /* loop over columns */

	datasize = (size_t) (outcolstart[ii + 1] - outcolstart[ii]);

	/* allocate memory for the compressed data */
	compressed_data = malloc(datasize);
	if (!compressed_data) {
            ffpmsg("data memory allocation error");
	    return(-1);
	}

	/* gzip compress the data */
	compress2mem_from_mem(buffer + outcolstart[ii], datasize,
	    &compressed_data,  &datasize, realloc, 
	    &dlen, status);        

	/* write the compressed data to the output column */
	fits_set_tscale(outfptr, ii + 1, 1.0, 0.0, status);  /* turn off any data scaling, first */
	fits_write_col(outfptr, TBYTE, ii + 1, 1, 1, dlen, compressed_data, status);

        cratio[ii] = (float) datasize / (float) dlen;	
	free(compressed_data);   /* don't need the compressed data any more */
	
	sprintf(results[ii]," %3d %10.10s %4d %c  %5.2f", ii+1, colname[ii], (int) inrepeat[ii],colcode[ii],cratio[ii]);
	trans_ratio[ii] = cratio[ii];
    }

    /* save the original PCOUNT value */
    fits_read_key(infptr, TLONG, "PCOUNT", &pcount, comm, status);
    fits_write_key(outfptr, TLONG, "ZPCOUNT", &pcount, comm, status);

    fits_write_key(outfptr, TLONGLONG, "ZNAXIS1", &naxis1, "original rows width",
	status);
    fits_write_key(outfptr, TLONGLONG, "ZNAXIS2", &nrows, "original number of rows",
	status);

    fits_write_key(outfptr, TLOGICAL, "TVIRTUAL", &ltrue, 
        "this is a virtual table", status);
    fits_write_key(outfptr, TSTRING, "ZMETHOD", "TRANSPOSED_SHUFFLED_GZIP",
	   "table compression method", status);

    fits_set_hdustruc(outfptr, status);

    /* copy the heap from input to output file */
    fits_gzip_heap(infptr, outfptr, status);
       	
    free(buffer);
    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_compress_table_rice(fitsfile *infptr, fitsfile *outfptr, int *status)

/*
  Transpose the elements in the input table columns from row-major order into
  column-major order, and write to the output table (which may be the same as
  the input table).   For example, a table with 10000 rows and 2 '1I' columns
  will be transformed into a 1 row table with 2 '10000I' columns.
  
  Integer columns are then compressed with Rice; all other columns compressed
  with GZIP.  In addition, the bytes in the floating point numeric data values 
  (columns with TFORM =  E, and D) are shuffled so that the most significant
  byte of every element occurs first in the array, followed by the next most
  significant byte, and so on to the least significant byte.   Thus, if you
  have 3 4-byte numeric values, the bytes 012301230123 get shuffled to
  000111222333
*/
{ 
    LONGLONG nrows, incolwidth[999], inrepeat[999], outcolstart[1000], outbytespan[999];
    LONGLONG headstart, datastart, dataend, startbyte, jj, kk, naxis1;
    long repeat, width, pcount;
    int ii, ncols, coltype, hdutype, ltrue = 1;
    char *buffer, *cptr, keyname[9], tform[40], colcode[999], tempstring[20];
    char comm[FLEN_COMMENT], *compressed_data;
    float cratio[999];
    
    size_t dlen, datasize;

    if (*status > 0)
        return(*status);
    
    fits_get_hdu_type(infptr, &hdutype, status);
    if (hdutype != BINARY_TBL) {
        *status = NOT_BTABLE;
        return(*status);
    }
        
    fits_get_num_rowsll(infptr, &nrows, status);
    fits_get_num_cols(infptr, &ncols, status);
    fits_read_key(infptr, TLONGLONG, "NAXIS1", &naxis1, NULL, status);
    if (*status > 0)
        return(*status);

    if (nrows < 1  || ncols < 1) {
	/* just copy the HDU if the table has 0 columns or rows */
	if (infptr != outfptr) {  /* copy input header to the output */
		fits_copy_hdu (infptr, outfptr, 0, status);
	}
	return(*status);
    }
   
    /* allocate space for the transposed table */
    buffer = calloc((size_t) naxis1, (size_t) nrows);
    if (!buffer) {
        ffpmsg("Could not allocate buffer for transformed table");
        *status = MEMORY_ALLOCATION;
        return(*status);
    }

    if (infptr != outfptr) {  /* copy input header to the output */
	fits_copy_header(infptr, outfptr, status);
    }

    outcolstart[0] = 0;
    for (ii = 0; ii < ncols; ii++) {

	/* get the column type, repeat count, and unit width */
	fits_make_keyn("TFORM", ii+1, keyname, status);
	fits_read_key(outfptr, TSTRING, keyname, tform, comm, status);

	/* preserve the original TFORM value and comment string */
	keyname[0] = 'Z';
	fits_write_key(outfptr, TSTRING, keyname, tform, comm, status);
	keyname[0] = 'T';
 
        fits_binary_tform(tform, &coltype, &repeat, &width, status);

   	/* BIT columns are a difficult case */
	 /* round up to a multiple of 8 bits */
/*
	if (coltype == TBIT) {  
	    repeat = (repeat + 7) / 8 * 8; 
	}
*/
	cptr = tform;
	while(isdigit(*cptr)) cptr++;
	colcode[ii] = *cptr; /* save the column type code */

/* all columns are now VLAs */
	fits_modify_key_str(outfptr, keyname, "1PB", "&", status);

	if (coltype == TBIT) {
	    repeat = (repeat + 7) / 8;  /* convert from bits to bytes */
	} else if (coltype == TSTRING) {
	    width = 1;  /* ignore the optional 'w' in 'rAw' format */
	} else if (coltype < 0) {  /* pointer to variable length array */
	    width = 8;
	    if (colcode[ii] == 'Q') width = 16;  /* this is a 'Q' not a 'P' column */
	    repeat = 1;
	}

	inrepeat[ii] = repeat;
	
	/* width (in bytes) of each element and field in the INPUT row-major table */
	incolwidth[ii] = repeat * width;
	
	/* starting offset of each field in the OUTPUT column-major table */
	outcolstart[ii + 1] = outcolstart[ii] + incolwidth[ii] * nrows;

	/* length of each sequence of bytes, after sorting them in signicant order */
	outbytespan[ii] = (incolwidth[ii] * nrows) / width;
    }

    /* the transformed table has only 1 row */
    /* output table width 8 bytes per column */

    fits_modify_key_lng(outfptr, "NAXIS2", 1, "&", status);
    fits_modify_key_lng(outfptr, "NAXIS1", ncols * 8, "&", status);

    /* move to the start of the input table */
    fits_get_hduaddrll(infptr, &headstart, &datastart, &dataend, status);
    ffmbyt(infptr, datastart, 0, status);

    for (jj = 0; jj < nrows; jj++)   {    /* loop over rows */
      for (ii = 0; ii < ncols; ii++) {  /* loop over columns */
      
	kk = 0;	

	switch (colcode[ii]) {
	/* separate the byte planes for the 2-byte, 4-byte, and 8-byte numeric columns */

	case 'E':
	  while(kk < incolwidth[ii]) {
	    cptr = buffer + (outcolstart[ii] + (jj * inrepeat[ii]) + kk/4);  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    cptr += outbytespan[ii];  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    cptr += outbytespan[ii];  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    cptr += outbytespan[ii];  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    kk += 4;
	  }
	  break;

	case 'D':
	case 'K':
	  while(kk < incolwidth[ii]) {
	    cptr = buffer + (outcolstart[ii] + (jj * inrepeat[ii]) + kk/8);  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    cptr += outbytespan[ii];  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    cptr += outbytespan[ii];  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    cptr += outbytespan[ii];  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    cptr += outbytespan[ii];  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    cptr += outbytespan[ii];  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    cptr += outbytespan[ii];  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    cptr += outbytespan[ii];  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    kk += 8;
	  }
	  break;

	default: /* don't bother separating the bytes for other column types */
	    cptr = buffer + (outcolstart[ii] + (jj * incolwidth[ii]));   /* addr to copy to */
	    startbyte = (infptr->Fptr)->bytepos;  /* save the starting byte location */

	    ffgbyt(infptr, incolwidth[ii], cptr, status);  /* copy all the bytes */

	    if (incolwidth[ii] >= MINDIRECT) { /* have to explicitly move to next byte */
		ffmbyt(infptr, startbyte + incolwidth[ii], 0, status);
	    }
	}
      }
    }

    fits_set_hdustruc(outfptr, status);
    
    for (ii = 0; ii < ncols; ii++) {  /* loop over columns */

	datasize = (size_t) (outcolstart[ii + 1] - outcolstart[ii]);
	/* allocate memory for the compressed data */
	compressed_data = malloc(datasize*2);
	if (!compressed_data) {
            ffpmsg("data memory allocation error");
	    return(-1);
	}


	switch (colcode[ii]) {


	case 'I':
#if BYTESWAPPED
                ffswap2((short *) (buffer + outcolstart[ii]),  datasize / 2); 
#endif
  	        dlen = fits_rcomp_short ((short *)(buffer + outcolstart[ii]), datasize / 2, (unsigned char *) compressed_data,
                       datasize * 2, 32);

	  	fits_make_keyn("ZCTYP", ii+1, keyname, status);
	  	fits_write_key(outfptr, TSTRING, keyname, "RICE_1",
	     	     "compression algorithm for column", status);

	  break;

	case 'J':

#if BYTESWAPPED
                ffswap4((int *) (buffer + outcolstart[ii]),  datasize / 4); 
#endif
   	        dlen = fits_rcomp ((int *)(buffer + outcolstart[ii]), datasize / 4, (unsigned char *) compressed_data,
                       datasize * 2, 32);

	  	fits_make_keyn("ZCTYP", ii+1, keyname, status);
	  	fits_write_key(outfptr, TSTRING, keyname, "RICE_1",
	     	     "compression algorithm for column", status);
	  break;

	case 'B':
  	        dlen = fits_rcomp_byte ((signed char *)(buffer + outcolstart[ii]), datasize, (unsigned char *) compressed_data,
                       datasize * 2, 32);

	  	fits_make_keyn("ZCTYP", ii+1, keyname, status);
	  	fits_write_key(outfptr, TSTRING, keyname, "RICE_1",
	     	     "compression algorithm for column", status);
	  break;

	default: 
		/* gzip compress the data */
		compress2mem_from_mem(buffer + outcolstart[ii], datasize,
	    		&compressed_data,  &datasize, realloc, &dlen, status);        


	  	fits_make_keyn("ZCTYP", ii+1, keyname, status);
	  	fits_write_key(outfptr, TSTRING, keyname, "GZIP_2",
	     	     "compression algorithm for column", status);

	} /* end of switch block */

	if (dlen != 0)
	cratio[ii] = (float) datasize / (float) dlen;  /* compression ratio of the column */

	/* write the compressed data to the output column */
	fits_set_tscale(outfptr, ii + 1, 1.0, 0.0, status);  /* turn off any data scaling, first */
	fits_write_col(outfptr, TBYTE, ii + 1, 1, 1, dlen, compressed_data, status);
	
	free(compressed_data);   /* don't need the compressed data any more */
/*	printf("   %c  %5.2f\n",colcode[ii],cratio[ii]); */

	    sprintf(tempstring,"  %5.2f\n",cratio[ii]);
/*
        if (colcode[ii] == 'I' || colcode[ii] == 'J' || colcode[ii] == 'B') 
	    sprintf(tempstring,"  %5.2f\n",cratio[ii]);
	else
	    sprintf(tempstring," \n");	
*/
	strcat(results[ii],tempstring);
    }  /* end of loop over ncols */

    printf("                       Trans   Shuf   Rice\n");
    for (ii = 0; ii < ncols; ii++) {  /* loop over columns */
      printf("%s", results[ii]);
    }
    
    /* save the original PCOUNT value */
    fits_read_key(infptr, TLONG, "PCOUNT", &pcount, comm, status);
    fits_write_key(outfptr, TLONG, "ZPCOUNT", &pcount, comm, status);

	
    fits_write_key(outfptr, TLONGLONG, "ZNAXIS1", &naxis1, "original rows width",
	status);
    fits_write_key(outfptr, TLONGLONG, "ZNAXIS2", &nrows, "original number of rows",
	status);

    fits_write_key(outfptr, TLOGICAL, "ZTABLE", &ltrue, 
        "this is a compressed table", status);

    free(buffer);

    fits_gzip_heap(infptr, outfptr, status);
    fits_set_hdustruc(outfptr, status);
       	
    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_compress_table_fast(fitsfile *infptr, fitsfile *outfptr, int *status)

/*
  Compress the input FITS binary table using the 'fast' method, which consists
  of (a) transposing the rows and columns, and (b) shuffling the bytes for
  the I, J, K, E, and D columns  so that the most significant byte of every
  element occurs first in the array, followed by the next most significant byte,
  and so on to the least significant byte.   Thus, if you have 3 4-byte numeric
  values, the bytes 012301230123 get shuffled to 000111222333

  Finally, (c) compress each column of bytes with gzip and copy to the output table.
  
*/
{ 
    LONGLONG nrows, incolwidth[999], inrepeat[999], outcolstart[1000], outbytespan[999];
    LONGLONG headstart, datastart, dataend, startbyte, jj, kk, naxis1;
    long repeat, width, pcount;
    int ii, ncols, coltype, hdutype, ltrue = 1;
    char *buffer, *cptr, keyname[9], tform[40], colcode[999];
    char comm[FLEN_COMMENT], *compressed_data, tempstring[20];
    size_t dlen, datasize;
    float cratio[999];
    
    if (*status > 0)
        return(*status);
    
    fits_get_hdu_type(infptr, &hdutype, status);
    if (hdutype != BINARY_TBL) {
        *status = NOT_BTABLE;
        return(*status);
    }
        
    fits_get_num_rowsll(infptr, &nrows, status);
    fits_get_num_cols(infptr, &ncols, status);
    fits_read_key(infptr, TLONGLONG, "NAXIS1", &naxis1, NULL, status);
    if (*status > 0)
        return(*status);

    if (nrows < 1  || ncols < 1) {
	/* just copy the HDU if the table has 0 columns or rows */
	if (infptr != outfptr) {  /* copy input header to the output */
		fits_copy_hdu (infptr, outfptr, 0, status);
	}
	return(*status);
    }
 
    /* allocate space for the transposed table */
    buffer = calloc((size_t) naxis1, (size_t) nrows);
    if (!buffer) {
        ffpmsg("Could not allocate buffer for transformed table");
        *status = MEMORY_ALLOCATION;
        return(*status);
    }

    if (infptr != outfptr) {  /* copy input header to the output */
	fits_copy_header(infptr, outfptr, status);
    }

    fits_write_key_log(outfptr, "ZTABLE", 1, 
          "extension contains compressed binary table", status);

    fits_write_key(outfptr, TLONGLONG, "ZTILELEN", &nrows,
          "number of rows in each tile", status);

    fits_write_key(outfptr, TLONGLONG, "ZNAXIS1", &naxis1, "original rows width",
	status);
    fits_write_key(outfptr, TLONGLONG, "ZNAXIS2", &nrows, "original number of rows",
	status);

    /* save the original PCOUNT value */
    fits_read_key(infptr, TLONG, "PCOUNT", &pcount, comm, status);
    fits_write_key(outfptr, TLONG, "ZPCOUNT", &pcount, comm, status);

    /* reset the PCOUNT keyword to zero */
    pcount = 0;
    fits_modify_key_lng(outfptr, "PCOUNT", pcount, NULL, status);

    outcolstart[0] = 0;
    for (ii = 0; ii < ncols; ii++) {

	/* get the column type, repeat count, and unit width */
	fits_make_keyn("TFORM", ii+1, keyname, status);
	fits_read_key(outfptr, TSTRING, keyname, tform, comm, status);

	/* preserve the original TFORM value and comment string */
	keyname[0] = 'Z';
	fits_write_key(outfptr, TSTRING, keyname, tform, comm, status);
	keyname[0] = 'T';
 
       /* all columns are now VLAs */
	fits_modify_key_str(outfptr, keyname, "1PB", "&", status);

 	fits_make_keyn("ZCTYP", ii+1, keyname, status);
	fits_write_key(outfptr, TSTRING, keyname, "GZIP_2",
	     "compression algorithm for column", status);

        fits_binary_tform(tform, &coltype, &repeat, &width, status);

	cptr = tform;
	while(isdigit(*cptr)) cptr++;
	colcode[ii] = *cptr; /* save the column type code */

	if (coltype == TBIT) {
	    repeat = (repeat + 7) / 8;  /* convert from bits to bytes */
	} else if (coltype == TSTRING) {
	    width = 1;  /* ignore the optional 'w' in 'rAw' format */
	} else if (coltype < 0) {  /* pointer to variable length array */
	    width = 8;
	    if (colcode[ii] == 'Q') width = 16;  /* this is a 'Q' not a 'P' column */
	    repeat = 1;
	}

	inrepeat[ii] = repeat;
	
	/* width (in bytes) of each element and field in the INPUT row-major table */
	incolwidth[ii] = repeat * width;
	
	/* starting offset of each field in the OUTPUT column-major table */
	outcolstart[ii + 1] = outcolstart[ii] + incolwidth[ii] * nrows;

	/* length of each sequence of bytes, after sorting them in signicant order */
	outbytespan[ii] = (incolwidth[ii] * nrows) / width;

    }

    /* the transformed table has only 1 row */
    /* output table width 8 bytes per column */

    fits_modify_key_lng(outfptr, "NAXIS2", 1, "&", status);
    fits_modify_key_lng(outfptr, "NAXIS1", ncols * 8, "&", status);

    /* move to the start of the input table */
    fits_get_hduaddrll(infptr, &headstart, &datastart, &dataend, status);
    ffmbyt(infptr, datastart, 0, status);

    for (jj = 0; jj < nrows; jj++)   {    /* loop over rows */
      for (ii = 0; ii < ncols; ii++) {  /* loop over columns */
      
	kk = 0;	

	switch (colcode[ii]) {
	/* separate the byte planes for the 2-byte, 4-byte, and 8-byte numeric columns */
	case 'I':
	  while(kk < incolwidth[ii]) {

	    cptr = buffer + (outcolstart[ii] + (jj * inrepeat[ii]) + kk/2);  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1st byte */
	    cptr += outbytespan[ii];  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 2nd byte */
	    kk += 2;
	  }

	  break;
	
	case 'J':
	case 'E':
	  while(kk < incolwidth[ii]) {
	    cptr = buffer + (outcolstart[ii] + (jj * inrepeat[ii]) + kk/4);  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    cptr += outbytespan[ii];  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    cptr += outbytespan[ii];  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    cptr += outbytespan[ii];  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    kk += 4;
	  }
	  break;

	case 'D':
	case 'K':
	  while(kk < incolwidth[ii]) {
	    cptr = buffer + (outcolstart[ii] + (jj * inrepeat[ii]) + kk/8);  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    cptr += outbytespan[ii];  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    cptr += outbytespan[ii];  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    cptr += outbytespan[ii];  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    cptr += outbytespan[ii];  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    cptr += outbytespan[ii];  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    cptr += outbytespan[ii];  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    cptr += outbytespan[ii];  
	    ffgbyt(infptr, 1, cptr, status);  /* copy 1 byte */
	    kk += 8;
	  }
	  break;

	default: /* don't bother separating the bytes for other column types */
	    cptr = buffer + (outcolstart[ii] + (jj * incolwidth[ii]));   /* addr to copy to */

	    startbyte = (infptr->Fptr)->bytepos;  /* save the starting byte location */

	    ffgbyt(infptr, incolwidth[ii], cptr, status);  /* copy all the bytes */

	    if (incolwidth[ii] >= MINDIRECT) { /* have to explicitly move to next byte */
		ffmbyt(infptr, startbyte + incolwidth[ii], 0, status);
	    }
	}
      }
    }

    fits_set_hdustruc(outfptr, status);
    
    /* now compress each column */
    for (ii = 0; ii < ncols; ii++) {  /* loop over columns */

	/* write the compression type code for this column */
	switch (colcode[ii]) {
	case 'I':
	case 'J':
	case 'K':
	case 'E':
	case 'D':
	  fits_make_keyn("ZCTYP", ii+1, keyname, status);
	  fits_write_key(outfptr, TSTRING, keyname, "GZIP_2",
	     "compression algorithm for column", status);
	  break;
	default:
	  fits_make_keyn("ZCTYP", ii+1, keyname, status);
	  fits_write_key(outfptr, TSTRING, keyname, "GZIP_1",
	     "compression algorithm for column", status);
        }

	datasize = (size_t) (outcolstart[ii + 1] - outcolstart[ii]);

	/* allocate memory for the compressed data */
	compressed_data = malloc(datasize);
	if (!compressed_data) {
            ffpmsg("data memory allocation error");
	    return(-1);
	}

	/* gzip compress the data */
	compress2mem_from_mem(buffer + outcolstart[ii], datasize,
	    &compressed_data,  &datasize, realloc, 
	    &dlen, status);        

	/* write the compressed data to the output column */
	fits_set_tscale(outfptr, ii + 1, 1.0, 0.0, status);  /* turn off any data scaling, first */
	fits_write_col(outfptr, TBYTE, ii + 1, 1, 1, dlen, compressed_data, status);

        cratio[ii] = (float) datasize / (float) dlen;	
	free(compressed_data);   /* don't need the compressed data any more */

	sprintf(tempstring,"  %5.2f",cratio[ii]);

	strcat(results[ii],tempstring);
    }

    free(buffer);

    /* shuffle and compress the input heap and append to the output file */

    fits_gzip_heap(infptr, outfptr, status);
    fits_set_hdustruc(outfptr, status);

    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_compress_table_best(fitsfile *infptr, fitsfile *outfptr, int *status)

/*
  Compress the input FITS binary table using the 'best' compression method, i.e,
  whichever method produces the highest compression for each column.
  
  First, transpose the rows and columns in the table, then, depending on the 
  data type of the column, try the different compression methods to see
  which one produces the highest amount of compression.
  
*/
{ 
    LONGLONG nrows, incolwidth[999], inrepeat[999], outcolstart[1000], outbytespan[999];
    LONGLONG headstart, datastart, dataend, startbyte, jj, naxis1;
    long repeat, width, pcount;
    int ii, ncols, coltype, hdutype, ltrue = 1;
    char *buffer, *cptr, keyname[9], tform[40], colcode[999];
    char comm[FLEN_COMMENT];
    char *gzip1_data = 0, *gzip2_data = 0, *rice_data = 0;
    size_t gzip1_len, gzip2_len, rice_len, datasize, buffsize;
    
    if (*status > 0)
        return(*status);
    
    fits_get_hdu_type(infptr, &hdutype, status);
    if (hdutype != BINARY_TBL) {
        *status = NOT_BTABLE;
        return(*status);
    }
        
    fits_get_num_rowsll(infptr, &nrows, status);
    fits_get_num_cols(infptr, &ncols, status);
    fits_read_key(infptr, TLONGLONG, "NAXIS1", &naxis1, NULL, status);
    if (*status > 0)
        return(*status);

    if (nrows < 1  || ncols < 1) {
	/* just copy the HDU if the table has 0 columns or rows */
	if (infptr != outfptr) {  /* copy input header to the output */
		fits_copy_hdu (infptr, outfptr, 0, status);
	}
	return(*status);
    }
 
    /* allocate space for the transposed table */
    buffer = calloc((size_t) naxis1, (size_t) nrows);
    if (!buffer) {
        ffpmsg("Could not allocate buffer for transformed table");
        *status = MEMORY_ALLOCATION;
        return(*status);
    }

    if (infptr != outfptr) {  /* copy input header to the output */
	fits_copy_header(infptr, outfptr, status);
    }

    fits_write_key_log(outfptr, "ZTABLE", 1, 
          "extension contains compressed binary table", status);

    fits_write_key(outfptr, TLONGLONG, "ZTILELEN", &nrows,
          "number of rows in each tile", status);

    fits_write_key(outfptr, TLONGLONG, "ZNAXIS1", &naxis1, "original rows width",
	status);
    fits_write_key(outfptr, TLONGLONG, "ZNAXIS2", &nrows, "original number of rows",
	status);

    /* save the original PCOUNT value */
    fits_read_key(infptr, TLONG, "PCOUNT", &pcount, comm, status);
    fits_write_key(outfptr, TLONG, "ZPCOUNT", &pcount, comm, status);
    /* reset the PCOUNT keyword to zero */
    pcount = 0;
    fits_modify_key_lng(outfptr, "PCOUNT", pcount, NULL, status);

    /* Modify the TFORMn keywords; all columns become variable-length arrays. */
    /* Save the original TFORMn values in the corresponding ZFORMn keyword.   */
    outcolstart[0] = 0;
    for (ii = 0; ii < ncols; ii++) {

	/* get the column type, repeat count, and unit width */
	fits_make_keyn("TFORM", ii+1, keyname, status);
	fits_read_key(outfptr, TSTRING, keyname, tform, comm, status);

	/* preserve the original TFORM value and comment string */
	keyname[0] = 'Z';
	fits_write_key(outfptr, TSTRING, keyname, tform, comm, status);
	keyname[0] = 'T';
 
       /* all columns are now VLAs */
	fits_modify_key_str(outfptr, keyname, "1PB", "&", status);

        fits_binary_tform(tform, &coltype, &repeat, &width, status);

	cptr = tform;
	while(isdigit(*cptr)) cptr++;
	colcode[ii] = *cptr; /* save the column type code */

        /* deal with special cases */
	if (coltype == TBIT) {
	    repeat = (repeat + 7) / 8;  /* convert from bits to bytes */
	} else if (coltype == TSTRING) {
	    width = 1;  /* ignore the optional 'w' in 'rAw' format */
	} else if (coltype < 0) {  /* pointer to variable length array */
	    repeat = 1;

	    if (colcode[ii] == 'Q')
	       width = 16;  /* this is a 'Q' column */
	    else
	       width = 8;  /* this is a 'P' column */
	}

	inrepeat[ii] = repeat;
	
	/* width (in bytes) of each element and field in the INPUT row-major table */
	incolwidth[ii] = repeat * width;
	
	/* starting offset of each field in the OUTPUT column-major table */
	outcolstart[ii + 1] = outcolstart[ii] + incolwidth[ii] * nrows;

	/* length of each sequence of bytes, after sorting them in signicant order */
	outbytespan[ii] = (incolwidth[ii] * nrows) / width;
    }

    /* the transformed table has only 1 row */
    /* output table width 8 bytes per column */

    fits_modify_key_lng(outfptr, "NAXIS2", 1, "&", status);
    fits_modify_key_lng(outfptr, "NAXIS1", ncols * 8, "&", status);

    /* move to the start of the input table */
    fits_get_hduaddrll(infptr, &headstart, &datastart, &dataend, status);
    ffmbyt(infptr, datastart, 0, status);

    /* now transpose the rows and columns in the table into an array in memory */
    for (jj = 0; jj < nrows; jj++)   {    /* loop over rows */
        for (ii = 0; ii < ncols; ii++) {  /* loop over columns */
      
	    cptr = buffer + (outcolstart[ii] + (jj * incolwidth[ii])); /* output address */
	    startbyte = (infptr->Fptr)->bytepos;  /* save the starting byte location */
	    ffgbyt(infptr, incolwidth[ii], cptr, status);  /* copy the column element */

	    if (incolwidth[ii] >= MINDIRECT) { /* have to explicitly move to next byte */
		ffmbyt(infptr, startbyte + incolwidth[ii], 0, status);
	    }
        }
    }

    fits_set_hdustruc(outfptr, status);  /* reinitialize internal pointers */
    
    /* Now compress each column.  Depending on the column data type, try */
    /* all the various available compression algorithms, then choose the one */
    /* that gives the most compression. */
    for (ii = 0; ii < ncols; ii++) {  /* loop over columns */

	datasize = (size_t) (outcolstart[ii + 1] - outcolstart[ii]);

	/* allocate memory for the gzip compressed data */
	gzip1_data = malloc(datasize);
	if (!gzip1_data) {
            ffpmsg("data memory allocation error");
	    return(-1);
	}
	buffsize = datasize;

	/*  First, simply compress the bytes with gzip (GZIP_1 algorithm code). */
	/*  This algorithm can be applied to every type of column. */
	compress2mem_from_mem(buffer + outcolstart[ii], datasize,
	    &gzip1_data,  &buffsize, realloc, &gzip1_len, status);        
	
	/* depending on the data type, try other compression methods */
	switch (colcode[ii]) {

	case 'I':   /* 2-byte Integer columns */

		/************* first, try rice compression *****************/
		rice_data = malloc(datasize * 2);  /* memory for the compressed bytes */
		if (!rice_data) {
                    ffpmsg("data memory allocation error");
		    return(-1);
		}

#if BYTESWAPPED
		/* have to swap the bytes on little endian machines */
                ffswap2((short *) (buffer + outcolstart[ii]),  datasize / 2); 
#endif
  	        rice_len = fits_rcomp_short ((short *)(buffer + outcolstart[ii]), datasize / 2, 
		   (unsigned char *) rice_data, datasize * 2, 32);
	
#if BYTESWAPPED
		/* un-swap the bytes, to restore the original order */
                ffswap2((short *) (buffer + outcolstart[ii]),  datasize / 2); 
#endif
	   
		/************* Second, try shuffled gzip compression *****************/
		fits_shuffle_2bytes(buffer + outcolstart[ii], datasize / 2, status);

		/* allocate memory for the shuffled gzip compressed data */
		gzip2_data = malloc(datasize);
		if (!gzip2_data) {
                    ffpmsg("data memory allocation error");
		    return(-1);
		}
		buffsize = datasize;

		compress2mem_from_mem(buffer + outcolstart[ii], datasize,
		    &gzip2_data,  &buffsize, realloc, &gzip2_len, status);        
           break;

	case 'J':   /* 4-byte Integer columns */

		/************* first, try rice compression *****************/
		rice_data = malloc(datasize * 2);  /* memory for the compressed bytes */
		if (!rice_data) {
                    ffpmsg("data memory allocation error");
		    return(-1);
		}
#if BYTESWAPPED
		/* have to swap the bytes on little endian machines */
                ffswap4((int *) (buffer + outcolstart[ii]),  datasize / 4); 
#endif
  	        rice_len = fits_rcomp ((int *)(buffer + outcolstart[ii]), datasize / 4, 
		   (unsigned char *) rice_data, datasize * 2, 32);
	
#if BYTESWAPPED
		/* un-swap the bytes, to restore the original order */
                ffswap4((int *) (buffer + outcolstart[ii]),  datasize / 4); 
#endif
	   
		/************* Second, try shuffled gzip compression *****************/
		fits_shuffle_4bytes(buffer + outcolstart[ii], datasize / 4, status);

		/* allocate memory for the shuffled gzip compressed data */
		gzip2_data = malloc(datasize);
		if (!gzip2_data) {
                    ffpmsg("data memory allocation error");
		    return(-1);
		}
		buffsize = datasize;

		compress2mem_from_mem(buffer + outcolstart[ii], datasize,
		    &gzip2_data,  &buffsize, realloc, &gzip2_len, status);        
           break;

	case 'E':   /* 4-byte floating-point */

		/************* try shuffled gzip compression *****************/
		fits_shuffle_4bytes(buffer + outcolstart[ii], datasize / 4, status);

		/* allocate memory for the gzip compressed data */
		gzip2_data = malloc(datasize);
		if (!gzip2_data) {
                    ffpmsg("data memory allocation error");
		    return(-1);
		}
		buffsize = datasize;
		
		compress2mem_from_mem(buffer + outcolstart[ii], datasize,
		    &gzip2_data,  &buffsize, realloc, &gzip2_len, status);        

		rice_len = 100 * datasize;  /* rice is not applicable to R*4 data */

	   break;

	case 'K':
	case 'D':  /* 8-byte floating-point or integers */

		/************* try shuffled gzip compression *****************/
		fits_shuffle_8bytes(buffer + outcolstart[ii], datasize / 8, status);

		/* allocate memory for the gzip compressed data */
		gzip2_data = malloc(datasize);
		if (!gzip2_data) {
                    ffpmsg("data memory allocation error");
		    return(-1);
		}
		buffsize = datasize;

		compress2mem_from_mem(buffer + outcolstart[ii], datasize,
		    &gzip2_data,  &buffsize, realloc, &gzip2_len, status);        

		rice_len = 100 * datasize;  /* rice is not applicable to R*8 or I*8 data */

	   break;

	default:  /* L, X, B, A, C, M, P, Q type columns: no other compression options */
		rice_len = 100 * datasize;   /* rice is not applicable */
		gzip2_len = 100 * datasize;  /* shuffled-gzip is not applicable */

	}  /* end of switch block */

	/* now write the compressed bytes from the best algorithm */
	fits_set_tscale(outfptr, ii + 1, 1.0, 0.0, status);  /* turn off any data scaling, first */
	if (gzip1_len <= gzip2_len && gzip1_len <= rice_len) {

	    fits_write_col(outfptr, TBYTE, ii + 1, 1, 1, gzip1_len, gzip1_data, status);
	    fits_make_keyn("ZCTYP", ii+1, keyname, status);
	    fits_write_key(outfptr, TSTRING, keyname, "GZIP_1",
	         "compression algorithm for column", status);
	} else if (gzip2_len <= gzip1_len && gzip2_len <= rice_len) {
	    fits_write_col(outfptr, TBYTE, ii + 1, 1, 1, gzip2_len, gzip2_data, status);
	    fits_make_keyn("ZCTYP", ii+1, keyname, status);
	    fits_write_key(outfptr, TSTRING, keyname, "GZIP_2",
	         "compression algorithm for column", status);
	} else {
	    fits_write_col(outfptr, TBYTE, ii + 1, 1, 1, rice_len, rice_data, status);
	    fits_make_keyn("ZCTYP", ii+1, keyname, status);
	    fits_write_key(outfptr, TSTRING, keyname, "RICE_1",
	         "compression algorithm for column", status);
	}

	/* free the temporary memory */
	if (gzip1_data) free(gzip1_data);   
	if (gzip2_data) free(gzip2_data);   
	gzip1_data = 0;
	gzip2_data = 0;
    }

    free(buffer);

    /* shuffle and compress the input heap and append to the output file */

    fits_gzip_heap(infptr, outfptr, status);
    fits_set_hdustruc(outfptr, status);

    return(*status);
}
/*--------------------------------------------------------------------------*/
int fits_uncompress_table(fitsfile *infptr, fitsfile *outfptr, int *status)

/*
  Uncompress the table that was compressed with fits_compress_table_fast or
  fits_compress_table_best.
*/
{ 
    LONGLONG nrows, rmajor_colwidth[999], rmajor_colstart[1000], cmajor_colstart[1000];
    LONGLONG cmajor_repeat[999], rmajor_repeat[999], cmajor_bytespan[999], kk;
    LONGLONG headstart, datastart, dataend;
    long repeat, width, vla_repeat;
    int  ncols, coltype, hdutype, anynull, tstatus, zctype[999];
    char *buffer, *transbuffer, *cptr, keyname[9], tform[40], colcode[999];
    long  pcount, zheapptr, naxis1, naxis2, ii, jj;
    char *ptr, comm[FLEN_COMMENT], zvalue[FLEN_VALUE];
    size_t dlen, fullsize;
    
    /**** do initial sanity checks *****/
    if (*status > 0)
        return(*status);
     
    fits_get_hdu_type(infptr, &hdutype, status);
    if (hdutype != BINARY_TBL) {
        *status = NOT_BTABLE;
        return(*status);
    }
 
    fits_get_num_rowsll(infptr, &nrows, status);
    fits_get_num_cols(infptr, &ncols, status);

    if (nrows != 1  || (ncols < 1)) {
	/* just copy the HDU if the table does not have 1 row and 
	   more than 0 columns */
	if (infptr != outfptr) { 
		fits_copy_hdu (infptr, outfptr, 0, status);
	}
	return(*status);
    }

    /**** get size of the uncompressed table */
    fits_read_key(infptr, TLONG, "ZNAXIS1", &naxis1, comm, status);
    if (*status > 0) {
        ffpmsg("Could not find the required ZNAXIS1 keyword");
        *status = 1;
        return(*status);
    }

    fits_read_key(infptr, TLONG, "ZNAXIS2", &naxis2, comm, status);
    if (*status > 0) {
        ffpmsg("Could not find the required ZNAXIS2 keyword");
        *status = 1;
        return(*status);
    }

    fits_read_key(infptr, TLONG, "ZPCOUNT", &pcount, comm, status);
    if (*status > 0) {
        ffpmsg("Could not find the required ZPCOUNT keyword");
        *status = 1;
        return(*status);
    }

    tstatus = 0;
    fits_read_key(infptr, TLONG, "ZHEAPPTR", &zheapptr, comm, &tstatus);
    if (tstatus > 0) {
        zheapptr = 0;  /* uncompressed table has no heap */
    }

    /**** recreate the uncompressed table header keywords ****/
    fits_copy_header(infptr, outfptr, status);

    /* reset the NAXISn keywords to what they were in the original uncompressed table */
    fits_modify_key_lng(outfptr, "NAXIS1", naxis1, "&", status);
    fits_modify_key_lng(outfptr, "NAXIS2", naxis2, "&", status);
    fits_modify_key_lng(outfptr, "PCOUNT", pcount, "&", status);

    fits_delete_key(outfptr, "ZTABLE", status);
    fits_delete_key(outfptr, "ZNAXIS1", status);
    fits_delete_key(outfptr, "ZNAXIS2", status);
    fits_delete_key(outfptr, "ZPCOUNT", status);
    fits_delete_key(outfptr, "ZTILELEN", status);
    tstatus = 0;
    fits_delete_key(outfptr, "ZHEAPPTR", &tstatus);

    /**** get the compression method that was used for each column ****/
    for (ii = 0; ii < ncols; ii++) {

	/* construct the ZCTYPn keyword name then read the keyword */
	fits_make_keyn("ZCTYP", ii+1, keyname, status);
	tstatus = 0;
        fits_read_key(infptr, TSTRING, keyname, zvalue, NULL, &tstatus);
	if (tstatus) {
           zctype[ii] = GZIP_2;
	} else {
	   if (!strcmp(zvalue, "GZIP_2")) {
               zctype[ii] = GZIP_2;
	   } else if (!strcmp(zvalue, "GZIP_1")) {
               zctype[ii] = GZIP_1;
	   } else if (!strcmp(zvalue, "RICE_1")) {
               zctype[ii] = RICE_1;
	   } else {
	       ffpmsg("Unrecognized ZCTYPn keyword compression code:");
	       ffpmsg(zvalue);
	       *status = DATA_DECOMPRESSION_ERR;
	       return(*status);
	   }
	   
	   /* delete this keyword from the uncompressed header */
	   fits_delete_key(outfptr, keyname, status);
	}
    }

    /**** allocate space for the full transposed and untransposed table ****/
    fullsize = naxis1 * naxis2;
    transbuffer = malloc(fullsize);
    if (!transbuffer) {
        ffpmsg("Could not allocate buffer for shuffled table");
        *status = MEMORY_ALLOCATION;
        return(*status);
    }

    buffer = malloc(fullsize);
    if (!buffer) {
        ffpmsg("Could not allocate buffer for unshuffled table");
        *status = MEMORY_ALLOCATION;
        return(*status);
    }

    /*** loop over each column: read and uncompress the bytes ****/
    rmajor_colstart[0] = 0;
    cmajor_colstart[0] = 0;
    for (ii = 0; ii < ncols; ii++) {

	/* get the original column type, repeat count, and unit width */
	fits_make_keyn("ZFORM", ii+1, keyname, status);
	fits_read_key(infptr, TSTRING, keyname, tform, comm, status);

	/* restore the original TFORM value and comment */
        keyname[0] = 'T';
	fits_modify_key_str(outfptr, keyname, tform, comm, status);

	/* now delete the ZFORM keyword */
        keyname[0] = 'Z';
	fits_delete_key(outfptr, keyname, status);

	cptr = tform;
	while(isdigit(*cptr)) cptr++;
	colcode[ii] = *cptr; /* save the column type code */

        fits_binary_tform(tform, &coltype, &repeat, &width, status);

	/* deal with special cases */
	if (coltype == TBIT) { 
	    repeat = (repeat + 7) / 8 ;   /* convert from bits to bytes */
	} else if (coltype == TSTRING) {
	    width = 1;
	} else if (coltype < 0) {  /* pointer to variable length array */
	    if (colcode[ii] == 'P')
	       width = 8;  /* this is a 'P' column */
	    else
	       width = 16;  /* this is a 'Q' not a 'P' column */
	}

	rmajor_repeat[ii] = repeat;
	cmajor_repeat[ii] = repeat * naxis2;

	/* width (in bytes) of each field in the row-major table */
	rmajor_colwidth[ii] = rmajor_repeat[ii] * width;

	/* starting offset of each field in the column-major table */
	cmajor_colstart[ii + 1] = cmajor_colstart[ii] + rmajor_colwidth[ii] * naxis2;

	/* length of each sequence of bytes, after sorting them in signicant order */
	cmajor_bytespan[ii] = (rmajor_colwidth[ii] * naxis2) / width;

	/* starting offset of each field in the  row-major table */
	rmajor_colstart[ii + 1] = rmajor_colstart[ii] + rmajor_colwidth[ii];

	/* read compressed bytes from input table */
	fits_read_descript(infptr, ii + 1, 1, &vla_repeat, NULL, status);
	
	/* allocate memory and read in the compressed bytes */
	ptr = malloc(vla_repeat);
	if (!ptr) {
            ffpmsg("Could not allocate buffer for compressed bytes");
            *status = MEMORY_ALLOCATION;
            return(*status);
	}

	fits_set_tscale(infptr, ii + 1, 1.0, 0.0, status);  /* turn off any data scaling, first */
	fits_read_col_byt(infptr, ii + 1, 1, 1, vla_repeat, 0, (unsigned char *) ptr, &anynull, status);

        cptr = transbuffer + cmajor_colstart[ii];
	
	fullsize = (size_t) (cmajor_colstart[ii+1] - cmajor_colstart[ii]);

	switch (colcode[ii]) {
	/* separate the byte planes for the 2-byte, 4-byte, and 8-byte numeric columns */


	case 'I':

	    if (zctype[ii] == RICE_1) {
   	        dlen = fits_rdecomp_short((unsigned char *)ptr, vla_repeat, (unsigned short *)cptr, 
		       fullsize / 2, 32);
#if BYTESWAPPED
                ffswap2((short *) cptr, fullsize / 2); 
#endif
	    } else { /* gunzip the data into the correct location */
	        uncompress2mem_from_mem(ptr, vla_repeat, &cptr, &fullsize, realloc, &dlen, status);        
	    }
	  break;

	case 'J':

	    if (zctype[ii] == RICE_1) {
   	        dlen = fits_rdecomp ((unsigned char *) ptr, vla_repeat, (unsigned int *)cptr, 
		     fullsize / 4, 32);
#if BYTESWAPPED
                ffswap4((int *) cptr,  fullsize / 4); 
#endif
	    } else { /* gunzip the data into the correct location */
	        uncompress2mem_from_mem(ptr, vla_repeat, &cptr, &fullsize, realloc, &dlen, status);        
	    }
	  break;

	case 'B':

	    if (zctype[ii] == RICE_1) {
   	        dlen = fits_rdecomp_byte ((unsigned char *) ptr, vla_repeat, (unsigned char *)cptr, 
		     fullsize, 32);
	    } else { /* gunzip the data into the correct location */
	        uncompress2mem_from_mem(ptr, vla_repeat, &cptr, &fullsize, realloc, &dlen, status);        
	    }
	  break;

	default: 
	    /* gunzip the data into the correct location in the full table buffer */
	    uncompress2mem_from_mem(ptr, vla_repeat,
	        &cptr,  &fullsize, realloc, &dlen, status);              

	} /* end of switch block */

	free(ptr);
    }

    /* now transpose the rows and columns (from transbuffer to buffer) */
    ptr = transbuffer;
    for (ii = 0; ii < ncols; ii++) {  /* loop over columns */

      if ((zctype[ii] == GZIP_2)) {  /*  need to unshuffle the bytes */

	switch (colcode[ii]) {
	
	/* recombine the byte planes for the 2-byte, 4-byte, and 8-byte numeric columns */

	case 'I':
	  /* get the 1st byte of each I*2 value */
          for (jj = 0; jj < naxis2; jj++) {  /* loop over number of rows in the output table */
	    cptr = buffer + (rmajor_colstart[ii] + (jj * rmajor_colstart[ncols]));  
	    for (kk = 0; kk < rmajor_repeat[ii]; kk++) {
	      *cptr = *ptr;  /* copy 1 byte */
	      ptr++;
	      cptr += 2;  
	    }
	  }

	  /* get the 2nd byte of each I*2 value */
          for (jj = 0; jj < naxis2; jj++) {  /* loop over number of rows in the output table */
	    cptr = buffer + (rmajor_colstart[ii] + (jj * rmajor_colstart[ncols]) + 1);  
	    for (kk = 0; kk < rmajor_repeat[ii]; kk++) {
	      *cptr = *ptr;  /* copy 1 byte */
	      ptr++;
	      cptr += 2;  
	    }
	  }

	  break;

	case 'J':
	case 'E':

	  /* get the 1st byte of each 4-byte value */
          for (jj = 0; jj < naxis2; jj++) {  /* loop over number of rows in the output table */
	    cptr = buffer + (rmajor_colstart[ii] + (jj * rmajor_colstart[ncols]));  
	    for (kk = 0; kk < rmajor_repeat[ii]; kk++) {
	      *cptr = *ptr;  /* copy 1 byte */
	      ptr++;
	      cptr += 4;  
	    }
	  }

	  /* get the 2nd byte  */
          for (jj = 0; jj < naxis2; jj++) {  /* loop over number of rows in the output table */
	    cptr = buffer + (rmajor_colstart[ii] + (jj * rmajor_colstart[ncols]) + 1);  
	    for (kk = 0; kk < rmajor_repeat[ii]; kk++) {
	      *cptr = *ptr;  /* copy 1 byte */
	      ptr++;
	      cptr += 4;  
	    }
	  }

	  /* get the 3rd byte  */
          for (jj = 0; jj < naxis2; jj++) {  /* loop over number of rows in the output table */
	    cptr = buffer + (rmajor_colstart[ii] + (jj * rmajor_colstart[ncols]) + 2);  
	    for (kk = 0; kk < rmajor_repeat[ii]; kk++) {
	      *cptr = *ptr;  /* copy 1 byte */
	      ptr++;
	      cptr += 4;  
	    }
	  }

	  /* get the 4th byte  */
          for (jj = 0; jj < naxis2; jj++) {  /* loop over number of rows in the output table */
	    cptr = buffer + (rmajor_colstart[ii] + (jj * rmajor_colstart[ncols]) + 3);  
	    for (kk = 0; kk < rmajor_repeat[ii]; kk++) {
	      *cptr = *ptr;  /* copy 1 byte */
	      ptr++;
	      cptr += 4;  
	    }
	  }

	  break;

	case 'D':
	case 'K':

	  /* get the 1st byte of each 8-byte value */
          for (jj = 0; jj < naxis2; jj++) {  /* loop over number of rows in the output table */
	    cptr = buffer + (rmajor_colstart[ii] + (jj * rmajor_colstart[ncols]));  
	    for (kk = 0; kk < rmajor_repeat[ii]; kk++) {
	      *cptr = *ptr;  /* copy 1 byte */
	      ptr++;
	      cptr += 8;  
	    }
	  }

	  /* get the 2nd byte  */
          for (jj = 0; jj < naxis2; jj++) {  /* loop over number of rows in the output table */
	    cptr = buffer + (rmajor_colstart[ii] + (jj * rmajor_colstart[ncols]) + 1);  
	    for (kk = 0; kk < rmajor_repeat[ii]; kk++) {
	      *cptr = *ptr;  /* copy 1 byte */
	      ptr++;
	      cptr += 8;  
	    }
	  }

	  /* get the 3rd byte  */
          for (jj = 0; jj < naxis2; jj++) {  /* loop over number of rows in the output table */
	    cptr = buffer + (rmajor_colstart[ii] + (jj * rmajor_colstart[ncols]) + 2);  
	    for (kk = 0; kk < rmajor_repeat[ii]; kk++) {
	      *cptr = *ptr;  /* copy 1 byte */
	      ptr++;
	      cptr += 8;  
	    }
	  }

	  /* get the 4th byte  */
          for (jj = 0; jj < naxis2; jj++) {  /* loop over number of rows in the output table */
	    cptr = buffer + (rmajor_colstart[ii] + (jj * rmajor_colstart[ncols]) + 3);  
	    for (kk = 0; kk < rmajor_repeat[ii]; kk++) {
	      *cptr = *ptr;  /* copy 1 byte */
	      ptr++;
	      cptr += 8;  
	    }
	  }

	  /* get the 5th byte */
          for (jj = 0; jj < naxis2; jj++) {  /* loop over number of rows in the output table */
	    cptr = buffer + (rmajor_colstart[ii] + (jj * rmajor_colstart[ncols]) + 4);  
	    for (kk = 0; kk < rmajor_repeat[ii]; kk++) {
	      *cptr = *ptr;  /* copy 1 byte */
	      ptr++;
	      cptr += 8;  
	    }
	  }

	  /* get the 6th byte  */
          for (jj = 0; jj < naxis2; jj++) {  /* loop over number of rows in the output table */
	    cptr = buffer + (rmajor_colstart[ii] + (jj * rmajor_colstart[ncols]) + 5);  
	    for (kk = 0; kk < rmajor_repeat[ii]; kk++) {
	      *cptr = *ptr;  /* copy 1 byte */
	      ptr++;
	      cptr += 8;  
	    }
	  }

	  /* get the 7th byte  */
          for (jj = 0; jj < naxis2; jj++) {  /* loop over number of rows in the output table */
	    cptr = buffer + (rmajor_colstart[ii] + (jj * rmajor_colstart[ncols]) + 6);  
	    for (kk = 0; kk < rmajor_repeat[ii]; kk++) {
	      *cptr = *ptr;  /* copy 1 byte */
	      ptr++;
	      cptr += 8;  
	    }
	  }

	  /* get the 8th byte  */
          for (jj = 0; jj < naxis2; jj++) {  /* loop over number of rows in the output table */
	    cptr = buffer + (rmajor_colstart[ii] + (jj * rmajor_colstart[ncols]) + 7);  
	    for (kk = 0; kk < rmajor_repeat[ii]; kk++) {
	      *cptr = *ptr;  /* copy 1 byte */
	      ptr++;
	      cptr += 8;  
	    }
	  }

	  break;
	default: /*  should never get here */
            ffpmsg("Error: unexpected use of GZIP_2 to compress a column");
	    *status = DATA_DECOMPRESSION_ERR;
            return(*status);

        }  /* end of switch */

      } else {  /* not GZIP_2, so just transpose the bytes */

          for (jj = 0; jj < naxis2; jj++) {  /* loop over number of rows in the output table */
	    cptr = buffer + (rmajor_colstart[ii] + jj * rmajor_colstart[ncols]);   /* addr to copy to */
	    memcpy(cptr, ptr, (size_t) rmajor_colwidth[ii]);
	 
	    ptr += (rmajor_colwidth[ii]);
	  }
      }

    }  /* end of ncols loop */

    /* copy the buffer of data to the output data unit */
    fits_get_hduaddrll(outfptr, &headstart, &datastart, &dataend, status);        
    ffmbyt(outfptr, datastart, 1, status);
    ffpbyt(outfptr, naxis1 * naxis2, buffer, status);
    free(buffer);
    free(transbuffer);
	
    /* reset internal table structure parameters */
    fits_set_hdustruc(outfptr, status);

    /* unshuffle the heap, if it exists */
    fits_gunzip_heap(infptr, outfptr, status);

    return(*status);
}

/*--------------------------------------------------------------------------*/
int fits_gzip_datablocks(fitsfile *fptr, size_t *size, int *status)
/*
  GZIP compress all the data blocks in the binary table HDU.
  Store the size of the compressed byte stream in the PCOUNT keyword.  
  Save the original PCOUNT value in the ZPCOUNT keyword.  
*/
{ 
    long headstart, datastart, dataend;
    char *ptr, *cptr, *iptr;
    size_t dlen, datasize, ii;

    /* allocate memory for the data and the compressed data */
    fits_get_hduaddr(fptr, &headstart, &datastart, &dataend, status); 
    datasize = dataend - datastart;
    ptr = malloc(datasize);
    cptr = malloc(datasize);
    if (!ptr || !cptr) {
        ffpmsg("data memory allocation error in fits_gzip_datablocks\n");
	return(-1);
    }

    /* copy the data into memory */
    ffmbyt(fptr,datastart, REPORT_EOF, status);
    iptr = ptr;
    for (ii = 0; ii < datasize; ii+= 2880) {
	ffgbyt(fptr, 2880, iptr, status);
	iptr += 2880;
    }
	
    /* gzip compress the data */
    compress2mem_from_mem(ptr, datasize,
	&cptr,  &datasize, realloc, 
	&dlen, status);        

    *size = dlen;

    free(cptr);   /* don't need the compressed data any more */
    free(ptr);  /* don't need the original data any more */
   
    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fits_gzip_heap(fitsfile *infptr, fitsfile *outfptr, int *status)

/*
  Compress the binary table heap in the input file and write it to the output file.
  First, shuffle the bytes for the numeric arrays in the heap, so that
  the bytes are sorted in order of decreasing significance.  Then gzip
  the entire heap as a single block of data.  Then append this compressed heap
  to the end of any existing data in the output file heap.
*/
{ 
    LONGLONG datastart, dataend, nrows, naxis1, heapsize, length, offset, pcount, jj;
    int coltype, ncols, ii;
    char *heap, *compheap, card[FLEN_CARD];
    size_t theapsize, compsize;
    
    if (*status > 0)
        return(*status);

    /* insert a set of COMMENT keyword to indicate that this is a compressed table */
    fits_read_card(outfptr, "TFIELDS", card, status);
    fits_insert_card(outfptr, "COMMENT [FPACK] This is a compressed binary table generated by fpack.", status);
    fits_insert_card(outfptr, "COMMENT [FPACK] It can be uncompressed using funpack.", status);
    fits_insert_card(outfptr, "COMMENT [FPACK] fpack and funpack are available from the HEASARC Web site.", status);     
    
    /* get the size of the heap (value of PCOUNT keyword) */
    fits_read_key(infptr, TLONGLONG, "PCOUNT", &heapsize, NULL, status);

    /* return if there is no heap */
    if (*status != 0 || heapsize == 0)
        return(*status);

    /* allocate memory for the heap and compressed heap */
         
    heap = malloc((size_t) heapsize);
    if (!heap) {
        ffpmsg("Could not allocate buffer for the heap (fits_gzip_heap");
        *status = MEMORY_ALLOCATION;
        return(*status);
    }

    compheap = malloc((size_t) heapsize);
    if (!compheap) {
        ffpmsg("Could not allocate buffer for compressed heap (fits_gzip_heap");
 	free(heap);
        *status = MEMORY_ALLOCATION;
        return(*status);
    }

    fits_get_hduaddrll(infptr, NULL, &datastart, NULL, status); 
    fits_get_num_rowsll(infptr, &nrows, status);
    fits_get_num_cols(infptr, &ncols, status);
    fits_read_key(infptr, TLONGLONG, "NAXIS1", &naxis1, NULL, status);

    /* move to start of the heap and copy the heap into memory */
    ffmbyt(infptr, datastart + (nrows * naxis1), REPORT_EOF, status);
    ffgbyt(infptr, heapsize, heap, status);
    
    /* shuffle the bytes for the numeric columns */
    for (ii = 1; ii <= ncols; ii++) {

        fits_get_coltype(infptr, ii, &coltype, NULL, NULL, status);

	if (coltype >= 0) continue;   /* only interested in variable length columns */
	
	coltype = coltype * (-1);
	
	switch (coltype) {
	/* shuffle the bytes for the 2-byte, 4-byte, and 8-byte numeric columns */
	case TSHORT:

	  for (jj = 1; jj <= nrows; jj++) {
	    fits_read_descriptll(infptr, ii, jj, &length, &offset, status);
	    fits_shuffle_2bytes(heap + offset, length, status);    
	  }
	  break;
	
	case TLONG:
	case TFLOAT:
	  for (jj = 1; jj <= nrows; jj++) {
	    fits_read_descriptll(infptr, ii, jj, &length, &offset, status);
	    fits_shuffle_4bytes(heap + offset, length, status);    
	  }
	  break;

	case TDOUBLE:
	case TLONGLONG:
	  for (jj = 1; jj <= nrows; jj++) {
	    fits_read_descriptll(infptr, ii, jj, &length, &offset, status);
	    fits_shuffle_8bytes(heap + offset, length,  status);    
	  }
	  break;

	default: /* don't have to do anything for other column types */
	  break;

	}   /* end of switch block */
    }

    /* gzip compress the shuffled heap */
    theapsize = (size_t) heapsize;
    compress2mem_from_mem(heap, (size_t) heapsize, &compheap,  &theapsize, 
                          realloc, &compsize, status);        
    free(heap);  /* don't need the uncompresse heap any more */
    
    /* update the internal pointers */
    fits_set_hdustruc(outfptr, status);
    
    /* save offset to the start of the compressed heap, relative to the
       start of the main data table in the ZHEAPPTR keyword, and
       update PCOUNT to the new extended heap size */
       
    fits_read_key(outfptr, TLONGLONG, "PCOUNT", &pcount, NULL, status);
    fits_get_num_rowsll(outfptr, &nrows, status);
    fits_read_key(outfptr, TLONGLONG, "NAXIS1", &naxis1, NULL, status);

    fits_write_key_lng(outfptr, "ZHEAPPTR", (LONGLONG) ((nrows * naxis1) + pcount), 
                   "byte offset to compressed heap", status);
    fits_modify_key_lng(outfptr, "PCOUNT", pcount + compsize, NULL, status);

    /* now append the compressed heap to the heap in the output file */
    dataend = (outfptr->Fptr)->datastart + (outfptr->Fptr)->heapstart + 
                    (outfptr->Fptr)->heapsize;

    ffmbyt(outfptr, dataend, IGNORE_EOF, status);
    ffpbyt(outfptr, compsize, compheap, status);
    free(compheap);   
    
    /* also update the internal pointer to the heap size */
    (outfptr->Fptr)->heapsize = (outfptr->Fptr)->heapsize + compsize;

    /* update the internal pointers again */
    fits_set_hdustruc(outfptr, status);
 
    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fits_shuffle_2bytes(char *heap, LONGLONG length, int *status)

/* shuffle the bytes in an array of 2-byte integers in the heap */

{
    LONGLONG ii;
    char *ptr, *cptr, *heapptr;
    
    ptr = malloc((size_t) (length * 2));
    heapptr = heap;
    cptr = ptr;
    
    for (ii = 0; ii < length; ii++) {
       *cptr = *heapptr;
       heapptr++;
       *(cptr + length) = *heapptr;
       heapptr++;
       cptr++;
    }
         
    memcpy(heap, ptr, (size_t) (length * 2));
    free(ptr);
    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fits_shuffle_4bytes(char *heap, LONGLONG length, int *status)

/* shuffle the bytes in an array of 4-byte integers or floats  */

{
    LONGLONG ii;
    char *ptr, *cptr, *heapptr;
    
    ptr = malloc((size_t) (length * 4));
    if (!ptr) {
      ffpmsg("malloc failed\n");
      return(*status);
    }

    heapptr = heap;
    cptr = ptr;
 
    for (ii = 0; ii < length; ii++) {
       *cptr = *heapptr;
       heapptr++;
       *(cptr + length) = *heapptr;
       heapptr++;
       *(cptr + (length * 2)) = *heapptr;
       heapptr++;
       *(cptr + (length * 3)) = *heapptr;
       heapptr++;
       cptr++;
    }
        
    memcpy(heap, ptr, (size_t) (length * 4));
    free(ptr);

    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fits_shuffle_8bytes(char *heap, LONGLONG length, int *status)

/* shuffle the bytes in an array of 8-byte integers or doubles in the heap */

{
    LONGLONG ii;
    char *ptr, *cptr, *heapptr;
    
    ptr = calloc(1, (size_t) (length * 8));
    heapptr = heap;
    
/* for some bizarre reason this loop fails to compile under OpenSolaris using
   the proprietary SunStudioExpress C compiler;  use the following equivalent
   loop instead.
   
    cptr = ptr;

    for (ii = 0; ii < length; ii++) {
       *cptr = *heapptr;
       heapptr++;
       *(cptr + length) = *heapptr;
       heapptr++;
       *(cptr + (length * 2)) = *heapptr;
       heapptr++;
       *(cptr + (length * 3)) = *heapptr;
       heapptr++;
       *(cptr + (length * 4)) = *heapptr;
       heapptr++;
       *(cptr + (length * 5)) = *heapptr;
       heapptr++;
       *(cptr + (length * 6)) = *heapptr;
       heapptr++;
       *(cptr + (length * 7)) = *heapptr;
       heapptr++;
       cptr++;
     }
*/
     for (ii = 0; ii < length; ii++) {
        cptr = ptr + ii;

        *cptr = *heapptr;

        heapptr++;
        cptr += length;
        *cptr = *heapptr;

        heapptr++;
        cptr += length;
        *cptr = *heapptr;

        heapptr++;
        cptr += length;
        *cptr = *heapptr;

        heapptr++;
        cptr += length;
        *cptr = *heapptr;

        heapptr++;
        cptr += length;
        *cptr = *heapptr;

        heapptr++;
        cptr += length;
        *cptr = *heapptr;

        heapptr++;
        cptr += length;
        *cptr = *heapptr;

        heapptr++;
     }
        
    memcpy(heap, ptr, (size_t) (length * 8));
    free(ptr);
 
    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fits_gunzip_heap(fitsfile *infptr, fitsfile *outfptr, int *status)

/*
   inverse of the fits_gzip_heap function: uncompress and unshuffle the heap
   in the input file and write it to the output file
*/
{ 
    LONGLONG datastart, nrows, naxis1, length, offset, pcount, jj;
    LONGLONG zpcount, zheapptr, cheapsize;
    int coltype, ncols, ii;
    char *heap, *compheap;
    size_t arraysize, theapsize;

    if (*status > 0)
        return(*status);

    /* first, delete any COMMENT keywords written by fits_gzip_heap */
    while (*status == 0) {
        fits_delete_str(outfptr, "COMMENT [FPACK]", status);
    }
    if (*status == KEY_NO_EXIST) *status = 0;

    /* ZPCOUNT = size of original uncompressed heap */
    fits_read_key(infptr, TLONGLONG, "ZPCOUNT", &zpcount, NULL, status);

    /* just return if there is no heap */
    if (*status != 0 || zpcount == 0)
        return(*status);

    fits_get_num_rowsll(infptr, &nrows, status);
    fits_read_key(infptr, TLONGLONG, "NAXIS1", &naxis1, NULL, status);

    /* ZHEAPPTR = offset to the start of the compressed heap */
    fits_read_key(infptr, TLONGLONG, "ZHEAPPTR", &zheapptr, NULL, status);

    /* PCOUNT = total size of the compressed 2D table plus the compressed heap */
    fits_read_key(infptr, TLONGLONG, "PCOUNT", &pcount, NULL, status);

    /* size of the compressed heap */
    cheapsize = pcount - (zheapptr - (naxis1 * nrows));

    /* allocate memory for the heap and uncompressed heap */
    arraysize = (size_t) zpcount;
    heap = malloc(arraysize);
    if (!heap) {
        ffpmsg("Could not allocate buffer for the heap (fits_gunzip_heap");
        *status = MEMORY_ALLOCATION;
        return(*status);
    }

    compheap = malloc((size_t) cheapsize);
    if (!compheap) {
        ffpmsg("Could not allocate buffer for compressed heap (fits_gunzip_heap");
 	free(heap);
        *status = MEMORY_ALLOCATION;
        return(*status);
    }

    fits_get_hduaddrll(infptr, NULL, &datastart, NULL, status); 

    /* read the compressed heap into memory */
    ffmbyt(infptr, datastart + zheapptr, REPORT_EOF, status);
    ffgbyt(infptr, cheapsize, compheap, status);

    /* uncompress the heap */
    theapsize = (size_t) zpcount;
    uncompress2mem_from_mem(compheap, (size_t) cheapsize, &heap, &arraysize, 
        realloc, &theapsize, status);        

    free(compheap);   /* don't need the compressed heap any more */

    if (theapsize != zpcount) {
       /* something is wrong */
       ffpmsg("uncompressed heap size != to ZPCOUNT");
       free(heap);
       *status = MEMORY_ALLOCATION;
       return(*status);
    }

    /* get dimensions of the uncompressed table */
    fits_get_num_rowsll(outfptr, &nrows, status);
    fits_read_key(outfptr, TLONGLONG, "NAXIS1", &naxis1, NULL, status);
    fits_get_num_cols(outfptr, &ncols, status);

    for (ii = ncols; ii > 0; ii--) {

        fits_get_coltype(outfptr, ii, &coltype, NULL, NULL, status);

	if (coltype >= 0) continue;   /* only interested in variable length columns */
	
	coltype = coltype * (-1);
	
	switch (coltype) {
	/* recombine the byte planes for the 2-byte, 4-byte, and 8-byte numeric columns */
	case TSHORT:

	  for (jj = nrows; jj > 0; jj--) {
	    fits_read_descriptll(outfptr, ii, jj, &length, &offset, status);
	    fits_unshuffle_2bytes(heap + offset, length, status);    
	  }
	  break;
	
	case TLONG:
	case TFLOAT:
	  for (jj = nrows; jj > 0; jj--) {
	    fits_read_descriptll(outfptr, ii, jj, &length, &offset, status);
	    fits_unshuffle_4bytes(heap + offset, length, status);    
	  }
	  break;

	case TDOUBLE:
	case TLONGLONG:
	  for (jj = nrows; jj > 0; jj--) {
	    fits_read_descriptll(outfptr, ii, jj, &length, &offset, status);
	    fits_unshuffle_8bytes(heap + offset, length,  status);    
	  }
	  break;

	default: /* don't need to recombine bytes for other column types */
	  break;

	}   /* end of switch block */
    }

    /* copy the unshuffled heap back to the output file */
    fits_get_hduaddrll(outfptr, NULL, &datastart, NULL, status); 

    ffmbyt(outfptr, datastart + (nrows * naxis1), IGNORE_EOF, status);
    ffpbyt(outfptr, zpcount, heap, status);

    free(heap);   

    /* also update the internal pointer to the heap size */
    (outfptr->Fptr)->heapsize = zpcount;

    /* update the internal pointers again */
    fits_set_hdustruc(outfptr, status);
 
    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fits_unshuffle_2bytes(char *heap, LONGLONG length, int *status)

/* unshuffle the bytes in an array of 2-byte integers */

{
    LONGLONG ii;
    char *ptr, *cptr, *heapptr;
    
    ptr = malloc((size_t) (length * 2));
    heapptr = heap + (2 * length) - 1;
    cptr = ptr + (2 * length) - 1;
    
    for (ii = 0; ii < length; ii++) {
       *cptr = *heapptr;
       cptr--;
       *cptr = *(heapptr - length);
       cptr--;
       heapptr--;
    }
         
    memcpy(heap, ptr, (size_t) (length * 2));
    free(ptr);
    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fits_unshuffle_4bytes(char *heap, LONGLONG length, int *status)

/* unshuffle the bytes in an array of 4-byte integers or floats */

{
    LONGLONG ii;
    char *ptr, *cptr, *heapptr;
    
    ptr = malloc((size_t) (length * 4));
    heapptr = heap + (4 * length) -1;
    cptr = ptr + (4 * length) -1;
 
    for (ii = 0; ii < length; ii++) {
       *cptr = *heapptr;
       cptr--;
       *cptr = *(heapptr - length);
       cptr--;
       *cptr = *(heapptr - (2 * length));
       cptr--;
       *cptr = *(heapptr - (3 * length));
       cptr--;
       heapptr--;
    }
        
    memcpy(heap, ptr, (size_t) (length * 4));
    free(ptr);
    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fits_unshuffle_8bytes(char *heap, LONGLONG length, int *status)

/* unshuffle the bytes in an array of 8-byte integers or doubles */

{
    LONGLONG ii;
    char *ptr, *cptr, *heapptr;
    
    ptr = malloc((size_t) (length * 8));
    heapptr = heap + (8 * length) - 1;
    cptr = ptr + (8 * length)  -1;
    
    for (ii = 0; ii < length; ii++) {
       *cptr = *heapptr;
       cptr--;
       *cptr = *(heapptr - length);
       cptr--;
       *cptr = *(heapptr - (2 * length));
       cptr--;
       *cptr = *(heapptr - (3 * length));
       cptr--;
       *cptr = *(heapptr - (4 * length));
       cptr--;
       *cptr = *(heapptr - (5 * length));
       cptr--;
       *cptr = *(heapptr - (6 * length));
       cptr--;
       *cptr = *(heapptr - (7 * length));
       cptr--;
       heapptr--;
    }
       
    memcpy(heap, ptr, (size_t) (length * 8));
    free(ptr);
    return(*status);
}
