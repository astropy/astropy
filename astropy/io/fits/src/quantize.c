/* $Id$
*/

/*****************************************************************************/
/*                                                                           */
/* This file, quantize.c, contains the code required to quantize a floating  */
/* point value.                                                              */
/*                                                                           */
/* Copyright (C) 2004 Association of Universities for Research in Astronomy  */
/* (AURA)                                                                    */
/*                                                                           */
/* Redistribution and use in source and binary forms, with or without        */
/* modification, are permitted provided that the following conditions are    */
/* met:                                                                      */
/*                                                                           */
/*    1. Redistributions of source code must retain the above copyright      */
/*      notice, this list of conditions and the following disclaimer.        */
/*                                                                           */
/*    2. Redistributions in binary form must reproduce the above             */
/*      copyright notice, this list of conditions and the following          */
/*      disclaimer in the documentation and/or other materials provided      */
/*      with the distribution.                                               */
/*                                                                           */
/*    3. The name of AURA and its representatives may not be used to         */
/*      endorse or promote products derived from this software without       */
/*      specific prior written permission.                                   */
/*                                                                           */
/* THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED    */
/* WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF      */
/* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE                  */
/* DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,    */
/* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,      */
/* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS     */
/* OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND    */
/* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR     */
/* TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE    */
/* USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH          */
/* DAMAGE.                                                                   */
/*                                                                           */
/* This file was copied intact from the FITSIO software that was written by  */
/* William Pence at the High Energy Astrophysic Science Archive Research     */
/* Center (HEASARC) at the NASA Goddard Space Flight Center.  That software  */
/* contained the following copyright and warranty notices:                   */
/*                                                                           */
/* Copyright (Unpublished--all rights reserved under the copyright laws of   */
/* the United States), U.S. Government as represented by the Administrator   */
/* of the National Aeronautics and Space Administration.  No copyright is    */
/* claimed in the United States under Title 17, U.S. Code.                   */
/*                                                                           */
/* Permission to freely use, copy, modify, and distribute this software      */
/* and its documentation without fee is hereby granted, provided that this   */
/* copyright notice and disclaimer of warranty appears in all copies.        */
/*                                                                           */
/* DISCLAIMER:                                                               */
/*                                                                           */
/* THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND,        */
/* EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO,   */
/* ANY WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY        */
/* IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR           */
/* PURPOSE, AND FREEDOM FROM INFRINGEMENT, AND ANY WARRANTY THAT THE         */
/* DOCUMENTATION WILL CONFORM TO THE SOFTWARE, OR ANY WARRANTY THAT THE      */
/* SOFTWARE WILL BE ERROR FREE.  IN NO EVENT SHALL NASA BE LIABLE FOR ANY    */
/* DAMAGES, INCLUDING, BUT NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR      */
/* CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN ANY WAY      */
/* CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY,         */
/* CONTRACT, TORT , OR OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY     */
/* PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS SUSTAINED   */
/* FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR          */
/* SERVICES PROVIDED HEREUNDER."                                             */
/*                                                                           */
/*****************************************************************************/

/*
  The following code is based on algorithms written by Richard White at STScI and made
  available for use in CFITSIO in July 1999 and updated in January 2008. 
*/

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <limits.h>
# include <float.h>

#include "fitsio.h"

/* nearest integer function */
# define NINT(x)  ((x >= 0.) ? (int) (x + 0.5) : (int) (x - 0.5))

#define NULL_VALUE -2147483647 /* value used to represent undefined pixels */
#define N_RESERVED_VALUES 1   /* number of reserved values, starting with */
                               /* and including NULL_VALUE.  These values */
                               /* may not be used to represent the quantized */
                               /* and scaled floating point pixel values */

/* more than this many standard deviations from the mean is an outlier */
# define SIGMA_CLIP     5.
# define NITER          3	/* number of sigma-clipping iterations */

static int FnMeanSigma_short(short *array, long npix, int nullcheck, 
  short nullvalue, long *ngoodpix, double *mean, double *sigma, int *status);       
static int FnMeanSigma_int(int *array, long npix, int nullcheck,
  int nullvalue, long *ngoodpix, double *mean, double *sigma, int *status);       
static int FnMeanSigma_float(float *array, long npix, int nullcheck,
  float nullvalue, long *ngoodpix, double *mean, double *sigma, int *status);       
static int FnMeanSigma_double(double *array, long npix, int nullcheck,
  double nullvalue, long *ngoodpix, double *mean, double *sigma, int *status);       

static int FnNoise3_short(short *array, long nx, long ny, int nullcheck, 
   short nullvalue, long *ngood, short *minval, short *maxval, double *noise, int *status);       
static int FnNoise3_int(int *array, long nx, long ny, int nullcheck, 
   int nullvalue, long *ngood, int *minval, int *maxval, double *noise, int *status);          
static int FnNoise3_float(float *array, long nx, long ny, int nullcheck, 
   float nullvalue, long *ngood, float *minval, float *maxval, double *noise, int *status);        
static int FnNoise3_double(double *array, long nx, long ny, int nullcheck, 
   double nullvalue, long *ngood, double *minval, double *maxval, double *noise, int *status);        

static int FnNoise1_short(short *array, long nx, long ny, 
   int nullcheck, short nullvalue, double *noise, int *status);       
static int FnNoise1_int(int *array, long nx, long ny, 
   int nullcheck, int nullvalue, double *noise, int *status);       
static int FnNoise1_float(float *array, long nx, long ny, 
   int nullcheck, float nullvalue, double *noise, int *status);       
static int FnNoise1_double(double *array, long nx, long ny, 
   int nullcheck, double nullvalue, double *noise, int *status);       

static int FnCompare_double (const void *, const void *);
static float quick_select_float(float arr[], int n);
static short quick_select_short(short arr[], int n);
static int quick_select_int(int arr[], int n);
static double quick_select_double(double arr[], int n);
/*---------------------------------------------------------------------------*/
int _pyfits_fits_quantize_float (float fdata[], long nxpix, long nypix,
        int nullcheck, 
	float in_null_value, float qlevel, int idata[], double *bscale,
	double *bzero, int *iminval, int *imaxval) {

/* arguments:
float fdata[]       i: array of image pixels to be compressed
long nxpix          i: number of pixels in each row of fdata
long nypix          i: number of rows in fdata
nullcheck           i: check for nullvalues in fdata?
float in_null_value i: value used to represent undefined pixels in fdata
int noise_bits      i: quantization level (number of bits)
int idata[]         o: values of fdata after applying bzero and bscale
double bscale       o: scale factor
double bzero        o: zero offset
int iminval         o: minimum quantized value that is returned
int imaxval         o: maximum quantized value that is returned

The function value will be one if the input fdata were copied to idata;
in this case the parameters bscale and bzero can be used to convert back to
nearly the original floating point values:  fdata ~= idata * bscale + bzero.
If the function value is zero, the data were not copied to idata.
*/

	int status, intflag, nshift, itemp, anynulls = 0;
	long i, nx, ngood = 0;
	double stdev;	/* mean and RMS of differences */
	float minval = 0., maxval = 0.;  /* min & max of fdata */
	double delta;		/* bscale, 1 in idata = delta in fdata */
	double zeropt;	        /* bzero */
	double temp;

	nx = nxpix * nypix;
	if (nx <= 1) {
	    *bscale = 1.;
	    *bzero  = 0.;
	    return (0);
	}

        *iminval = INT32_MAX;
        *imaxval = INT32_MIN;

	/* Check to see if data are "floating point integer." */
        /* This also catches the case where all the pixels are null */

        /* Since idata and fdata may point to the same memory location, */
	/* we cannot write to idata unless we are sure we don't need   */
	/* the corresponding float value any more */
	
	intflag = 1;		/* initial value */
	for (i = 0;  i < nx;  i++) {
            if (nullcheck && fdata[i] == in_null_value) {
                anynulls = 1;
            }
	    else if (fdata[i] > INT32_MAX || 
                     fdata[i] < NULL_VALUE + N_RESERVED_VALUES) {
		intflag = 0;	/* not integer */
		break;
	    }
            else {
  	        itemp = (int)(fdata[i] + 0.5f);

	        if (itemp != fdata[i]) {
		    intflag = 0;	/* not integer */
		    break;
                }
	    }
	}

        if (intflag) { /* data are "floating point integer" */
	  for (i = 0;  i < nx;  i++) {
            if (nullcheck && fdata[i] == in_null_value) {
                idata[i] = NULL_VALUE;
                anynulls = 1;
            }
            else {
  	        idata[i] = (int)(fdata[i] + 0.5);
                *iminval = minvalue(idata[i], *iminval);
                *imaxval = maxvalue(idata[i], *imaxval);
	    }
	  }
	}

	if (intflag) {  /* data are "floating point integer" */
            if (anynulls) {
                /* Shift the range of values so they lie close to NULL_VALUE. */
                /* This will make the compression more efficient.             */
                /* Maximum allowed shift is 2^31 - 1 = 2147483646 */
                /* Can't use 2147483647 because OSF says this is not a legal number */

                if (*iminval >= 0) {
		   nshift = -(NULL_VALUE + 1) - N_RESERVED_VALUES;
		} else {
                  nshift = *iminval - NULL_VALUE - N_RESERVED_VALUES;
                }

                for (i = 0;  i < nx;  i++) {
                    if (idata[i] != NULL_VALUE) {
                        idata[i] -= nshift;
                    }
                }
                *iminval = *iminval - nshift;
                *imaxval = *imaxval - nshift;
  	        *bscale = 1.;
	        *bzero = (double) nshift;
            }
            else {
                /* there were no null values, so no need to shift the range */
  	        *bscale = 1.;
	        *bzero = 0.;
            }
	    return (1);
	}

	/* ************************************************************ */
        /* data are not "floating point integer"; need to quantize them */

        if (qlevel >= 0.) {

	    /* estimate background noise using 3rd order absolute pixel differences */
	    FnNoise3_float(fdata, nxpix, nypix, nullcheck, in_null_value, &ngood,
	        &minval, &maxval, &stdev, &status);      

	    /* substitute sigma-clipping if median is zero */
	    if (stdev == 0.0) {

		FnNoise1_float(fdata, nxpix, nypix, nullcheck, in_null_value, 
		&stdev, &status);       
	    }

	    if (qlevel == 0.)
	        delta = stdev / 16.;  /* default quantization */
	    else
	        delta = stdev / qlevel;

	    if (delta == 0)
	        return (0);   /* Zero variance in differences!  Don't quantize. */

	} else {
	    /* negative value represents the absolute quantization level */
	    delta = -qlevel;

	    /* only nned to calculate the min and max values */
	    FnNoise3_float(fdata, nxpix, nypix, nullcheck, in_null_value, 0,
	        &minval, &maxval, 0, &status);      
 	}

        /* check that the range of quantized levels is not > range of int */
	if ((maxval - minval) / delta > 2. * 2147483647. - N_RESERVED_VALUES )
	    return (0);			/* don't quantize */

        if (ngood == nx) {   /* don't have to check for nulls */
            /* return all positive values, if possible since some */
            /* compression algorithms either only work for positive integers, */
            /* or are more efficient.  */
            if ((maxval - minval) / delta < 2147483647. - N_RESERVED_VALUES )
            {
                zeropt = minval;
            }
            else
            {
                /* center the quantized levels around zero */
                zeropt = (minval + maxval) / 2.;
            }

       	    for (i = 0;  i < nx;  i++) {
	        idata[i] = NINT ((fdata[i] - zeropt) / delta);
            }
        }
        else {
            /* data contains null values; shift the range to be */
            /* close to the value used to represent null values */
            zeropt = minval - delta * (NULL_VALUE + N_RESERVED_VALUES);

	    for (i = 0;  i < nx;  i++) {
                if (fdata[i] != in_null_value) {
	            idata[i] = NINT ((fdata[i] - zeropt) / delta);
                }
                else
                    idata[i] = NULL_VALUE;
            }
	}

        /* calc min and max values */
        temp = (minval - zeropt) / delta;
        *iminval =  NINT (temp);
        temp = (maxval - zeropt) / delta;
        *imaxval =  NINT (temp);

	*bscale = delta;
	*bzero = zeropt;

	return (1);			/* yes, data have been quantized */
}
/*---------------------------------------------------------------------------*/
int _pyfits_fits_quantize_double (double fdata[], long nxpix, long nypix,
        int nullcheck, 
	double in_null_value, float qlevel, int idata[], double *bscale,
	double *bzero, int *iminval, int *imaxval) {

/* arguments:
double fdata[]       i: array of image pixels to be compressed
long nxpix          i: number of pixels in each row of fdata
long nypix          i: number of rows in fdata
nullcheck           i: check for nullvalues in fdata?
double in_null_value i: value used to represent undefined pixels in fdata
int noise_bits      i: quantization level (number of bits)
int idata[]         o: values of fdata after applying bzero and bscale
double bscale       o: scale factor
double bzero        o: zero offset
int iminval         o: minimum quantized value that is returned
int imaxval         o: maximum quantized value that is returned

The function value will be one if the input fdata were copied to idata;
in this case the parameters bscale and bzero can be used to convert back to
nearly the original floating point values:  fdata ~= idata * bscale + bzero.
If the function value is zero, the data were not copied to idata.
*/

	int status, intflag, nshift, itemp, anynulls = 0;
	long i, nx, ngood = 0;
	double stdev;	/* mean and RMS of differences */
	double minval = 0., maxval = 0.;  /* min & max of fdata */
	double delta;		/* bscale, 1 in idata = delta in fdata */
	double zeropt;	        /* bzero */
	double temp;

	nx = nxpix * nypix;
	if (nx <= 1) {
	    *bscale = 1.;
	    *bzero  = 0.;
	    return (0);
	}

        *iminval = INT32_MAX;
        *imaxval = INT32_MIN;

	/* Check to see if data are "floating point integer." */
        /* This also catches the case where all the pixels are null */

        /* Since idata and fdata may point to the same memory location, */
	/* we cannot write to idata unless we are sure we don't need   */
	/* the corresponding float value any more */
	
	intflag = 1;		/* initial value */
	for (i = 0;  i < nx;  i++) {
            if (nullcheck && fdata[i] == in_null_value) {
                anynulls = 1;
            }
	    else if (fdata[i] > INT32_MAX || 
                     fdata[i] < NULL_VALUE + N_RESERVED_VALUES) {
		intflag = 0;	/* not integer */
		break;
	    }
            else {
  	        itemp = (int)(fdata[i] + 0.5);

	        if (itemp != fdata[i]) {
		    intflag = 0;	/* not integer */
		    break;
                }
	    }
	}

        if (intflag) { /* data are "floating point integer" */
	  for (i = 0;  i < nx;  i++) {
            if (nullcheck && fdata[i] == in_null_value) {
                idata[i] = NULL_VALUE;
                anynulls = 1;
            }
            else {
  	        idata[i] = (int)(fdata[i] + 0.5);
                *iminval = minvalue(idata[i], *iminval);
                *imaxval = maxvalue(idata[i], *imaxval);
	    }
	  }
	}

	if (intflag) {  /* data are "floating point integer" */
            if (anynulls) {
                /* Shift the range of values so they lie close to NULL_VALUE. */
                /* This will make the compression more efficient.             */
                /* Maximum allowed shift is 2^31 - 1 = 2147483646 */
                /* Can't use 2147483647 because OSF says this is not a legal number */

                if (*iminval >= 0) {
		   nshift = -(NULL_VALUE + 1) - N_RESERVED_VALUES;
		} else {
                  nshift = *iminval - NULL_VALUE - N_RESERVED_VALUES;
                }

                for (i = 0;  i < nx;  i++) {
                    if (idata[i] != NULL_VALUE) {
                        idata[i] -= nshift;
                    }
                }
                *iminval = *iminval - nshift;
                *imaxval = *imaxval - nshift;
  	        *bscale = 1.;
	        *bzero = (double) nshift;
            }
            else {
                /* there were no null values, so no need to shift the range */
  	        *bscale = 1.;
	        *bzero = 0.;
            }
	    return (1);
	}

	/* ************************************************************ */
        /* data are not "floating point integer"; need to quantize them */

        if (qlevel >= 0.) {

	    /* estimate background noise using 3rd order absolute pixel differences */
	    FnNoise3_double(fdata, nxpix, nypix, nullcheck, in_null_value, &ngood,
	        &minval, &maxval, &stdev, &status);      

	    /* substitute sigma-clipping if median is zero */
	    if (stdev == 0.0) {

		FnNoise1_double(fdata, nxpix, nypix, nullcheck, in_null_value, 
		&stdev, &status);       
	    }

	    if (qlevel == 0.)
	        delta = stdev / 16.;  /* default quantization */
	    else
	        delta = stdev / qlevel;

	    if (delta == 0)
	        return (0);   /* Zero variance in differences!  Don't quantize. */

	} else {
	    /* negative value represents the absolute quantization level */
	    delta = -qlevel;

	    /* only nned to calculate the min and max values */
	    FnNoise3_double(fdata, nxpix, nypix, nullcheck, in_null_value, 0,
	        &minval, &maxval, 0, &status);      
 	}

        /* check that the range of quantized levels is not > range of int */
	if ((maxval - minval) / delta > 2. * 2147483647. - N_RESERVED_VALUES )
	    return (0);			/* don't quantize */

        if (ngood == nx) {   /* don't have to check for nulls */
            /* return all positive values, if possible since some */
            /* compression algorithms either only work for positive integers, */
            /* or are more efficient.  */
            if ((maxval - minval) / delta < 2147483647. - N_RESERVED_VALUES )
            {
                zeropt = minval;
            }
            else
            {
                /* center the quantized levels around zero */
                zeropt = (minval + maxval) / 2.;
            }

       	    for (i = 0;  i < nx;  i++) {
	        temp = (fdata[i] - zeropt) / delta;
	        idata[i] = NINT (temp);
            }
        }
        else {
            /* data contains null values; shift the range to be */
            /* close to the value used to represent null values */
            zeropt = minval - delta * (NULL_VALUE + N_RESERVED_VALUES);

	    for (i = 0;  i < nx;  i++) {
                if (fdata[i] != in_null_value) {
	            temp = (fdata[i] - zeropt) / delta;
	            idata[i] = NINT (temp);
                }
                else
                    idata[i] = NULL_VALUE;
            }
	}

        /* calc min and max values */
        temp = (minval - zeropt) / delta;
        *iminval =  NINT (temp);
        temp = (maxval - zeropt) / delta;
        *imaxval =  NINT (temp);

	*bscale = delta;
	*bzero = zeropt;

	return (1);			/* yes, data have been quantized */
}
/*--------------------------------------------------------------------------*/
static int fits_img_stats_short(
        short *array, /*  2 dimensional array of image pixels */
        long nx,            /* number of pixels in each row of the image */
	long ny,            /* number of rows in the image */
	                    /* (if this is a 3D image, then ny should be the */
			    /* product of the no. of rows times the no. of planes) */
	int nullcheck,      /* check for null values, if true */
	short nullvalue,    /* value of null pixels, if nullcheck is true */

   /* returned parameters (if the pointer is not null)  */
	long *ngoodpix,     /* number of non-null pixels in the image */
	short *minvalue,    /* returned minimum non-null value in the array */
	short *maxvalue,    /* returned maximum non-null value in the array */
	double *mean,       /* returned mean value of all non-null pixels */
	double *sigma,      /* returned R.M.S. value of all non-null pixels */
	double *noise1,     /* 1st order estimate of noise in image background level */
	double *noise3,     /* 3rd order estimate of noise in image background level */
	int *status)        /* error status */

/*
    Compute statistics of the input short integer image.
*/
{
	long ngood = 0;
	short minval = 0, maxval = 0;
	double xmean = 0., xsigma = 0., xnoise = 0;

	/* need to calculate mean and/or sigma and/or limits? */
	if (mean || sigma ) {
		FnMeanSigma_short(array, nx * ny, nullcheck, nullvalue, 
			&ngood, &xmean, &xsigma, status);

	    if (ngoodpix) *ngoodpix = ngood;
	    if (mean)     *mean = xmean;
	    if (sigma)    *sigma = xsigma;
	}

	if (noise1) {
		FnNoise1_short(array, nx, ny, nullcheck, nullvalue, 
		  &xnoise, status);

		*noise1  = xnoise;
	}

	if (minvalue || maxvalue || noise3) {
		FnNoise3_short(array, nx, ny, nullcheck, nullvalue, 
			&ngood, &minval, &maxval, &xnoise, status);

		if (ngoodpix) *ngoodpix = ngood;
		if (minvalue) *minvalue= minval;
		if (maxvalue) *maxvalue = maxval;
		*noise3  = xnoise;
	}
	return(*status);
}
/*--------------------------------------------------------------------------*/
int _pyfits_fits_img_stats_int(
        int *array, /*  2 dimensional array of image pixels */
        long nx,            /* number of pixels in each row of the image */
	long ny,            /* number of rows in the image */
	                    /* (if this is a 3D image, then ny should be the */
			    /* product of the no. of rows times the no. of planes) */
	int nullcheck,      /* check for null values, if true */
	int nullvalue,    /* value of null pixels, if nullcheck is true */

   /* returned parameters (if the pointer is not null)  */
	long *ngoodpix,     /* number of non-null pixels in the image */
	int *minvalue,    /* returned minimum non-null value in the array */
	int *maxvalue,    /* returned maximum non-null value in the array */
	double *mean,       /* returned mean value of all non-null pixels */
	double *sigma,      /* returned R.M.S. value of all non-null pixels */
	double *noise1,     /* 1st order estimate of noise in image background level */
	double *noise3,     /* 3rd order estimate of noise in image background level */
	int *status)        /* error status */

/*
    Compute statistics of the input short integer image.
*/
{
	long ngood = 0;
	int minval = 0, maxval = 0;
	double xmean = 0., xsigma = 0., xnoise = 0;

	/* need to calculate mean and/or sigma and/or limits? */
	if (mean || sigma ) {
		FnMeanSigma_int(array, nx * ny, nullcheck, nullvalue, 
			&ngood, &xmean, &xsigma, status);

	    if (ngoodpix) *ngoodpix = ngood;
	    if (mean)     *mean = xmean;
	    if (sigma)    *sigma = xsigma;
	}

	if (noise1) {
		FnNoise1_int(array, nx, ny, nullcheck, nullvalue, 
		  &xnoise, status);

		*noise1  = xnoise;
	}

	if (minvalue || maxvalue || noise3) {
		FnNoise3_int(array, nx, ny, nullcheck, nullvalue, 
			&ngood, &minval, &maxval, &xnoise, status);

		if (ngoodpix) *ngoodpix = ngood;
		if (minvalue) *minvalue= minval;
		if (maxvalue) *maxvalue = maxval;
		*noise3  = xnoise;
	}
	return(*status);
}
/*--------------------------------------------------------------------------*/
static int fits_img_stats_float(
        float *array, /*  2 dimensional array of image pixels */
        long nx,            /* number of pixels in each row of the image */
	long ny,            /* number of rows in the image */
	                    /* (if this is a 3D image, then ny should be the */
			    /* product of the no. of rows times the no. of planes) */
	int nullcheck,      /* check for null values, if true */
	float nullvalue,    /* value of null pixels, if nullcheck is true */

   /* returned parameters (if the pointer is not null)  */
	long *ngoodpix,     /* number of non-null pixels in the image */
	float *minvalue,    /* returned minimum non-null value in the array */
	float *maxvalue,    /* returned maximum non-null value in the array */
	double *mean,       /* returned mean value of all non-null pixels */
	double *sigma,      /* returned R.M.S. value of all non-null pixels */
	double *noise1,     /* 1st order estimate of noise in image background level */
	double *noise3,     /* 3rd order estimate of noise in image background level */
	int *status)        /* error status */

/*
    Compute statistics of the input short integer image.
*/
{
	long ngood;
	float minval, maxval;
	double xmean = 0., xsigma = 0., xnoise = 0;

	/* need to calculate mean and/or sigma and/or limits? */
	if (mean || sigma ) {
		FnMeanSigma_float(array, nx * ny, nullcheck, nullvalue, 
			&ngood, &xmean, &xsigma, status);

	    if (ngoodpix) *ngoodpix = ngood;
	    if (mean)     *mean = xmean;
	    if (sigma)    *sigma = xsigma;
	}

	if (noise1) {
		FnNoise1_float(array, nx, ny, nullcheck, nullvalue, 
		  &xnoise, status);

		*noise1  = xnoise;
	}

	if (minvalue || maxvalue || noise3) {
		FnNoise3_float(array, nx, ny, nullcheck, nullvalue, 
			&ngood, &minval, &maxval, &xnoise, status);

		if (ngoodpix) *ngoodpix = ngood;
		if (minvalue) *minvalue= minval;
		if (maxvalue) *maxvalue = maxval;
		*noise3  = xnoise;
	}
	return(*status);
}
/*--------------------------------------------------------------------------*/
static int FnMeanSigma_short
       (short *array,       /*  2 dimensional array of image pixels */
        long npix,          /* number of pixels in the image */
	int nullcheck,      /* check for null values, if true */
	short nullvalue,    /* value of null pixels, if nullcheck is true */

   /* returned parameters */
   
	long *ngoodpix,     /* number of non-null pixels in the image */
	double *mean,       /* returned mean value of all non-null pixels */
	double *sigma,      /* returned R.M.S. value of all non-null pixels */
	int *status)        /* error status */

/*
Compute mean and RMS sigma of the non-null pixels in the input array.
*/
{
	long ii, ngood = 0;
	short *value;
	double sum = 0., sum2 = 0., xtemp;

	value = array;
	    
	if (nullcheck) {
	        for (ii = 0; ii < npix; ii++, value++) {
		    if (*value != nullvalue) {
		        ngood++;
		        xtemp = (double) *value;
		        sum += xtemp;
		        sum2 += (xtemp * xtemp);
		    }
		}
	} else {
	        ngood = npix;
	        for (ii = 0; ii < npix; ii++, value++) {
		        xtemp = (double) *value;
		        sum += xtemp;
		        sum2 += (xtemp * xtemp);
		}
	}

	if (ngood > 1) {
		if (ngoodpix) *ngoodpix = ngood;
		xtemp = sum / ngood;
		if (mean)     *mean = xtemp;
		if (sigma)    *sigma = sqrt((sum2 / ngood) - (xtemp * xtemp));
	} else if (ngood == 1){
		if (ngoodpix) *ngoodpix = 1;
		if (mean)     *mean = sum;
		if (sigma)    *sigma = 0.0;
	} else {
		if (ngoodpix) *ngoodpix = 0;
	        if (mean)     *mean = 0.;
		if (sigma)    *sigma = 0.;
	}	    
	return(*status);
}
/*--------------------------------------------------------------------------*/
static int FnMeanSigma_int
       (int *array,       /*  2 dimensional array of image pixels */
        long npix,          /* number of pixels in the image */
	int nullcheck,      /* check for null values, if true */
	int nullvalue,    /* value of null pixels, if nullcheck is true */

   /* returned parameters */
   
	long *ngoodpix,     /* number of non-null pixels in the image */
	double *mean,       /* returned mean value of all non-null pixels */
	double *sigma,      /* returned R.M.S. value of all non-null pixels */
	int *status)        /* error status */

/*
Compute mean and RMS sigma of the non-null pixels in the input array.
*/
{
	long ii, ngood = 0;
	int *value;
	double sum = 0., sum2 = 0., xtemp;

	value = array;
	    
	if (nullcheck) {
	        for (ii = 0; ii < npix; ii++, value++) {
		    if (*value != nullvalue) {
		        ngood++;
		        xtemp = (double) *value;
		        sum += xtemp;
		        sum2 += (xtemp * xtemp);
		    }
		}
	} else {
	        ngood = npix;
	        for (ii = 0; ii < npix; ii++, value++) {
		        xtemp = (double) *value;
		        sum += xtemp;
		        sum2 += (xtemp * xtemp);
		}
	}

	if (ngood > 1) {
		if (ngoodpix) *ngoodpix = ngood;
		xtemp = sum / ngood;
		if (mean)     *mean = xtemp;
		if (sigma)    *sigma = sqrt((sum2 / ngood) - (xtemp * xtemp));
	} else if (ngood == 1){
		if (ngoodpix) *ngoodpix = 1;
		if (mean)     *mean = sum;
		if (sigma)    *sigma = 0.0;
	} else {
		if (ngoodpix) *ngoodpix = 0;
	        if (mean)     *mean = 0.;
		if (sigma)    *sigma = 0.;
	}	    
	return(*status);
}
/*--------------------------------------------------------------------------*/
static int FnMeanSigma_float
       (float *array,       /*  2 dimensional array of image pixels */
        long npix,          /* number of pixels in the image */
	int nullcheck,      /* check for null values, if true */
	float nullvalue,    /* value of null pixels, if nullcheck is true */

   /* returned parameters */
   
	long *ngoodpix,     /* number of non-null pixels in the image */
	double *mean,       /* returned mean value of all non-null pixels */
	double *sigma,      /* returned R.M.S. value of all non-null pixels */
	int *status)        /* error status */

/*
Compute mean and RMS sigma of the non-null pixels in the input array.
*/
{
	long ii, ngood = 0;
	float *value;
	double sum = 0., sum2 = 0., xtemp;

	value = array;
	    
	if (nullcheck) {
	        for (ii = 0; ii < npix; ii++, value++) {
		    if (*value != nullvalue) {
		        ngood++;
		        xtemp = (double) *value;
		        sum += xtemp;
		        sum2 += (xtemp * xtemp);
		    }
		}
	} else {
	        ngood = npix;
	        for (ii = 0; ii < npix; ii++, value++) {
		        xtemp = (double) *value;
		        sum += xtemp;
		        sum2 += (xtemp * xtemp);
		}
	}

	if (ngood > 1) {
		if (ngoodpix) *ngoodpix = ngood;
		xtemp = sum / ngood;
		if (mean)     *mean = xtemp;
		if (sigma)    *sigma = sqrt((sum2 / ngood) - (xtemp * xtemp));
	} else if (ngood == 1){
		if (ngoodpix) *ngoodpix = 1;
		if (mean)     *mean = sum;
		if (sigma)    *sigma = 0.0;
	} else {
		if (ngoodpix) *ngoodpix = 0;
	        if (mean)     *mean = 0.;
		if (sigma)    *sigma = 0.;
	}	    
	return(*status);
}
/*--------------------------------------------------------------------------*/
static int FnMeanSigma_double
       (double *array,       /*  2 dimensional array of image pixels */
        long npix,          /* number of pixels in the image */
	int nullcheck,      /* check for null values, if true */
	double nullvalue,    /* value of null pixels, if nullcheck is true */

   /* returned parameters */
   
	long *ngoodpix,     /* number of non-null pixels in the image */
	double *mean,       /* returned mean value of all non-null pixels */
	double *sigma,      /* returned R.M.S. value of all non-null pixels */
	int *status)        /* error status */

/*
Compute mean and RMS sigma of the non-null pixels in the input array.
*/
{
	long ii, ngood = 0;
	double *value;
	double sum = 0., sum2 = 0., xtemp;

	value = array;
	    
	if (nullcheck) {
	        for (ii = 0; ii < npix; ii++, value++) {
		    if (*value != nullvalue) {
		        ngood++;
		        xtemp = *value;
		        sum += xtemp;
		        sum2 += (xtemp * xtemp);
		    }
		}
	} else {
	        ngood = npix;
	        for (ii = 0; ii < npix; ii++, value++) {
		        xtemp = *value;
		        sum += xtemp;
		        sum2 += (xtemp * xtemp);
		}
	}

	if (ngood > 1) {
		if (ngoodpix) *ngoodpix = ngood;
		xtemp = sum / ngood;
		if (mean)     *mean = xtemp;
		if (sigma)    *sigma = sqrt((sum2 / ngood) - (xtemp * xtemp));
	} else if (ngood == 1){
		if (ngoodpix) *ngoodpix = 1;
		if (mean)     *mean = sum;
		if (sigma)    *sigma = 0.0;
	} else {
		if (ngoodpix) *ngoodpix = 0;
	        if (mean)     *mean = 0.;
		if (sigma)    *sigma = 0.;
	}	    
	return(*status);
}
/*--------------------------------------------------------------------------*/
static int FnNoise3_short
       (short *array,       /*  2 dimensional array of image pixels */
        long nx,            /* number of pixels in each row of the image */
        long ny,            /* number of rows in the image */
	int nullcheck,      /* check for null values, if true */
	short nullvalue,    /* value of null pixels, if nullcheck is true */
   /* returned parameters */   
	long *ngood,        /* number of good, non-null pixels? */
	short *minval,    /* minimum non-null value */
	short *maxval,    /* maximum non-null value */
	double *noise,      /* returned R.M.S. value of all non-null pixels */
	int *status)        /* error status */

/*
Estimate the median and background noise in the input image using 3rd order differences.

The noise in the background of the image is calculated using the 3rd order algorithm 
developed for deriving the signal to noise ratio in spectra
(see issue #42 of the ST-ECF newsletter, http://www.stecf.org/documents/newsletter/)

  noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))

The returned estimates are the median of the values that are computed for each 
row of the image.
*/
{
	long ii, jj, nrows = 0, nvals, ngoodpix = 0;
	short *differences, *rowpix, v1, v2, v3, v4, v5;
	short xminval = SHRT_MAX, xmaxval = SHRT_MIN, do_range = 0;
	double *diffs, xnoise = 0, sigma;
	
	if (nx < 5) {
		/* treat entire array as an image with a single row */
		nx = nx * ny;
		ny = 1;
	}

	/* rows must have at least 5 pixels */
	if (nx < 5) {

		for (ii = 0; ii < nx; ii++) {
		    if (nullcheck && array[ii] == nullvalue)
		        continue;
		    else {
			if (array[ii] < xminval) xminval = array[ii];
			if (array[ii] > xmaxval) xmaxval = array[ii];
			ngoodpix++;
		    }
		}
		if (minval) *minval = xminval;
		if (maxval) *maxval = xmaxval;
		if (ngood) *ngood = ngoodpix;
		if (noise) *noise = 0.;
		return(*status);
	}

	/* do we need to compute the min and max value? */
	if (minval || maxval) do_range = 1;
	
        /* allocate arrays used to compute the median and noise estimates */
	differences = calloc(nx, sizeof(short));
	if (!differences) {
        	*status = MEMORY_ALLOCATION;
		return(*status);
	}

	diffs = calloc(ny, sizeof(double));
	if (!diffs) {
		free(differences);
        	*status = MEMORY_ALLOCATION;
		return(*status);
	}

	/* loop over each row of the image */
	for (jj=0; jj < ny; jj++) {

                rowpix = array + (jj * nx); /* point to first pixel in the row */

		/***** find the first valid pixel in row */
		ii = 0;
		if (nullcheck)
		    while (ii < nx && rowpix[ii] == nullvalue) ii++;

		if (ii == nx) continue;  /* hit end of row */
		v1 = rowpix[ii];  /* store the good pixel value */

		if (do_range) {
			if (v1 < xminval) xminval = v1;
			if (v1 > xmaxval) xmaxval = v1;
		}

		/***** find the 2nd valid pixel in row (which we will skip over) */
		ii++;
		if (nullcheck)
		    while (ii < nx && rowpix[ii] == nullvalue) ii++;

		if (ii == nx) continue;  /* hit end of row */
		v2 = rowpix[ii];  /* store the good pixel value */
		
		if (do_range) {
			if (v2 < xminval) xminval = v2;
			if (v2 > xmaxval) xmaxval = v2;
		}

		/***** find the 3rd valid pixel in row */
		ii++;
		if (nullcheck)
		    while (ii < nx && rowpix[ii] == nullvalue) ii++;

		if (ii == nx) continue;  /* hit end of row */
		v3 = rowpix[ii];  /* store the good pixel value */

		if (do_range) {
			if (v3 < xminval) xminval = v3;
			if (v3 > xmaxval) xmaxval = v3;
		}
				
		/* find the 4nd valid pixel in row (to be skipped) */
		ii++;
		if (nullcheck)
		    while (ii < nx && rowpix[ii] == nullvalue) ii++;

		if (ii == nx) continue;  /* hit end of row */
		v4 = rowpix[ii];  /* store the good pixel value */

		if (do_range) {
			if (v4 < xminval) xminval = v4;
			if (v4 > xmaxval) xmaxval = v4;
		}
		
		/* now populate the differences arrays */
		/* for the remaining pixels in the row */
		nvals = 0;
		for (ii++; ii < nx; ii++) {

		    /* find the next valid pixel in row */
                    if (nullcheck)
		        while (ii < nx && rowpix[ii] == nullvalue) ii++;
		     
		    if (ii == nx) break;  /* hit end of row */
		    v5 = rowpix[ii];  /* store the good pixel value */

		    if (do_range) {
			if (v5 < xminval) xminval = v5;
			if (v5 > xmaxval) xmaxval = v5;
		    }

		    /* construct array of 3rd order absolute differences */
		    if (!(v1 == v2 && v2 == v3 && v3 == v4 && v4 == v5)) {
		        differences[nvals] = abs((2 * v3) - v1 - v5);
		        nvals++;  
		    } else {
		        /* ignore constant background regions */
			ngoodpix++;
		    }


		    /* shift over 1 pixel */
		    v1 = v2;
		    v2 = v3;
		    v3 = v4;
		    v4 = v5;
	        }  /* end of loop over pixels in the row */

		/* compute the 3rd order diffs */
		/* Note that there are 4 more pixel values than there are diffs values. */
		ngoodpix += (nvals + 4);

		if (nvals == 0) {
		    continue;  /* cannot compute medians on this row */
		} else if (nvals == 1) {
		    diffs[nrows] = differences[0];
		} else {
                    /* quick_select returns the median MUCH faster than using qsort */
                    diffs[nrows] = quick_select_short(differences, nvals);
		}

		nrows++;
	}  /* end of loop over rows */

	    /* compute median of the values for each row */
	if (nrows == 0) { 
	       xnoise = 0;
	} else if (nrows == 1) {
	       xnoise = diffs[0];
	} else {	    


	       qsort(diffs, nrows, sizeof(double), FnCompare_double);
	       xnoise =  (diffs[(nrows - 1)/2] + diffs[nrows/2]) / 2.;

              FnMeanSigma_double(diffs, nrows, 0, 0.0, 0, &xnoise, &sigma, status); 

	      /* do a 4.5 sigma rejection of outliers */
	      jj = 0;
	      sigma = 4.5 * sigma;
	      for (ii = 0; ii < nrows; ii++) {
		if ( fabs(diffs[ii] - xnoise) <= sigma)	 {
		   if (jj != ii)
		       diffs[jj] = diffs[ii];
		   jj++;
	        } 
	      }
	      if (ii != jj)
                FnMeanSigma_double(diffs, jj, 0, 0.0, 0, &xnoise, &sigma, status); 
	}

	if (ngood)  *ngood  = ngoodpix;
	if (minval) *minval = xminval;
	if (maxval) *maxval = xmaxval;
	if (noise)  *noise  = 0.6052697 * xnoise;

	free(diffs);
	free(differences);

	return(*status);
}
/*--------------------------------------------------------------------------*/
static int FnNoise3_int
       (int *array,       /*  2 dimensional array of image pixels */
        long nx,            /* number of pixels in each row of the image */
        long ny,            /* number of rows in the image */
	int nullcheck,      /* check for null values, if true */
	int nullvalue,    /* value of null pixels, if nullcheck is true */
   /* returned parameters */   
	long *ngood,        /* number of good, non-null pixels? */
	int *minval,    /* minimum non-null value */
	int *maxval,    /* maximum non-null value */
	double *noise,      /* returned R.M.S. value of all non-null pixels */
	int *status)        /* error status */

/*
Estimate the background noise in the input image using 3rd order differences.

The noise in the background of the image is calculated using the 3rd order algorithm 
developed for deriving the signal to noise ratio in spectra
(see issue #42 of the ST-ECF newsletter, http://www.stecf.org/documents/newsletter/)

  noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))

The returned estimates are the median of the values that are computed for each 
row of the image.
*/
{
	long ii, jj, nrows = 0, nvals, ngoodpix = 0;
	int *differences, *rowpix, v1, v2, v3, v4, v5;
	int xminval = INT_MAX, xmaxval = INT_MIN, do_range = 0;
	double *diffs, xnoise = 0, sigma;
	
	if (nx < 5) {
		/* treat entire array as an image with a single row */
		nx = nx * ny;
		ny = 1;
	}

	/* rows must have at least 5 pixels */
	if (nx < 5) {

		for (ii = 0; ii < nx; ii++) {
		    if (nullcheck && array[ii] == nullvalue)
		        continue;
		    else {
			if (array[ii] < xminval) xminval = array[ii];
			if (array[ii] > xmaxval) xmaxval = array[ii];
			ngoodpix++;
		    }
		}
		if (minval) *minval = xminval;
		if (maxval) *maxval = xmaxval;
		if (ngood) *ngood = ngoodpix;
		if (noise) *noise = 0.;
		return(*status);
	}

	/* do we need to compute the min and max value? */
	if (minval || maxval) do_range = 1;
	
        /* allocate arrays used to compute the median and noise estimates */
	differences = calloc(nx, sizeof(int));
	if (!differences) {
        	*status = MEMORY_ALLOCATION;
		return(*status);
	}

	diffs = calloc(ny, sizeof(double));
	if (!diffs) {
		free(differences);
        	*status = MEMORY_ALLOCATION;
		return(*status);
	}

	/* loop over each row of the image */
	for (jj=0; jj < ny; jj++) {

                rowpix = array + (jj * nx); /* point to first pixel in the row */

		/***** find the first valid pixel in row */
		ii = 0;
		if (nullcheck)
		    while (ii < nx && rowpix[ii] == nullvalue) ii++;

		if (ii == nx) continue;  /* hit end of row */
		v1 = rowpix[ii];  /* store the good pixel value */

		if (do_range) {
			if (v1 < xminval) xminval = v1;
			if (v1 > xmaxval) xmaxval = v1;
		}

		/***** find the 2nd valid pixel in row (which we will skip over) */
		ii++;
		if (nullcheck)
		    while (ii < nx && rowpix[ii] == nullvalue) ii++;

		if (ii == nx) continue;  /* hit end of row */
		v2 = rowpix[ii];  /* store the good pixel value */
		
		if (do_range) {
			if (v2 < xminval) xminval = v2;
			if (v2 > xmaxval) xmaxval = v2;
		}

		/***** find the 3rd valid pixel in row */
		ii++;
		if (nullcheck)
		    while (ii < nx && rowpix[ii] == nullvalue) ii++;

		if (ii == nx) continue;  /* hit end of row */
		v3 = rowpix[ii];  /* store the good pixel value */

		if (do_range) {
			if (v3 < xminval) xminval = v3;
			if (v3 > xmaxval) xmaxval = v3;
		}
				
		/* find the 4nd valid pixel in row (to be skipped) */
		ii++;
		if (nullcheck)
		    while (ii < nx && rowpix[ii] == nullvalue) ii++;

		if (ii == nx) continue;  /* hit end of row */
		v4 = rowpix[ii];  /* store the good pixel value */

		if (do_range) {
			if (v4 < xminval) xminval = v4;
			if (v4 > xmaxval) xmaxval = v4;
		}
		
		/* now populate the differences arrays */
		/* for the remaining pixels in the row */
		nvals = 0;
		for (ii++; ii < nx; ii++) {

		    /* find the next valid pixel in row */
                    if (nullcheck)
		        while (ii < nx && rowpix[ii] == nullvalue) ii++;
		     
		    if (ii == nx) break;  /* hit end of row */
		    v5 = rowpix[ii];  /* store the good pixel value */

		    if (do_range) {
			if (v5 < xminval) xminval = v5;
			if (v5 > xmaxval) xmaxval = v5;
		    }

		    /* construct array of 3rd order absolute differences */
		    if (!(v1 == v2 && v2 == v3 && v3 == v4 && v4 == v5)) {
		        differences[nvals] = abs((2 * v3) - v1 - v5);
		        nvals++;  
		    } else {
		        /* ignore constant background regions */
			ngoodpix++;
		    }

		    /* shift over 1 pixel */
		    v1 = v2;
		    v2 = v3;
		    v3 = v4;
		    v4 = v5;
	        }  /* end of loop over pixels in the row */

		/* compute the 3rd order diffs */
		/* Note that there are 4 more pixel values than there are diffs values. */
		ngoodpix += (nvals + 4);

		if (nvals == 0) {
		    continue;  /* cannot compute medians on this row */
		} else if (nvals == 1) {
		    diffs[nrows] = differences[0];
		} else {
                    /* quick_select returns the median MUCH faster than using qsort */
                    diffs[nrows] = quick_select_int(differences, nvals);
		}

		nrows++;
	}  /* end of loop over rows */

	    /* compute median of the values for each row */
	if (nrows == 0) { 
	       xnoise = 0;
	} else if (nrows == 1) {
	       xnoise = diffs[0];
	} else {	    

	       qsort(diffs, nrows, sizeof(double), FnCompare_double);
	       xnoise =  (diffs[(nrows - 1)/2] + diffs[nrows/2]) / 2.;

              FnMeanSigma_double(diffs, nrows, 0, 0.0, 0, &xnoise, &sigma, status); 

	      /* do a 4.5 sigma rejection of outliers */
	      jj = 0;
	      sigma = 4.5 * sigma;
	      for (ii = 0; ii < nrows; ii++) {
		if ( fabs(diffs[ii] - xnoise) <= sigma)	 {
		   if (jj != ii)
		       diffs[jj] = diffs[ii];
		   jj++;
	        }
	      }
	      if (ii != jj)
                FnMeanSigma_double(diffs, jj, 0, 0.0, 0, &xnoise, &sigma, status); 
	}

	if (ngood)  *ngood  = ngoodpix;
	if (minval) *minval = xminval;
	if (maxval) *maxval = xmaxval;
	if (noise)  *noise  = 0.6052697 * xnoise;

	free(diffs);
	free(differences);

	return(*status);
}
/*--------------------------------------------------------------------------*/
static int FnNoise3_float
       (float *array,       /*  2 dimensional array of image pixels */
        long nx,            /* number of pixels in each row of the image */
        long ny,            /* number of rows in the image */
	int nullcheck,      /* check for null values, if true */
	float nullvalue,    /* value of null pixels, if nullcheck is true */
   /* returned parameters */   
	long *ngood,        /* number of good, non-null pixels? */
	float *minval,    /* minimum non-null value */
	float *maxval,    /* maximum non-null value */
	double *noise,      /* returned R.M.S. value of all non-null pixels */
	int *status)        /* error status */

/*
Estimate the median and background noise in the input image using 3rd order differences.

The noise in the background of the image is calculated using the 3rd order algorithm 
developed for deriving the signal to noise ratio in spectra
(see issue #42 of the ST-ECF newsletter, http://www.stecf.org/documents/newsletter/)

  noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))

The returned estimates are the median of the values that are computed for each 
row of the image.
*/
{
	long ii, jj, nrows = 0, nvals, ngoodpix = 0;
	float *differences, *rowpix, v1, v2, v3, v4, v5;
	float xminval = FLT_MAX, xmaxval = FLT_MIN;
	int do_range = 0;
	double *diffs, xnoise = 0;
	
	if (nx < 5) {
		/* treat entire array as an image with a single row */
		nx = nx * ny;
		ny = 1;
	}

	/* rows must have at least 5 pixels to calc noise, so just calc min, max, ngood */
	if (nx < 5) {

		for (ii = 0; ii < nx; ii++) {
		    if (nullcheck && array[ii] == nullvalue)
		        continue;
		    else {
			if (array[ii] < xminval) xminval = array[ii];
			if (array[ii] > xmaxval) xmaxval = array[ii];
			ngoodpix++;
		    }
		}
		if (minval) *minval = xminval;
		if (maxval) *maxval = xmaxval;
		if (ngood) *ngood = ngoodpix;
		if (noise) *noise = 0.;
		return(*status);
	}

	/* do we need to compute the min and max value? */
	if (minval || maxval) do_range = 1;
	
        /* allocate arrays used to compute the median and noise estimates */
        differences = 0;
        diffs = 0;
	if (noise) {
	    differences = calloc(nx, sizeof(float));
	    if (!differences) {
        	*status = MEMORY_ALLOCATION;
		return(*status);
	    }

	    diffs = calloc(ny, sizeof(double));
	    if (!diffs) {
		free(differences);
        	*status = MEMORY_ALLOCATION;
		return(*status);
	    }
	}

	/* loop over each row of the image */
	for (jj=0; jj < ny; jj++) {

                rowpix = array + (jj * nx); /* point to first pixel in the row */

		/***** find the first valid pixel in row */
		ii = 0;
		if (nullcheck)
		    while (ii < nx && rowpix[ii] == nullvalue) ii++;

		if (ii == nx) continue;  /* hit end of row */
		v1 = rowpix[ii];  /* store the good pixel value */

		if (do_range) {
			if (v1 < xminval) xminval = v1;
			if (v1 > xmaxval) xmaxval = v1;
		}

		/***** find the 2nd valid pixel in row (which we will skip over) */
		ii++;
		if (nullcheck)
		    while (ii < nx && rowpix[ii] == nullvalue) ii++;

		if (ii == nx) continue;  /* hit end of row */
		v2 = rowpix[ii];  /* store the good pixel value */
		
		if (do_range) {
			if (v2 < xminval) xminval = v2;
			if (v2 > xmaxval) xmaxval = v2;
		}

		/***** find the 3rd valid pixel in row */
		ii++;
		if (nullcheck)
		    while (ii < nx && rowpix[ii] == nullvalue) ii++;

		if (ii == nx) continue;  /* hit end of row */
		v3 = rowpix[ii];  /* store the good pixel value */

		if (do_range) {
			if (v3 < xminval) xminval = v3;
			if (v3 > xmaxval) xmaxval = v3;
		}
				
		/* find the 4nd valid pixel in row (to be skipped) */
		ii++;
		if (nullcheck)
		    while (ii < nx && rowpix[ii] == nullvalue) ii++;

		if (ii == nx) continue;  /* hit end of row */
		v4 = rowpix[ii];  /* store the good pixel value */

		if (do_range) {
			if (v4 < xminval) xminval = v4;
			if (v4 > xmaxval) xmaxval = v4;
		}
		
		/* now populate the differences arrays */
		/* for the remaining pixels in the row */
		nvals = 0;
		for (ii++; ii < nx; ii++) {

		    /* find the next valid pixel in row */
                    if (nullcheck)
		        while (ii < nx && rowpix[ii] == nullvalue) ii++;
		     
		    if (ii == nx) break;  /* hit end of row */
		    v5 = rowpix[ii];  /* store the good pixel value */

		    if (do_range) {
			if (v5 < xminval) xminval = v5;
			if (v5 > xmaxval) xmaxval = v5;
		    }

		    /* construct array of 3rd order absolute differences */
		    if (noise) {
		        if (!(v1 == v2 && v2 == v3 && v3 == v4 && v4 == v5)) {

		            differences[nvals] = fabs((2. * v3) - v1 - v5);
		            nvals++;  
			}
		    } else {
		        /* ignore constant background regions */
			ngoodpix++;
		    }

		    /* shift over 1 pixel */
		    v1 = v2;
		    v2 = v3;
		    v3 = v4;
		    v4 = v5;
	        }  /* end of loop over pixels in the row */

		/* compute the 3rd order diffs */
		/* Note that there are 4 more pixel values than there are diffs values. */
		ngoodpix += (nvals + 4);

		if (noise) {
		    if (nvals == 0) {
		        continue;  /* cannot compute medians on this row */
		    } else if (nvals == 1) {
		        diffs[nrows] = differences[0];
		    } else {
                        /* quick_select returns the median MUCH faster than using qsort */
                        diffs[nrows] = quick_select_float(differences, nvals);
		    }
		}
		nrows++;
	}  /* end of loop over rows */

	    /* compute median of the values for each row */
	if (noise) {
	    if (nrows == 0) { 
	       xnoise = 0;
	    } else if (nrows == 1) {
	       xnoise = diffs[0];
	    } else {	    
	       qsort(diffs, nrows, sizeof(double), FnCompare_double);
	       xnoise =  (diffs[(nrows - 1)/2] + diffs[nrows/2]) / 2.;
	    }
	}

	if (ngood)  *ngood  = ngoodpix;
	if (minval) *minval = xminval;
	if (maxval) *maxval = xmaxval;
	if (noise) {
		*noise  = 0.6052697 * xnoise;
		free(diffs);
		free(differences);
	}

	return(*status);
}
/*--------------------------------------------------------------------------*/
static int FnNoise3_double
       (double *array,       /*  2 dimensional array of image pixels */
        long nx,            /* number of pixels in each row of the image */
        long ny,            /* number of rows in the image */
	int nullcheck,      /* check for null values, if true */
	double nullvalue,    /* value of null pixels, if nullcheck is true */
   /* returned parameters */   
	long *ngood,        /* number of good, non-null pixels? */
	double *minval,    /* minimum non-null value */
	double *maxval,    /* maximum non-null value */
	double *noise,      /* returned R.M.S. value of all non-null pixels */
	int *status)        /* error status */

/*
Estimate the median and background noise in the input image using 3rd order differences.

The noise in the background of the image is calculated using the 3rd order algorithm 
developed for deriving the signal to noise ratio in spectra
(see issue #42 of the ST-ECF newsletter, http://www.stecf.org/documents/newsletter/)

  noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))

The returned estimates are the median of the values that are computed for each 
row of the image.
*/
{
	long ii, jj, nrows = 0, nvals, ngoodpix = 0;
	double *differences, *rowpix, v1, v2, v3, v4, v5;
	double xminval = DBL_MAX, xmaxval = DBL_MIN;
	int do_range = 0;
	double *diffs, xnoise = 0;
	
	if (nx < 5) {
		/* treat entire array as an image with a single row */
		nx = nx * ny;
		ny = 1;
	}

	/* rows must have at least 5 pixels */
	if (nx < 5) {

		for (ii = 0; ii < nx; ii++) {
		    if (nullcheck && array[ii] == nullvalue)
		        continue;
		    else {
			if (array[ii] < xminval) xminval = array[ii];
			if (array[ii] > xmaxval) xmaxval = array[ii];
			ngoodpix++;
		    }
		}
		if (minval) *minval = xminval;
		if (maxval) *maxval = xmaxval;
		if (ngood) *ngood = ngoodpix;
		if (noise) *noise = 0.;
		return(*status);
	}

	/* do we need to compute the min and max value? */
	if (minval || maxval) do_range = 1;
	
        /* allocate arrays used to compute the median and noise estimates */
        differences = 0;
        diffs = 0;
	if (noise) {
	    differences = calloc(nx, sizeof(double));
	    if (!differences) {
        	*status = MEMORY_ALLOCATION;
		return(*status);
	    }

	    diffs = calloc(ny, sizeof(double));
	    if (!diffs) {
		free(differences);
        	*status = MEMORY_ALLOCATION;
		return(*status);
	    }
	}

	/* loop over each row of the image */
	for (jj=0; jj < ny; jj++) {

                rowpix = array + (jj * nx); /* point to first pixel in the row */

		/***** find the first valid pixel in row */
		ii = 0;
		if (nullcheck)
		    while (ii < nx && rowpix[ii] == nullvalue) ii++;

		if (ii == nx) continue;  /* hit end of row */
		v1 = rowpix[ii];  /* store the good pixel value */

		if (do_range) {
			if (v1 < xminval) xminval = v1;
			if (v1 > xmaxval) xmaxval = v1;
		}

		/***** find the 2nd valid pixel in row (which we will skip over) */
		ii++;
		if (nullcheck)
		    while (ii < nx && rowpix[ii] == nullvalue) ii++;

		if (ii == nx) continue;  /* hit end of row */
		v2 = rowpix[ii];  /* store the good pixel value */
		
		if (do_range) {
			if (v2 < xminval) xminval = v2;
			if (v2 > xmaxval) xmaxval = v2;
		}

		/***** find the 3rd valid pixel in row */
		ii++;
		if (nullcheck)
		    while (ii < nx && rowpix[ii] == nullvalue) ii++;

		if (ii == nx) continue;  /* hit end of row */
		v3 = rowpix[ii];  /* store the good pixel value */

		if (do_range) {
			if (v3 < xminval) xminval = v3;
			if (v3 > xmaxval) xmaxval = v3;
		}
				
		/* find the 4nd valid pixel in row (to be skipped) */
		ii++;
		if (nullcheck)
		    while (ii < nx && rowpix[ii] == nullvalue) ii++;

		if (ii == nx) continue;  /* hit end of row */
		v4 = rowpix[ii];  /* store the good pixel value */

		if (do_range) {
			if (v4 < xminval) xminval = v4;
			if (v4 > xmaxval) xmaxval = v4;
		}
		
		/* now populate the differences arrays */
		/* for the remaining pixels in the row */
		nvals = 0;
		for (ii++; ii < nx; ii++) {

		    /* find the next valid pixel in row */
                    if (nullcheck)
		        while (ii < nx && rowpix[ii] == nullvalue) ii++;
		     
		    if (ii == nx) break;  /* hit end of row */
		    v5 = rowpix[ii];  /* store the good pixel value */

		    if (do_range) {
			if (v5 < xminval) xminval = v5;
			if (v5 > xmaxval) xmaxval = v5;
		    }

		    /* construct array of 3rd order absolute differences */
		    if (noise) {
		        if (!(v1 == v2 && v2 == v3 && v3 == v4 && v4 == v5)) {

		            differences[nvals] = fabs((2. * v3) - v1 - v5);
		            nvals++;  
			}
		    } else {
		        /* ignore constant background regions */
			ngoodpix++;
		    }

		    /* shift over 1 pixel */
		    v1 = v2;
		    v2 = v3;
		    v3 = v4;
		    v4 = v5;
	        }  /* end of loop over pixels in the row */

		/* compute the 3rd order diffs */
		/* Note that there are 4 more pixel values than there are diffs values. */
		ngoodpix += (nvals + 4);

		if (noise) {
		    if (nvals == 0) {
		        continue;  /* cannot compute medians on this row */
		    } else if (nvals == 1) {
		        diffs[nrows] = differences[0];
		    } else {
                        /* quick_select returns the median MUCH faster than using qsort */
                        diffs[nrows] = quick_select_double(differences, nvals);
		    }
		}
		nrows++;
	}  /* end of loop over rows */

	    /* compute median of the values for each row */
	if (noise) {
	    if (nrows == 0) { 
	       xnoise = 0;
	    } else if (nrows == 1) {
	       xnoise = diffs[0];
	    } else {	    
	       qsort(diffs, nrows, sizeof(double), FnCompare_double);
	       xnoise =  (diffs[(nrows - 1)/2] + diffs[nrows/2]) / 2.;
	    }
	}

	if (ngood)  *ngood  = ngoodpix;
	if (minval) *minval = xminval;
	if (maxval) *maxval = xmaxval;
	if (noise) {
		*noise  = 0.6052697 * xnoise;
		free(diffs);
		free(differences);
	}

	return(*status);
}
/*--------------------------------------------------------------------------*/
static int FnNoise1_short
       (short *array,       /*  2 dimensional array of image pixels */
        long nx,            /* number of pixels in each row of the image */
        long ny,            /* number of rows in the image */
	int nullcheck,      /* check for null values, if true */
	short nullvalue,    /* value of null pixels, if nullcheck is true */
   /* returned parameters */   
	double *noise,      /* returned R.M.S. value of all non-null pixels */
	int *status)        /* error status */
/*
Estimate the background noise in the input image using sigma of 1st order differences.

  noise = 1.0 / sqrt(2) * rms of (flux[i] - flux[i-1])

The returned estimate is the median of the values that are computed for each 
row of the image.
*/
{
	int iter;
	long ii, jj, kk, nrows = 0, nvals;
	short *differences, *rowpix, v1;
	double  *diffs, xnoise, mean, stdev;

	/* rows must have at least 3 pixels to estimate noise */
	if (nx < 3) {
		*noise = 0;
		return(*status);
	}
	
        /* allocate arrays used to compute the median and noise estimates */
	differences = calloc(nx, sizeof(short));
	if (!differences) {
        	*status = MEMORY_ALLOCATION;
		return(*status);
	}

	diffs = calloc(ny, sizeof(double));
	if (!diffs) {
		free(differences);
        	*status = MEMORY_ALLOCATION;
		return(*status);
	}

	/* loop over each row of the image */
	for (jj=0; jj < ny; jj++) {

                rowpix = array + (jj * nx); /* point to first pixel in the row */

		/***** find the first valid pixel in row */
		ii = 0;
		if (nullcheck)
		    while (ii < nx && rowpix[ii] == nullvalue) ii++;

		if (ii == nx) continue;  /* hit end of row */
		v1 = rowpix[ii];  /* store the good pixel value */

		/* now continue populating the differences arrays */
		/* for the remaining pixels in the row */
		nvals = 0;
		for (ii++; ii < nx; ii++) {

		    /* find the next valid pixel in row */
                    if (nullcheck)
		        while (ii < nx && rowpix[ii] == nullvalue) ii++;
		     
		    if (ii == nx) break;  /* hit end of row */
		
		    /* construct array of 1st order differences */
		    differences[nvals] = v1 - rowpix[ii];

		    nvals++;  
		    /* shift over 1 pixel */
		    v1 = rowpix[ii];
	        }  /* end of loop over pixels in the row */

		if (nvals < 2)
		   continue;
		else {

		    FnMeanSigma_short(differences, nvals, 0, 0, 0, &mean, &stdev, status);

		    if (stdev > 0.) {
		        for (iter = 0;  iter < NITER;  iter++) {
		            kk = 0;
		            for (ii = 0;  ii < nvals;  ii++) {
		                if (fabs (differences[ii] - mean) < SIGMA_CLIP * stdev) {
			            if (kk < ii)
			                differences[kk] = differences[ii];
			            kk++;
		                }
		            }
		            if (kk == nvals) break;

		            nvals = kk;
		            FnMeanSigma_short(differences, nvals, 0, 0, 0, &mean, &stdev, status);
	              }
		   }

		   diffs[nrows] = stdev;
		   nrows++;
		}
	}  /* end of loop over rows */

	/* compute median of the values for each row */
	if (nrows == 0) { 
	       xnoise = 0;
	} else if (nrows == 1) {
	       xnoise = diffs[0];
	} else {
	       qsort(diffs, nrows, sizeof(double), FnCompare_double);
	       xnoise =  (diffs[(nrows - 1)/2] + diffs[nrows/2]) / 2.;
	}

	*noise = .70710678 * xnoise;

	free(diffs);
	free(differences);

	return(*status);
}
/*--------------------------------------------------------------------------*/
static int FnNoise1_int
       (int *array,       /*  2 dimensional array of image pixels */
        long nx,            /* number of pixels in each row of the image */
        long ny,            /* number of rows in the image */
	int nullcheck,      /* check for null values, if true */
	int nullvalue,    /* value of null pixels, if nullcheck is true */
   /* returned parameters */   
	double *noise,      /* returned R.M.S. value of all non-null pixels */
	int *status)        /* error status */
/*
Estimate the background noise in the input image using sigma of 1st order differences.

  noise = 1.0 / sqrt(2) * rms of (flux[i] - flux[i-1])

The returned estimate is the median of the values that are computed for each 
row of the image.
*/
{
	int iter;
	long ii, jj, kk, nrows = 0, nvals;
	int *differences, *rowpix, v1;
	double  *diffs, xnoise, mean, stdev;

	/* rows must have at least 3 pixels to estimate noise */
	if (nx < 3) {
		*noise = 0;
		return(*status);
	}
	
        /* allocate arrays used to compute the median and noise estimates */
	differences = calloc(nx, sizeof(int));
	if (!differences) {
        	*status = MEMORY_ALLOCATION;
		return(*status);
	}

	diffs = calloc(ny, sizeof(double));
	if (!diffs) {
		free(differences);
        	*status = MEMORY_ALLOCATION;
		return(*status);
	}

	/* loop over each row of the image */
	for (jj=0; jj < ny; jj++) {

                rowpix = array + (jj * nx); /* point to first pixel in the row */

		/***** find the first valid pixel in row */
		ii = 0;
		if (nullcheck)
		    while (ii < nx && rowpix[ii] == nullvalue) ii++;

		if (ii == nx) continue;  /* hit end of row */
		v1 = rowpix[ii];  /* store the good pixel value */

		/* now continue populating the differences arrays */
		/* for the remaining pixels in the row */
		nvals = 0;
		for (ii++; ii < nx; ii++) {

		    /* find the next valid pixel in row */
                    if (nullcheck)
		        while (ii < nx && rowpix[ii] == nullvalue) ii++;
		     
		    if (ii == nx) break;  /* hit end of row */
		
		    /* construct array of 1st order differences */
		    differences[nvals] = v1 - rowpix[ii];

		    nvals++;  
		    /* shift over 1 pixel */
		    v1 = rowpix[ii];
	        }  /* end of loop over pixels in the row */

		if (nvals < 2)
		   continue;
		else {

		    FnMeanSigma_int(differences, nvals, 0, 0, 0, &mean, &stdev, status);

		    if (stdev > 0.) {
		        for (iter = 0;  iter < NITER;  iter++) {
		            kk = 0;
		            for (ii = 0;  ii < nvals;  ii++) {
		                if (fabs (differences[ii] - mean) < SIGMA_CLIP * stdev) {
			            if (kk < ii)
			                differences[kk] = differences[ii];
			            kk++;
		                }
		            }
		            if (kk == nvals) break;

		            nvals = kk;
		            FnMeanSigma_int(differences, nvals, 0, 0, 0, &mean, &stdev, status);
	              }
		   }

		   diffs[nrows] = stdev;
		   nrows++;
		}
	}  /* end of loop over rows */

	/* compute median of the values for each row */
	if (nrows == 0) { 
	       xnoise = 0;
	} else if (nrows == 1) {
	       xnoise = diffs[0];
	} else {
	       qsort(diffs, nrows, sizeof(double), FnCompare_double);
	       xnoise =  (diffs[(nrows - 1)/2] + diffs[nrows/2]) / 2.;
	}

	*noise = .70710678 * xnoise;

	free(diffs);
	free(differences);

	return(*status);
}
/*--------------------------------------------------------------------------*/
static int FnNoise1_float
       (float *array,       /*  2 dimensional array of image pixels */
        long nx,            /* number of pixels in each row of the image */
        long ny,            /* number of rows in the image */
	int nullcheck,      /* check for null values, if true */
	float nullvalue,    /* value of null pixels, if nullcheck is true */
   /* returned parameters */   
	double *noise,      /* returned R.M.S. value of all non-null pixels */
	int *status)        /* error status */
/*
Estimate the background noise in the input image using sigma of 1st order differences.

  noise = 1.0 / sqrt(2) * rms of (flux[i] - flux[i-1])

The returned estimate is the median of the values that are computed for each 
row of the image.
*/
{
	int iter;
	long ii, jj, kk, nrows = 0, nvals;
	float *differences, *rowpix, v1;
	double  *diffs, xnoise, mean, stdev;

	/* rows must have at least 3 pixels to estimate noise */
	if (nx < 3) {
		*noise = 0;
		return(*status);
	}
	
        /* allocate arrays used to compute the median and noise estimates */
	differences = calloc(nx, sizeof(float));
	if (!differences) {
        	*status = MEMORY_ALLOCATION;
		return(*status);
	}

	diffs = calloc(ny, sizeof(double));
	if (!diffs) {
		free(differences);
        	*status = MEMORY_ALLOCATION;
		return(*status);
	}

	/* loop over each row of the image */
	for (jj=0; jj < ny; jj++) {

                rowpix = array + (jj * nx); /* point to first pixel in the row */

		/***** find the first valid pixel in row */
		ii = 0;
		if (nullcheck)
		    while (ii < nx && rowpix[ii] == nullvalue) ii++;

		if (ii == nx) continue;  /* hit end of row */
		v1 = rowpix[ii];  /* store the good pixel value */

		/* now continue populating the differences arrays */
		/* for the remaining pixels in the row */
		nvals = 0;
		for (ii++; ii < nx; ii++) {

		    /* find the next valid pixel in row */
                    if (nullcheck)
		        while (ii < nx && rowpix[ii] == nullvalue) ii++;
		     
		    if (ii == nx) break;  /* hit end of row */
		
		    /* construct array of 1st order differences */
		    differences[nvals] = v1 - rowpix[ii];

		    nvals++;  
		    /* shift over 1 pixel */
		    v1 = rowpix[ii];
	        }  /* end of loop over pixels in the row */

		if (nvals < 2)
		   continue;
		else {

		    FnMeanSigma_float(differences, nvals, 0, 0, 0, &mean, &stdev, status);

		    if (stdev > 0.) {
		        for (iter = 0;  iter < NITER;  iter++) {
		            kk = 0;
		            for (ii = 0;  ii < nvals;  ii++) {
		                if (fabs (differences[ii] - mean) < SIGMA_CLIP * stdev) {
			            if (kk < ii)
			                differences[kk] = differences[ii];
			            kk++;
		                }
		            }
		            if (kk == nvals) break;

		            nvals = kk;
		            FnMeanSigma_float(differences, nvals, 0, 0, 0, &mean, &stdev, status);
	              }
		   }

		   diffs[nrows] = stdev;
		   nrows++;
		}
	}  /* end of loop over rows */

	/* compute median of the values for each row */
	if (nrows == 0) { 
	       xnoise = 0;
	} else if (nrows == 1) {
	       xnoise = diffs[0];
	} else {
	       qsort(diffs, nrows, sizeof(double), FnCompare_double);
	       xnoise =  (diffs[(nrows - 1)/2] + diffs[nrows/2]) / 2.;
	}

	*noise = .70710678 * xnoise;

	free(diffs);
	free(differences);

	return(*status);
}
/*--------------------------------------------------------------------------*/
static int FnNoise1_double
       (double *array,       /*  2 dimensional array of image pixels */
        long nx,            /* number of pixels in each row of the image */
        long ny,            /* number of rows in the image */
	int nullcheck,      /* check for null values, if true */
	double nullvalue,    /* value of null pixels, if nullcheck is true */
   /* returned parameters */   
	double *noise,      /* returned R.M.S. value of all non-null pixels */
	int *status)        /* error status */
/*
Estimate the background noise in the input image using sigma of 1st order differences.

  noise = 1.0 / sqrt(2) * rms of (flux[i] - flux[i-1])

The returned estimate is the median of the values that are computed for each 
row of the image.
*/
{
	int iter;
	long ii, jj, kk, nrows = 0, nvals;
	double *differences, *rowpix, v1;
	double  *diffs, xnoise, mean, stdev;

	/* rows must have at least 3 pixels to estimate noise */
	if (nx < 3) {
		*noise = 0;
		return(*status);
	}
	
        /* allocate arrays used to compute the median and noise estimates */
	differences = calloc(nx, sizeof(double));
	if (!differences) {
        	*status = MEMORY_ALLOCATION;
		return(*status);
	}

	diffs = calloc(ny, sizeof(double));
	if (!diffs) {
		free(differences);
        	*status = MEMORY_ALLOCATION;
		return(*status);
	}

	/* loop over each row of the image */
	for (jj=0; jj < ny; jj++) {

                rowpix = array + (jj * nx); /* point to first pixel in the row */

		/***** find the first valid pixel in row */
		ii = 0;
		if (nullcheck)
		    while (ii < nx && rowpix[ii] == nullvalue) ii++;

		if (ii == nx) continue;  /* hit end of row */
		v1 = rowpix[ii];  /* store the good pixel value */

		/* now continue populating the differences arrays */
		/* for the remaining pixels in the row */
		nvals = 0;
		for (ii++; ii < nx; ii++) {

		    /* find the next valid pixel in row */
                    if (nullcheck)
		        while (ii < nx && rowpix[ii] == nullvalue) ii++;
		     
		    if (ii == nx) break;  /* hit end of row */
		
		    /* construct array of 1st order differences */
		    differences[nvals] = v1 - rowpix[ii];

		    nvals++;  
		    /* shift over 1 pixel */
		    v1 = rowpix[ii];
	        }  /* end of loop over pixels in the row */

		if (nvals < 2)
		   continue;
		else {

		    FnMeanSigma_double(differences, nvals, 0, 0, 0, &mean, &stdev, status);

		    if (stdev > 0.) {
		        for (iter = 0;  iter < NITER;  iter++) {
		            kk = 0;
		            for (ii = 0;  ii < nvals;  ii++) {
		                if (fabs (differences[ii] - mean) < SIGMA_CLIP * stdev) {
			            if (kk < ii)
			                differences[kk] = differences[ii];
			            kk++;
		                }
		            }
		            if (kk == nvals) break;

		            nvals = kk;
		            FnMeanSigma_double(differences, nvals, 0, 0, 0, &mean, &stdev, status);
	              }
		   }

		   diffs[nrows] = stdev;
		   nrows++;
		}
	}  /* end of loop over rows */

	/* compute median of the values for each row */
	if (nrows == 0) { 
	       xnoise = 0;
	} else if (nrows == 1) {
	       xnoise = diffs[0];
	} else {
	       qsort(diffs, nrows, sizeof(double), FnCompare_double);
	       xnoise =  (diffs[(nrows - 1)/2] + diffs[nrows/2]) / 2.;
	}

	*noise = .70710678 * xnoise;

	free(diffs);
	free(differences);

	return(*status);
}
/*--------------------------------------------------------------------------*/
static int FnCompare_double(const void *v1, const void *v2)
{
   const double *i1 = v1;
   const double *i2 = v2;
   
   if (*i1 < *i2)
     return(-1);
   else if (*i1 > *i2)
     return(1);
   else
     return(0);
}
/*--------------------------------------------------------------------------*/

/*
 *  These Quickselect routines are based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 */

/*--------------------------------------------------------------------------*/

#define ELEM_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }

static float quick_select_float(float arr[], int n) 
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        ELEM_SWAP(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
}

#undef ELEM_SWAP

/*--------------------------------------------------------------------------*/

#define ELEM_SWAP(a,b) { register short t=(a);(a)=(b);(b)=t; }

static short quick_select_short(short arr[], int n) 
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        ELEM_SWAP(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
}

#undef ELEM_SWAP

/*--------------------------------------------------------------------------*/

#define ELEM_SWAP(a,b) { register int t=(a);(a)=(b);(b)=t; }

static int quick_select_int(int arr[], int n) 
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        ELEM_SWAP(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
}

#undef ELEM_SWAP

/*--------------------------------------------------------------------------*/

#define ELEM_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }

static double quick_select_double(double arr[], int n) 
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        ELEM_SWAP(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
}

#undef ELEM_SWAP


