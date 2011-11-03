/* $Id$
*/

/*****************************************************************************/
/*                                                                           */
/* This file, fitsio.c, contains the core of a set of FITSIO routines that   */
/* are used to compress and uncompress image data in FITS binary tables.     */
/* The code was copied and modified from the FITSIO software written at      */
/* HEASRC.  The goal for the pyfitsComp module was to take this code nearly  */
/* intact.  In FITSIO, interaction with the FITS file is accomplished        */
/* directly within the FITSIO code.  With pyfitsComp, interaction with the   */
/* FITS file is accomplished from within pyfits.  This may make some of the  */
/* constructs within the FISTIO code seem confusing when viewed from the     */
/* perspective of pyfitsComp.  It should be noted that the FITSfile          */
/* structure acts as the interface to the file in both cases.  In FITSIO it  */
/* contains the file handle in order to access the file, and in pyfitsComp   */
/* it holds the file data, either compressed or uncompressed.                */
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
/* This file contains source code that was copied and modified from the      */
/* FITSIO software that was written by William Pence at the High Energy      */
/* Astrophysic Science Archive Research Center (HEASARC) at the NASA Goddard */
/* Space Flight Center.  That software contained the following copyright and */
/* warranty notices:                                                         */
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

#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include "fitsio.h"

/*****************************************************************************/
/*                                                                           */
/* The following code was copied and modified from the FITSIO source code    */
/* file fitscore.c.                                                          */
/*                                                                           */
/*****************************************************************************/

#define GetMesg    4 /* pop and return the message */
#define PutMesg    5 /* put a new message */

static int imcomp_nullfloats(float *fdata, long tilelen, int *idata,
    int nullcheck, float nullflagval, int nullval, int *status);
static int imcomp_nullscalefloats(float *fdata, long tilelen, int *idata,
    double scale, double zero, int nullcheck, float nullflagval, int nullval,
    int *status);
static int imcomp_nulldoubles(double *fdata, long tilelen, int *idata,
    int nullcheck, double nullflagval, int nullval, int *status);
static int imcomp_nullscaledoubles(double *fdata, long tilelen, int *idata,
    double scale, double zero, int nullcheck, double nullflagval, int nullval,
    int *status);
static int fits_write_compressed_pixels(fitsfile *fptr,
            int  datatype, LONGLONG  fpixel, LONGLONG npixels,
            int nullcheck,  void *array, void *nulval,
            int  *status);
static int fits_write_compressed_img_plane(fitsfile *fptr, int  datatype,
      int  bytesperpixel,  long   nplane, long *firstcoord, long *lastcoord,
      long *naxes,  int  nullcheck,
      void *array,  void *nullval, long *nread, int  *status);
static int fits_read_compressed_pixels(fitsfile *fptr,
            int  datatype, LONGLONG  fpixel, LONGLONG npixels,
            int nullcheck, void *nulval,  void *array, char *nullarray,
            int  *anynul, int  *status);
static int fits_read_compressed_img_plane(fitsfile *fptr, int  datatype,
      int  bytesperpixel,  long   nplane, LONGLONG *firstcoord,
      LONGLONG *lastcoord, long *inc, long *naxes, int  nullcheck,
      void *nullval, void *array, char *nullarray, int  *anynul, long *nread,
      int  *status);
static int imcomp_decompress_tile (fitsfile *infptr,
          int nrow, int tilesize, int datatype, int nullcheck,
          void *nulval, void *buffer, char *bnullarray, int *anynul,
          int *status);
static int imcomp_copy_overlap (char *tile, int pixlen, int ndim,
         long *tfpixel, long *tlpixel, char *bnullarray, char *image,
         long *fpixel, long *lpixel, long *inc, int nullcheck, char *nullarray,
         int *status);
static int imcomp_merge_overlap (char *tile, int pixlen, int ndim,
         long *tfpixel, long *tlpixel, char *bnullarray, char *image,
         long *fpixel, long *lpixel, int nullcheck, int *status);

/*--------------------------------------------------------------------------*/
static void ffxmsg( int action,
                    char *errmsg)
/*
  general routine to get or put an error message to the error buffer.

  Action Code:
   GetMesg    4  pop and return the message
   PutMesg    5  add a new message to the buffer

*/
{
    static char errbuff[81];  /* error buffer */
    static int firstTime = 1;

    if (firstTime)
    {
        firstTime = 0;
        errbuff[0] = '\0';
    }

    if (action == GetMesg)  /* pop and return oldest message from stack */
    {                            /* ignoring markers */
        strcpy(errmsg, errbuff);   /* copy message to output */
    }
    else if (action == PutMesg)  /* add new message to stack */
    {
        strncpy(errbuff, errmsg, 80);
    }

    return;
}
/*--------------------------------------------------------------------------*/
void _pyfits_ffpmsg(const char *err_message)
/*
  put message in error buffer
*/
{
    ffxmsg(PutMesg, (char *)err_message);
    return;
}
/*--------------------------------------------------------------------------*/
int _pyfits_ffgmsg(char *err_message)
/*
  get the error message from the error buffer
*/
{
    ffxmsg(GetMesg, err_message);
    return(*err_message);
}
/*--------------------------------------------------------------------------*/
static int ffpxsz(int datatype)
/*
   return the number of bytes per pixel associated with the datatype
*/
{
    if (datatype == TBYTE)
       return(sizeof(char));
    else if (datatype == TUSHORT)
       return(sizeof(short));
    else if (datatype == TSHORT)
       return(sizeof(short));
    else if (datatype == TULONG)
       return(sizeof(long));
    else if (datatype == TLONG)
       return(sizeof(long));
    else if (datatype == TINT)
       return(sizeof(int));
    else if (datatype == TUINT)
       return(sizeof(int));
    else if (datatype == TFLOAT)
       return(sizeof(float));
    else if (datatype == TDOUBLE)
       return(sizeof(double));
    else if (datatype == TLOGICAL)
       return(sizeof(char));
    else
       return(0);
}

/*****************************************************************************/
/*                                                                           */
/* The following code was copied and modified from the FITSIO source code    */
/* file swapproc.c.                                                          */
/*                                                                           */
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
static void ffswap2(short *svalues,  /* IO - pointer to shorts to be swapped */
                    long nvals)      /* I  - number of shorts to be swapped  */
/*
  swap the bytes in the input short integers: ( 0 1 -> 1 0 )
*/
{
    register char *cvalues;
    register long ii;

    union u_tag {
        char cvals[2];   /* equivalence an array of 4 bytes with */
        short sval;      /* a short */
    } u;

    cvalues = (char *) svalues;      /* copy the initial pointer value */

    for (ii = 0; ii < nvals;)
    {
        u.sval = svalues[ii++];  /* copy next short to temporary buffer */

        *cvalues++ = u.cvals[1]; /* copy the 2 bytes to output in turn */
        *cvalues++ = u.cvals[0];
    }
    return;
}
/*--------------------------------------------------------------------------*/
static void ffswap4(INT32BIT *ivalues, /* IO - pointer to floats to be swapped*/
                    long nvals)        /* I  - number of floats to be swapped */
/*
  swap the bytes in the input 4-byte integer: ( 0 1 2 3 -> 3 2 1 0 )
*/
{
    register char *cvalues;
    register long ii;

    union u_tag {
        char cvals[4];      /* equivalence an array of 4 bytes with */
        INT32BIT ival;      /* a float */
    } u;

    cvalues = (char *) ivalues;   /* copy the initial pointer value */

    for (ii = 0; ii < nvals;)
    {
        u.ival = ivalues[ii++];  /* copy next float to buffer */

        *cvalues++ = u.cvals[3]; /* copy the 4 bytes in turn */
        *cvalues++ = u.cvals[2];
        *cvalues++ = u.cvals[1];
        *cvalues++ = u.cvals[0];
    }
    return;
}

/*****************************************************************************/
/*                                                                           */
/* The following code was copied and modified from the FITSIO source code    */
/* file getcolb.c.                                                           */
/*                                                                           */
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
static int ffgpvb(
            fitsfile *fptr,          /* I - FITS file pointer                */
            long  group,      /* I - group to read (1 = 1st group)           */
            LONGLONG  firstelem,  /* I - first vector element to read        */
                                  /* (1 = 1st)                               */
            LONGLONG  nelem,      /* I - number of values to read            */
            unsigned char nulval, /* I - value for undefined pixels          */
            unsigned char *array, /* O - array of values that are returned   */
            int  *anynul,     /* O - set to 1 if any values are null; else 0 */
            int  *status)     /* IO - error status                           */
/*
  Read an array of values from the primary array. Data conversion
  and scaling will be performed if necessary (e.g, if the datatype of
  the FITS array is not the same as the array being read).
  Undefined elements will be set equal to NULVAL, unless NULVAL=0
  in which case no checking for undefined values will be performed.
  ANYNUL is returned with a value of .true. if any pixels are undefined.
*/
{
    int nullcheck = 1;
    unsigned char nullvalue;

    nullvalue = nulval;  /* set local variable */

    fits_read_compressed_pixels(fptr, TBYTE, firstelem, nelem,
            nullcheck, &nullvalue, array, NULL, anynul, status);
    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fffi1i1(
            unsigned char *input, /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            unsigned char tnull,  /* I - value of FITS TNULLn keyword if any */
            unsigned char nullval,/* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            unsigned char *output,/* O - array of converted pixels           */
            int *status)          /* IO - error status                       */
/*
  Copy input to output following reading of the input from a FITS file.
  Check for null values and do datatype conversion and scaling if required.
  The nullcheck code value determines how any null values in the input array
  are treated.  A null value is an input pixel that is equal to tnull.  If
  nullcheck = 0, then no checking for nulls is performed and any null values
  will be transformed just like any other pixel.  If nullcheck = 1, then the
  output pixel will be set = nullval if the corresponding input pixel is null.
  If nullcheck = 2, then if the pixel is null then the corresponding value of
  nullarray will be set to 1; the value of nullarray for non-null pixels
  will = 0.  The anynull parameter will be set = 1 if any of the returned
  pixels are null, otherwise anynull will be returned with a value = 0;
*/
{
    long ii;
    double dvalue;

    if (nullcheck == 0)     /* no null checking required */
    {
        if (scale == 1. && zero == 0.)      /* no scaling */
        {              /* this routine is normally not called in this case */
           memcpy(output, input, ntodo );
        }
        else             /* must scale the data */
        {
            for (ii = 0; ii < ntodo; ii++)
            {
                dvalue = input[ii] * scale + zero;

                if (dvalue < DUCHAR_MIN)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = 0;
                }
                else if (dvalue > DUCHAR_MAX)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = UCHAR_MAX;
                }
                else
                    output[ii] = (unsigned char) dvalue;
            }
        }
    }
    else        /* must check for null values */
    {
        if (scale == 1. && zero == 0.)  /* no scaling */
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
                    output[ii] = input[ii];
            }
        }
        else                  /* must scale the data */
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
                    dvalue = input[ii] * scale + zero;

                    if (dvalue < DUCHAR_MIN)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = 0;
                    }
                    else if (dvalue > DUCHAR_MAX)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = UCHAR_MAX;
                    }
                    else
                        output[ii] = (unsigned char) dvalue;
                }
            }
        }
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fffi2i1(
            short *input,         /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            short tnull,          /* I - value of FITS TNULLn keyword if any */
            unsigned char nullval,/* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            unsigned char *output,/* O - array of converted pixels           */
            int *status)          /* IO - error status                       */
/*
  Copy input to output following reading of the input from a FITS file.
  Check for null values and do datatype conversion and scaling if required.
  The nullcheck code value determines how any null values in the input array
  are treated.  A null value is an input pixel that is equal to tnull.  If
  nullcheck = 0, then no checking for nulls is performed and any null values
  will be transformed just like any other pixel.  If nullcheck = 1, then the
  output pixel will be set = nullval if the corresponding input pixel is null.
  If nullcheck = 2, then if the pixel is null then the corresponding value of
  nullarray will be set to 1; the value of nullarray for non-null pixels
  will = 0.  The anynull parameter will be set = 1 if any of the returned
  pixels are null, otherwise anynull will be returned with a value = 0;
*/
{
    long ii;
    double dvalue;

    if (nullcheck == 0)     /* no null checking required */
    {
        if (scale == 1. && zero == 0.)      /* no scaling */
        {
            for (ii = 0; ii < ntodo; ii++)
            {
                if (input[ii] < 0)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = 0;
                }
                else if (input[ii] > UCHAR_MAX)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = UCHAR_MAX;
                }
                else
                    output[ii] = (unsigned char) input[ii];
            }
        }
        else             /* must scale the data */
        {
            for (ii = 0; ii < ntodo; ii++)
            {
                dvalue = input[ii] * scale + zero;

                if (dvalue < DUCHAR_MIN)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = 0;
                }
                else if (dvalue > DUCHAR_MAX)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = UCHAR_MAX;
                }
                else
                    output[ii] = (unsigned char) dvalue;
            }
        }
    }
    else        /* must check for null values */
    {
        if (scale == 1. && zero == 0.)  /* no scaling */
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
                    if (input[ii] < 0)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = 0;
                    }
                    else if (input[ii] > UCHAR_MAX)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = UCHAR_MAX;
                    }
                    else
                        output[ii] = (unsigned char) input[ii];
                }
            }
        }
        else                  /* must scale the data */
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
                    dvalue = input[ii] * scale + zero;

                    if (dvalue < DUCHAR_MIN)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = 0;
                    }
                    else if (dvalue > DUCHAR_MAX)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = UCHAR_MAX;
                    }
                    else
                        output[ii] = (unsigned char) dvalue;
                }
            }
        }
    }
    return(*status);
}
/*---------------------------------------------------------------------------*/
static int fffi4i1(
            INT32BIT *input,      /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            INT32BIT tnull,       /* I - value of FITS TNULLn keyword if any */
            unsigned char nullval,/* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            unsigned char *output,/* O - array of converted pixels           */
            int *status)          /* IO - error status                       */
/*
  Copy input to output following reading of the input from a FITS file.
  Check for null values and do datatype conversion and scaling if required.
  The nullcheck code value determines how any null values in the input array
  are treated.  A null value is an input pixel that is equal to tnull.  If
  nullcheck = 0, then no checking for nulls is performed and any null values
  will be transformed just like any other pixel.  If nullcheck = 1, then the
  output pixel will be set = nullval if the corresponding input pixel is null.
  If nullcheck = 2, then if the pixel is null then the corresponding value of
  nullarray will be set to 1; the value of nullarray for non-null pixels
  will = 0.  The anynull parameter will be set = 1 if any of the returned
  pixels are null, otherwise anynull will be returned with a value = 0;
*/
{
    long ii;
    double dvalue;

    if (nullcheck == 0)     /* no null checking required */
    {
        if (scale == 1. && zero == 0.)      /* no scaling */
        {
            for (ii = 0; ii < ntodo; ii++)
            {
                if (input[ii] < 0)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = 0;
                }
                else if (input[ii] > UCHAR_MAX)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = UCHAR_MAX;
                }
                else
                    output[ii] = (unsigned char) input[ii];
            }
        }
        else             /* must scale the data */
        {
            for (ii = 0; ii < ntodo; ii++)
            {
                dvalue = input[ii] * scale + zero;

                if (dvalue < DUCHAR_MIN)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = 0;
                }
                else if (dvalue > DUCHAR_MAX)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = UCHAR_MAX;
                }
                else
                    output[ii] = (unsigned char) dvalue;
            }
        }
    }
    else        /* must check for null values */
    {
        if (scale == 1. && zero == 0.)  /* no scaling */
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
                    if (input[ii] < 0)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = 0;
                    }
                    else if (input[ii] > UCHAR_MAX)
                    {
                        output[ii] = UCHAR_MAX;
                    }
                    else
                        output[ii] = (unsigned char) input[ii];
                }
            }
        }
        else                  /* must scale the data */
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
                    dvalue = input[ii] * scale + zero;

                    if (dvalue < DUCHAR_MIN)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = 0;
                    }
                    else if (dvalue > DUCHAR_MAX)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = UCHAR_MAX;
                    }
                    else
                        output[ii] = (unsigned char) dvalue;
                }
            }
        }
    }
    return(*status);
}

/*****************************************************************************/
/*                                                                           */
/* The following code was copied and modified from the FITSIO source code    */
/* file getcold.c.                                                           */
/*                                                                           */
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
static int ffgpvd( 
            fitsfile *fptr,        /* I - FITS file pointer                  */
            long  group,      /* I - group to read (1 = 1st group)           */
            LONGLONG  firstelem,  /* I - first vector element to read        */
                                  /* (1 = 1st)                               */
            LONGLONG  nelem,      /* I - number of values to read            */
            double nulval,    /* I - value for undefined pixels              */
            double *array,    /* O - array of values that are returned       */
            int  *anynul,     /* O - set to 1 if any values are null; else 0 */
            int  *status)     /* IO - error status                           */
/*
  Read an array of values from the primary array. Data conversion
  and scaling will be performed if necessary (e.g, if the datatype of
  the FITS array is not the same as the array being read).
  Undefined elements will be set equal to NULVAL, unless NULVAL=0
  in which case no checking for undefined values will be performed.
  ANYNUL is returned with a value of .true. if any pixels are undefined.
*/
{
    int nullcheck = 1;
    double nullvalue;

    nullvalue = nulval;  /* set local variable */

    fits_read_compressed_pixels(fptr, TDOUBLE, firstelem, nelem,
            nullcheck, &nullvalue, array, NULL, anynul, status);
    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fffi1r8(
            unsigned char *input, /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            unsigned char tnull,  /* I - value of FITS TNULLn keyword if any */
            double nullval,       /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            double *output,       /* O - array of converted pixels           */
            int *status)          /* IO - error status                       */
/*
  Copy input to output following reading of the input from a FITS file.
  Check for null values and do datatype conversion and scaling if required.
  The nullcheck code value determines how any null values in the input array
  are treated.  A null value is an input pixel that is equal to tnull.  If
  nullcheck = 0, then no checking for nulls is performed and any null values
  will be transformed just like any other pixel.  If nullcheck = 1, then the
  output pixel will be set = nullval if the corresponding input pixel is null.
  If nullcheck = 2, then if the pixel is null then the corresponding value of
  nullarray will be set to 1; the value of nullarray for non-null pixels
  will = 0.  The anynull parameter will be set = 1 if any of the returned
  pixels are null, otherwise anynull will be returned with a value = 0;
*/
{
    long ii;

    if (nullcheck == 0)     /* no null checking required */
    {
        if (scale == 1. && zero == 0.)      /* no scaling */
        {
            for (ii = 0; ii < ntodo; ii++)
                output[ii] = (double) input[ii]; /* copy input to output */
        }
        else             /* must scale the data */
        {
            for (ii = 0; ii < ntodo; ii++)
            {
                output[ii] = input[ii] * scale + zero;
            }
        }
    }
    else        /* must check for null values */
    {
        if (scale == 1. && zero == 0.)  /* no scaling */
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
                    output[ii] = (double) input[ii];
            }
        }
        else                  /* must scale the data */
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
                    output[ii] = input[ii] * scale + zero;
                }
            }
        }
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fffi2r8(
            short *input,         /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            short tnull,          /* I - value of FITS TNULLn keyword if any */
            double nullval,       /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            double *output,       /* O - array of converted pixels           */
            int *status)          /* IO - error status                       */
/*
  Copy input to output following reading of the input from a FITS file.
  Check for null values and do datatype conversion and scaling if required.
  The nullcheck code value determines how any null values in the input array
  are treated.  A null value is an input pixel that is equal to tnull.  If
  nullcheck = 0, then no checking for nulls is performed and any null values
  will be transformed just like any other pixel.  If nullcheck = 1, then the
  output pixel will be set = nullval if the corresponding input pixel is null.
  If nullcheck = 2, then if the pixel is null then the corresponding value of
  nullarray will be set to 1; the value of nullarray for non-null pixels
  will = 0.  The anynull parameter will be set = 1 if any of the returned
  pixels are null, otherwise anynull will be returned with a value = 0;
*/
{
    long ii;

    if (nullcheck == 0)     /* no null checking required */
    {
        if (scale == 1. && zero == 0.)      /* no scaling */
        {
            for (ii = 0; ii < ntodo; ii++)
                output[ii] = (double) input[ii]; /* copy input to output */
        }
        else             /* must scale the data */
        {
            for (ii = 0; ii < ntodo; ii++)
            {
                output[ii] = input[ii] * scale + zero;
            }
        }
    }
    else        /* must check for null values */
    {
        if (scale == 1. && zero == 0.)  /* no scaling */
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
                    output[ii] = (double) input[ii];
            }
        }
        else                  /* must scale the data */
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
                    output[ii] = input[ii] * scale + zero;
                }
            }
        }
    }
    return(*status);
}
/*---------------------------------------------------------------------------*/
static int fffi4r8(
            INT32BIT *input,      /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            INT32BIT tnull,       /* I - value of FITS TNULLn keyword if any */
            double nullval,       /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            double *output,       /* O - array of converted pixels           */
            int *status)          /* IO - error status                       */
/*
  Copy input to output following reading of the input from a FITS file.
  Check for null values and do datatype conversion and scaling if required.
  The nullcheck code value determines how any null values in the input array
  are treated.  A null value is an input pixel that is equal to tnull.  If
  nullcheck = 0, then no checking for nulls is performed and any null values
  will be transformed just like any other pixel.  If nullcheck = 1, then the
  output pixel will be set = nullval if the corresponding input pixel is null.
  If nullcheck = 2, then if the pixel is null then the corresponding value of
  nullarray will be set to 1; the value of nullarray for non-null pixels
  will = 0.  The anynull parameter will be set = 1 if any of the returned
  pixels are null, otherwise anynull will be returned with a value = 0;
*/
{
    long ii;

    if (nullcheck == 0)     /* no null checking required */
    {
        if (scale == 1. && zero == 0.)      /* no scaling */
        {
            for (ii = 0; ii < ntodo; ii++)
                output[ii] = (double) input[ii]; /* copy input to output */
        }
        else             /* must scale the data */
        {
            for (ii = 0; ii < ntodo; ii++)
            {
                output[ii] = input[ii] * scale + zero;
            }
        }
    }
    else        /* must check for null values */
    {
        if (scale == 1. && zero == 0.)  /* no scaling */
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
                    output[ii] = (double) input[ii];
            }
        }
        else                  /* must scale the data */
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
                    output[ii] = input[ii] * scale + zero;
                }
            }
        }
    }
    return(*status);
}

/*****************************************************************************/
/*                                                                           */
/* The following code was copied and modified from the FITSIO source code    */
/* file getcole.c.                                                           */
/*                                                                           */
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
static int ffgpve(
            fitsfile *fptr,   /* I - FITS file pointer                       */
            long  group,      /* I - group to read (1 = 1st group)           */
            LONGLONG  firstelem,  /* I - first vector element to read        */
                                  /* (1 = 1st)                               */
            LONGLONG  nelem,      /* I - number of values to read            */
            float nulval,     /* I - value for undefined pixels              */
            float *array,     /* O - array of values that are returned       */
            int  *anynul,     /* O - set to 1 if any values are null; else 0 */
            int  *status)     /* IO - error status                           */
/*
  Read an array of values from the primary array. Data conversion
  and scaling will be performed if necessary (e.g, if the datatype of
  the FITS array is not the same as the array being read).
  Undefined elements will be set equal to NULVAL, unless NULVAL=0
  in which case no checking for undefined values will be performed.
  ANYNUL is returned with a value of .true. if any pixels are undefined.
*/
{
    int nullcheck = 1;
    float nullvalue;

    nullvalue = nulval;  /* set local variable */

    fits_read_compressed_pixels(fptr, TFLOAT, firstelem, nelem,
            nullcheck, &nullvalue, array, NULL, anynul, status);
    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fffi1r4(
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
  Copy input to output following reading of the input from a FITS file.
  Check for null values and do datatype conversion and scaling if required.
  The nullcheck code value determines how any null values in the input array
  are treated.  A null value is an input pixel that is equal to tnull.  If
  nullcheck = 0, then no checking for nulls is performed and any null values
  will be transformed just like any other pixel.  If nullcheck = 1, then the
  output pixel will be set = nullval if the corresponding input pixel is null.
  If nullcheck = 2, then if the pixel is null then the corresponding value of
  nullarray will be set to 1; the value of nullarray for non-null pixels
  will = 0.  The anynull parameter will be set = 1 if any of the returned
  pixels are null, otherwise anynull will be returned with a value = 0;
*/
{
    long ii;

    if (nullcheck == 0)     /* no null checking required */
    {
        if (scale == 1. && zero == 0.)      /* no scaling */
        {
            for (ii = 0; ii < ntodo; ii++)
                output[ii] = (float) input[ii];  /* copy input to output */
        }
        else             /* must scale the data */
        {
            for (ii = 0; ii < ntodo; ii++)
            {
                output[ii] = (float) (( (double) input[ii] ) * scale + zero);
            }
        }
    }
    else        /* must check for null values */
    {
        if (scale == 1. && zero == 0.)  /* no scaling */
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
                    output[ii] = (float) input[ii];
            }
        }
        else                  /* must scale the data */
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
                    output[ii] = (float) (( (double) input[ii] ) * scale + zero);
                }
            }
        }
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fffi2r4(
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
  Copy input to output following reading of the input from a FITS file.
  Check for null values and do datatype conversion and scaling if required.
  The nullcheck code value determines how any null values in the input array
  are treated.  A null value is an input pixel that is equal to tnull.  If
  nullcheck = 0, then no checking for nulls is performed and any null values
  will be transformed just like any other pixel.  If nullcheck = 1, then the
  output pixel will be set = nullval if the corresponding input pixel is null.
  If nullcheck = 2, then if the pixel is null then the corresponding value of
  nullarray will be set to 1; the value of nullarray for non-null pixels
  will = 0.  The anynull parameter will be set = 1 if any of the returned
  pixels are null, otherwise anynull will be returned with a value = 0;
*/
{
    long ii;

    if (nullcheck == 0)     /* no null checking required */
    {
        if (scale == 1. && zero == 0.)      /* no scaling */
        {
            for (ii = 0; ii < ntodo; ii++)
                output[ii] = (float) input[ii];  /* copy input to output */
        }
        else             /* must scale the data */
        {
            for (ii = 0; ii < ntodo; ii++)
            {
                output[ii] = (float) (input[ii] * scale + zero);
            }
        }
    }
    else        /* must check for null values */
    {
        if (scale == 1. && zero == 0.)  /* no scaling */
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
                    output[ii] = (float) input[ii];
            }
        }
        else                  /* must scale the data */
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
                    output[ii] = (float) (input[ii] * scale + zero);
                }
            }
        }
    }
    return(*status);
}
/*---------------------------------------------------------------------------*/
static int fffi4r4(
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
  Copy input to output following reading of the input from a FITS file.
  Check for null values and do datatype conversion and scaling if required.
  The nullcheck code value determines how any null values in the input array
  are treated.  A null value is an input pixel that is equal to tnull.  If
  nullcheck = 0, then no checking for nulls is performed and any null values
  will be transformed just like any other pixel.  If nullcheck = 1, then the
  output pixel will be set = nullval if the corresponding input pixel is null.
  If nullcheck = 2, then if the pixel is null then the corresponding value of
  nullarray will be set to 1; the value of nullarray for non-null pixels
  will = 0.  The anynull parameter will be set = 1 if any of the returned
  pixels are null, otherwise anynull will be returned with a value = 0;
*/
{
    long ii;

    if (nullcheck == 0)     /* no null checking required */
    {
        if (scale == 1. && zero == 0.)      /* no scaling */
        {
            for (ii = 0; ii < ntodo; ii++)
                output[ii] = (float) input[ii];  /* copy input to output */
        }
        else             /* must scale the data */
        {
            for (ii = 0; ii < ntodo; ii++)
            {
                output[ii] = (float) (input[ii] * scale + zero);
            }
        }
    }
    else        /* must check for null values */
    {
        if (scale == 1. && zero == 0.)  /* no scaling */
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
                    output[ii] = (float) input[ii];
            }
        }
        else                  /* must scale the data */
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
                    output[ii] = (float) (input[ii] * scale + zero);
                }
            }
        }
    }
    return(*status);
}

/*****************************************************************************/
/*                                                                           */
/* The following code was copied and modified from the FITSIO source code    */
/* file getcoli.c.                                                           */
/*                                                                           */
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
static int ffgpvi(
            fitsfile *fptr,   /* I - FITS file pointer                       */
            long  group,      /* I - group to read (1 = 1st group)           */
            LONGLONG  firstelem,  /* I - first vector element to read        */
                                  /* (1 = 1st)                               */
            LONGLONG  nelem,      /* I - number of values to read            */
            short nulval,     /* I - value for undefined pixels              */
            short *array,     /* O - array of values that are returned       */
            int  *anynul,     /* O - set to 1 if any values are null; else 0 */
            int  *status)     /* IO - error status                           */
/*
  Read an array of values from the primary array. Data conversion
  and scaling will be performed if necessary (e.g, if the datatype of
  the FITS array is not the same as the array being read).
  Undefined elements will be set equal to NULVAL, unless NULVAL=0
  in which case no checking for undefined values will be performed.
  ANYNUL is returned with a value of .true. if any pixels are undefined.
*/
{
    int nullcheck = 1;
    short nullvalue;

    nullvalue = nulval;  /* set local variable */
    
    fits_read_compressed_pixels(fptr, TSHORT, firstelem, nelem,
            nullcheck, &nullvalue, array, NULL, anynul, status);
    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fffi1i2(
            unsigned char *input, /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            unsigned char tnull,  /* I - value of FITS TNULLn keyword if any */
            short nullval,        /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            short *output,        /* O - array of converted pixels           */
            int *status)          /* IO - error status                       */
/*
  Copy input to output following reading of the input from a FITS file.
  Check for null values and do datatype conversion and scaling if required.
  The nullcheck code value determines how any null values in the input array
  are treated.  A null value is an input pixel that is equal to tnull.  If
  nullcheck = 0, then no checking for nulls is performed and any null values
  will be transformed just like any other pixel.  If nullcheck = 1, then the
  output pixel will be set = nullval if the corresponding input pixel is null.
  If nullcheck = 2, then if the pixel is null then the corresponding value of
  nullarray will be set to 1; the value of nullarray for non-null pixels
  will = 0.  The anynull parameter will be set = 1 if any of the returned
  pixels are null, otherwise anynull will be returned with a value = 0;
*/
{
    long ii;
    double dvalue;

    if (nullcheck == 0)     /* no null checking required */
    {
        if (scale == 1. && zero == 0.)      /* no scaling */
        {
            for (ii = 0; ii < ntodo; ii++)
                output[ii] = (short) input[ii];  /* copy input to output */
        }
        else             /* must scale the data */
        {
            for (ii = 0; ii < ntodo; ii++)
            {
                dvalue = input[ii] * scale + zero;

                if (dvalue < DSHRT_MIN)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = SHRT_MIN;
                }
                else if (dvalue > DSHRT_MAX)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = SHRT_MAX;
                }
                else
                    output[ii] = (short) dvalue;
            }
        }
    }
    else        /* must check for null values */
    {
        if (scale == 1. && zero == 0.)  /* no scaling */
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
                    output[ii] = (short) input[ii];
            }
        }
        else                  /* must scale the data */
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
                    dvalue = input[ii] * scale + zero;

                    if (dvalue < DSHRT_MIN)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = SHRT_MIN;
                    }
                    else if (dvalue > DSHRT_MAX)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = SHRT_MAX;
                    }
                    else
                        output[ii] = (short) dvalue;
                }
            }
        }
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fffi2i2(
            short *input,         /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            short tnull,          /* I - value of FITS TNULLn keyword if any */
            short nullval,        /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            short *output,        /* O - array of converted pixels           */
            int *status)          /* IO - error status                       */
/*
  Copy input to output following reading of the input from a FITS file.
  Check for null values and do datatype conversion and scaling if required.
  The nullcheck code value determines how any null values in the input array
  are treated.  A null value is an input pixel that is equal to tnull.  If
  nullcheck = 0, then no checking for nulls is performed and any null values
  will be transformed just like any other pixel.  If nullcheck = 1, then the
  output pixel will be set = nullval if the corresponding input pixel is null.
  If nullcheck = 2, then if the pixel is null then the corresponding value of
  nullarray will be set to 1; the value of nullarray for non-null pixels
  will = 0.  The anynull parameter will be set = 1 if any of the returned
  pixels are null, otherwise anynull will be returned with a value = 0;
*/
{
    long ii;
    double dvalue;

    if (nullcheck == 0)     /* no null checking required */
    {
        if (scale == 1. && zero == 0.)      /* no scaling */
        {
            memcpy(output, input, ntodo * sizeof(short) );
        }
        else             /* must scale the data */
        {
            for (ii = 0; ii < ntodo; ii++)
            {
                dvalue = input[ii] * scale + zero;

                if (dvalue < DSHRT_MIN)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = SHRT_MIN;
                }
                else if (dvalue > DSHRT_MAX)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = SHRT_MAX;
                }
                else
                    output[ii] = (short) dvalue;
            }
        }
    }
    else        /* must check for null values */
    {
        if (scale == 1. && zero == 0.)  /* no scaling */
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
                    output[ii] = input[ii];
            }
        }
        else                  /* must scale the data */
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
                    dvalue = input[ii] * scale + zero;

                    if (dvalue < DSHRT_MIN)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = SHRT_MIN;
                    }
                    else if (dvalue > DSHRT_MAX)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = SHRT_MAX;
                    }
                    else
                        output[ii] = (short) dvalue;
                }
            }
        }
    }
    return(*status);
}
/*---------------------------------------------------------------------------*/
static int fffi4i2(
            INT32BIT *input,      /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            INT32BIT tnull,       /* I - value of FITS TNULLn keyword if any */
            short nullval,        /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            short *output,        /* O - array of converted pixels           */
            int *status)          /* IO - error status                       */
/*
  Copy input to output following reading of the input from a FITS file.
  Check for null values and do datatype conversion and scaling if required.
  The nullcheck code value determines how any null values in the input array
  are treated.  A null value is an input pixel that is equal to tnull.  If
  nullcheck = 0, then no checking for nulls is performed and any null values
  will be transformed just like any other pixel.  If nullcheck = 1, then the
  output pixel will be set = nullval if the corresponding input pixel is null.
  If nullcheck = 2, then if the pixel is null then the corresponding value of
  nullarray will be set to 1; the value of nullarray for non-null pixels
  will = 0.  The anynull parameter will be set = 1 if any of the returned
  pixels are null, otherwise anynull will be returned with a value = 0;
*/
{
    long ii;
    double dvalue;

    if (nullcheck == 0)     /* no null checking required */
    {
        if (scale == 1. && zero == 0.)      /* no scaling */
        {
            for (ii = 0; ii < ntodo; ii++)
            {
                if (input[ii] < SHRT_MIN)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = SHRT_MIN;
                }
                else if (input[ii] > SHRT_MAX)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = SHRT_MAX;
                }
                else
                    output[ii] = (short) input[ii];
            }
        }
        else             /* must scale the data */
        {
            for (ii = 0; ii < ntodo; ii++)
            {
                dvalue = input[ii] * scale + zero;

                if (dvalue < DSHRT_MIN)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = SHRT_MIN;
                }
                else if (dvalue > DSHRT_MAX)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = SHRT_MAX;
                }
                else
                    output[ii] = (short) dvalue;
            }
        }
    }
    else        /* must check for null values */
    {
        if (scale == 1. && zero == 0.)  /* no scaling */
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
                    if (input[ii] < SHRT_MIN)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = SHRT_MIN;
                    }
                    else if (input[ii] > SHRT_MAX)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = SHRT_MAX;
                    }
                    else
                        output[ii] = (short) input[ii];
                }
            }
        }
        else                  /* must scale the data */
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
                    dvalue = input[ii] * scale + zero;

                    if (dvalue < DSHRT_MIN)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = SHRT_MIN;
                    }
                    else if (dvalue > DSHRT_MAX)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = SHRT_MAX;
                    }
                    else
                        output[ii] = (short) dvalue;
                }
            }
        }
    }
    return(*status);
}

/*****************************************************************************/
/*                                                                           */
/* The following code was copied and modified from the FITSIO source code    */
/* file getcolk.c.                                                           */
/*                                                                           */
/*****************************************************************************/

/*--------------------------------------------------------------------------*/
static int ffgpvk(
            fitsfile *fptr,   /* I - FITS file pointer                       */
            long  group,      /* I - group to read (1 = 1st group)           */
            LONGLONG  firstelem,  /* I - first vector element to read (1 = 1st)  */
            LONGLONG  nelem,      /* I - number of values to read                */
            int   nulval,     /* I - value for undefined pixels              */
            int   *array,     /* O - array of values that are returned       */
            int  *anynul,     /* O - set to 1 if any values are null; else 0 */
            int  *status)     /* IO - error status                           */
/*
  Read an array of values from the primary array. Data conversion
  and scaling will be performed if necessary (e.g, if the datatype of
  the FITS array is not the same as the array being read).
  Undefined elements will be set equal to NULVAL, unless NULVAL=0
  in which case no checking for undefined values will be performed.
  ANYNUL is returned with a value of .true. if any pixels are undefined.
*/
{
    int nullcheck = 1;
    int nullvalue;

    nullvalue = nulval;  /* set local variable */

    fits_read_compressed_pixels(fptr, TINT, firstelem, nelem,
            nullcheck, &nullvalue, array, NULL, anynul, status);
    return(*status);
}

/*--------------------------------------------------------------------------*/
static int fffi1int(
            unsigned char *input,/* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            unsigned char tnull,  /* I - value of FITS TNULLn keyword if any */
            int  nullval,         /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            int  *output,         /* O - array of converted pixels           */
            int *status)          /* IO - error status                       */
/*
  Copy input to output following reading of the input from a FITS file.
  Check for null values and do datatype conversion and scaling if required.
  The nullcheck code value determines how any null values in the input array
  are treated.  A null value is an input pixel that is equal to tnull.  If
  nullcheck = 0, then no checking for nulls is performed and any null values
  will be transformed just like any other pixel.  If nullcheck = 1, then the
  output pixel will be set = nullval if the corresponding input pixel is null.
  If nullcheck = 2, then if the pixel is null then the corresponding value of
  nullarray will be set to 1; the value of nullarray for non-null pixels
  will = 0.  The anynull parameter will be set = 1 if any of the returned
  pixels are null, otherwise anynull will be returned with a value = 0;
*/
{
    long ii;
    double dvalue;

    if (nullcheck == 0)     /* no null checking required */
    {
        if (scale == 1. && zero == 0.)      /* no scaling */
        {
            for (ii = 0; ii < ntodo; ii++)
                output[ii] = (int) input[ii];  /* copy input to output */
        }
        else             /* must scale the data */
        {
            for (ii = 0; ii < ntodo; ii++)
            {
                dvalue = input[ii] * scale + zero;

                if (dvalue < DINT_MIN)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = INT_MIN;
                }
                else if (dvalue > DINT_MAX)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = INT_MAX;
                }
                else
                    output[ii] = (int) dvalue;
            }
        }
    }
    else        /* must check for null values */
    {
        if (scale == 1. && zero == 0.)  /* no scaling */
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
                    output[ii] = (int) input[ii];
            }
        }
        else                  /* must scale the data */
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
                    dvalue = input[ii] * scale + zero;

                    if (dvalue < DINT_MIN)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = INT_MIN;
                    }
                    else if (dvalue > DINT_MAX)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = INT_MAX;
                    }
                    else
                        output[ii] = (int) dvalue;
                }
            }
        }
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fffi2int(
            short *input,        /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            short tnull,          /* I - value of FITS TNULLn keyword if any */
            int  nullval,         /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            int  *output,         /* O - array of converted pixels           */
            int *status)          /* IO - error status                       */
/*
  Copy input to output following reading of the input from a FITS file.
  Check for null values and do datatype conversion and scaling if required.
  The nullcheck code value determines how any null values in the input array
  are treated.  A null value is an input pixel that is equal to tnull.  If
  nullcheck = 0, then no checking for nulls is performed and any null values
  will be transformed just like any other pixel.  If nullcheck = 1, then the
  output pixel will be set = nullval if the corresponding input pixel is null.
  If nullcheck = 2, then if the pixel is null then the corresponding value of
  nullarray will be set to 1; the value of nullarray for non-null pixels
  will = 0.  The anynull parameter will be set = 1 if any of the returned
  pixels are null, otherwise anynull will be returned with a value = 0;
*/
{
    long ii;
    double dvalue;

    if (nullcheck == 0)     /* no null checking required */
    {
        if (scale == 1. && zero == 0.)      /* no scaling */
        {
            for (ii = 0; ii < ntodo; ii++)
                output[ii] = (int) input[ii];   /* copy input to output */
        }
        else             /* must scale the data */
        {
            for (ii = 0; ii < ntodo; ii++)
            {
                dvalue = input[ii] * scale + zero;

                if (dvalue < DINT_MIN)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = INT_MIN;
                }
                else if (dvalue > DINT_MAX)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = INT_MAX;
                }
                else
                    output[ii] = (int) dvalue;
            }
        }
    }
    else        /* must check for null values */
    {
        if (scale == 1. && zero == 0.)  /* no scaling */
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
                    output[ii] = (int) input[ii];
            }
        }
        else                  /* must scale the data */
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
                    dvalue = input[ii] * scale + zero;

                    if (dvalue < DINT_MIN)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = INT_MIN;
                    }
                    else if (dvalue > DINT_MAX)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = INT_MAX;
                    }
                    else
                        output[ii] = (int) dvalue;
                }
            }
        }
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fffi4int(
            INT32BIT *input,     /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            INT32BIT tnull,       /* I - value of FITS TNULLn keyword if any */
            int  nullval,         /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            int  *output,         /* O - array of converted pixels           */
            int *status)          /* IO - error status                       */
/*
  Copy input to output following reading of the input from a FITS file.
  Check for null values and do datatype conversion and scaling if required.
  The nullcheck code value determines how any null values in the input array
  are treated.  A null value is an input pixel that is equal to tnull.  If
  nullcheck = 0, then no checking for nulls is performed and any null values
  will be transformed just like any other pixel.  If nullcheck = 1, then the
  output pixel will be set = nullval if the corresponding input pixel is null.
  If nullcheck = 2, then if the pixel is null then the corresponding value of
  nullarray will be set to 1; the value of nullarray for non-null pixels
  will = 0.  The anynull parameter will be set = 1 if any of the returned
  pixels are null, otherwise anynull will be returned with a value = 0;
*/
{
    long ii;
    double dvalue;
    if (nullcheck == 0)     /* no null checking required */
    {
        if (scale == 1. && zero == 0.)      /* no scaling */
        {
            for (ii = 0; ii < ntodo; ii++)
                output[ii] = (int) input[ii];   /* copy input to output */
        }
        else             /* must scale the data */
        {
            for (ii = 0; ii < ntodo; ii++)
            {
                dvalue = input[ii] * scale + zero;

                if (dvalue < DINT_MIN)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = INT_MIN;
                }
                else if (dvalue > DINT_MAX)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = INT_MAX;
                }
                else
                    output[ii] = (int) dvalue;
            }
        }
    }
    else        /* must check for null values */
    {
        if (scale == 1. && zero == 0.)  /* no scaling */
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
                    output[ii] = (int) input[ii];

            }
        }
        else                  /* must scale the data */
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
                    dvalue = input[ii] * scale + zero;

                    if (dvalue < DINT_MIN)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = INT_MIN;
                    }
                    else if (dvalue > DINT_MAX)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = INT_MAX;
                    }
                    else
                        output[ii] = (int) dvalue;
                }
            }
        }
    }
    return(*status);
}

/*****************************************************************************/
/*                                                                           */
/* The following code was copied and modified from the FITSIO source code    */
/* file getcolj.c.                                                           */
/*                                                                           */
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
static int ffgpvj(
            fitsfile *fptr,   /* I - FITS file pointer                       */
            long  group,      /* I - group to read (1 = 1st group)           */
            LONGLONG  firstelem,  /* I - first vector element to read        */
                                  /* (1 = 1st)                               */
            LONGLONG  nelem,      /* I - number of values to read            */
            long  nulval,     /* I - value for undefined pixels              */
            long  *array,     /* O - array of values that are returned       */
            int  *anynul,     /* O - set to 1 if any values are null; else 0 */
            int  *status)     /* IO - error status                           */
/*
  Read an array of values from the primary array. Data conversion
  and scaling will be performed if necessary (e.g, if the datatype of
  the FITS array is not the same as the array being read).
  Undefined elements will be set equal to NULVAL, unless NULVAL=0
  in which case no checking for undefined values will be performed.
  ANYNUL is returned with a value of .true. if any pixels are undefined.
*/
{
    int nullcheck = 1;
    long nullvalue;

    nullvalue = nulval;  /* set local variable */

    fits_read_compressed_pixels(fptr, TLONG, firstelem, nelem,
            nullcheck, &nullvalue, array, NULL, anynul, status);
    return(*status);
}

/*--------------------------------------------------------------------------*/
static int fffi1i4(
            unsigned char *input, /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            unsigned char tnull,  /* I - value of FITS TNULLn keyword if any */
            long nullval,         /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            long *output,         /* O - array of converted pixels           */
            int *status)          /* IO - error status                       */
/*
  Copy input to output following reading of the input from a FITS file.
  Check for null values and do datatype conversion and scaling if required.
  The nullcheck code value determines how any null values in the input array
  are treated.  A null value is an input pixel that is equal to tnull.  If
  nullcheck = 0, then no checking for nulls is performed and any null values
  will be transformed just like any other pixel.  If nullcheck = 1, then the
  output pixel will be set = nullval if the corresponding input pixel is null.
  If nullcheck = 2, then if the pixel is null then the corresponding value of
  nullarray will be set to 1; the value of nullarray for non-null pixels
  will = 0.  The anynull parameter will be set = 1 if any of the returned
  pixels are null, otherwise anynull will be returned with a value = 0;
*/
{
    long ii;
    double dvalue;

    if (nullcheck == 0)     /* no null checking required */
    {
        if (scale == 1. && zero == 0.)      /* no scaling */
        {
            for (ii = 0; ii < ntodo; ii++)
                output[ii] = (long) input[ii];  /* copy input to output */
        }
        else             /* must scale the data */
        {
            for (ii = 0; ii < ntodo; ii++)
            {
                dvalue = input[ii] * scale + zero;

                if (dvalue < DLONG_MIN)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MIN;
                }
                else if (dvalue > DLONG_MAX)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MAX;
                }
                else
                    output[ii] = (long) dvalue;
            }
        }
    }
    else        /* must check for null values */
    {
        if (scale == 1. && zero == 0.)  /* no scaling */
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
                    output[ii] = (long) input[ii];
            }
        }
        else                  /* must scale the data */
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
                    dvalue = input[ii] * scale + zero;

                    if (dvalue < DLONG_MIN)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONG_MIN;
                    }
                    else if (dvalue > DLONG_MAX)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONG_MAX;
                    }
                    else
                        output[ii] = (long) dvalue;
                }
            }
        }
    }
    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fffi2i4(
            short *input,         /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            short tnull,          /* I - value of FITS TNULLn keyword if any */
            long nullval,         /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            long *output,         /* O - array of converted pixels           */
            int *status)          /* IO - error status                       */
/*
  Copy input to output following reading of the input from a FITS file.
  Check for null values and do datatype conversion and scaling if required.
  The nullcheck code value determines how any null values in the input array
  are treated.  A null value is an input pixel that is equal to tnull.  If
  nullcheck = 0, then no checking for nulls is performed and any null values
  will be transformed just like any other pixel.  If nullcheck = 1, then the
  output pixel will be set = nullval if the corresponding input pixel is null.
  If nullcheck = 2, then if the pixel is null then the corresponding value of
  nullarray will be set to 1; the value of nullarray for non-null pixels
  will = 0.  The anynull parameter will be set = 1 if any of the returned
  pixels are null, otherwise anynull will be returned with a value = 0;
*/
{
    long ii;
    double dvalue;

    if (nullcheck == 0)     /* no null checking required */
    {
        if (scale == 1. && zero == 0.)      /* no scaling */
        {
            for (ii = 0; ii < ntodo; ii++)
                output[ii] = (long) input[ii];   /* copy input to output */
        }
        else             /* must scale the data */
        {
            for (ii = 0; ii < ntodo; ii++)
            {
                dvalue = input[ii] * scale + zero;

                if (dvalue < DLONG_MIN)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MIN;
                }
                else if (dvalue > DLONG_MAX)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MAX;
                }
                else
                    output[ii] = (long) dvalue;
            }
        }
    }
    else        /* must check for null values */
    {
        if (scale == 1. && zero == 0.)  /* no scaling */
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
                    output[ii] = (long) input[ii];
            }
        }
        else                  /* must scale the data */
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
                    dvalue = input[ii] * scale + zero;

                    if (dvalue < DLONG_MIN)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONG_MIN;
                    }
                    else if (dvalue > DLONG_MAX)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONG_MAX;
                    }
                    else
                        output[ii] = (long) dvalue;
                }
            }
        }
    }
    return(*status);
}
/*---------------------------------------------------------------------------*/
static int fffi4i4(
            INT32BIT *input,      /* I - array of values to be converted     */
            long ntodo,           /* I - number of elements in the array     */
            double scale,         /* I - FITS TSCALn or BSCALE value         */
            double zero,          /* I - FITS TZEROn or BZERO  value         */
            int nullcheck,        /* I - null checking code; 0 = don't check */
                                  /*     1:set null pixels = nullval         */
                                  /*     2: if null pixel, set nullarray = 1 */
            INT32BIT tnull,       /* I - value of FITS TNULLn keyword if any */
            long nullval,         /* I - set null pixels, if nullcheck = 1   */
            char *nullarray,      /* I - bad pixel array, if nullcheck = 2   */
            int  *anynull,        /* O - set to 1 if any pixels are null     */
            long *output,         /* O - array of converted pixels           */
            int *status)          /* IO - error status                       */
/*
  Copy input to output following reading of the input from a FITS file.
  Check for null values and do datatype conversion and scaling if required.
  The nullcheck code value determines how any null values in the input array
  are treated.  A null value is an input pixel that is equal to tnull.  If
  nullcheck = 0, then no checking for nulls is performed and any null values
  will be transformed just like any other pixel.  If nullcheck = 1, then the
  output pixel will be set = nullval if the corresponding input pixel is null.
  If nullcheck = 2, then if the pixel is null then the corresponding value of
  nullarray will be set to 1; the value of nullarray for non-null pixels
  will = 0.  The anynull parameter will be set = 1 if any of the returned
  pixels are null, otherwise anynull will be returned with a value = 0;

  Process the array of data in reverse order, to handle the case where
  the input data is 4-bytes and the output is  8-bytes and the conversion
  is being done in place in the same array.
*/
{
    long ii;
    double dvalue;

    if (nullcheck == 0)     /* no null checking required */
    {
        if (scale == 1. && zero == 0.)      /* no scaling */
        {
            for (ii = ntodo - 1; ii >= 0; ii--)
                output[ii] = (long) input[ii];   /* copy input to output */
        }
        else             /* must scale the data */
        {
            for (ii = ntodo - 1; ii >= 0; ii--)
            {
                dvalue = input[ii] * scale + zero;

                if (dvalue < DLONG_MIN)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MIN;
                }
                else if (dvalue > DLONG_MAX)
                {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MAX;
                }
                else
                    output[ii] = (long) dvalue;
            }
        }
    }
    else        /* must check for null values */
    {
        if (scale == 1. && zero == 0.)  /* no scaling */
        {
            for (ii = ntodo - 1; ii >= 0; ii--)
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
                    output[ii] = input[ii];

            }
        }
        else                  /* must scale the data */
        {
            for (ii = ntodo - 1; ii >= 0; ii--)
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
                    dvalue = input[ii] * scale + zero;

                    if (dvalue < DLONG_MIN)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONG_MIN;
                    }
                    else if (dvalue > DLONG_MAX)
                    {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONG_MAX;
                    }
                    else
                        output[ii] = (long) dvalue;
                }
            }
        }
    }
    return(*status);
}
/*---------------------------------------------------------------------------*/
static int ffgpvjj(
            fitsfile *fptr,   /* I - FITS file pointer                       */
            long  group,      /* I - group to read (1 = 1st group)           */
            LONGLONG  firstelem,  /* I - first vector element to read        */
                                  /* (1 = 1st)                               */
            LONGLONG  nelem,      /* I - number of values to read            */
            LONGLONG  nulval, /* I - value for undefined pixels              */
            LONGLONG  *array, /* O - array of values that are returned       */
            int  *anynul,     /* O - set to 1 if any values are null; else 0 */
            int  *status)     /* IO - error status                           */
/*
  Read an array of values from the primary array. Data conversion
  and scaling will be performed if necessary (e.g, if the datatype of
  the FITS array is not the same as the array being read).
  Undefined elements will be set equal to NULVAL, unless NULVAL=0
  in which case no checking for undefined values will be performed.
  ANYNUL is returned with a value of .true. if any pixels are undefined.
*/
{
    int nullcheck = 1;
    LONGLONG nullvalue;

    nullvalue = nulval;  /* set local variable */

    fits_read_compressed_pixels(fptr, TLONGLONG, firstelem, nelem,
            nullcheck, &nullvalue, array, NULL, anynul, status);
    return(*status);
}

/*****************************************************************************/
/*                                                                           */
/* The following code was copied and modified from the FITSIO source code    */
/* file getcol.c.                                                            */
/*                                                                           */
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
int _pyfits_ffgpv(
            fitsfile *fptr,   /* I - FITS file pointer                       */
            int  datatype,    /* I - datatype of the value                   */
            LONGLONG firstelem,   /* I - first vector element to read        */
                                  /* (1 = 1st)                               */
            LONGLONG nelem,       /* I - number of values to read            */
            void *nulval,     /* I - value for undefined pixels              */
            void *array,      /* O - array of values that are returned       */
            int  *anynul,     /* O - set to 1 if any values are null; else 0 */
            int  *status)     /* IO - error status                           */
/*
  Read an array of values from the primary array. The datatype of the
  input array is defined by the 2nd argument.  Data conversion
  and scaling will be performed if necessary (e.g, if the datatype of
  the FITS array is not the same as the array being read).
  Undefined elements will be set equal to NULVAL, unless NULVAL=0
  in which case no checking for undefined values will be performed.
  ANYNUL is returned with a value of .true. if any pixels are undefined.
*/
{

    if (*status > 0 || nelem == 0)   /* inherit input status value if > 0 */
        return(*status);

    /*
      the primary array is represented as a binary table:
      each group of the primary array is a row in the table,
      where the first column contains the group parameters
      and the second column contains the image itself.
    */

    if (datatype == TBYTE)
    {
      if (nulval == 0)
        ffgpvb(fptr, 1, firstelem, nelem, 0,
               (unsigned char *) array, anynul, status);
      else
        ffgpvb(fptr, 1, firstelem, nelem, *(unsigned char *) nulval,
               (unsigned char *) array, anynul, status);
    }
    else if (datatype == TSHORT)
    {
      if (nulval == 0)
        ffgpvi(fptr, 1, firstelem, nelem, 0,
               (short *) array, anynul, status);
      else
        ffgpvi(fptr, 1, firstelem, nelem, *(short *) nulval,
               (short *) array, anynul, status);
    }
    else if (datatype == TINT)
    {
      if (nulval == 0)
        ffgpvk(fptr, 1, firstelem, nelem, 0,
               (int *) array, anynul, status);
      else
        ffgpvk(fptr, 1, firstelem, nelem, *(int *) nulval,
               (int *) array, anynul, status);
    }
    else if (datatype == TLONG)
    {
      if (nulval == 0)
        ffgpvj(fptr, 1, firstelem, nelem, 0,
               (long *) array, anynul, status);
      else
        ffgpvj(fptr, 1, firstelem, nelem, *(long *) nulval,
               (long *) array, anynul, status);
    }
    else if (datatype == TLONGLONG)
    {
      if (nulval == 0)
        ffgpvjj(fptr, 1, firstelem, nelem, 0,
               (LONGLONG *) array, anynul, status);
      else
        ffgpvjj(fptr, 1, firstelem, nelem, *(LONGLONG *) nulval,
               (LONGLONG *) array, anynul, status);
    }
    else if (datatype == TFLOAT)
    {
      if (nulval == 0)
        ffgpve(fptr, 1, firstelem, nelem, 0,
               (float *) array, anynul, status);
      else
        ffgpve(fptr, 1, firstelem, nelem, *(float *) nulval,
               (float *) array, anynul, status);
    }
    else if (datatype == TDOUBLE)
    {
      if (nulval == 0)
        ffgpvd(fptr, 1, firstelem, nelem, 0,
               (double *) array, anynul, status);
      else
      {
        ffgpvd(fptr, 1, firstelem, nelem, *(double *) nulval,
               (double *) array, anynul, status);
      }
    }
    else
      *status = BAD_DATATYPE;

    return(*status);
}

/*****************************************************************************/
/*                                                                           */
/* The following code was copied and modified from the FITSIO source code    */
/* file putcolb.c.                                                           */
/*                                                                           */
/*****************************************************************************/

/*--------------------------------------------------------------------------*/
static int ffpprb(
            fitsfile *fptr,  /* I - FITS file pointer                       */
            long  group,     /* I - group to write(1 = 1st group)           */
            LONGLONG  firstelem, /* I - first vector element to write       */
                                 /* (1 = 1st)                               */
            LONGLONG  nelem,     /* I - number of values to write           */
            unsigned char *array, /* I - array of values that are written   */
            int  *status)    /* IO - error status                           */
/*
  Write an array of values to the primary array. Data conversion
  and scaling will be performed if necessary (e.g, if the datatype of
  the FITS array is not the same as the array being written).
*/
{
    unsigned char nullvalue;

    /*
      the primary array is represented as a binary table:
      each group of the primary array is a row in the table,
      where the first column contains the group parameters
      and the second column contains the image itself.
    */

    fits_write_compressed_pixels(fptr, TBYTE, firstelem, nelem,
            0, array, &nullvalue, status);
    return(*status);
}

/*****************************************************************************/
/*                                                                           */
/* The following code was copied and modified from the FITSIO source code    */
/* file putcold.c.                                                           */
/*                                                                           */
/*****************************************************************************/

/*--------------------------------------------------------------------------*/
static int ffpprd(
            fitsfile *fptr,  /* I - FITS file pointer                       */
            long  group,     /* I - group to write(1 = 1st group)           */
            LONGLONG  firstelem, /* I - first vector element to write       */
                                 /* (1 = 1st)                               */
            LONGLONG  nelem,     /* I - number of values to write           */
            double *array,   /* I - array of values that are written        */
            int  *status)    /* IO - error status                           */
/*
  Write an array of values to the primary array. Data conversion
  and scaling will be performed if necessary (e.g, if the datatype of
  the FITS array is not the same as the array being written).
*/
{
    double nullvalue;

    /*
      the primary array is represented as a binary table:
      each group of the primary array is a row in the table,
      where the first column contains the group parameters
      and the second column contains the image itself.
    */

    fits_write_compressed_pixels(fptr, TDOUBLE, firstelem, nelem,
            0, array, &nullvalue, status);
    return(*status);
}

/*****************************************************************************/
/*                                                                           */
/* The following code was copied and modified from the FITSIO source code    */
/* file putcole.c.                                                           */
/*                                                                           */
/*****************************************************************************/

/*--------------------------------------------------------------------------*/
static int ffppre(
            fitsfile *fptr,  /* I - FITS file pointer                       */
            long  group,     /* I - group to write(1 = 1st group)           */
            LONGLONG firstelem, /* I - first vector element to write        */
                                /* (1 = 1st)                                */
            LONGLONG nelem,     /* I - number of values to write            */
            float *array,    /* I - array of values that are written        */
            int  *status)    /* IO - error status                           */
/*
  Write an array of values to the primary array. Data conversion
  and scaling will be performed if necessary (e.g, if the datatype of
  the FITS array is not the same as the array being written).

  This routine cannot be called directly by users to write to large
  arrays with > 2**31 pixels (although CFITSIO can do so by passing
  the firstelem thru a LONGLONG sized global variable)
*/
{
    float nullvalue;

    /*
      the primary array is represented as a binary table:
      each group of the primary array is a row in the table,
      where the first column contains the group parameters
      and the second column contains the image itself.
    */

    fits_write_compressed_pixels(fptr, TFLOAT, firstelem, nelem,
            0, array, &nullvalue, status);
    return(*status);
}

/*****************************************************************************/
/*                                                                           */
/* The following code was copied and modified from the FITSIO source code    */
/* file putcoli.c.                                                            */
/*                                                                           */
/*****************************************************************************/

/*--------------------------------------------------------------------------*/
static int ffppri(
            fitsfile *fptr,  /* I - FITS file pointer                       */
            long  group,     /* I - group to write (1 = 1st group)          */
            LONGLONG  firstelem, /* I - first vector element to write       */
                                 /* (1 = 1st)                               */
            LONGLONG  nelem,     /* I - number of values to write           */
            short *array,    /* I - array of values that are written        */
            int  *status)    /* IO - error status                           */
/*
  Write an array of values to the primary array. Data conversion
  and scaling will be performed if necessary (e.g, if the datatype of
  the FITS array is not the same as the array being written).
*/
{
    short nullvalue;

    /*
      the primary array is represented as a binary table:
      each group of the primary array is a row in the table,
      where the first column contains the group parameters
      and the second column contains the image itself.
    */

    fits_write_compressed_pixels(fptr, TSHORT, firstelem, nelem,
            0, array, &nullvalue, status);
    return(*status);
}

/*****************************************************************************/
/*                                                                           */
/* The following code was copied and modified from the FITSIO source code    */
/* file putcolk.c.                                                            */
/*                                                                           */
/*****************************************************************************/

/*--------------------------------------------------------------------------*/
static int ffpprk(
            fitsfile *fptr,  /* I - FITS file pointer                       */
            long  group,     /* I - group to write(1 = 1st group)           */
            LONGLONG  firstelem, /* I - first vector element to write       */
                                 /* (1 = 1st)                               */
            LONGLONG  nelem,     /* I - number of values to write           */
            int   *array,    /* I - array of values that are written        */
            int  *status)    /* IO - error status                           */
/*
  Write an array of values to the primary array. Data conversion
  and scaling will be performed if necessary (e.g, if the datatype of
  the FITS array is not the same as the array being written).
*/
{
    int nullvalue;

    /*
      the primary array is represented as a binary table:
      each group of the primary array is a row in the table,
      where the first column contains the group parameters
      and the second column contains the image itself.
    */

    fits_write_compressed_pixels(fptr, TINT, firstelem, nelem,
            0, array, &nullvalue, status);
        return(*status);
    return(*status);
}

/*****************************************************************************/
/*                                                                           */
/* The following code was copied and modified from the FITSIO source code    */
/* file putcolj.c.                                                            */
/*                                                                           */
/*****************************************************************************/

/*--------------------------------------------------------------------------*/
static int ffpprj(
            fitsfile *fptr,  /* I - FITS file pointer                       */
            long  group,     /* I - group to write(1 = 1st group)           */
            LONGLONG  firstelem, /* I - first vector element to write       */
                                 /* (1 = 1st)                               */
            LONGLONG  nelem,     /* I - number of values to write           */
            long  *array,    /* I - array of values that are written        */
            int  *status)    /* IO - error status                           */
/*
  Write an array of values to the primary array. Data conversion
  and scaling will be performed if necessary (e.g, if the datatype of
  the FITS array is not the same as the array being written).
*/
{
    long nullvalue;

    /*
      the primary array is represented as a binary table:
      each group of the primary array is a row in the table,
      where the first column contains the group parameters
      and the second column contains the image itself.
    */

    fits_write_compressed_pixels(fptr, TLONG, firstelem, nelem,
            0, array, &nullvalue, status);
    return(*status);
}
/* ======================================================================== */
/*      the following routines support the 'long long' data type            */
/* ======================================================================== */

/*--------------------------------------------------------------------------*/
static int ffpprjj(
            fitsfile *fptr,  /* I - FITS file pointer                       */
            long  group,     /* I - group to write(1 = 1st group)           */
            LONGLONG  firstelem, /* I - first vector element to write       */
                                 /* (1 = 1st)                               */
            LONGLONG  nelem,     /* I - number of values to write           */
            LONGLONG  *array, /* I - array of values that are written       */
            int  *status)    /* IO - error status                           */
/*
  Write an array of values to the primary array. Data conversion
  and scaling will be performed if necessary (e.g, if the datatype of
  the FITS array is not the same as the array being written).
*/
{
    /*
      the primary array is represented as a binary table:
      each group of the primary array is a row in the table,
      where the first column contains the group parameters
      and the second column contains the image itself.
    */

    _pyfits_ffpmsg("writing to compressed image is not supported");

    return(*status = DATA_COMPRESSION_ERR);
}

/*****************************************************************************/
/*                                                                           */
/* The following code was copied and modified from the FITSIO source code    */
/* file putcol.c.                                                            */
/*                                                                           */
/*****************************************************************************/

/*--------------------------------------------------------------------------*/
int _pyfits_ffppr(
            fitsfile *fptr,  /* I - FITS file pointer                       */
            int  datatype,   /* I - datatype of the value                   */
            LONGLONG  firstelem, /* I - first vector element to write       */
                                 /* (1 = 1st)                               */
            LONGLONG  nelem,     /* I - number of values to write           */
            void  *array,    /* I - array of values that are written        */
            int  *status)    /* IO - error status                           */
/*
  Write an array of values to the primary array.  The datatype of the
  input array is defined by the 2nd argument. Data conversion
  and scaling will be performed if necessary (e.g, if the datatype of
  the FITS array is not the same as the array being written).

*/
{
    long group = 1;

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    if (datatype == TBYTE)
    {
      ffpprb(fptr, group, firstelem, nelem, (unsigned char *) array, status);
    }
    else if (datatype == TSHORT)
    {
      ffppri(fptr, group, firstelem, nelem, (short *) array, status);
    }
    else if (datatype == TINT)
    {
      ffpprk(fptr, group, firstelem, nelem, (int *) array, status);
    }
    else if (datatype == TLONG)
    {
      ffpprj(fptr, group, firstelem, nelem, (long *) array, status);
    }
    else if (datatype == TLONGLONG)
    {
      ffpprjj(fptr, group, firstelem, nelem, (LONGLONG *) array, status);
    }
    else if (datatype == TFLOAT)
    {
      ffppre(fptr, group, firstelem, nelem, (float *) array, status);
    }
    else if (datatype == TDOUBLE)
    {
      ffpprd(fptr, group, firstelem, nelem, (double *) array, status);
    }
    else
    {
      *status = BAD_DATATYPE;
    }

    return(*status);
}

/*****************************************************************************/
/*                                                                           */
/* The following code was copied and modified from the FITSIO source code    */
/* file imcompress.c.                                                        */
/*                                                                           */
/*****************************************************************************/

#define NULL_VALUE -2147483647 /* value used to represent undefined pixels */
#define N_RESERVED_VALUES 1   /* number of reserved values, starting with */
                               /* and including NULL_VALUE.  These values */
                               /* may not be used to represent the quantized */
                               /* and scaled floating point pixel values */

/* nearest integer function */
# define NINT(x)  ((x >= 0.) ? (int) (x + 0.5) : (int) (x - 0.5))

/* ######################################################################## */
/* ###                 Image Compression Routines                       ### */
/* ######################################################################## */

/*--------------------------------------------------------------------------*/
int _pyfits_imcomp_calc_max_elem (int comptype, int nx, int zbitpix,
                                  int blocksize)

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
    else if (comptype == GZIP_1)
    {
        /* gzip usually compressed by at least a factor of 2 for I*4 images */
        /* and somewhat less for I*2 images */
        /* If this size turns out to be too small, then the gzip */
        /* compression routine will allocate more space as required */

        if (zbitpix == 16 || zbitpix == 8)
            return(nx * sizeof(short) / 1.3);
        else
            return(nx * sizeof(int) / 2);
    }
    else if (comptype == HCOMPRESS_1)
    {
        /* Imperical evidence suggests in the worst case,
           the compressed stream could be up to 10% larger than the original
           image.  Add 26 byte overhead, only significant for very small tiles

           Possible improvement: may need to allow a larger size for 32-bit
                                 images */

        if (zbitpix == 16 || zbitpix == 8)

            return( (int) (nx * 2.2 + 26));   /* will be compressing 16-bit */
                                              /* int  array                 */
        else
            return( (int) (nx * 4.4 + 26));   /* will be compressing 32-bit */
                                              /* int array                  */
    }
    else
        return(nx * sizeof(int));
}
/*--------------------------------------------------------------------------*/
static int imcomp_compress_tile (
    fitsfile *outfptr,
    long row,
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
   are written to the output table.

   This input array may be modified by this routine.  If the array is of type TINT
   or TFLOAT, and the compression type is HCOMPRESS, then it must have been
   allocated to be twice as large (8 bytes per pixel) to provide scratch space.

  Note that this routine does not fully support the implicit datatype conversion that
  is supported when writing to normal FITS images.  The datatype of the input array
  must have the same datatype (either signed or unsigned) as the output (compressed)
  FITS image in most cases.
*/
{
    int *idata, *itemp;         /* quantized integer data */
    short *cbuf;        /* compressed data */
    short *sbuff;
    int clen;           /* size of cbuf */
    size_t gzip_clen;
    int flag = 1; /* true by default; only = 0 if float data couldn't be */
                  /* quantized                                           */
    int iminval = 0, imaxval = 0;  /* min and max quantized integers */
    double bscale[1] = {1.}, bzero[1] = {0.};   /* scaling parameters */
    double scale, zero;
    int  nelem = 0;             /* number of bytes */
    size_t gzip_nelem = 0;
    long ii, hcomp_len;
    LONGLONG *lldata;
    unsigned char *usbbuff;
    int ihcompscale, cn_zblank, zbitpix, nullval, flagval, gotnulls;
    int intlength;  /* size of integers to be compressed */
    float floatnull, hcompscale;
    float fminval, fmaxval, delta, zeropt, *fdata, *ftemp;
    double doublenull, noise3;

    if (*status > 0)
        return(*status);

    idata = (int *) tiledata;
    hcompscale = (outfptr->Fptr)->hcomp_scale;
    zbitpix = (outfptr->Fptr)->zbitpix;

    /* if the tile/image has an integer datatype, see if a null value has */
    /* been defined (with the BLANK keyword in a normal FITS image).  */
    /* If so, and if the input tile array also contains null pixels, */
    /* (represented by pixels that have a value = nullflagval) then  */
    /* any pixels whose value = nullflagval, must be set to the      */
    /* value = nullval before the pixel array is compressed.  These  */
    /* null pixel values must not be inverse scaled by the           */
    /* BSCALE/BZERO values, if present.                              */

    cn_zblank = (outfptr->Fptr)->cn_zblank;
    nullval = (outfptr->Fptr)->zblank;

    if (zbitpix > 0 && cn_zblank != -1)  /* If the integer image has no */
                                         /* defined null                */
        nullcheck = 0;    /* value, then don't bother checking input array */
                          /* for nulls.                                    */

    /* if the BSCALE and BZERO keywords exist, then the input values must */
    /* be inverse scaled by this factor, before the values are compressed. */
    /* (The program may have turned off scaling, which over rides the      */
    /* keywords)                                                           */

    scale = (outfptr->Fptr)->cn_bscale;
    zero  = (outfptr->Fptr)->cn_bzero;

    /* =================================================================== */
    /*  Convert input tile array in place to 4 or 8-byte ints for          */
    /*  compression, if needed.  Do null value substitution if needed      */
    /*  Note that the calling routine must have allocated the array big    */
    /*  enough to be able to do this.                                      */

    if (datatype == TSHORT)
    {
       /* datatype of input array is TSHORT.  We only support writing this 
          datatype to a FITS image with BITPIX = 16 and with BZERO = 0 and 
          BSCALE = 1.  */

       if (zbitpix != SHORT_IMG || scale != 1.0 || zero != 0.0) {
           _pyfits_ffpmsg("Datatype conversion/scaling is not supported when writing to compressed images");
           return(*status = DATA_COMPRESSION_ERR);
       }

       sbuff = (short *) tiledata;

       if (((outfptr->Fptr)->compress_type == RICE_1 || 
            (outfptr->Fptr)->compress_type == GZIP_1))
       {
           /* don't have to convert to int if using gzip or Rice compression */
           intlength = 2;

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
       } else {
           /* have to convert to int if using HCOMPRESS or PLIO */
           intlength = 4;

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
    else if (datatype == TINT || (datatype == TLONG && sizeof(long) == 4))
    {
       /* datatype of input array is int.  We only support writing this datatype
          to a FITS image with BITPIX = 32 and with BZERO = 0 and BSCALE = 1.
       */

       if (zbitpix != LONG_IMG || scale != 1.0 || zero != 0.) {
           _pyfits_ffpmsg("Implicit datatype conversion is not supported when writing to compressed images");
           return(*status = DATA_COMPRESSION_ERR);
       }

       intlength = 4;

       if (nullcheck == 1) {
               /* no datatype conversion is required for any of the compression
                  algorithms, except possibly for HCOMPRESS (to I*8), which is
                  handled later.  Just reset pixels equal to flagval to the
                  FITS null value */
               flagval = *(int *) (nullflagval);
               if (flagval != nullval) {
                  for (ii = tilelen - 1; ii >= 0; ii--) {
                    if (idata[ii] == flagval)
                       idata[ii] = nullval;
                  }
               }
       }
    }
    else if (datatype == TBYTE)
    {
       /* datatype of input array is unsigned byte.  We only support writing
          this datatype to a FITS image with BITPIX = 8 and with BZERO = 0 and
          BSCALE = 1.  */

       if (zbitpix != BYTE_IMG || scale != 1.0 || zero != 0.) {
           _pyfits_ffpmsg("Implicit datatype conversion is not supported when writing to compressed images");
           return(*status = DATA_COMPRESSION_ERR);
       }

       usbbuff = (unsigned char *) tiledata;

       if (((outfptr->Fptr)->compress_type == RICE_1 ||
            (outfptr->Fptr)->compress_type == GZIP_1))
       {
           /* don't have to convert to int if using gzip or Rice compression */
           intlength = 1;

           if (nullcheck == 1) {
               /* reset pixels equal to flagval to the FITS null value, prior
                  to compression */
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
           intlength = 4;

           if (nullcheck == 1) {
               /* reset pixels equal to flagval to the FITS null value, prior
                  to compression */
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
    }
    else if (datatype == TLONG && sizeof(long) == 8)
    {
           _pyfits_ffpmsg("Integer*8 Long datatype is not supported when writing to compressed images");
           return(*status = DATA_COMPRESSION_ERR);
    }
    else if (datatype == TFLOAT)
    {
           intlength = 4;

          /* if the tile-compressed table contains zscale and zzero columns */
          /* then scale and quantize the input floating point data.    */
          /* Otherwise, just truncate the floats to (scaled) integers.     */
          if ((outfptr->Fptr)->cn_zscale > 0) {
            if (nullcheck == 1)
              floatnull = *(float *) (nullflagval);
            else
              floatnull = FLOATNULLVALUE;  /* NaNs are represented by this, by
                                              default */

            if ((outfptr->Fptr)->quantize_level < 0)  {

              /* negative value represents the absolute quantization level. */
              /* We don't have to calculate the noise in the image, so do */
              /* this simple linear scaling in line (here) for efficiency, */
              /* instead of calling fits_quantize_float */

              delta = ((outfptr->Fptr)->quantize_level) * -1.;

              fdata = tiledata;
              gotnulls = 0;

              /* set min and max value = first valid pixel value */
              ftemp = fdata;
              fminval = 0;
              fmaxval = 0;
              for (ii = 0; ii < tilelen; ftemp++, ii++) {
                  if (*fdata != floatnull) {
                      fminval = *ftemp;
                      fmaxval = *ftemp;
                      break;
                  }
              }

              /* find min and max values */
              ftemp = fdata;
              for (ii = 0; ii < tilelen; ftemp++, ii++) {
                  if (*ftemp == floatnull) {
                      gotnulls = 1;
                  } else if (*ftemp < fminval) {
                      fminval = *ftemp;
                  } else if (*ftemp > fmaxval) {
                      fmaxval = *ftemp;
                  }
              }

              /* check that the range of quantized levels is not > range of int
              */
              if ((fmaxval - fminval) / delta > 2. * 2147483647. -
                  N_RESERVED_VALUES ) {
                  flag = 0;                     /* don't quantize */
              } else {

                  flag = 1;
                  if (!gotnulls) {   /* don't have to check for nulls */
                  /* return all positive values, if possible since some */
                  /* compression algorithms either only work for positive */
                  /* integers, or are more efficient.  */
                      if ((fmaxval - fminval) / delta < 2147483647. - 
                          N_RESERVED_VALUES ) {
                          zeropt = fminval;
                          ftemp = fdata;
                          itemp = idata;
                          for (ii = 0;  ii < tilelen;  ftemp++, itemp++, ii++) {
                              *itemp = (int) (((*ftemp - zeropt) / delta) + 
                                              0.5f);
                          }
                       } else {
                          /* center the quantized levels around zero */
                          zeropt = (fminval + fmaxval) / 2.;
                          for (ii = 0;  ii < tilelen;  ii++) {
                              idata[ii] = NINT((fdata[ii] - zeropt) / delta);
                          }
                      }
                  } else {
                      /* data contains null values; shift the range to be */
                      /* close to the value used to represent null values */
                     zeropt = fminval - delta*(NULL_VALUE + N_RESERVED_VALUES);

                      for (ii = 0;  ii < tilelen;  ii++) {
                          if (fdata[ii] != floatnull) {
                              idata[ii] = NINT ((fdata[ii] - zeropt) / delta);
                          } else  {
                              idata[ii] = NULL_VALUE;
                          }
                      }
                 }

                 /* calc min and max values of the integer array */

                 bscale[0] = delta;
                 bzero[0]  = zeropt;
              }
            } else {
                /* quantize level is positive, so we have to calculate the */
                /* noise quantize the float values into integers */
                flag = _pyfits_fits_quantize_float ((float *) tiledata, tilenx,
                   tileny, nullcheck, floatnull,
                   (outfptr->Fptr)->quantize_level, idata,
                   bscale, bzero, &iminval, &imaxval);
            }
          }
          else  /* input float data is implicitly converted (truncated) to
                   integers */
          {
            if ((scale != 1. || zero != 0.))  /* must scale the values */
               imcomp_nullscalefloats((float *) tiledata, tilelen, idata,
                   scale, zero, nullcheck, *(float *) (nullflagval), nullval,
                   status);
             else
               imcomp_nullfloats((float *) tiledata, tilelen, idata,
                   nullcheck, *(float *) (nullflagval), nullval,  status);
          }
    }
    else if (datatype == TDOUBLE)
    {
           intlength = 4;

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
            flag = _pyfits_fits_quantize_double ((double *) tiledata, tilenx,
               tileny, nullcheck, doublenull, (outfptr->Fptr)->quantize_level,
               idata, bscale, bzero, &iminval, &imaxval);
          }
          else  /* input double data is implicitly converted (truncated) to
                   integers */
          {
             if ((scale != 1. || zero != 0.))  /* must scale the values */
               imcomp_nullscaledoubles((double *) tiledata, tilelen, idata,
                   scale, zero, nullcheck, *(double *) (nullflagval), nullval,
                   status);
             else
               imcomp_nulldoubles((double *) tiledata, tilelen, idata,
                   nullcheck, *(double *) (nullflagval), nullval,  status);
          }
    }
    else
    {
          _pyfits_ffpmsg("unsupported datatype (imcomp_compress_tile)");
          return(*status = BAD_DATATYPE);
    }

    /* ==================================================================== */

    if (flag)   /* we can now compress the int array */
    {
        /* allocate buffer for the compressed tile bytes */
        clen = (outfptr->Fptr)->maxelem;
        cbuf = (short *) calloc (clen, sizeof (unsigned char));

        if (cbuf == NULL)
        {
            _pyfits_ffpmsg("Out of memory. (imcomp_compress_tile)");
            return (*status = MEMORY_ALLOCATION);
        }

    /* =================================================================== */

        /* Compress the integer data, then write the compressed bytes */
        if ( (outfptr->Fptr)->compress_type == RICE_1)
        {
            if (intlength == 2) {
                nelem = _pyfits_fits_rcomp_short ((short *)idata, tilelen,
                       (unsigned char *) cbuf,
                       clen, (outfptr->Fptr)->rice_blocksize);
            } else if (intlength == 1) {
                nelem = _pyfits_fits_rcomp_byte ((signed char *)idata, tilelen,
                       (unsigned char *) cbuf,
                       clen, (outfptr->Fptr)->rice_blocksize);
            } else {
                nelem = _pyfits_fits_rcomp (idata, tilelen,
                       (unsigned char *) cbuf,
                       clen, (outfptr->Fptr)->rice_blocksize);
            }

                /* Write the compressed byte stream. */

                (outfptr->Fptr)->dataLen[row-1] = nelem;
                (outfptr->Fptr)->data[row-1] = (unsigned char *)cbuf;

        }
        else if ( (outfptr->Fptr)->compress_type == PLIO_1)
        {
              for (ii = 0; ii < tilelen; ii++)  {
                if (idata[ii] < 0 || idata[ii] > 16777215)
                {
                   /* plio algorithn only supports positive 24 bit ints */
                   _pyfits_ffpmsg("data out of range for PLIO compression (0 - 2**24)");
                   return(*status = DATA_COMPRESSION_ERR);
                }
              }

                nelem = _pyfits_pl_p2li (idata, 1, cbuf, tilelen);

                /* Write the compressed byte stream. */

                (outfptr->Fptr)->dataLen[row-1] = nelem*2;
                (outfptr->Fptr)->data[row-1] = (unsigned char *)cbuf;
        }
        else if ( (outfptr->Fptr)->compress_type == GZIP_1)
        {

#if BYTESWAPPED
           if (intlength == 2)
               ffswap2((short *) idata, tilelen);
           else if (intlength == 4)
               ffswap4(idata, tilelen);
#endif

           gzip_clen = clen;
           if (intlength == 2) {
                 _pyfits_compress2mem_from_mem(
                 (char *) idata, tilelen * sizeof(short),
                 (char **) &cbuf, (size_t *) &gzip_clen, realloc,
                 &gzip_nelem, status);
           } else if (intlength == 1) {
                _pyfits_compress2mem_from_mem((char *) idata, 
                 tilelen * sizeof(unsigned char),
                 (char **) &cbuf, (size_t *) &gzip_clen, realloc,
                 &gzip_nelem, status);
           } else {
                _pyfits_compress2mem_from_mem(
                 (char *) idata, tilelen * sizeof(int),
                 (char **) &cbuf, (size_t *) &gzip_clen, realloc,
                 &gzip_nelem, status);
           }

                /* Write the compressed byte stream. */

                (outfptr->Fptr)->dataLen[row-1] = gzip_nelem;
                (outfptr->Fptr)->data[row-1] = (unsigned char *)cbuf;
        }
        else if ( (outfptr->Fptr)->compress_type == HCOMPRESS_1)
        {
            /*
              if hcompscale is positive, then we have to multiply
              the value by the RMS background noise to get the
              absolute scale value.  If negative, then it gives the
              absolute scale value directly.
            */
            hcompscale = (outfptr->Fptr)->hcomp_scale;

            if (hcompscale > 0.) {
               _pyfits_fits_img_stats_int(idata, tilenx, tileny, nullcheck,
                        nullval, 0,0,0,0,0,0,&noise3,status);

                hcompscale = hcompscale * noise3;

            } else if (hcompscale < 0.) {

                hcompscale = hcompscale * -1.0;
            }

            ihcompscale = (int) (hcompscale + 0.5);

            hcomp_len = clen;  /* allocated size of the buffer */

            if (zbitpix == BYTE_IMG || zbitpix == SHORT_IMG) {
                _pyfits_fits_hcompress(idata, tilenx, tileny,
                  ihcompscale, (char *) cbuf, &hcomp_len, status);

            } else {
                 /* have to convert idata to an I*8 array, in place */
                 /* idata must have been allocated large enough to do this */
                lldata = (LONGLONG *) idata;

                for (ii = tilelen - 1; ii >= 0; ii--) {
                    lldata[ii] = idata[ii];;
                }

                _pyfits_fits_hcompress64(lldata, tilenx, tileny,
                  ihcompscale, (char *) cbuf, &hcomp_len, status);
            }

            /* Write the compressed byte stream. */

            (outfptr->Fptr)->dataLen[row-1] = hcomp_len;
            (outfptr->Fptr)->data[row-1] = (unsigned char *)cbuf;
        }

        if (nelem < 0)  /* error condition */
        {
            free (cbuf); cbuf = 0;
            _pyfits_ffpmsg
                ("error compressing row of the image (imcomp_compress_tile)");
            return (*status = DATA_COMPRESSION_ERR);
        }

        if ((outfptr->Fptr)->cn_zscale > 0)
        {
              /* write the linear scaling parameters */
              (outfptr->Fptr)->c_zscale[row-1] = bscale[0];
              (outfptr->Fptr)->c_zzero[row-1] = bzero[0];
        }
    }
    else     /* floating point data couldn't be quantized */
    {
         /* Write the original floating point data. */

         if ((outfptr->Fptr)->cn_uncompressed > 0)
         {
             (outfptr->Fptr)->ucDataLen[row-1] = tilelen;
             (outfptr->Fptr)->ucData[row-1] = (void*)malloc(
                                                       tilelen*sizeof(double));

             if (datatype == TFLOAT)
             {
                 for(ii=0; ii < tilelen; ii++)
                 {
                    ((double**)((outfptr->Fptr)->ucData))[row-1][ii] =
                                                        ((float*)tiledata)[ii];
                 }
             }
             else if (datatype == TDOUBLE)
             {
                 for(ii=0; ii < tilelen; ii++)
                 {
                    ((double**)((outfptr->Fptr)->ucData))[row-1][ii] =
                                                       ((double*)tiledata)[ii];
                 }
             }
         }
         else
         {
             _pyfits_ffpmsg("There is no UNCOMPRESSED_DATA column in the table.");
             return(*status = BAD_COL_NUM);
         }

    }

    return (*status);
}
/*---------------------------------------------------------------------------*/
static int imcomp_nullscale(
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
static int imcomp_nullvalues(
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
static int imcomp_scalevalues(
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
static int imcomp_nullfloats(
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
static int imcomp_nullscalefloats(
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
static int imcomp_nulldoubles(
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
static int imcomp_nullscaledoubles(
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
static int fits_write_compressed_img(
            fitsfile *fptr,   /* I - FITS file pointer     */
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
    int ii, i5, i4, i3, i2, i1, i0, ndim, irow, pixlen;
    int buffpixsiz;
    void *buffer;
    char *bnullarray = 0; //, card[FLEN_CARD];

    if (*status > 0)
        return(*status);

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
        _pyfits_ffpmsg("unsupported datatype for compressing image");
        return(*status = BAD_DATATYPE);
    }

    /* allocate scratch space for processing one tile of the image */
    buffpixsiz = pixlen;  /* this is the minimum pixel size */

    if ( (fptr->Fptr)->compress_type == HCOMPRESS_1) { 
        /* need 4 or 8 bytes per pixel */
        if ((fptr->Fptr)->zbitpix == BYTE_IMG ||
            (fptr->Fptr)->zbitpix == SHORT_IMG )
                buffpixsiz = maxvalue(buffpixsiz, 4);
        else
                buffpixsiz = 8;
    }
    else if ( (fptr->Fptr)->compress_type == PLIO_1) {
                /* need 4 bytes per pixel */
                buffpixsiz = maxvalue(buffpixsiz, 4);
    }
    else if ( (fptr->Fptr)->compress_type == RICE_1  ||
              (fptr->Fptr)->compress_type == GZIP_1) { 
         /* need 1, 2, or 4 bytes per pixel */
        if ((fptr->Fptr)->zbitpix == BYTE_IMG)
            buffpixsiz = maxvalue(buffpixsiz, 1);
        else if ((fptr->Fptr)->zbitpix == SHORT_IMG)
            buffpixsiz = maxvalue(buffpixsiz, 2);
        else
            buffpixsiz = maxvalue(buffpixsiz, 4);
    }
    else
    {
        _pyfits_ffpmsg("unsupported image compression algorithm");
        return(*status = BAD_DATATYPE);
    }

    /* cast to double to force alignment on 8-byte addresses */
    buffer = (double *) calloc ((fptr->Fptr)->maxtilelen, buffpixsiz);

    if (buffer == NULL)
    {
            _pyfits_ffpmsg("Out of memory (fits_write_compress_img)");
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
        fpixel[ii] = infpixel[ii];
        lpixel[ii] = inlpixel[ii];

        /* calc number of tiles in each dimension, and tile containing */
        /* the first and last pixel we want to read in each dimension  */
        naxis[ii] = (fptr->Fptr)->znaxis[ii];
        if (fpixel[ii] < 1)
        {
            free(buffer); buffer = 0;
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

              memset(buffer, 0, pixlen * thistilesize[0]);

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

    free(buffer); buffer = 0;

    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fits_write_compressed_pixels(
            fitsfile *fptr, /* I - FITS file pointer   */
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
    naxis = (fptr->Fptr)->zndim;

    for ( ii = 0; ii < (fptr->Fptr)->zndim; ii++)
    {
        naxes[ii] = ((fptr->Fptr)->znaxis)[ii];
    }
    /* calc the cumulative number of pixels in each successive dimension */

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
        _pyfits_ffpmsg("only 1D, 2D, or 3D images are currently supported");
        return(*status = DATA_COMPRESSION_ERR);
    }

    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fits_write_compressed_img_plane(
            fitsfile *fptr, /* I - FITS file    */
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

/*---------------------------------------------------------------------------*/
static int fits_read_compressed_img(
            fitsfile *fptr,   /* I - FITS file pointer      */
            int  datatype,  /* I - datatype of the array to be returned      */
            LONGLONG  *infpixel, /* I - 'bottom left corner' of the          */
                                 /* subsection                               */
            LONGLONG  *inlpixel, /* I - 'top right corner' of the subsection */
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
        _pyfits_ffpmsg("unsupported datatype for uncompressing image");
        return(*status = BAD_DATATYPE);
    }

    /* If nullcheck ==1 and nullval == 0, then this means that the */
    /* calling routine does not want to check for null pixels in the array */
    if (nullcheck == 1 && testnullval == 0.)
        nullcheck = 0;

    if (buffer == NULL)
    {
        _pyfits_ffpmsg("Out of memory (fits_read_compress_img)");
        return (*status = MEMORY_ALLOCATION);
    }

    /* allocate memory for a null flag array, if needed */
    if (nullcheck == 2)
    {
        bnullarray = calloc ((fptr->Fptr)->maxtilelen, sizeof (char));

        if (bnullarray == NULL)
        {
            _pyfits_ffpmsg("Out of memory (fits_read_compress_img)");
            free(buffer); buffer = 0;
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
                free(bnullarray); bnullarray = 0;
            }
            free(buffer); buffer = 0;
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

              /* read and uncompress this row (tile) of the table */
              /* also do type conversion and undefined pixel substitution */
              /* at this point */

              imcomp_decompress_tile(fptr, irow, thistilesize[0],
                    datatype, nullcheck, nullval, buffer, bnullarray, &tilenul,
                    status);

              if (tilenul && anynul)
                  *anynul = 1;  /* there are null pixels */
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
       free(bnullarray); bnullarray = 0;
    }
    free(buffer); buffer = 0;

    return(*status);
}
/*--------------------------------------------------------------------------*/
static int fits_read_compressed_pixels(
            fitsfile *fptr, /* I - FITS file pointer    */
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
    naxis = (fptr->Fptr)->zndim;

    for ( ii = 0; ii < (fptr->Fptr)->zndim; ii++)
    {
        naxes[ii] = ((fptr->Fptr)->znaxis)[ii];
    }
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
        _pyfits_ffpmsg("only 1D, 2D, or 3D images are currently supported");
        return(*status = DATA_DECOMPRESSION_ERR);
    }

    return(*status);
}
///*--------------------------------------------------------------------------*/
static int fits_read_compressed_img_plane(
            fitsfile *fptr, /* I - FITS file   */
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
static int imcomp_decompress_tile (
          fitsfile *infptr,
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
    static int *idata = 0;     /* this variable must persist */
    int tiledatatype, pixlen;          /* uncompressed integer data */
    LONGLONG *lldata = 0;
    size_t idatalen, tilebytesize;
    int ii, tnull = 0;        /* value in the data which represents nulls */
    short *sbuf;
    int blocksize;
    double bscale, bzero, dummy = 0;    /* scaling parameters */
    long nelem = 0;      /* number of bytes */
    int smooth, nx, ny, scale;  /* hcompress parameters */

    if (*status > 0)
       return(*status);

    nelem = (infptr->Fptr)->dataLen[nrow-1];

    if (nelem == 0)  /* tile was not compressed; read uncompressed data */
    {
        if ((infptr->Fptr)->cn_uncompressed < 1 )
        {
            return (*status = NO_COMPRESSED_TILE);
        }

        /* no compressed data, so simply read the uncompressed data */
        /* directly from the UNCOMPRESSED_DATA column, then return */

        for (ii = 0; ii < ((infptr->Fptr)->ucDataLen)[nrow-1]; ii++)
        {
            if (datatype == TFLOAT)
            {
                ((float*)buffer)[ii] =
                              ((float**)((infptr->Fptr)->ucData))[nrow-1][ii];
            }
            else if (datatype == TDOUBLE)
            {
                ((double*)buffer)[ii] =
                             ((double**)((infptr->Fptr)->ucData))[nrow-1][ii];
            }
        }

        ((infptr->Fptr)->dataLen)[nrow-1] = ((infptr->Fptr)->ucDataLen)[nrow-1];

        return(*status);
    }

    /* **************************************************************** */

    if (nullcheck == 2)
    {
        for (ii = 0; ii < tilelen; ii++)  /* initialize the null array */
            bnullarray[ii] = 0;
    }

    if (anynul)
       *anynul = 0;

    /* get linear scaling and offset values, if they exist */

    if ((infptr->Fptr)->cn_zscale == 0)
    {
         /* set default scaling, if scaling is not defined */
         bscale = 1.;
         bzero = 0.;
    }
    else if ((infptr->Fptr)->cn_zscale == -1)
    {
        bscale = (infptr->Fptr)->zscale;
        bzero  = (infptr->Fptr)->zzero;
    }
    else
    {
        /* read the linear scale and offset values for this row */
        bscale = ((infptr->Fptr)->bscale)[nrow-1];
        bzero = ((infptr->Fptr)->bzero)[nrow-1];

        /* test if floating-point FITS image also has non-default BSCALE and */
        /* BZERO keywords.  If so, we have to combine the 2 linear scaling   */
        /* factors.                                                          */

        if ( ((infptr->Fptr)->zbitpix == FLOAT_IMG ||
              (infptr->Fptr)->zbitpix == DOUBLE_IMG )
            &&
              ((infptr->Fptr)->cn_bscale != 1.0 ||
               (infptr->Fptr)->cn_bzero  != 0.0 )    )
            {
               bscale = bscale * (infptr->Fptr)->cn_bscale;
               bzero  = bzero  * (infptr->Fptr)->cn_bscale + 
                        (infptr->Fptr)->cn_bzero;
            }
    }

    if (bscale == 1.0 && bzero == 0.0 )
    {
      /* if no other scaling has been specified, try using the values
         given by the BSCALE and BZERO keywords, if any */

        bscale = (infptr->Fptr)->cn_bscale;
        bzero  = (infptr->Fptr)->cn_bzero;
    }

    /* ************************************************************* */
    /* get the value used to represent nulls in the int array */
    if ((infptr->Fptr)->cn_zblank == 0)
    {
        nullcheck = 0;  /* no null value; don't check for nulls */
    }
    else if ((infptr->Fptr)->cn_zblank == -1)
    {
        tnull = (infptr->Fptr)->zblank;  /* use the the ZBLANK keyword */
    }
    else
    {
        /* read the null value for this row */
        tnull = ((infptr->Fptr)->blank)[nrow-1];
    }

    /* ************************************************************* */

    /* allocate memory for the uncompressed array of tile integers */

    if ((infptr->Fptr)->compress_type == HCOMPRESS_1 &&
          ((infptr->Fptr)->zbitpix != BYTE_IMG &&
           (infptr->Fptr)->zbitpix != SHORT_IMG) ) {

           /*  must allocate 8 bytes per pixel of scratch space */
           lldata = (LONGLONG*) malloc (tilelen * sizeof (LONGLONG));
           idata = (int *) lldata;
    } else if ( (infptr->Fptr)->compress_type == RICE_1 &&
               (infptr->Fptr)->zbitpix == BYTE_IMG &&
               (infptr->Fptr)->rice_bytepix == 1) {

           /*  must allocate 1 byte per pixel of scratch space */
           idatalen = tilelen;
           idata = (int *) malloc (idatalen);
    } else if ( (infptr->Fptr)->compress_type == RICE_1 &&
               (infptr->Fptr)->zbitpix == SHORT_IMG &&
               (infptr->Fptr)->rice_bytepix == 2) {

           /*  must allocate 2 bytes per pixel of scratch space */
           idatalen = tilelen * sizeof(short);
           idata = (int *) malloc (idatalen);
     } else if ( (infptr->Fptr)->compress_type == GZIP_1 &&
               (infptr->Fptr)->zbitpix == SHORT_IMG ) {

           /*  must allocate 2 bytes per pixel of scratch space */
           idatalen = tilelen * sizeof(short);
           idata = (int *) malloc (idatalen);
    } else if ((infptr->Fptr)->compress_type == GZIP_1 &&
               (infptr->Fptr)->zbitpix == BYTE_IMG ) {

           /*  must allocate 1 byte per pixel of scratch space */
           idatalen = tilelen * sizeof(char);
           idata = (int *) malloc (idatalen);
    } else {
           /* all other cases have int pixels */
           idatalen = tilelen * sizeof(int);
           idata = (int*) malloc (idatalen);
    }

    if (idata == NULL)
    {
            _pyfits_ffpmsg("Out of memory for idata. (imcomp_decompress_tile)");
            return (*status = MEMORY_ALLOCATION);
    }

    /* ************************************************************* */
    /*    call the algorithm-specific code to uncompress the tile */

    /* default uncomopressed pixels have int data type */
    tiledatatype = TINT;

    if ((infptr->Fptr)->compress_type == RICE_1)
    {

        /* uncompress the data */
        blocksize = (infptr->Fptr)->rice_blocksize;

         if ((infptr->Fptr)->rice_bytepix == 1 ) {
            if ((*status = _pyfits_fits_rdecomp_byte (
                ((infptr->Fptr)->data)[nrow-1],
                ((infptr->Fptr)->dataLen)[nrow-1],
                (unsigned char *)idata,
                tilelen, blocksize)))
            {
                free(idata);
                return (*status);
            }
            tiledatatype = TBYTE;
        } else if ((infptr->Fptr)->rice_bytepix == 2 ) {
            if ((*status = _pyfits_fits_rdecomp_short (
                ((infptr->Fptr)->data)[nrow-1],
                ((infptr->Fptr)->dataLen)[nrow-1],
                (unsigned short *)idata,
                tilelen, blocksize)))
            {
                free(idata);
                return (*status);
            }
            tiledatatype = TSHORT;
        } else {
            if ((*status = _pyfits_fits_rdecomp (
                ((infptr->Fptr)->data)[nrow-1],
                ((infptr->Fptr)->dataLen)[nrow-1],
                (unsigned int *)idata,
                tilelen, blocksize)))
            {
                free(idata);
                return (*status);
            }
            tiledatatype = TINT;
        }

    }

    /* ************************************************************* */

    else if ((infptr->Fptr)->compress_type == HCOMPRESS_1)
    {
        /* uncompress the data */
        smooth = (infptr->Fptr)->hcomp_smooth;

        if ( ((infptr->Fptr)->zbitpix == BYTE_IMG ||
               (infptr->Fptr)->zbitpix == SHORT_IMG) )  {

            if ((*status = _pyfits_fits_hdecompress(
                                          ((infptr->Fptr)->data)[nrow-1],
                                            smooth, idata, &nx, &ny,
                                            &scale, status)))
            {
                free(idata); idata = 0;
                return (*status);
            }

        } else {

            /* idata must have been allocated twice as large for this to work */
            if ((*status = _pyfits_fits_hdecompress64(
                                            ((infptr->Fptr)->data)[nrow-1],
                                              smooth, lldata, &nx, &ny,
                &scale, status)))
            {
                free(idata); idata = 0;
                return (*status);
            }
        }

    }

    /* ************************************************************* */

    else if ((infptr->Fptr)->compress_type == PLIO_1)
    {
        nelem = ((infptr->Fptr)->dataLen)[nrow-1];
        nelem = nelem/2;

        sbuf = ((short *)(((infptr->Fptr)->data)[nrow-1]));

#if BYTESWAPPED
        ffswap2(sbuf, nelem); /* reverse order of bytes */
#endif

        _pyfits_pl_l2pi (sbuf, 1, idata, tilelen);  /* uncompress the data */

    }

    /* ************************************************************* */

    else if ((infptr->Fptr)->compress_type == GZIP_1)
    {
        /* uncompress the data */

        if (_pyfits_uncompress2mem_from_mem
                                    ((char *)((infptr->Fptr)->data)[nrow-1],
                                     ((infptr->Fptr)->dataLen)[nrow-1],
             (char **) &idata, &idatalen, realloc, &tilebytesize, status))
        {
            _pyfits_ffpmsg("_pyfits_uncompress2mem_from_mem returned with an error");
            free(idata); idata = 0;
            return (*status);
        }

        if (tilebytesize == tilelen * 2) {
            /* this is a short I*2 array */
            tiledatatype = TSHORT;

#if BYTESWAPPED
            ffswap2((short *) idata, tilelen);
#endif

        } else if (tilebytesize == tilelen * 4) {
            /* this is a int I*4 array */
            tiledatatype = TINT;

#if BYTESWAPPED
         ffswap4(idata, tilelen); /* reverse order of bytes */
#endif

        } else if (tilebytesize == tilelen) {

            /* this is an unsigned char I*1 array */
            tiledatatype = TBYTE;

        } else {
            _pyfits_ffpmsg("error: uncompressed tile has wrong size");
            free(idata);
            return (*status = DATA_DECOMPRESSION_ERR);
        }
    }

    /* ************************************************************* */
    else
    {
        _pyfits_ffpmsg("unknown compression algorithm");
        free(idata); idata = 0;
        return (*status = DATA_DECOMPRESSION_ERR);
    }

    /* ************************************************************* */
    /* copy the uncompressed tile data to the output buffer, doing */
    /* null checking, datatype conversion and linear scaling, if necessary */

    if (nulval == 0)
         nulval = &dummy;  /* set address to dummy value */

    if (datatype == TSHORT)
    {
        pixlen = sizeof(short);

        if (tiledatatype == TINT)
          fffi4i2(idata, tilelen, bscale, bzero, nullcheck, tnull,
           *(short *) nulval, bnullarray, anynul,
          (short *) buffer, status);
        else if (tiledatatype == TSHORT)
        {
          fffi2i2((short *)idata, tilelen, bscale, bzero, nullcheck, tnull,
           *(short *) nulval, bnullarray, anynul,
          (short *) buffer, status);
        }
        else if (tiledatatype == TBYTE)
          fffi1i2((unsigned char *)idata, tilelen, bscale, bzero, nullcheck,
           tnull, *(short *) nulval, bnullarray, anynul,
           (short *) buffer, status);
    }
    else if (datatype == TINT)
    {
        pixlen = sizeof(int);
        if (tiledatatype == TINT)
          fffi4int(idata, (long) tilelen, bscale, bzero, nullcheck, tnull,
           *(int *) nulval, bnullarray, anynul,
           (int *) buffer, status);
        else if (tiledatatype == TSHORT)
          fffi2int((short *)idata, tilelen, bscale, bzero, nullcheck, tnull,
           *(int *) nulval, bnullarray, anynul,
           (int *) buffer, status);
        else if (tiledatatype == TBYTE)
          fffi1int((unsigned char *)idata, tilelen, bscale, bzero, nullcheck,
           tnull, *(int *) nulval, bnullarray, anynul,
           (int *) buffer, status);
    }
    else if (datatype == TLONG)
    {
        pixlen = sizeof(long);
        if (tiledatatype == TINT)
          fffi4i4(idata, tilelen, bscale, bzero, nullcheck, tnull,
           *(long *) nulval, bnullarray, anynul,
            (long *) buffer, status);
        else if (tiledatatype == TSHORT)
          fffi2i4((short *)idata, tilelen, bscale, bzero, nullcheck, tnull,
           *(long *) nulval, bnullarray, anynul,
            (long *) buffer, status);
        else if (tiledatatype == TBYTE)
          fffi1i4((unsigned char *)idata, tilelen, bscale, bzero, nullcheck,
            tnull, *(long *) nulval, bnullarray, anynul,
            (long *) buffer, status);
    }
    else if (datatype == TFLOAT)
    {
        pixlen = sizeof(float);
        if (tiledatatype == TINT)
          fffi4r4(idata, tilelen, bscale, bzero, nullcheck, tnull,
           *(float *) nulval, bnullarray, anynul,
            (float *) buffer, status);
        else if (tiledatatype == TSHORT)
          fffi2r4((short *)idata, tilelen, bscale, bzero, nullcheck, tnull,
           *(float *) nulval, bnullarray, anynul,
            (float *) buffer, status);
        else if (tiledatatype == TBYTE)
          fffi1r4((unsigned char *)idata, tilelen, bscale, bzero, nullcheck,
            tnull, *(float *) nulval, bnullarray, anynul,
            (float *) buffer, status);
    }
    else if (datatype == TDOUBLE)
    {
        pixlen = sizeof(double);
        if (tiledatatype == TINT)
          fffi4r8(idata, tilelen, bscale, bzero, nullcheck, tnull,
           *(double *) nulval, bnullarray, anynul,
            (double *) buffer, status);
        else if (tiledatatype == TSHORT)
          fffi2r8((short *)idata, tilelen, bscale, bzero, nullcheck, tnull,
           *(double *) nulval, bnullarray, anynul,
            (double *) buffer, status);
        else if (tiledatatype == TBYTE)
          fffi1r8((unsigned char *)idata, tilelen, bscale, bzero, nullcheck,
            tnull, *(double *) nulval, bnullarray, anynul,
            (double *) buffer, status);
    }
    else if (datatype == TBYTE)
    {
        pixlen = sizeof(char);
        if (tiledatatype == TINT)
          fffi4i1(idata, tilelen, bscale, bzero, nullcheck, tnull,
           *(unsigned char *) nulval, bnullarray, anynul,
            (unsigned char *) buffer, status);
        else if (tiledatatype == TSHORT)
          fffi2i1((short *)idata, tilelen, bscale, bzero, nullcheck, tnull,
           *(unsigned char *) nulval, bnullarray, anynul,
            (unsigned char *) buffer, status);
        else if (tiledatatype == TBYTE)
          fffi1i1((unsigned char *)idata, tilelen, bscale, bzero, nullcheck,
            tnull, *(unsigned char *) nulval, bnullarray, anynul,
            (unsigned char *) buffer, status);
    }
    else
    {
         *status = BAD_DATATYPE;
    }

    free(idata); idata = 0;

    return (*status);
}
/*--------------------------------------------------------------------------*/
static int imcomp_copy_overlap (
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
             /* offset to byte within the row */
             if (inc[0] > 0)
                imgpix = imgfpix[0] + im1;
             else
                imgpix = imgdim[0] - 1 - imgfpix[0] + im1;
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
static int imcomp_merge_overlap (
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
             /* offset to byte within the row */
             if (inc[0] > 0)
                imgpix = imgfpix[0] + im1;
             else
                imgpix = imgdim[0] - 1 - imgfpix[0] + im1;
             /* loop over pixels along one row of the image */
             for (ipos = imgfpix[0]; ipos <= imglpix[0]; ipos += overlap_flags)
             {
               /* convert from image pixel to byte offset */
               tilepixbyte = tilepix * pixlen;
               imgpixbyte  = imgpix  * pixlen;
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





