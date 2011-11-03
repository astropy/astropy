/* $Id$
*/

/*****************************************************************************/
/*                                                                           */
/* This file, fitsio.h, contains the header information for the core of a    */
/* set of FITSIO routines that are used to compress and uncompress image     */
/* data in FITS binary tables.  The code was copied and modified from the    */
/* FITSIO software written at HEASRC.  The goal for the pyfitsComp module    */
/* was to take this code nearly intact.  In FITSIO, interaction with the     */
/* FITS file is accomplished directly within the FITSIO code.  With          */
/* pyfitsComp, interaction with the FITS file is accomplished from within    */
/* pyfits.  This may make some of the constructs within the FISTIO code seem */
/* confusing when viewed from the perspective of pyfitsComp.  It should be   */
/* noted that the FITSfile structure acts as the interface to the file in    */
/* both cases.  In FITSIO it contains the file handle in order to access the */
/* file, and in pyfitsComp it holds the file data, either compressed or      */
/* uncompressed.                                                             */
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

#ifndef _FITSIO_H
#define _FITSIO_H

#include <stdio.h>

/* the following was provided by Michael Greason (GSFC) to fix a */
/*  C/Fortran compatibility problem on an SGI Altix system running */
/*  SGI ProPack 4 [this is a Novell SuSE Enterprise 9 derivative]  */
/*  and using the Intel C++ and Fortran compilers (version 9.1)  */
#if defined(__INTEL_COMPILER) && defined(__itanium__)
#  define mipsFortran 1
#  define _MIPS_SZLONG 64
#endif

#if defined(linux) || defined(__APPLE__) || defined(__sgi)
#  include <sys/types.h>  /* apparently needed on debian linux systems */
#endif                    /* to define off_t                           */

#include <stdlib.h>  /* apparently needed to define size_t with gcc 2.8.1 */
#include <limits.h>  /* needed for LLONG_MAX and INT64_MAX definitions */

/* Define the datatype for variables which store file offset values. */
/* The newer 'off_t' datatype should be used for this purpose, but some */
/* older compilers do not recognize this type, in which case we use 'long' */
/* instead.  Note that _OFF_T is defined (or not) in stdio.h depending */
/* on whether _LARGEFILE_SOURCE is defined in sys/feature_tests.h  */
/* (at least on Solaris platforms using cc)  */

/*  Debian systems require the 2nd test, below,         */
/*  i.e, "(defined(linux) && defined(__off_t_defined))" */
#if defined(_OFF_T) || (defined(linux) && defined(__off_t_defined)) || defined(_MIPS_SZLONG) || defined(__APPLE__) || defined(_AIX)
#    define OFF_T off_t
#else
#    define OFF_T long
#endif

/* this block determines if the the string function name is 
    strtol or strtoll, and whether to use %ld or %lld in printf statements */

/* 
   The following 2 cases for that Athon64 were removed on 4 Jan 2006;  
   they appear to be incorrect now that LONGLONG is always typedef'ed 
   to 'long long'
    ||  defined(__ia64__)   \
    ||  defined(__x86_64__) \
*/
#if (defined(__alpha) && ( defined(__unix__) || defined(__NetBSD__) )) \
    ||  defined(__sparcv9)  \
    ||  defined(__powerpc64__) || defined(__64BIT__) \
    ||  (defined(_MIPS_SZLONG) &&  _MIPS_SZLONG == 64) \
    ||  defined( _MSC_VER)|| defined(__BORLANDC__)
    
#   define USE_LL_SUFFIX 0
#else
#   define USE_LL_SUFFIX 1
#endif

/* 
   Determine what 8-byte integer data type is available.
  'long long' is now supported by most compilers, but
  older MS Visual C++ compilers before V7.0 use '__int64' instead.
*/

#ifndef LONGLONG_TYPE   /* this may have been previously defined */
#if defined(_MSC_VER)   /* Microsoft Visual C++ */

#if (_MSC_VER < 1300)   /* versions earlier than V7.0 do not have 'long long' */
    typedef __int64 LONGLONG;
#else                   /* newer versions do support 'long long' */
    typedef long long LONGLONG; 
#endif

#elif defined( __BORLANDC__)  /* for the Borland 5.5 compiler, in particular */
    typedef __int64 LONGLONG;
#else
    typedef long long LONGLONG; 
#endif

#define LONGLONG_TYPE
#endif  

#define TBIT          1  /* codes for FITS table data types */
#define TBYTE        11
#define TSBYTE       12
#define TLOGICAL     14
#define TUSHORT      20
#define TSHORT       21
#define TUINT        30
#define TINT         31
#define TULONG       40
#define TLONG        41
#define TFLOAT       42
#define TLONGLONG    81
#define TDOUBLE      82

#define INT32BIT int  /* 32-bit integer datatype.  Currently this       */
                      /* datatype is an 'int' on all useful platforms   */
                      /* however, it is possible that that are cases    */
                      /* where 'int' is a 2-byte integer, in which case */
                      /* INT32BIT would need to be defined as 'long'.   */

#define BYTE_IMG      8  /* BITPIX code values for FITS image types */
#define SHORT_IMG    16
#define LONG_IMG     32
#define LONGLONG_IMG 64
#define FLOAT_IMG   -32
#define DOUBLE_IMG  -64

/* adopt a hopefully obscure number to use as a null value flag */
/* could be problems if the FITS files contain data with these values */
#define FLOATNULLVALUE -9.11912E-36F
#define DOUBLENULLVALUE -9.1191291391491E-36
 
/* Image compression algorithm types */
#define MAX_COMPRESS_DIM     6
#define RICE_1      11
#define GZIP_1      21
#define PLIO_1      31
#define HCOMPRESS_1 41
#define NOCOMPRESS  0

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif


typedef struct      /* structure used to store basic FITS file information */
{
    unsigned char** data;  /* Compressed data */
    int* dataLen;          /* Length of compressed data */
    double* c_zscale;        /* Scale factors for zscale column */
    double* c_zzero;         /* Zero values for zzero column */
    void** ucData;           /* Uncompressed data for failure to quantize */
    int* ucDataLen;          /* Length of uncompressed data */
    int compress_type;      /* type of compression algorithm */
    int zbitpix;            /* FITS data type of image (BITPIX) */
    int zndim;              /* dimension of image */
    long znaxis[MAX_COMPRESS_DIM];  /* length of each axis */
    long tilesize[MAX_COMPRESS_DIM]; /* size of compression tiles */
    long maxtilelen;        /* max number of pixels in each image tile */
    long maxelem;		/* maximum length of variable length arrays */
    int cn_uncompressed;    /* column number for UNCOMPRESSED_DATA column */
    int cn_zscale;	    /* column number for ZSCALE column */
    int cn_zzero;	    /* column number for ZZERO column */
    int cn_zblank;          /* column number for the ZBLANK column */
    double* bscale;        /* BSCALE value for each tile in the image */
    double* bzero;         /* BZERO value for each tile in the image */
    double zscale;          /* scaling value, if same for all tiles */
    double zzero;           /* zero pt, if same for all tiles */
    double cn_bscale;       /* value of the BSCALE keyword in header */
    double cn_bzero;        /* value of the BZERO keyword in header */
    long* blank;          /* value for null pixels for each tile in the image */
    int zblank;             /* value for null pixels, if not a column */
    int rice_blocksize;     /* first compression parameter */
    int rice_bytepix;       /* 2nd compression parameter: bytes/pixel */
    float quantize_level;   /* floating point quantization level */
    float hcomp_scale;      /* 1st hcompress compression parameter */
    int hcomp_smooth;       /* 2nd hcompress compression parameter */
} FITSfile;

typedef struct         /* structure used to store basic HDU information */
{
    FITSfile *Fptr;   /* pointer to FITS file structure */
}fitsfile;

/* error status codes */

#define OVERFLOW_ERR      -11  /* overflow during datatype conversion */
#define MEMORY_ALLOCATION 113  /* Could not allocate memory */
#define BAD_COL_NUM       302  /* column number < 1 or > tfields */
#define BAD_PIX_NUM       321  /* first pixel number greater than last pixel */
#define NEG_AXIS          323  /* illegal axis length < 1 */
#define BAD_DATATYPE      410  /* bad keyword datatype code */
# define DATA_COMPRESSION_ERR 413  /* error in imcompress routines */
# define DATA_DECOMPRESSION_ERR 414 /* error in imcompress routines */
# define NO_COMPRESSED_TILE  415 /* compressed tile doesn't exist */

/* ======================================================================= */
/* The following logic is used to determine the type machine,              */
/*  whether the bytes are swapped, and the number of bits in a long value  */
/* ======================================================================= */

/*   The following platforms have sizeof(long) == 8               */
/*   This block of code should match a similar block in fitsio.h  */
/*   and the block of code at the beginning of f77_wrap.h         */

#if defined(__alpha) && ( defined(__unix__) || defined(__NetBSD__) )
                                  /* old Dec Alpha platforms running OSF */
#define BYTESWAPPED TRUE
#define LONGSIZE 64

#elif defined(__sparcv9)
                               /*  SUN Solaris7 in 64-bit mode */
#define BYTESWAPPED FALSE
#define MACHINE NATIVE
#define LONGSIZE 64

#elif defined(__ia64__)  || defined(__x86_64__)
                  /*  Intel itanium 64-bit PC, or AMD opteron 64-bit PC */
#define BYTESWAPPED TRUE
#define LONGSIZE 64

#elif defined(_SX)             /* Nec SuperUx */

#define BYTESWAPPED FALSE
#define MACHINE NATIVE
#define LONGSIZE 64

#elif defined(__powerpc64__) || defined(__64BIT__) /* IBM 64-bit AIX powerpc*/
                              /* could also test for __ppc64__ or __PPC64 */
#define BYTESWAPPED FALSE
#define MACHINE NATIVE
#define LONGSIZE 64

#elif defined(_MIPS_SZLONG)

#  if defined(MIPSEL)
#    define BYTESWAPPED TRUE
#  else
#    define BYTESWAPPED FALSE
#    define MACHINE NATIVE
#  endif

#  if _MIPS_SZLONG == 32
#    define LONGSIZE 32
#  elif _MIPS_SZLONG == 64
#    define LONGSIZE 64
#  else
#    error "can't handle long size given by _MIPS_SZLONG"
#  endif

/* ============================================================== */
/*  the following are all 32-bit byteswapped platforms            */

#elif defined(vax) && defined(VMS)

#define MACHINE VAXVMS
#define BYTESWAPPED TRUE

#elif defined(__alpha) && defined(__VMS)

#if (__D_FLOAT == TRUE)

/* this float option is the same as for VAX/VMS machines. */
#define MACHINE VAXVMS
#define BYTESWAPPED TRUE

#elif  (__G_FLOAT == TRUE)

/*  G_FLOAT is the default for ALPHA VMS systems */
#define MACHINE ALPHAVMS
#define BYTESWAPPED TRUE
#define FLOATTYPE GFLOAT

#elif  (__IEEE_FLOAT == TRUE)

#define MACHINE ALPHAVMS
#define BYTESWAPPED TRUE
#define FLOATTYPE IEEEFLOAT

#endif  /* end of alpha VMS case */

#elif defined(ultrix) && defined(unix)
 /* old Dec ultrix machines */
#define BYTESWAPPED TRUE

#elif defined(__i386) || defined(__i386__) || defined(__i486__) || defined(__i586__) \
  || defined(_MSC_VER) || defined(__BORLANDC__) || defined(__TURBOC__) \
  || defined(_NI_mswin_) || defined(__EMX__)

/*  generic 32-bit IBM PC */
#define MACHINE IBMPC
#define BYTESWAPPED TRUE

#elif defined(__arm__)

/* This assumes all ARM are little endian.  In the future, it might be  */
/* necessary to use  "if defined(__ARMEL__)"  to distinguish little from big. */
/* (__ARMEL__ would be defined on little-endian, but not on big-endian). */

#define BYTESWAPPED TRUE

#else

/*  assume all other machine uses the same IEEE formats as used in FITS files */
/*  e.g., Macs fall into this category  */

#define MACHINE NATIVE
#define BYTESWAPPED FALSE

#endif

#ifndef MACHINE
#define MACHINE  OTHERTYPE
#endif

/*  assume longs are 4 bytes long, unless previously set otherwise */
#ifndef LONGSIZE
#define LONGSIZE 32
#endif

/*       end of block that determine long size and byte swapping        */
/* ==================================================================== */

#define maxvalue(A,B) ((A) > (B) ? (A) : (B))
#define minvalue(A,B) ((A) < (B) ? (A) : (B))

#define DSCHAR_MAX  127.49 /* max double value that fits in an signed char */
#define DSCHAR_MIN -128.49 /* min double value that fits in an signed char */
#define DUCHAR_MAX  255.49 /* max double value that fits in an unsigned char */
#define DUCHAR_MIN -0.49   /* min double value that fits in an unsigned char */
#define DUSHRT_MAX  65535.49 /* max double value that fits in a unsigned short*/
#define DUSHRT_MIN -0.49   /* min double value that fits in an unsigned short */
#define DSHRT_MAX  32767.49 /* max double value that fits in a short */
#define DSHRT_MIN -32768.49 /* min double value that fits in a short */

#if LONGSIZE == 32
#  define DLONG_MAX  2147483647.49 /* max double value that fits in a long */
#  define DLONG_MIN -2147483648.49 /* min double value that fits in a long */
#  define DULONG_MAX 4294967295.49 /* max double that fits in a unsigned long */
#else
#  define DLONG_MAX   9.2233720368547752E18 /* max double value  long */
#  define DLONG_MIN  -9.2233720368547752E18 /* min double value  long */
#  define DULONG_MAX 1.84467440737095504E19 /* max double value  ulong */
#endif

#define DULONG_MIN -0.49   /* min double value that fits in an unsigned long */

#define DUINT_MAX 4294967295.49 /* max dbl that fits in a unsigned 4-byte int */
#define DUINT_MIN -0.49   /* min dbl that fits in an unsigned 4-byte int */
#define DINT_MAX  2147483647.49 /* max double value that fits in a 4-byte int */
#define DINT_MIN -2147483648.49 /* min double value that fits in a 4-byte int */

#ifndef INT32_MAX
#define INT32_MAX  2147483647 /* max 32-bit integer */
#endif
#ifndef INT32_MIN
#define INT32_MIN (-INT32_MAX -1) /* min 32-bit integer */
#endif

/*---------------- utility routines -------------*/

void _pyfits_ffpmsg(const char *err_message);
int  _pyfits_ffgmsg(char *err_message);
int  _pyfits_ffgpv(fitsfile *fptr, int  datatype, LONGLONG firstelem,
          LONGLONG nelem, void *nulval, void *array, int *anynul, int  *status);
int  _pyfits_ffppr(fitsfile *fptr, int datatype, LONGLONG  firstelem,
          LONGLONG nelem, void *array, int *status);

/*--------------------- group template parser routines ------------------*/

int _pyfits_fits_img_stats_int(int *array,long nx, long ny, int nullcheck,
    int nullvalue,long *ngoodpix, int *minvalue, int *maxvalue, double *mean,
    double *sigma, double *noise1, double *noise3, int *status);

/* H-compress routines */
int _pyfits_fits_hcompress(int *a, int nx, int ny, int scale, char *output, 
    long *nbytes, int *status);
int _pyfits_fits_hcompress64(LONGLONG *a, int nx, int ny, int scale,
    char *output, long *nbytes, int *status);
int _pyfits_fits_hdecompress(unsigned char *input, int smooth, int *a, int *nx, 
    int *ny, int *scale, int *status);
int _pyfits_fits_hdecompress64(unsigned char *input, int smooth, LONGLONG *a,
    int *nx, int *ny, int *scale, int *status);


/*  image compression routines */

int _pyfits_imcomp_calc_max_elem (int comptype, int nx, int zbitpix,
                                  int blocksize);


/*  image decompression routines */
int _pyfits_fits_quantize_float (float fdata[], long nx, long ny, int nullcheck,
         float in_null_value,
         float quantize_level, int idata[], double *bscale, double *bzero,
         int *iminval, int *imaxval);
int _pyfits_fits_quantize_double (double fdata[], long nx, long ny,
         int nullcheck, double in_null_value,
         float quantize_level, int idata[], double *bscale, double *bzero,
         int *iminval, int *imaxval);
int _pyfits_fits_rcomp(int a[], int nx, unsigned char *c, int clen,int nblock);
int _pyfits_fits_rcomp_short(short a[], int nx, unsigned char *c, int clen,
                             int nblock);
int _pyfits_fits_rcomp_byte(signed char a[], int nx, unsigned char *c, int clen,
                            int nblock);
int _pyfits_fits_rdecomp (unsigned char *c, int clen, unsigned int array[],
                          int nx, int nblock);
int _pyfits_fits_rdecomp_short (unsigned char *c, int clen,
                                unsigned short array[], int nx, int nblock);
int _pyfits_fits_rdecomp_byte (unsigned char *c, int clen,
                               unsigned char array[], int nx, int nblock);
int _pyfits_pl_p2li (int *pxsrc, int xs, short *lldst, int npix);
int _pyfits_pl_l2pi (short *ll_src, int xs, int *px_dst, int npix);

int _pyfits_uncompress2mem_from_mem(
             char *inmemptr,
             size_t inmemsize,
             char **buffptr,
             size_t *buffsize,
             void *(*mem_realloc)(void *p, size_t newsize),
             size_t *filesize,
             int *status);

int _pyfits_compress2mem_from_mem(
             char *inmemptr,
             size_t inmemsize,
             char **buffptr,
             size_t *buffsize,
             void *(*mem_realloc)(void *p, size_t newsize),
             size_t *filesize,
             int *status);

/* Translate the long names for some routines to their actual short names */

#define _pyfits_fits_read_img  _pyfits_ffgpv
#define _pyfits_fits_write_img _pyfits_ffppr
#define fits_write_img_usht    ffpprui

#endif

