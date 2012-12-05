/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */
/*

Copyright (Unpublished--all rights reserved under the copyright laws of
the United States), U.S. Government as represented by the Administrator
of the National Aeronautics and Space Administration.  No copyright is
claimed in the United States under Title 17, U.S. Code.

Permission to freely use, copy, modify, and distribute this software
and its documentation without fee is hereby granted, provided that this
copyright notice and disclaimer of warranty appears in all copies.

DISCLAIMER:

THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND,
EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO,
ANY WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY
IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE, AND FREEDOM FROM INFRINGEMENT, AND ANY WARRANTY THAT THE
DOCUMENTATION WILL CONFORM TO THE SOFTWARE, OR ANY WARRANTY THAT THE
SOFTWARE WILL BE ERROR FREE.  IN NO EVENT SHALL NASA BE LIABLE FOR ANY
DAMAGES, INCLUDING, BUT NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR
CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN ANY WAY
CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY,
CONTRACT, TORT , OR OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY
PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS SUSTAINED
FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR
SERVICES PROVIDED HEREUNDER."

*/

#ifndef _FITSIO_H
#define _FITSIO_H

#define CFITSIO_VERSION 3.30
#define CFITSIO_MINOR 30
#define CFITSIO_MAJOR 3

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
#elif defined(_MSC_VER) &&  (_MSC_VER>= 1400)
#    define OFF_T long long
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
    ||  defined(__sparcv9) || (defined(__sparc__) && defined(__arch64__))  \
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

#ifndef LONGLONG_MAX

#ifdef LLONG_MAX
/* Linux and Solaris definition */
#define LONGLONG_MAX LLONG_MAX
#define LONGLONG_MIN LLONG_MIN

#elif defined(LONG_LONG_MAX)
#define LONGLONG_MAX LONG_LONG_MAX
#define LONGLONG_MIN LONG_LONG_MIN

#elif defined(__LONG_LONG_MAX__)
/* Mac OS X & CYGWIN defintion */
#define LONGLONG_MAX __LONG_LONG_MAX__
#define LONGLONG_MIN (-LONGLONG_MAX -1LL)

#elif defined(INT64_MAX)
/* windows definition */
#define LONGLONG_MAX INT64_MAX
#define LONGLONG_MIN INT64_MIN

#elif defined(_I64_MAX)
/* windows definition */
#define LONGLONG_MAX _I64_MAX
#define LONGLONG_MIN _I64_MIN

#elif (defined(__alpha) && ( defined(__unix__) || defined(__NetBSD__) )) \
    ||  defined(__sparcv9)  \
    ||  defined(__ia64__)   \
    ||  defined(__x86_64__) \
    ||  defined(_SX)        \
    ||  defined(__powerpc64__) || defined(__64BIT__) \
    ||  (defined(_MIPS_SZLONG) &&  _MIPS_SZLONG == 64)
/* sizeof(long) = 64 */
#define LONGLONG_MAX  9223372036854775807L /* max 64-bit integer */
#define LONGLONG_MIN (-LONGLONG_MAX -1L)   /* min 64-bit integer */

#else
/*  define a default value, even if it is never used */
#define LONGLONG_MAX  9223372036854775807LL /* max 64-bit integer */
#define LONGLONG_MIN (-LONGLONG_MAX -1LL)   /* min 64-bit integer */

#endif
#endif  /* end of ndef LONGLONG_MAX section */


/* ================================================================= */


/*  The following exclusion if __CINT__ is defined is needed for ROOT */
#ifndef __CINT__
#include "longnam.h"
#endif
 
#define NIOBUF  40  /* number of IO buffers to create (default = 40) */
          /* !! Significantly increasing NIOBUF may degrade performance !! */

#define IOBUFLEN 2880    /* size in bytes of each IO buffer (DONT CHANGE!) */

/* global variables */
 
#define FLEN_FILENAME 1025 /* max length of a filename  */
#define FLEN_KEYWORD   72  /* max length of a keyword (HIERARCH convention) */
#define FLEN_CARD      81  /* length of a FITS header card */
#define FLEN_VALUE     71  /* max length of a keyword value string */
#define FLEN_COMMENT   73  /* max length of a keyword comment string */
#define FLEN_ERRMSG    81  /* max length of a FITSIO error message */
#define FLEN_STATUS    31  /* max length of a FITSIO status text string */
 
#define TBIT          1  /* codes for FITS table data types */
#define TBYTE        11
#define TSBYTE       12
#define TLOGICAL     14
#define TSTRING      16
#define TUSHORT      20
#define TSHORT       21
#define TUINT        30
#define TINT         31
#define TULONG       40
#define TLONG        41
#define TINT32BIT    41  /* used when returning datatype of a column */
#define TFLOAT       42
#define TLONGLONG    81
#define TDOUBLE      82
#define TCOMPLEX     83
#define TDBLCOMPLEX 163

#define TYP_STRUC_KEY 10
#define TYP_CMPRS_KEY 20
#define TYP_SCAL_KEY  30
#define TYP_NULL_KEY  40
#define TYP_DIM_KEY   50
#define TYP_RANG_KEY  60
#define TYP_UNIT_KEY  70
#define TYP_DISP_KEY  80
#define TYP_HDUID_KEY 90
#define TYP_CKSUM_KEY 100
#define TYP_WCS_KEY   110
#define TYP_REFSYS_KEY 120
#define TYP_COMM_KEY  130
#define TYP_CONT_KEY  140
#define TYP_USER_KEY  150


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
                         /* The following 2 codes are not true FITS         */
                         /* datatypes; these codes are only used internally */
                         /* within cfitsio to make it easier for users      */
                         /* to deal with unsigned integers.                 */
#define SBYTE_IMG    10
#define USHORT_IMG   20
#define ULONG_IMG    40

#define IMAGE_HDU  0  /* Primary Array or IMAGE HDU */
#define ASCII_TBL  1  /* ASCII table HDU  */
#define BINARY_TBL 2  /* Binary table HDU */
#define ANY_HDU   -1  /* matches any HDU type */

#define READONLY  0    /* options when opening a file */
#define READWRITE 1

/* adopt a hopefully obscure number to use as a null value flag */
/* could be problems if the FITS files contain data with these values */
#define FLOATNULLVALUE -9.11912E-36F
#define DOUBLENULLVALUE -9.1191291391491E-36
 
/* compression algorithm type codes */
#define SUBTRACTIVE_DITHER_1 1
#define MAX_COMPRESS_DIM     6
#define RICE_1      11
#define GZIP_1      21
#define GZIP_2      22
#define PLIO_1      31
#define HCOMPRESS_1 41
#define BZIP2_1     51  /* not publicly supported; only for test purposes */
#define NOCOMPRESS  0

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#define CASESEN   1   /* do case-sensitive string match */
#define CASEINSEN 0   /* do case-insensitive string match */
 
#define GT_ID_ALL_URI  0   /* hierarchical grouping parameters */
#define GT_ID_REF      1
#define GT_ID_POS      2
#define GT_ID_ALL      3
#define GT_ID_REF_URI 11
#define GT_ID_POS_URI 12

#define OPT_RM_GPT      0
#define OPT_RM_ENTRY    1
#define OPT_RM_MBR      2
#define OPT_RM_ALL      3

#define OPT_GCP_GPT     0
#define OPT_GCP_MBR     1
#define OPT_GCP_ALL     2

#define OPT_MCP_ADD     0
#define OPT_MCP_NADD    1
#define OPT_MCP_REPL    2
#define OPT_MCP_MOV     3

#define OPT_MRG_COPY    0
#define OPT_MRG_MOV     1

#define OPT_CMT_MBR      1
#define OPT_CMT_MBR_DEL 11

typedef struct        /* structure used to store table column information */
{
    char ttype[70];   /* column name = FITS TTYPEn keyword; */
    LONGLONG tbcol;       /* offset in row to first byte of each column */
    int  tdatatype;   /* datatype code of each column */
    LONGLONG trepeat;    /* repeat count of column; number of elements */
    double tscale;    /* FITS TSCALn linear scaling factor */
    double tzero;     /* FITS TZEROn linear scaling zero point */
    LONGLONG tnull;   /* FITS null value for int image or binary table cols */
    char strnull[20]; /* FITS null value string for ASCII table columns */
    char tform[10];   /* FITS tform keyword value  */
    long  twidth;     /* width of each ASCII table column */
}tcolumn;

#define VALIDSTRUC 555  /* magic value used to identify if structure is valid */

typedef struct      /* structure used to store basic FITS file information */
{
    int filehandle;   /* handle returned by the file open function */
    int driver;       /* defines which set of I/O drivers should be used */
    int open_count;   /* number of opened 'fitsfiles' using this structure */
    char *filename;   /* file name */
    int validcode;    /* magic value used to verify that structure is valid */
    int only_one;     /* flag meaning only copy the specified extension */
    LONGLONG filesize; /* current size of the physical disk file in bytes */
    LONGLONG logfilesize; /* logical size of file, including unflushed buffers */
    int lasthdu;      /* is this the last HDU in the file? 0 = no, else yes */
    LONGLONG bytepos; /* current logical I/O pointer position in file */
    LONGLONG io_pos;  /* current I/O pointer position in the physical file */
    int curbuf;       /* number of I/O buffer currently in use */ 
    int curhdu;       /* current HDU number; 0 = primary array */
    int hdutype;      /* 0 = primary array, 1 = ASCII table, 2 = binary table */
    int writemode;    /* 0 = readonly, 1 = readwrite */
    int maxhdu;       /* highest numbered HDU known to exist in the file */
    int MAXHDU;       /* dynamically allocated dimension of headstart array */
    LONGLONG *headstart; /* byte offset in file to start of each HDU */
    LONGLONG headend;   /* byte offest in file to end of the current HDU header */
    LONGLONG ENDpos;    /* byte offest to where the END keyword was last written */
    LONGLONG nextkey;   /* byte offset in file to beginning of next keyword */
    LONGLONG datastart; /* byte offset in file to start of the current data unit */
    int imgdim;         /* dimension of image; cached for fast access */
    LONGLONG imgnaxis[99]; /* length of each axis; cached for fast access */
    int tfield;          /* number of fields in the table (primary array has 2 */
    LONGLONG origrows;   /* original number of rows (value of NAXIS2 keyword)  */
    LONGLONG numrows;    /* number of rows in the table (dynamically updated) */
    LONGLONG rowlength;  /* length of a table row or image size (bytes) */
    tcolumn *tableptr;   /* pointer to the table structure */
    LONGLONG heapstart;  /* heap start byte relative to start of data unit */
    LONGLONG heapsize;   /* size of the heap, in bytes */

         /* the following elements are related to compressed images */
    int request_compress_type;  /* requested image compression algorithm */
    long request_tilesize[MAX_COMPRESS_DIM]; /* requested tiling size */

    float request_hcomp_scale;    /* requested HCOMPRESS scale factor */
    int request_hcomp_smooth;     /* requested HCOMPRESS smooth parameter */
    int request_quantize_dither ; /* requested dithering mode when quantizing */
                                  /* floating point images to integer */
    int request_dither_offset;    /* starting offset into the array of random dithering */
    int request_lossy_int_compress;  /* lossy compress integer image as if float image? */

    int compressimg; /* 1 if HDU contains a compressed image, else 0 */
    int quantize_dither;   /* floating point pixel quantization algorithm */
    char zcmptype[12];      /* compression type string */
    int compress_type;      /* type of compression algorithm */
    int zbitpix;            /* FITS data type of image (BITPIX) */
    int zndim;              /* dimension of image */
    long znaxis[MAX_COMPRESS_DIM];  /* length of each axis */
    long tilesize[MAX_COMPRESS_DIM]; /* size of compression tiles */
    long maxtilelen;        /* max number of pixels in each image tile */
    long maxelem;	    /* maximum byte length of tile compressed arrays */

    int cn_compressed;	    /* column number for COMPRESSED_DATA column */
    int cn_uncompressed;    /* column number for UNCOMPRESSED_DATA column */
    int cn_gzip_data;       /* column number for GZIP2 lossless compressed data */
    int cn_zscale;	    /* column number for ZSCALE column */
    int cn_zzero;	    /* column number for ZZERO column */
    int cn_zblank;          /* column number for the ZBLANK column */

    double zscale;          /* scaling value, if same for all tiles */
    double zzero;           /* zero pt, if same for all tiles */
    double cn_bscale;       /* value of the BSCALE keyword in header */
    double cn_bzero;        /* value of the BZERO keyword (may be reset) */
    double cn_actual_bzero; /* actual value of the BZERO keyword  */
    int zblank;             /* value for null pixels, if not a column */

    int rice_blocksize;     /* first compression parameter: pixels/block */
    int rice_bytepix;       /* 2nd compression parameter: bytes/pixel */
    float quantize_level;   /* floating point quantization level */
    int dither_offset;      /* starting offset into the array of random dithering */
    float hcomp_scale;      /* 1st hcompress compression parameter */
    int hcomp_smooth;       /* 2nd hcompress compression parameter */

    int  tilerow;           /* row number of the uncompressed tiledata */
    long tiledatasize;       /* length of the tile data in bytes */
    int tiletype;           /* datatype of the tile (TINT, TSHORT, etc) */
    void *tiledata;         /* uncompressed tile of data, for row tilerow */
    char *tilenullarray;    /* optional array of null value flags */
    int tileanynull;        /* anynulls in this tile? */

    char *iobuffer;         /* pointer to FITS file I/O buffers */
    long bufrecnum[NIOBUF]; /* file record number of each of the buffers */
    int dirty[NIOBUF];     /* has the corresponding buffer been modified? */
    int ageindex[NIOBUF];  /* relative age of each buffer */  
} FITSfile;

typedef struct         /* structure used to store basic HDU information */
{
    int HDUposition;  /* HDU position in file; 0 = first HDU */
    FITSfile *Fptr;   /* pointer to FITS file structure */
}fitsfile;

typedef struct  /* structure for the iterator function column information */
{  
     /* elements required as input to fits_iterate_data: */

    fitsfile *fptr;     /* pointer to the HDU containing the column */
    int      colnum;    /* column number in the table (use name if < 1) */
    char     colname[70]; /* name (= TTYPEn value) of the column (optional) */
    int      datatype;  /* output datatype (converted if necessary  */
    int      iotype;    /* = InputCol, InputOutputCol, or OutputCol */

    /* output elements that may be useful for the work function: */

    void     *array;    /* pointer to the array (and the null value) */
    long     repeat;    /* binary table vector repeat value */
    long     tlmin;     /* legal minimum data value */
    long     tlmax;     /* legal maximum data value */
    char     tunit[70]; /* physical unit string */
    char     tdisp[70]; /* suggested display format */

} iteratorCol;

#define InputCol         0  /* flag for input only iterator column       */
#define InputOutputCol   1  /* flag for input and output iterator column */
#define OutputCol        2  /* flag for output only iterator column      */

/*=============================================================================
*
*       The following wtbarr typedef is used in the fits_read_wcstab() routine,
*       which is intended for use with the WCSLIB library written by Mark
*       Calabretta, http://www.atnf.csiro.au/~mcalabre/index.html
*
*       In order to maintain WCSLIB and CFITSIO as independent libraries it
*       was not permissible for any CFITSIO library code to include WCSLIB
*       header files, or vice versa.  However, the CFITSIO function
*       fits_read_wcstab() accepts an array of structs defined by wcs.h within
*       WCSLIB.  The problem then was to define this struct within fitsio.h
*       without including wcs.h, especially noting that wcs.h will often (but
*       not always) be included together with fitsio.h in an applications
*       program that uses fits_read_wcstab().
*
*       Of the various possibilities, the solution adopted was for WCSLIB to
*       define "struct wtbarr" while fitsio.h defines "typedef wtbarr", a
*       untagged struct with identical members.  This allows both wcs.h and
*       fitsio.h to define a wtbarr data type without conflict by virtue of
*       the fact that structure tags and typedef names share different
*       namespaces in C. Therefore, declarations within WCSLIB look like
*
*          struct wtbarr *w;
*
*       while within CFITSIO they are simply
*
*          wtbarr *w;
*
*       but as suggested by the commonality of the names, these are really the
*       same aggregate data type.  However, in passing a (struct wtbarr *) to
*       fits_read_wcstab() a cast to (wtbarr *) is formally required.
*===========================================================================*/

#ifndef WCSLIB_GETWCSTAB
#define WCSLIB_GETWCSTAB

typedef struct {
   int  i;                      /* Image axis number.                       */
   int  m;                      /* Array axis number for index vectors.     */
   int  kind;                   /* Array type, 'c' (coord) or 'i' (index).  */
   char extnam[72];             /* EXTNAME of binary table extension.       */
   int  extver;                 /* EXTVER  of binary table extension.       */
   int  extlev;                 /* EXTLEV  of binary table extension.       */
   char ttype[72];              /* TTYPEn of column containing the array.   */
   long row;                    /* Table row number.                        */
   int  ndim;                   /* Expected array dimensionality.           */
   int  *dimlen;                /* Where to write the array axis lengths.   */
   double **arrayp;             /* Where to write the address of the array  */
                                /* allocated to store the array.            */
} wtbarr;

int fits_read_wcstab(fitsfile *fptr, int nwtb, wtbarr *wtb, int *status);

#endif /* WCSLIB_GETWCSTAB */

/* error status codes */

#define CREATE_DISK_FILE -106 /* create disk file, without extended filename syntax */
#define OPEN_DISK_FILE   -105 /* open disk file, without extended filename syntax */
#define SKIP_TABLE       -104 /* move to 1st image when opening file */
#define SKIP_IMAGE       -103 /* move to 1st table when opening file */
#define SKIP_NULL_PRIMARY -102 /* skip null primary array when opening file */
#define USE_MEM_BUFF     -101  /* use memory buffer when opening file */
#define OVERFLOW_ERR      -11  /* overflow during datatype conversion */
#define PREPEND_PRIMARY    -9  /* used in ffiimg to insert new primary array */
#define SAME_FILE         101  /* input and output files are the same */
#define TOO_MANY_FILES    103  /* tried to open too many FITS files */
#define FILE_NOT_OPENED   104  /* could not open the named file */
#define FILE_NOT_CREATED  105  /* could not create the named file */
#define WRITE_ERROR       106  /* error writing to FITS file */
#define END_OF_FILE       107  /* tried to move past end of file */
#define READ_ERROR        108  /* error reading from FITS file */
#define FILE_NOT_CLOSED   110  /* could not close the file */
#define ARRAY_TOO_BIG     111  /* array dimensions exceed internal limit */
#define READONLY_FILE     112  /* Cannot write to readonly file */
#define MEMORY_ALLOCATION 113  /* Could not allocate memory */
#define BAD_FILEPTR       114  /* invalid fitsfile pointer */
#define NULL_INPUT_PTR    115  /* NULL input pointer to routine */
#define SEEK_ERROR        116  /* error seeking position in file */

#define BAD_URL_PREFIX    121  /* invalid URL prefix on file name */
#define TOO_MANY_DRIVERS  122  /* tried to register too many IO drivers */
#define DRIVER_INIT_FAILED 123  /* driver initialization failed */
#define NO_MATCHING_DRIVER 124  /* matching driver is not registered */
#define URL_PARSE_ERROR    125  /* failed to parse input file URL */
#define RANGE_PARSE_ERROR  126  /* failed to parse input file URL */

#define	SHARED_ERRBASE	(150)
#define	SHARED_BADARG	(SHARED_ERRBASE + 1)
#define	SHARED_NULPTR	(SHARED_ERRBASE + 2)
#define	SHARED_TABFULL	(SHARED_ERRBASE + 3)
#define	SHARED_NOTINIT	(SHARED_ERRBASE + 4)
#define	SHARED_IPCERR	(SHARED_ERRBASE + 5)
#define	SHARED_NOMEM	(SHARED_ERRBASE + 6)
#define	SHARED_AGAIN	(SHARED_ERRBASE + 7)
#define	SHARED_NOFILE	(SHARED_ERRBASE + 8)
#define	SHARED_NORESIZE	(SHARED_ERRBASE + 9)

#define HEADER_NOT_EMPTY  201  /* header already contains keywords */
#define KEY_NO_EXIST      202  /* keyword not found in header */
#define KEY_OUT_BOUNDS    203  /* keyword record number is out of bounds */
#define VALUE_UNDEFINED   204  /* keyword value field is blank */
#define NO_QUOTE          205  /* string is missing the closing quote */
#define BAD_INDEX_KEY     206  /* illegal indexed keyword name */
#define BAD_KEYCHAR       207  /* illegal character in keyword name or card */
#define BAD_ORDER         208  /* required keywords out of order */
#define NOT_POS_INT       209  /* keyword value is not a positive integer */
#define NO_END            210  /* couldn't find END keyword */
#define BAD_BITPIX        211  /* illegal BITPIX keyword value*/
#define BAD_NAXIS         212  /* illegal NAXIS keyword value */
#define BAD_NAXES         213  /* illegal NAXISn keyword value */
#define BAD_PCOUNT        214  /* illegal PCOUNT keyword value */
#define BAD_GCOUNT        215  /* illegal GCOUNT keyword value */
#define BAD_TFIELDS       216  /* illegal TFIELDS keyword value */
#define NEG_WIDTH         217  /* negative table row size */
#define NEG_ROWS          218  /* negative number of rows in table */
#define COL_NOT_FOUND     219  /* column with this name not found in table */
#define BAD_SIMPLE        220  /* illegal value of SIMPLE keyword  */
#define NO_SIMPLE         221  /* Primary array doesn't start with SIMPLE */
#define NO_BITPIX         222  /* Second keyword not BITPIX */
#define NO_NAXIS          223  /* Third keyword not NAXIS */
#define NO_NAXES          224  /* Couldn't find all the NAXISn keywords */
#define NO_XTENSION       225  /* HDU doesn't start with XTENSION keyword */
#define NOT_ATABLE        226  /* the CHDU is not an ASCII table extension */
#define NOT_BTABLE        227  /* the CHDU is not a binary table extension */
#define NO_PCOUNT         228  /* couldn't find PCOUNT keyword */
#define NO_GCOUNT         229  /* couldn't find GCOUNT keyword */
#define NO_TFIELDS        230  /* couldn't find TFIELDS keyword */
#define NO_TBCOL          231  /* couldn't find TBCOLn keyword */
#define NO_TFORM          232  /* couldn't find TFORMn keyword */
#define NOT_IMAGE         233  /* the CHDU is not an IMAGE extension */
#define BAD_TBCOL         234  /* TBCOLn keyword value < 0 or > rowlength */
#define NOT_TABLE         235  /* the CHDU is not a table */
#define COL_TOO_WIDE      236  /* column is too wide to fit in table */
#define COL_NOT_UNIQUE    237  /* more than 1 column name matches template */
#define BAD_ROW_WIDTH     241  /* sum of column widths not = NAXIS1 */
#define UNKNOWN_EXT       251  /* unrecognizable FITS extension type */
#define UNKNOWN_REC       252  /* unrecognizable FITS record */
#define END_JUNK          253  /* END keyword is not blank */
#define BAD_HEADER_FILL   254  /* Header fill area not blank */
#define BAD_DATA_FILL     255  /* Data fill area not blank or zero */
#define BAD_TFORM         261  /* illegal TFORM format code */
#define BAD_TFORM_DTYPE   262  /* unrecognizable TFORM datatype code */
#define BAD_TDIM          263  /* illegal TDIMn keyword value */
#define BAD_HEAP_PTR      264  /* invalid BINTABLE heap address */
 
#define BAD_HDU_NUM       301  /* HDU number < 1 or > MAXHDU */
#define BAD_COL_NUM       302  /* column number < 1 or > tfields */
#define NEG_FILE_POS      304  /* tried to move before beginning of file  */
#define NEG_BYTES         306  /* tried to read or write negative bytes */
#define BAD_ROW_NUM       307  /* illegal starting row number in table */
#define BAD_ELEM_NUM      308  /* illegal starting element number in vector */
#define NOT_ASCII_COL     309  /* this is not an ASCII string column */
#define NOT_LOGICAL_COL   310  /* this is not a logical datatype column */
#define BAD_ATABLE_FORMAT 311  /* ASCII table column has wrong format */
#define BAD_BTABLE_FORMAT 312  /* Binary table column has wrong format */
#define NO_NULL           314  /* null value has not been defined */
#define NOT_VARI_LEN      317  /* this is not a variable length column */
#define BAD_DIMEN         320  /* illegal number of dimensions in array */
#define BAD_PIX_NUM       321  /* first pixel number greater than last pixel */
#define ZERO_SCALE        322  /* illegal BSCALE or TSCALn keyword = 0 */
#define NEG_AXIS          323  /* illegal axis length < 1 */
 
#define NOT_GROUP_TABLE         340
#define HDU_ALREADY_MEMBER      341
#define MEMBER_NOT_FOUND        342
#define GROUP_NOT_FOUND         343
#define BAD_GROUP_ID            344
#define TOO_MANY_HDUS_TRACKED   345
#define HDU_ALREADY_TRACKED     346
#define BAD_OPTION              347
#define IDENTICAL_POINTERS      348
#define BAD_GROUP_ATTACH        349
#define BAD_GROUP_DETACH        350

#define BAD_I2C           401  /* bad int to formatted string conversion */
#define BAD_F2C           402  /* bad float to formatted string conversion */
#define BAD_INTKEY        403  /* can't interprete keyword value as integer */
#define BAD_LOGICALKEY    404  /* can't interprete keyword value as logical */
#define BAD_FLOATKEY      405  /* can't interprete keyword value as float */
#define BAD_DOUBLEKEY     406  /* can't interprete keyword value as double */
#define BAD_C2I           407  /* bad formatted string to int conversion */
#define BAD_C2F           408  /* bad formatted string to float conversion */
#define BAD_C2D           409  /* bad formatted string to double conversion */
#define BAD_DATATYPE      410  /* bad keyword datatype code */
#define BAD_DECIM         411  /* bad number of decimal places specified */
#define NUM_OVERFLOW      412  /* overflow during datatype conversion */

# define DATA_COMPRESSION_ERR 413  /* error in imcompress routines */
# define DATA_DECOMPRESSION_ERR 414 /* error in imcompress routines */
# define NO_COMPRESSED_TILE  415 /* compressed tile doesn't exist */

#define BAD_DATE          420  /* error in date or time conversion */

#define PARSE_SYNTAX_ERR  431  /* syntax error in parser expression */
#define PARSE_BAD_TYPE    432  /* expression did not evaluate to desired type */
#define PARSE_LRG_VECTOR  433  /* vector result too large to return in array */
#define PARSE_NO_OUTPUT   434  /* data parser failed not sent an out column */
#define PARSE_BAD_COL     435  /* bad data encounter while parsing column */
#define PARSE_BAD_OUTPUT  436  /* Output file not of proper type          */

#define ANGLE_TOO_BIG     501  /* celestial angle too large for projection */
#define BAD_WCS_VAL       502  /* bad celestial coordinate or pixel value */
#define WCS_ERROR         503  /* error in celestial coordinate calculation */
#define BAD_WCS_PROJ      504  /* unsupported type of celestial projection */
#define NO_WCS_KEY        505  /* celestial coordinate keywords not found */
#define APPROX_WCS_KEY    506  /* approximate WCS keywords were calculated */

#define NO_CLOSE_ERROR    999  /* special value used internally to switch off */
                               /* the error message from ffclos and ffchdu */

/*------- following error codes are used in the grparser.c file -----------*/
#define	NGP_ERRBASE		(360)			/* base chosen so not to interfere with CFITSIO */
#define	NGP_OK			(0)
#define	NGP_NO_MEMORY		(NGP_ERRBASE + 0)	/* malloc failed */
#define	NGP_READ_ERR		(NGP_ERRBASE + 1)	/* read error from file */
#define	NGP_NUL_PTR		(NGP_ERRBASE + 2)	/* null pointer passed as argument */
#define	NGP_EMPTY_CURLINE	(NGP_ERRBASE + 3)	/* line read seems to be empty */
#define	NGP_UNREAD_QUEUE_FULL	(NGP_ERRBASE + 4)	/* cannot unread more then 1 line (or single line twice) */
#define	NGP_INC_NESTING		(NGP_ERRBASE + 5)	/* too deep include file nesting (inf. loop ?) */
#define	NGP_ERR_FOPEN		(NGP_ERRBASE + 6)	/* fopen() failed, cannot open file */
#define	NGP_EOF			(NGP_ERRBASE + 7)	/* end of file encountered */
#define	NGP_BAD_ARG		(NGP_ERRBASE + 8)	/* bad arguments passed */
#define	NGP_TOKEN_NOT_EXPECT	(NGP_ERRBASE + 9)	/* token not expected here */

/*  The following exclusion if __CINT__ is defined is needed for ROOT */
#ifndef __CINT__
/*  the following 3 lines are needed to support C++ compilers */
#ifdef __cplusplus
extern "C" {
#endif
#endif

int CFITS2Unit( fitsfile *fptr );
fitsfile* CUnit2FITS(int unit);

/*----------------  FITS file URL parsing routines -------------*/
int fits_get_token(char **ptr, char *delimiter, char *token, int *isanumber);
char *fits_split_names(char *list);
int ffiurl(  char *url,  char *urltype, char *infile,
                    char *outfile, char *extspec, char *rowfilter,
                    char *binspec, char *colspec, int *status);
int ffifile (char *url,  char *urltype, char *infile,
                    char *outfile, char *extspec, char *rowfilter,
                    char *binspec, char *colspec, char *pixfilter, int *status);
int ffrtnm(char *url, char *rootname, int *status);
int ffexist(const char *infile, int *exists, int *status);
int ffexts(char *extspec, int *extnum,  char *extname, int *extvers,
          int *hdutype, char *colname, char *rowexpress, int *status);
int ffextn(char *url, int *extension_num, int *status);
int ffurlt(fitsfile *fptr, char *urlType, int *status);
int ffbins(char *binspec, int *imagetype, int *haxis, 
                      char colname[4][FLEN_VALUE], double *minin,
                      double *maxin, double *binsizein,
                      char minname[4][FLEN_VALUE], char maxname[4][FLEN_VALUE],
                      char binname[4][FLEN_VALUE], double *weight, char *wtname,
                      int *recip, int *status);
int ffbinr(char **binspec, char *colname, double *minin, 
                        double *maxin, double *binsizein, char *minname,
                        char *maxname, char *binname, int *status);
int fits_copy_cell2image(fitsfile *fptr, fitsfile *newptr, char *colname,
                      long rownum, int *status);
int fits_copy_image2cell(fitsfile *fptr, fitsfile *newptr, char *colname,
                      long rownum, int copykeyflag, int *status);
int fits_copy_pixlist2image(fitsfile *infptr, fitsfile *outfptr, int firstkey,       /* I - first HDU record number to start with */
           int naxis, int *colnum, int *status);
int ffimport_file( char *filename, char **contents, int *status );
int ffrwrg( char *rowlist, LONGLONG maxrows, int maxranges, int *numranges,
      long *minrow, long *maxrow, int *status);
int ffrwrgll( char *rowlist, LONGLONG maxrows, int maxranges, int *numranges,
      LONGLONG *minrow, LONGLONG *maxrow, int *status);
/*----------------  FITS file I/O routines -------------*/
int fits_init_cfitsio(void);
int ffomem(fitsfile **fptr, const char *name, int mode, void **buffptr,
           size_t *buffsize, size_t deltasize,
           void *(*mem_realloc)(void *p, size_t newsize),
           int *status);
int ffopen(fitsfile **fptr, const char *filename, int iomode, int *status);
int ffopentest(double version, fitsfile **fptr, const char *filename, int iomode, int *status);

int ffdopn(fitsfile **fptr, const char *filename, int iomode, int *status);
int fftopn(fitsfile **fptr, const char *filename, int iomode, int *status);
int ffiopn(fitsfile **fptr, const char *filename, int iomode, int *status);
int ffdkopn(fitsfile **fptr, const char *filename, int iomode, int *status);
int ffreopen(fitsfile *openfptr, fitsfile **newfptr, int *status); 
int ffinit(  fitsfile **fptr, const char *filename, int *status);
int ffdkinit(fitsfile **fptr, const char *filename, int *status);
int ffimem(fitsfile **fptr,  void **buffptr,
           size_t *buffsize, size_t deltasize,
           void *(*mem_realloc)(void *p, size_t newsize),
           int *status);
int fftplt(fitsfile **fptr, const char *filename, const char *tempname,
           int *status);
int ffflus(fitsfile *fptr, int *status);
int ffflsh(fitsfile *fptr, int clearbuf, int *status);
int ffclos(fitsfile *fptr, int *status);
int ffdelt(fitsfile *fptr, int *status);
int ffflnm(fitsfile *fptr, char *filename, int *status);
int ffflmd(fitsfile *fptr, int *filemode, int *status);
int fits_delete_iraf_file(char *filename, int *status);

/*---------------- utility routines -------------*/

float ffvers(float *version);
void ffupch(char *string);
void ffgerr(int status, char *errtext);
void ffpmsg(const char *err_message);
void ffpmrk(void);
int  ffgmsg(char *err_message);
void ffcmsg(void);
void ffcmrk(void);
void ffrprt(FILE *stream, int status);
void ffcmps(char *templt, char *colname, int  casesen, int *match,
           int *exact);
int fftkey(const char *keyword, int *status);
int fftrec(char *card, int *status);
int ffnchk(fitsfile *fptr, int *status);
int ffkeyn(const char *keyroot, int value, char *keyname, int *status);
int ffnkey(int value, char *keyroot, char *keyname, int *status);
int ffgkcl(char *card);
int ffdtyp(char *cval, char *dtype, int *status);
int ffinttyp(char *cval, int *datatype, int *negative, int *status);
int ffpsvc(char *card, char *value, char *comm, int *status);
int ffgknm(char *card, char *name, int *length, int *status);
int ffgthd(char *tmplt, char *card, int *hdtype, int *status);
int fits_translate_keyword(char *inrec, char *outrec, char *patterns[][2],
          int npat, int n_value, int n_offset, int n_range, int *pat_num,
          int *i, int *j,  int *m, int *n, int *status);
int fits_translate_keywords(fitsfile *infptr, fitsfile *outfptr,
          int firstkey, char *patterns[][2],
          int npat, int n_value, int n_offset, int n_range, int *status);    
int ffasfm(char *tform, int *datacode, long *width, int *decim, int *status);
int ffbnfm(char *tform, int *datacode, long *repeat, long *width, int *status);
int ffbnfmll(char *tform, int *datacode, LONGLONG *repeat, long *width, int *status);
int ffgabc(int tfields, char **tform, int space, long *rowlen, long *tbcol,
           int *status);
int fits_get_section_range(char **ptr,long *secmin,long *secmax,long *incre,
              int *status);
/* ffmbyt should not normally be used in application programs, but it is
   defined here as a publicly available routine because there are a few
   rare cases where it is needed
*/ 
int ffmbyt(fitsfile *fptr, LONGLONG bytpos, int ignore_err, int *status);
/*----------------- write single keywords --------------*/
int ffpky(fitsfile *fptr, int datatype, const char *keyname, void *value,
          const char *comm, int *status);
int ffprec(fitsfile *fptr, const char *card, int *status);
int ffpcom(fitsfile *fptr, const char *comm, int *status);
int ffpunt(fitsfile *fptr, const char *keyname, char *unit, int *status);
int ffphis(fitsfile *fptr, const char *history, int *status);
int ffpdat(fitsfile *fptr, int *status);
int ffverifydate(int year, int month, int day, int *status);
int ffgstm(char *timestr, int *timeref, int *status);
int ffgsdt(int *day, int *month, int *year, int *status);
int ffdt2s(int year, int month, int day, char *datestr, int *status);
int fftm2s(int year, int month, int day, int hour, int minute, double second,
          int decimals, char *datestr, int *status);
int ffs2dt(char *datestr, int *year, int *month, int *day, int *status);
int ffs2tm(char *datestr, int *year, int *month, int *day, int *hour,
          int *minute, double *second, int *status);
int ffpkyu(fitsfile *fptr, const char *keyname, const char *comm, int *status);
int ffpkys(fitsfile *fptr, const char *keyname, char *value, const char *comm,int *status);
int ffpkls(fitsfile *fptr, const char *keyname, const char *value, const char *comm,int *status);
int ffplsw(fitsfile *fptr, int *status);
int ffpkyl(fitsfile *fptr, const char *keyname, int  value, const char *comm, int *status);
int ffpkyj(fitsfile *fptr, const char *keyname, LONGLONG value, const char *comm, int *status);
int ffpkyf(fitsfile *fptr, const char *keyname, float value, int decim, const char *comm,
          int *status);
int ffpkye(fitsfile *fptr, const char *keyname, float  value, int decim, const char *comm,
          int *status);
int ffpkyg(fitsfile *fptr, const char *keyname, double value, int decim, const char *comm,
          int *status);
int ffpkyd(fitsfile *fptr, const char *keyname, double value, int decim, const char *comm,
          int *status);
int ffpkyc(fitsfile *fptr, const char *keyname, float *value, int decim, const char *comm,
          int *status);
int ffpkym(fitsfile *fptr, const char *keyname, double *value, int decim, const char *comm,
          int *status);
int ffpkfc(fitsfile *fptr, const char *keyname, float *value, int decim, const char *comm,
          int *status);
int ffpkfm(fitsfile *fptr, const char *keyname, double *value, int decim, const char *comm,
          int *status);
int ffpkyt(fitsfile *fptr, const char *keyname, long intval, double frac, const char *comm,
          int *status);
int ffptdm( fitsfile *fptr, int colnum, int naxis, long naxes[], int *status);
int ffptdmll( fitsfile *fptr, int colnum, int naxis, LONGLONG naxes[], int *status);

/*----------------- write array of keywords --------------*/
int ffpkns(fitsfile *fptr, const char *keyroot, int nstart, int nkey, char *value[],
           char *comm[], int *status);
int ffpknl(fitsfile *fptr, const char *keyroot, int nstart, int nkey, int *value,
           char *comm[], int *status);
int ffpknj(fitsfile *fptr, const char *keyroot, int nstart, int nkey, long *value,
           char *comm[], int *status);
int ffpknjj(fitsfile *fptr, const char *keyroot, int nstart, int nkey, LONGLONG *value,
           char *comm[], int *status);
int ffpknf(fitsfile *fptr, const char *keyroot, int nstart, int nkey, float *value,
           int decim, char *comm[], int *status);
int ffpkne(fitsfile *fptr, const char *keyroot, int nstart, int nkey, float *value,
           int decim, char *comm[], int *status);
int ffpkng(fitsfile *fptr, const char *keyroot, int nstart, int nkey, double *value,
           int decim, char *comm[], int *status);
int ffpknd(fitsfile *fptr, const char *keyroot, int nstart, int nkey, double *value,
           int decim, char *comm[], int *status);
int ffcpky(fitsfile *infptr,fitsfile *outfptr,int incol,int outcol,
           char *rootname, int *status); 

/*----------------- write required header keywords --------------*/
int ffphps( fitsfile *fptr, int bitpix, int naxis, long naxes[], int *status);
int ffphpsll( fitsfile *fptr, int bitpix, int naxis, LONGLONG naxes[], int *status);
int ffphpr( fitsfile *fptr, int simple, int bitpix, int naxis, long naxes[],
            LONGLONG pcount, LONGLONG gcount, int extend, int *status);
int ffphprll( fitsfile *fptr, int simple, int bitpix, int naxis, LONGLONG naxes[],
            LONGLONG pcount, LONGLONG gcount, int extend, int *status);
int ffphtb(fitsfile *fptr, LONGLONG naxis1, LONGLONG naxis2, int tfields, char **ttype,
          long *tbcol, char **tform, char **tunit, const char *extname, int *status);
int ffphbn(fitsfile *fptr, LONGLONG naxis2, int tfields, char **ttype,
          char **tform, char **tunit, const char *extname, LONGLONG pcount, int *status);
int ffphext( fitsfile *fptr, const char *xtension, int bitpix, int naxis, long naxes[],
            LONGLONG pcount, LONGLONG gcount, int *status);
/*----------------- write template keywords --------------*/
int ffpktp(fitsfile *fptr, const char *filename, int *status);

/*------------------ get header information --------------*/
int ffghsp(fitsfile *fptr, int *nexist, int *nmore, int *status);
int ffghps(fitsfile *fptr, int *nexist, int *position, int *status);
 
/*------------------ move position in header -------------*/
int ffmaky(fitsfile *fptr, int nrec, int *status);
int ffmrky(fitsfile *fptr, int nrec, int *status);
 
/*------------------ read single keywords -----------------*/
int ffgnxk(fitsfile *fptr, char **inclist, int ninc, char **exclist,
           int nexc, char *card, int  *status);
int ffgrec(fitsfile *fptr, int nrec,      char *card, int *status);
int ffgcrd(fitsfile *fptr, const char *keyname, char *card, int *status);
int ffgstr(fitsfile *fptr, const char *string, char *card, int *status);
int ffgunt(fitsfile *fptr, const char *keyname, char *unit, int  *status);
int ffgkyn(fitsfile *fptr, int nkey, char *keyname, char *keyval, char *comm,
           int *status);
int ffgkey(fitsfile *fptr, const char *keyname, char *keyval, char *comm,
           int *status);
 
int ffgky( fitsfile *fptr, int datatype, const char *keyname, void *value,
           char *comm, int *status);
int ffgkys(fitsfile *fptr, const char *keyname, char *value, char *comm, int *status);
int ffgkls(fitsfile *fptr, const char *keyname, char **value, char *comm, int *status);
int fffkls(char *value, int *status);
int ffgkyl(fitsfile *fptr, const char *keyname, int *value, char *comm, int *status);
int ffgkyj(fitsfile *fptr, const char *keyname, long *value, char *comm, int *status);
int ffgkyjj(fitsfile *fptr, const char *keyname, LONGLONG *value, char *comm, int *status);
int ffgkye(fitsfile *fptr, const char *keyname, float *value, char *comm,int *status);
int ffgkyd(fitsfile *fptr, const char *keyname, double *value,char *comm,int *status);
int ffgkyc(fitsfile *fptr, const char *keyname, float *value, char *comm,int *status);
int ffgkym(fitsfile *fptr, const char *keyname, double *value,char *comm,int *status);
int ffgkyt(fitsfile *fptr, const char *keyname, long *ivalue, double *dvalue,
           char *comm, int *status);
int ffgtdm(fitsfile *fptr, int colnum, int maxdim, int *naxis, long naxes[],
           int *status);
int ffgtdmll(fitsfile *fptr, int colnum, int maxdim, int *naxis, LONGLONG naxes[],
           int *status);
int ffdtdm(fitsfile *fptr, char *tdimstr, int colnum, int maxdim,
           int *naxis, long naxes[], int *status);
int ffdtdmll(fitsfile *fptr, char *tdimstr, int colnum, int maxdim,
           int *naxis, LONGLONG naxes[], int *status);

/*------------------ read array of keywords -----------------*/
int ffgkns(fitsfile *fptr, const char *keyname, int nstart, int nmax, char *value[],
           int *nfound,  int *status);
int ffgknl(fitsfile *fptr, const char *keyname, int nstart, int nmax, int *value,
           int *nfound, int *status);
int ffgknj(fitsfile *fptr, const char *keyname, int nstart, int nmax, long *value,
           int *nfound, int *status);
int ffgknjj(fitsfile *fptr, const char *keyname, int nstart, int nmax, LONGLONG *value,
           int *nfound, int *status);
int ffgkne(fitsfile *fptr, const char *keyname, int nstart, int nmax, float *value,
           int *nfound, int *status);
int ffgknd(fitsfile *fptr, const char *keyname, int nstart, int nmax, double *value,
           int *nfound, int *status);
int ffh2st(fitsfile *fptr, char **header, int  *status);
int ffhdr2str( fitsfile *fptr,  int exclude_comm, char **exclist,
   int nexc, char **header, int *nkeys, int  *status);
int ffcnvthdr2str( fitsfile *fptr,  int exclude_comm, char **exclist,
   int nexc, char **header, int *nkeys, int  *status);

/*----------------- read required header keywords --------------*/
int ffghpr(fitsfile *fptr, int maxdim, int *simple, int *bitpix, int *naxis,
          long naxes[], long *pcount, long *gcount, int *extend, int *status);
 
int ffghprll(fitsfile *fptr, int maxdim, int *simple, int *bitpix, int *naxis,
          LONGLONG naxes[], long *pcount, long *gcount, int *extend, int *status);

int ffghtb(fitsfile *fptr,int maxfield, long *naxis1, long *naxis2,
           int *tfields, char **ttype, long *tbcol, char **tform, char **tunit,
           char *extname,  int *status);

int ffghtbll(fitsfile *fptr,int maxfield, LONGLONG *naxis1, LONGLONG *naxis2,
           int *tfields, char **ttype, LONGLONG *tbcol, char **tform, char **tunit,
           char *extname,  int *status);
 
 
int ffghbn(fitsfile *fptr, int maxfield, long *naxis2, int *tfields,
           char **ttype, char **tform, char **tunit, char *extname,
           long *pcount, int *status);

int ffghbnll(fitsfile *fptr, int maxfield, LONGLONG *naxis2, int *tfields,
           char **ttype, char **tform, char **tunit, char *extname,
           LONGLONG *pcount, int *status);

/*--------------------- update keywords ---------------*/
int ffuky(fitsfile *fptr, int datatype, const char *keyname, void *value,
          char *comm, int *status);
int ffucrd(fitsfile *fptr, const char *keyname, char *card, int *status);
int ffukyu(fitsfile *fptr, const char *keyname, char *comm, int *status);
int ffukys(fitsfile *fptr, const char *keyname, char *value, char *comm, int *status);
int ffukls(fitsfile *fptr, const char *keyname, char *value, char *comm, int *status);
int ffukyl(fitsfile *fptr, const char *keyname, int value, char *comm, int *status);
int ffukyj(fitsfile *fptr, const char *keyname, LONGLONG value, char *comm, int *status);
int ffukyf(fitsfile *fptr, const char *keyname, float value, int decim, char *comm,
          int *status);
int ffukye(fitsfile *fptr, const char *keyname, float value, int decim, char *comm,
          int *status);
int ffukyg(fitsfile *fptr, const char *keyname, double value, int decim, char *comm,
          int *status);
int ffukyd(fitsfile *fptr, const char *keyname, double value, int decim, char *comm,
          int *status);
int ffukyc(fitsfile *fptr, const char *keyname, float *value, int decim, char *comm,
          int *status);
int ffukym(fitsfile *fptr, const char *keyname, double *value, int decim, char *comm,
          int *status);
int ffukfc(fitsfile *fptr, const char *keyname, float *value, int decim, char *comm,
          int *status);
int ffukfm(fitsfile *fptr, const char *keyname, double *value, int decim, char *comm,
          int *status);

/*--------------------- modify keywords ---------------*/
int ffmrec(fitsfile *fptr, int nkey, char *card, int *status);
int ffmcrd(fitsfile *fptr, const char *keyname, char *card, int *status);
int ffmnam(fitsfile *fptr, const char *oldname, const char *newname, int *status);
int ffmcom(fitsfile *fptr, const char *keyname, char *comm, int *status);
int ffmkyu(fitsfile *fptr, const char *keyname, char *comm, int *status);
int ffmkys(fitsfile *fptr, const char *keyname, char *value, char *comm,int *status);
int ffmkls(fitsfile *fptr, const char *keyname, char *value, char *comm,int *status);
int ffmkyl(fitsfile *fptr, const char *keyname, int value, char *comm, int *status);
int ffmkyj(fitsfile *fptr, const char *keyname, LONGLONG value, char *comm, int *status);
int ffmkyf(fitsfile *fptr, const char *keyname, float value, int decim, char *comm,
          int *status);
int ffmkye(fitsfile *fptr, const char *keyname, float value, int decim, char *comm,
          int *status);
int ffmkyg(fitsfile *fptr, const char *keyname, double value, int decim, char *comm,
          int *status);
int ffmkyd(fitsfile *fptr, const char *keyname, double value, int decim, char *comm,
          int *status);
int ffmkyc(fitsfile *fptr, const char *keyname, float *value, int decim, char *comm,
          int *status);
int ffmkym(fitsfile *fptr, const char *keyname, double *value, int decim, char *comm,
          int *status);
int ffmkfc(fitsfile *fptr, const char *keyname, float *value, int decim, char *comm,
          int *status);
int ffmkfm(fitsfile *fptr, const char *keyname, double *value, int decim, char *comm,
          int *status);
 
/*--------------------- insert keywords ---------------*/
int ffirec(fitsfile *fptr, int nkey, char *card, int *status);
int ffikey(fitsfile *fptr, char *card, int *status);
int ffikyu(fitsfile *fptr, const char *keyname, char *comm, int *status);
int ffikys(fitsfile *fptr, const char *keyname, char *value, char *comm,int *status);
int ffikls(fitsfile *fptr, const char *keyname, char *value, char *comm,int *status);
int ffikyl(fitsfile *fptr, const char *keyname, int value, char *comm, int *status);
int ffikyj(fitsfile *fptr, const char *keyname, LONGLONG value, char *comm, int *status);
int ffikyf(fitsfile *fptr, const char *keyname, float value, int decim, char *comm,
          int *status);
int ffikye(fitsfile *fptr, const char *keyname, float value, int decim, char *comm,
          int *status);
int ffikyg(fitsfile *fptr, const char *keyname, double value, int decim, char *comm,
          int *status);
int ffikyd(fitsfile *fptr, const char *keyname, double value, int decim, char *comm,
          int *status);
int ffikyc(fitsfile *fptr, const char *keyname, float *value, int decim, char *comm,
          int *status);
int ffikym(fitsfile *fptr, const char *keyname, double *value, int decim, char *comm,
          int *status);
int ffikfc(fitsfile *fptr, const char *keyname, float *value, int decim, char *comm,
          int *status);
int ffikfm(fitsfile *fptr, const char *keyname, double *value, int decim, char *comm,
          int *status);

/*--------------------- delete keywords ---------------*/
int ffdkey(fitsfile *fptr, const char *keyname, int *status);
int ffdstr(fitsfile *fptr, const char *string, int *status);
int ffdrec(fitsfile *fptr, int keypos, int *status);
 
/*--------------------- get HDU information -------------*/
int ffghdn(fitsfile *fptr, int *chdunum);
int ffghdt(fitsfile *fptr, int *exttype, int *status);
int ffghad(fitsfile *fptr, long *headstart, long *datastart, long *dataend,
           int *status);
int ffghadll(fitsfile *fptr, LONGLONG *headstart, LONGLONG *datastart,
           LONGLONG *dataend, int *status);
int ffghof(fitsfile *fptr, OFF_T *headstart, OFF_T *datastart, OFF_T *dataend,
           int *status);
int ffgipr(fitsfile *fptr, int maxaxis, int *imgtype, int *naxis,
           long *naxes, int *status);
int ffgiprll(fitsfile *fptr, int maxaxis, int *imgtype, int *naxis,
           LONGLONG *naxes, int *status);
int ffgidt(fitsfile *fptr, int *imgtype, int *status);
int ffgiet(fitsfile *fptr, int *imgtype, int *status);
int ffgidm(fitsfile *fptr, int *naxis,  int *status);
int ffgisz(fitsfile *fptr, int nlen, long *naxes, int *status);
int ffgiszll(fitsfile *fptr, int nlen, LONGLONG *naxes, int *status);

/*--------------------- HDU operations -------------*/
int ffmahd(fitsfile *fptr, int hdunum, int *exttype, int *status);
int ffmrhd(fitsfile *fptr, int hdumov, int *exttype, int *status);
int ffmnhd(fitsfile *fptr, int exttype, char *hduname, int hduvers,
           int *status);
int ffthdu(fitsfile *fptr, int *nhdu, int *status);
int ffcrhd(fitsfile *fptr, int *status);
int ffcrim(fitsfile *fptr, int bitpix, int naxis, long *naxes, int *status);
int ffcrimll(fitsfile *fptr, int bitpix, int naxis, LONGLONG *naxes, int *status);
int ffcrtb(fitsfile *fptr, int tbltype, LONGLONG naxis2, int tfields, char **ttype,
           char **tform, char **tunit, const char *extname, int *status);
int ffiimg(fitsfile *fptr, int bitpix, int naxis, long *naxes, int *status);
int ffiimgll(fitsfile *fptr, int bitpix, int naxis, LONGLONG *naxes, int *status);
int ffitab(fitsfile *fptr, LONGLONG naxis1, LONGLONG naxis2, int tfields, char **ttype,
           long *tbcol, char **tform, char **tunit, const char *extname, int *status);
int ffibin(fitsfile *fptr, LONGLONG naxis2, int tfields, char **ttype, char **tform,
           char **tunit, const char *extname, LONGLONG pcount, int *status);
int ffrsim(fitsfile *fptr, int bitpix, int naxis, long *naxes, int *status);
int ffrsimll(fitsfile *fptr, int bitpix, int naxis, LONGLONG *naxes, int *status);
int ffdhdu(fitsfile *fptr, int *hdutype, int *status);
int ffcopy(fitsfile *infptr, fitsfile *outfptr, int morekeys, int *status);
int ffcpfl(fitsfile *infptr, fitsfile *outfptr, int prev, int cur, int follow,
            int *status);
int ffcphd(fitsfile *infptr, fitsfile *outfptr, int *status);
int ffcpdt(fitsfile *infptr, fitsfile *outfptr, int *status);
int ffchfl(fitsfile *fptr, int *status);
int ffcdfl(fitsfile *fptr, int *status);
int ffwrhdu(fitsfile *fptr, FILE *outstream, int *status);

int ffrdef(fitsfile *fptr, int *status);
int ffhdef(fitsfile *fptr, int morekeys, int *status);
int ffpthp(fitsfile *fptr, long theap, int *status);
 
int ffcsum(fitsfile *fptr, long nrec, unsigned long *sum, int *status);
void ffesum(unsigned long sum, int complm, char *ascii);
unsigned long ffdsum(char *ascii, int complm, unsigned long *sum);
int ffpcks(fitsfile *fptr, int *status);
int ffupck(fitsfile *fptr, int *status);
int ffvcks(fitsfile *fptr, int *datastatus, int *hdustatus, int *status);
int ffgcks(fitsfile *fptr, unsigned long *datasum, unsigned long *hdusum,
    int *status);
 
/*--------------------- define scaling or null values -------------*/
int ffpscl(fitsfile *fptr, double scale, double zero, int *status);
int ffpnul(fitsfile *fptr, LONGLONG nulvalue, int *status);
int fftscl(fitsfile *fptr, int colnum, double scale, double zero, int *status);
int fftnul(fitsfile *fptr, int colnum, LONGLONG nulvalue, int *status);
int ffsnul(fitsfile *fptr, int colnum, char *nulstring, int *status);
 
/*--------------------- get column information -------------*/
int ffgcno(fitsfile *fptr, int casesen, char *templt, int  *colnum,
           int *status);
int ffgcnn(fitsfile *fptr, int casesen, char *templt, char *colname,
           int *colnum, int *status);
 
int ffgtcl(fitsfile *fptr, int colnum, int *typecode, long *repeat,
           long *width, int *status);
int ffgtclll(fitsfile *fptr, int colnum, int *typecode, LONGLONG *repeat,
           LONGLONG *width, int *status);
int ffeqty(fitsfile *fptr, int colnum, int *typecode, long *repeat,
           long *width, int *status);
int ffeqtyll(fitsfile *fptr, int colnum, int *typecode, LONGLONG *repeat,
           LONGLONG *width, int *status);
int ffgncl(fitsfile *fptr, int  *ncols, int *status);
int ffgnrw(fitsfile *fptr, long *nrows, int *status);
int ffgnrwll(fitsfile *fptr, LONGLONG *nrows, int *status);
int ffgacl(fitsfile *fptr, int colnum, char *ttype, long *tbcol,
           char *tunit, char *tform, double *tscal, double *tzero,
           char *tnull, char *tdisp, int *status);
int ffgbcl(fitsfile *fptr, int colnum, char *ttype, char *tunit,
           char *dtype, long *repeat, double *tscal, double *tzero,
           long *tnull, char *tdisp, int  *status);
int ffgbclll(fitsfile *fptr, int colnum, char *ttype, char *tunit,
           char *dtype, LONGLONG *repeat, double *tscal, double *tzero,
           LONGLONG *tnull, char *tdisp, int  *status);
int ffgrsz(fitsfile *fptr, long *nrows, int *status);
int ffgcdw(fitsfile *fptr, int colnum, int *width, int *status);

/*--------------------- read primary array or image elements -------------*/
int ffgpxv(fitsfile *fptr, int  datatype, long *firstpix, LONGLONG nelem,
          void *nulval, void *array, int *anynul, int *status);
int ffgpxvll(fitsfile *fptr, int  datatype, LONGLONG *firstpix, LONGLONG nelem,
          void *nulval, void *array, int *anynul, int *status);
int ffgpxf(fitsfile *fptr, int  datatype, long *firstpix, LONGLONG nelem,
           void *array, char *nullarray, int *anynul, int *status);
int ffgpxfll(fitsfile *fptr, int  datatype, LONGLONG *firstpix, LONGLONG nelem,
           void *array, char *nullarray, int *anynul, int *status);
int ffgsv(fitsfile *fptr, int datatype, long *blc, long *trc, long *inc,
          void *nulval, void *array, int *anynul, int  *status);

int ffgpv(fitsfile *fptr, int  datatype, LONGLONG firstelem, LONGLONG nelem,
          void *nulval, void *array, int *anynul, int  *status);
int ffgpf(fitsfile *fptr, int  datatype, LONGLONG firstelem, LONGLONG nelem,
          void *array, char *nullarray, int  *anynul, int  *status);
int ffgpvb(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem, unsigned
           char nulval, unsigned char *array, int *anynul, int *status);
int ffgpvsb(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem, signed
           char nulval, signed char *array, int *anynul, int *status);
int ffgpvui(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           unsigned short nulval, unsigned short *array, int *anynul, 
           int *status);
int ffgpvi(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           short nulval, short *array, int *anynul, int *status);
int ffgpvuj(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           unsigned long nulval, unsigned long *array, int *anynul, 
           int *status);
int ffgpvj(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           long nulval, long *array, int *anynul, int *status);
int ffgpvjj(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           LONGLONG nulval, LONGLONG *array, int *anynul, int *status);
int ffgpvuk(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           unsigned int nulval, unsigned int *array, int *anynul, int *status);
int ffgpvk(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           int nulval, int *array, int *anynul, int *status);
int ffgpve(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           float nulval, float *array, int *anynul, int *status);
int ffgpvd(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           double nulval, double *array, int *anynul, int *status);
 
int ffgpfb(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           unsigned char *array, char *nularray, int *anynul, int *status);
int ffgpfsb(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           signed char *array, char *nularray, int *anynul, int *status);
int ffgpfui(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           unsigned short *array, char *nularray, int *anynul, int *status);
int ffgpfi(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           short *array, char *nularray, int *anynul, int *status);
int ffgpfuj(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           unsigned long *array, char *nularray, int *anynul, int *status);
int ffgpfj(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           long *array, char *nularray, int *anynul, int *status);
int ffgpfjj(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           LONGLONG *array, char *nularray, int *anynul, int *status);
int ffgpfuk(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           unsigned int *array, char *nularray, int *anynul, int *status);
int ffgpfk(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           int *array, char *nularray, int *anynul, int *status);
int ffgpfe(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           float *array, char *nularray, int *anynul, int *status);
int ffgpfd(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           double *array, char *nularray, int *anynul, int *status);
 
int ffg2db(fitsfile *fptr, long group, unsigned char nulval, LONGLONG ncols,
           LONGLONG naxis1, LONGLONG naxis2, unsigned char *array,
           int *anynul, int *status);
int ffg2dsb(fitsfile *fptr, long group, signed char nulval, LONGLONG ncols,
           LONGLONG naxis1, LONGLONG naxis2, signed char *array,
           int *anynul, int *status);
int ffg2dui(fitsfile *fptr, long group, unsigned short nulval, LONGLONG ncols,
           LONGLONG naxis1, LONGLONG naxis2, unsigned short *array,
           int *anynul, int *status);
int ffg2di(fitsfile *fptr, long group, short nulval, LONGLONG ncols,
           LONGLONG naxis1, LONGLONG naxis2, short *array,
           int *anynul, int *status);
int ffg2duj(fitsfile *fptr, long group, unsigned long nulval, LONGLONG ncols,
           LONGLONG naxis1, LONGLONG naxis2, unsigned long *array,
           int *anynul, int *status);
int ffg2dj(fitsfile *fptr, long group, long nulval, LONGLONG ncols,
           LONGLONG naxis1, LONGLONG naxis2, long *array,
           int *anynul, int *status);
int ffg2djj(fitsfile *fptr, long group, LONGLONG nulval, LONGLONG ncols,
           LONGLONG naxis1, LONGLONG naxis2, LONGLONG *array,
           int *anynul, int *status);
int ffg2duk(fitsfile *fptr, long group, unsigned int nulval, LONGLONG ncols,
           LONGLONG naxis1, LONGLONG naxis2, unsigned int *array,
           int *anynul, int *status);
int ffg2dk(fitsfile *fptr, long group, int nulval, LONGLONG ncols,
           LONGLONG naxis1, LONGLONG naxis2, int *array,
           int *anynul, int *status);
int ffg2de(fitsfile *fptr, long group, float nulval, LONGLONG ncols,
           LONGLONG naxis1, LONGLONG naxis2, float *array,
           int *anynul, int *status);
int ffg2dd(fitsfile *fptr, long group, double nulval, LONGLONG ncols,
           LONGLONG naxis1, LONGLONG naxis2, double *array,
           int *anynul, int *status);
 
int ffg3db(fitsfile *fptr, long group, unsigned char nulval, LONGLONG ncols,
           LONGLONG nrows, LONGLONG naxis1, LONGLONG naxis2, LONGLONG naxis3,
           unsigned char *array, int *anynul, int *status);
int ffg3dsb(fitsfile *fptr, long group, signed char nulval, LONGLONG ncols,
           LONGLONG nrows, LONGLONG naxis1, LONGLONG naxis2, LONGLONG naxis3,
           signed char *array, int *anynul, int *status);
int ffg3dui(fitsfile *fptr, long group, unsigned short nulval, LONGLONG ncols,
           LONGLONG nrows, LONGLONG naxis1, LONGLONG naxis2, LONGLONG naxis3,
           unsigned short *array, int *anynul, int *status);
int ffg3di(fitsfile *fptr, long group, short nulval, LONGLONG ncols,
           LONGLONG nrows, LONGLONG naxis1, LONGLONG naxis2, LONGLONG naxis3,
           short *array, int *anynul, int *status);
int ffg3duj(fitsfile *fptr, long group, unsigned long nulval, LONGLONG ncols,
           LONGLONG nrows, LONGLONG naxis1, LONGLONG naxis2, LONGLONG naxis3,
           unsigned long *array, int *anynul, int *status);
int ffg3dj(fitsfile *fptr, long group, long nulval, LONGLONG ncols,
           LONGLONG nrows, LONGLONG naxis1, LONGLONG naxis2, LONGLONG naxis3,
           long *array, int *anynul, int *status);
int ffg3djj(fitsfile *fptr, long group, LONGLONG nulval, LONGLONG ncols,
           LONGLONG nrows, LONGLONG naxis1, LONGLONG naxis2, LONGLONG naxis3,
           LONGLONG *array, int *anynul, int *status);
int ffg3duk(fitsfile *fptr, long group, unsigned int nulval, LONGLONG ncols,
           LONGLONG nrows, LONGLONG naxis1, LONGLONG naxis2, LONGLONG naxis3,
           unsigned int *array, int *anynul, int *status);
int ffg3dk(fitsfile *fptr, long group, int nulval, LONGLONG ncols,
           LONGLONG nrows, LONGLONG naxis1, LONGLONG naxis2, LONGLONG naxis3,
           int *array, int *anynul, int *status);
int ffg3de(fitsfile *fptr, long group, float nulval, LONGLONG ncols,
           LONGLONG nrows, LONGLONG naxis1, LONGLONG naxis2, LONGLONG naxis3,
           float *array, int *anynul, int *status);
int ffg3dd(fitsfile *fptr, long group, double nulval, LONGLONG ncols,
           LONGLONG nrows, LONGLONG naxis1, LONGLONG naxis2, LONGLONG naxis3,
           double *array, int *anynul, int *status);
 
int ffgsvb(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, unsigned char nulval, unsigned char *array,
  int *anynul, int *status);
int ffgsvsb(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, signed char nulval, signed char *array,
  int *anynul, int *status);
int ffgsvui(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, unsigned short nulval, unsigned short *array, 
  int *anynul, int *status);
int ffgsvi(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, short nulval, short *array, int *anynul, int *status);
int ffgsvuj(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, unsigned long nulval, unsigned long *array, 
  int *anynul, int *status);
int ffgsvj(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, long nulval, long *array, int *anynul, int *status);
int ffgsvjj(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, LONGLONG nulval, LONGLONG *array, int *anynul,
  int *status);
int ffgsvuk(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, unsigned int nulval, unsigned int *array,
  int *anynul, int *status);
int ffgsvk(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, int nulval, int *array, int *anynul, int *status);
int ffgsve(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, float nulval, float *array, int *anynul, int *status);
int ffgsvd(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, double nulval, double *array, int *anynul,
  int *status);
 
int ffgsfb(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, unsigned char *array, char *flagval,
  int *anynul, int *status);
int ffgsfsb(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, signed char *array, char *flagval,
  int *anynul, int *status);
int ffgsfui(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, unsigned short *array, char *flagval, int *anynul, 
  int *status);
int ffgsfi(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, short *array, char *flagval, int *anynul, int *status);
int ffgsfuj(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long  *trc, long *inc, unsigned long *array, char *flagval, int *anynul,
  int *status);
int ffgsfj(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long  *trc, long *inc, long *array, char *flagval, int *anynul, int *status);
int ffgsfjj(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long  *trc, long *inc, LONGLONG *array, char *flagval, int *anynul,
  int *status);
int ffgsfuk(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long  *trc, long *inc, unsigned int *array, char *flagval, int *anynul,
  int *status);
int ffgsfk(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long  *trc, long *inc, int *array, char *flagval, int *anynul, int *status);
int ffgsfe(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, float *array, char *flagval, int *anynul, int *status);
int ffgsfd(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, double *array, char *flagval, int *anynul,
  int *status);
 
int ffggpb(fitsfile *fptr, long group, long firstelem, long nelem,
           unsigned char *array, int *status);
int ffggpsb(fitsfile *fptr, long group, long firstelem, long nelem,
           signed char *array, int *status);
int ffggpui(fitsfile *fptr, long group, long firstelem, long nelem,
           unsigned short *array, int *status);
int ffggpi(fitsfile *fptr, long group, long firstelem, long nelem,
           short *array, int *status);
int ffggpuj(fitsfile *fptr, long group, long firstelem, long nelem,
           unsigned long *array, int *status);
int ffggpj(fitsfile *fptr, long group, long firstelem, long nelem,
           long *array, int *status);
int ffggpjj(fitsfile *fptr, long group, long firstelem, long nelem,
           LONGLONG *array, int *status);
int ffggpuk(fitsfile *fptr, long group, long firstelem, long nelem,
           unsigned int *array, int *status);
int ffggpk(fitsfile *fptr, long group, long firstelem, long nelem,
           int *array, int *status);
int ffggpe(fitsfile *fptr, long group, long firstelem, long nelem,
           float *array, int *status);
int ffggpd(fitsfile *fptr, long group, long firstelem, long nelem,
           double *array, int *status);
 
/*--------------------- read column elements -------------*/
int ffgcv( fitsfile *fptr, int datatype, int colnum, LONGLONG firstrow,
           LONGLONG firstelem, LONGLONG nelem, void *nulval, void *array, int *anynul,
           int  *status);
int ffgcf( fitsfile *fptr, int datatype, int colnum, LONGLONG firstrow,
           LONGLONG firstelem, LONGLONG nelem, void *array, char *nullarray,
           int *anynul, int *status);
int ffgcvs(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, char *nulval, char **array, int *anynul, int *status);
int ffgcl (fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, char *array, int  *status);
int ffgcvl (fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, char nulval, char *array, int *anynul, int  *status);
int ffgcvb(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, unsigned char nulval, unsigned char *array,
           int *anynul, int *status);
int ffgcvsb(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, signed char nulval, signed char *array,
           int *anynul, int *status);
int ffgcvui(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, unsigned short nulval, unsigned short *array, 
           int *anynul, int *status);
int ffgcvi(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, short nulval, short *array, int *anynul, int *status);
int ffgcvuj(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, unsigned long nulval, unsigned long *array, int *anynul,
           int *status);
int ffgcvj(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, long nulval, long *array, int *anynul, int *status);
int ffgcvjj(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, LONGLONG nulval, LONGLONG *array, int *anynul,
           int *status);
int ffgcvuk(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, unsigned int nulval, unsigned int *array, int *anynul,
           int *status);
int ffgcvk(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, int nulval, int *array, int *anynul, int *status);
int ffgcve(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, float nulval, float *array, int *anynul, int *status);
int ffgcvd(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
         LONGLONG nelem, double nulval, double *array, int *anynul, int *status);
int ffgcvc(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, float nulval, float *array, int *anynul, int *status);
int ffgcvm(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
         LONGLONG nelem, double nulval, double *array, int *anynul, int *status);

int ffgcx(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstbit,
            LONGLONG nbits, char *larray, int *status);
int ffgcxui(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG nrows,
            long firstbit, int nbits, unsigned short *array, int *status);
int ffgcxuk(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG nrows,
            long firstbit, int nbits, unsigned int *array, int *status);

int ffgcfs(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem, 
      LONGLONG nelem, char **array, char *nularray, int *anynul, int *status);
int ffgcfl(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
      LONGLONG nelem, char *array, char *nularray, int *anynul, int *status);
int ffgcfb(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem, 
      LONGLONG nelem, unsigned char *array, char *nularray, int *anynul, int *status);
int ffgcfsb(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
      LONGLONG nelem, signed char *array, char *nularray, int *anynul, int *status);
int ffgcfui(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
      LONGLONG nelem, unsigned short *array, char *nularray, int *anynul, 
      int *status);
int ffgcfi(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
      LONGLONG nelem, short *array, char *nularray, int *anynul, int *status);
int ffgcfuj(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
      LONGLONG nelem, unsigned long *array, char *nularray, int *anynul,
      int *status);
int ffgcfj(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
      LONGLONG nelem, long *array, char *nularray, int *anynul, int *status);
int ffgcfjj(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
      LONGLONG nelem, LONGLONG *array, char *nularray, int *anynul, int *status);
int ffgcfuk(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
      LONGLONG nelem, unsigned int *array, char *nularray, int *anynul,
      int *status);
int ffgcfk(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
      LONGLONG nelem, int *array, char *nularray, int *anynul, int *status);
int ffgcfe(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
      LONGLONG nelem, float *array, char *nularray, int *anynul, int *status);
int ffgcfd(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
      LONGLONG nelem, double *array, char *nularray, int *anynul, int *status);
int ffgcfc(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
      LONGLONG nelem, float *array, char *nularray, int *anynul, int *status);
int ffgcfm(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
      LONGLONG nelem, double *array, char *nularray, int *anynul, int *status);
 
int ffgdes(fitsfile *fptr, int colnum, LONGLONG rownum, long *length,
           long *heapaddr, int *status);
int ffgdesll(fitsfile *fptr, int colnum, LONGLONG rownum, LONGLONG *length,
           LONGLONG *heapaddr, int *status);
int ffgdess(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG nrows, long *length,
           long *heapaddr, int *status);
int ffgdessll(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG nrows, LONGLONG *length,
           LONGLONG *heapaddr, int *status);
int ffpdes(fitsfile *fptr, int colnum, LONGLONG rownum, LONGLONG length,
           LONGLONG heapaddr, int *status);
int fftheap(fitsfile *fptr, LONGLONG *heapsize, LONGLONG *unused, LONGLONG *overlap,
            int *valid, int *status);
int ffcmph(fitsfile *fptr, int *status);

int ffgtbb(fitsfile *fptr, LONGLONG firstrow, LONGLONG firstchar, LONGLONG nchars,
           unsigned char *values, int *status);

int ffgextn(fitsfile *fptr, LONGLONG offset, LONGLONG nelem, void *array, int *status);
int ffpextn(fitsfile *fptr, LONGLONG offset, LONGLONG nelem, void *array, int *status);

/*------------ write primary array or image elements -------------*/
int ffppx(fitsfile *fptr, int datatype, long *firstpix, LONGLONG nelem,
          void *array, int *status);
int ffppxll(fitsfile *fptr, int datatype, LONGLONG *firstpix, LONGLONG nelem,
          void *array, int *status);
int ffppxn(fitsfile *fptr, int datatype, long *firstpix, LONGLONG nelem,
          void *array, void *nulval, int *status);
int ffppxnll(fitsfile *fptr, int datatype, LONGLONG *firstpix, LONGLONG nelem,
          void *array, void *nulval, int *status);
int ffppr(fitsfile *fptr, int datatype, LONGLONG  firstelem,
           LONGLONG nelem, void *array, int *status);
int ffpprb(fitsfile *fptr, long group, LONGLONG firstelem,
           LONGLONG nelem, unsigned char *array, int *status);
int ffpprsb(fitsfile *fptr, long group, LONGLONG firstelem,
           LONGLONG nelem, signed char *array, int *status);
int ffpprui(fitsfile *fptr, long group, LONGLONG firstelem,
           LONGLONG nelem, unsigned short *array, int *status);
int ffppri(fitsfile *fptr, long group, LONGLONG firstelem,
           LONGLONG nelem, short *array, int *status);
int ffppruj(fitsfile *fptr, long group, LONGLONG firstelem,
           LONGLONG nelem, unsigned long *array, int *status);
int ffpprj(fitsfile *fptr, long group, LONGLONG firstelem,
           LONGLONG nelem, long *array, int *status);
int ffppruk(fitsfile *fptr, long group, LONGLONG firstelem,
           LONGLONG nelem, unsigned int *array, int *status);
int ffpprk(fitsfile *fptr, long group, LONGLONG firstelem,
           LONGLONG nelem, int *array, int *status);
int ffppre(fitsfile *fptr, long group, LONGLONG firstelem,
           LONGLONG nelem, float *array, int *status);
int ffpprd(fitsfile *fptr, long group, LONGLONG firstelem,
           LONGLONG nelem, double *array, int *status);
int ffpprjj(fitsfile *fptr, long group, LONGLONG firstelem,
           LONGLONG nelem, LONGLONG *array, int *status);

int ffppru(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           int *status);
int ffpprn(fitsfile *fptr, LONGLONG firstelem, LONGLONG nelem, int *status);
 
int ffppn(fitsfile *fptr, int datatype, LONGLONG  firstelem, LONGLONG nelem,
          void  *array, void *nulval, int  *status);
int ffppnb(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           unsigned char *array, unsigned char nulval, int *status);
int ffppnsb(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           signed char *array, signed char nulval, int *status);
int ffppnui(fitsfile *fptr, long group, LONGLONG firstelem,
           LONGLONG nelem, unsigned short *array, unsigned short nulval,
           int *status);
int ffppni(fitsfile *fptr, long group, LONGLONG firstelem,
           LONGLONG nelem, short *array, short nulval, int *status);
int ffppnj(fitsfile *fptr, long group, LONGLONG firstelem,
           LONGLONG nelem, long *array, long nulval, int *status);
int ffppnuj(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           unsigned long *array, unsigned long nulval, int *status);
int ffppnuk(fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelem,
           unsigned int *array, unsigned int nulval, int *status);
int ffppnk(fitsfile *fptr, long group, LONGLONG firstelem,
           LONGLONG nelem, int *array, int nulval, int *status);
int ffppne(fitsfile *fptr, long group, LONGLONG firstelem,
           LONGLONG nelem, float *array, float nulval, int *status);
int ffppnd(fitsfile *fptr, long group, LONGLONG firstelem,
           LONGLONG nelem, double *array, double nulval, int *status);
int ffppnjj(fitsfile *fptr, long group, LONGLONG firstelem,
           LONGLONG nelem, LONGLONG *array, LONGLONG nulval, int *status);

int ffp2db(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG naxis1,
           LONGLONG naxis2, unsigned char *array, int *status);
int ffp2dsb(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG naxis1,
           LONGLONG naxis2, signed char *array, int *status);
int ffp2dui(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG naxis1,
           LONGLONG naxis2, unsigned short *array, int *status);
int ffp2di(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG naxis1,
           LONGLONG naxis2, short *array, int *status);
int ffp2duj(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG naxis1,
           LONGLONG naxis2, unsigned long *array, int *status);
int ffp2dj(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG naxis1,
           LONGLONG naxis2, long *array, int *status);
int ffp2duk(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG naxis1,
           LONGLONG naxis2, unsigned int *array, int *status);
int ffp2dk(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG naxis1,
           LONGLONG naxis2, int *array, int *status);
int ffp2de(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG naxis1,
           LONGLONG naxis2, float *array, int *status);
int ffp2dd(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG naxis1,
           LONGLONG naxis2, double *array, int *status);
int ffp2djj(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG naxis1,
           LONGLONG naxis2, LONGLONG *array, int *status);

int ffp3db(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG nrows, LONGLONG naxis1,
           LONGLONG naxis2, LONGLONG naxis3, unsigned char *array, int *status);
int ffp3dsb(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG nrows, LONGLONG naxis1,
           LONGLONG naxis2, LONGLONG naxis3, signed char *array, int *status);
int ffp3dui(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG nrows, LONGLONG naxis1,
           LONGLONG naxis2, LONGLONG naxis3, unsigned short *array, int *status);
int ffp3di(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG nrows, LONGLONG naxis1,
           LONGLONG naxis2, LONGLONG naxis3, short *array, int *status);
int ffp3duj(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG nrows, LONGLONG naxis1,
           LONGLONG naxis2, LONGLONG naxis3, unsigned long *array, int *status);
int ffp3dj(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG nrows, LONGLONG naxis1,
           LONGLONG naxis2, LONGLONG naxis3, long *array, int *status);
int ffp3duk(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG nrows, LONGLONG naxis1,
           LONGLONG naxis2, LONGLONG naxis3, unsigned int *array, int *status);
int ffp3dk(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG nrows, LONGLONG naxis1,
           LONGLONG naxis2, LONGLONG naxis3, int *array, int *status);
int ffp3de(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG nrows, LONGLONG naxis1,
           LONGLONG naxis2, LONGLONG naxis3, float *array, int *status);
int ffp3dd(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG nrows, LONGLONG naxis1,
           LONGLONG naxis2, LONGLONG naxis3, double *array, int *status);
int ffp3djj(fitsfile *fptr, long group, LONGLONG ncols, LONGLONG nrows, LONGLONG naxis1,
           LONGLONG naxis2, LONGLONG naxis3, LONGLONG *array, int *status);

int ffpss(fitsfile *fptr, int datatype,
           long *fpixel, long *lpixel, void *array, int *status);
int ffpssb(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, unsigned char *array, int *status);
int ffpsssb(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, signed char *array, int *status);
int ffpssui(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, unsigned short *array, int *status);
int ffpssi(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, short *array, int *status);
int ffpssuj(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, unsigned long *array, int *status);
int ffpssj(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, long *array, int *status);
int ffpssuk(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, unsigned int *array, int *status);
int ffpssk(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, int *array, int *status);
int ffpsse(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, float *array, int *status);
int ffpssd(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, double *array, int *status);
int ffpssjj(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, LONGLONG *array, int *status);

int ffpgpb(fitsfile *fptr, long group, long firstelem,
           long nelem, unsigned char *array, int *status);
int ffpgpsb(fitsfile *fptr, long group, long firstelem,
           long nelem, signed char *array, int *status);
int ffpgpui(fitsfile *fptr, long group, long firstelem,
           long nelem, unsigned short *array, int *status);
int ffpgpi(fitsfile *fptr, long group, long firstelem,
           long nelem, short *array, int *status);
int ffpgpuj(fitsfile *fptr, long group, long firstelem,
           long nelem, unsigned long *array, int *status);
int ffpgpj(fitsfile *fptr, long group, long firstelem,
           long nelem, long *array, int *status);
int ffpgpuk(fitsfile *fptr, long group, long firstelem,
           long nelem, unsigned int *array, int *status);
int ffpgpk(fitsfile *fptr, long group, long firstelem,
           long nelem, int *array, int *status);
int ffpgpe(fitsfile *fptr, long group, long firstelem,
           long nelem, float *array, int *status);
int ffpgpd(fitsfile *fptr, long group, long firstelem,
           long nelem, double *array, int *status);
int ffpgpjj(fitsfile *fptr, long group, long firstelem,
           long nelem, LONGLONG *array, int *status);

/*--------------------- iterator functions -------------*/
int fits_iter_set_by_name(iteratorCol *col, fitsfile *fptr, char *colname,
          int datatype,  int iotype);
int fits_iter_set_by_num(iteratorCol *col, fitsfile *fptr, int colnum,
          int datatype,  int iotype);
int fits_iter_set_file(iteratorCol *col, fitsfile *fptr);
int fits_iter_set_colname(iteratorCol *col, char *colname);
int fits_iter_set_colnum(iteratorCol *col, int colnum);
int fits_iter_set_datatype(iteratorCol *col, int datatype);
int fits_iter_set_iotype(iteratorCol *col, int iotype);

fitsfile * fits_iter_get_file(iteratorCol *col);
char * fits_iter_get_colname(iteratorCol *col);
int fits_iter_get_colnum(iteratorCol *col);
int fits_iter_get_datatype(iteratorCol *col);
int fits_iter_get_iotype(iteratorCol *col);
void * fits_iter_get_array(iteratorCol *col);
long fits_iter_get_tlmin(iteratorCol *col);
long fits_iter_get_tlmax(iteratorCol *col);
long fits_iter_get_repeat(iteratorCol *col);
char * fits_iter_get_tunit(iteratorCol *col);
char * fits_iter_get_tdisp(iteratorCol *col);

int ffiter(int ncols,  iteratorCol *data, long offset, long nPerLoop,
           int (*workFn)( long totaln, long offset, long firstn,
             long nvalues, int narrays, iteratorCol *data, void *userPointer),
           void *userPointer, int *status);

/*--------------------- write column elements -------------*/
int ffpcl(fitsfile *fptr, int datatype, int colnum, LONGLONG firstrow,
          LONGLONG firstelem, LONGLONG nelem, void *array, int *status);
int ffpcls(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, char **array, int *status);
int ffpcll(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, char *array, int *status);
int ffpclb(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, unsigned char *array, int *status);
int ffpclsb(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, signed char *array, int *status);
int ffpclui(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, unsigned short *array, int *status);
int ffpcli(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, short *array, int *status);
int ffpcluj(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, unsigned long *array, int *status);
int ffpclj(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, long *array, int *status);
int ffpcluk(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, unsigned int *array, int *status);
int ffpclk(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, int *array, int *status);
int ffpcle(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, float *array, int *status);
int ffpcld(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, double *array, int *status);
int ffpclc(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, float *array, int *status);
int ffpclm(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, double *array, int *status);
int ffpclu(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, int *status);
int ffprwu(fitsfile *fptr, LONGLONG firstrow, LONGLONG nrows, int *status);
int ffpcljj(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, LONGLONG *array, int *status);
int ffpclx(fitsfile *fptr, int colnum, LONGLONG frow, long fbit, long nbit,
            char *larray, int *status);

int ffpcn(fitsfile *fptr, int datatype, int colnum, LONGLONG firstrow, LONGLONG firstelem,
          LONGLONG nelem, void *array, void *nulval, int *status);
int ffpcns( fitsfile *fptr, int  colnum, LONGLONG  firstrow, LONGLONG  firstelem,
            LONGLONG  nelem, char **array, char  *nulvalue, int  *status);
int ffpcnl( fitsfile *fptr, int  colnum, LONGLONG  firstrow, LONGLONG  firstelem,
            LONGLONG  nelem, char *array, char  nulvalue,  int  *status);
int ffpcnb(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, unsigned char *array, unsigned char nulvalue,
           int *status);
int ffpcnsb(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, signed char *array, signed char nulvalue,
           int *status);
int ffpcnui(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, unsigned short *array, unsigned short nulvalue,
           int *status);
int ffpcni(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, short *array, short nulvalue, int *status);
int ffpcnuj(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, unsigned long *array, unsigned long nulvalue,
           int *status);
int ffpcnj(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, long *array, long nulvalue, int *status);
int ffpcnuk(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, unsigned int *array, unsigned int nulvalue,
           int *status);
int ffpcnk(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, int *array, int nulvalue, int *status);
int ffpcne(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, float *array, float nulvalue, int *status);
int ffpcnd(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, double *array, double nulvalue, int *status);
int ffpcnjj(fitsfile *fptr, int colnum, LONGLONG firstrow, LONGLONG firstelem,
           LONGLONG nelem, LONGLONG *array, LONGLONG nulvalue, int *status);
int ffptbb(fitsfile *fptr, LONGLONG firstrow, LONGLONG firstchar, LONGLONG nchars,
           unsigned char *values, int *status);
 
int ffirow(fitsfile *fptr, LONGLONG firstrow, LONGLONG nrows, int *status);
int ffdrow(fitsfile *fptr, LONGLONG firstrow, LONGLONG nrows, int *status);
int ffdrrg(fitsfile *fptr, char *ranges, int *status);
int ffdrws(fitsfile *fptr, long *rownum,  long nrows, int *status);
int ffdrwsll(fitsfile *fptr, LONGLONG *rownum,  LONGLONG nrows, int *status);
int fficol(fitsfile *fptr, int numcol, char *ttype, char *tform, int *status);
int fficls(fitsfile *fptr, int firstcol, int ncols, char **ttype,
           char **tform, int *status);
int ffmvec(fitsfile *fptr, int colnum, LONGLONG newveclen, int *status);
int ffdcol(fitsfile *fptr, int numcol, int *status);
int ffcpcl(fitsfile *infptr, fitsfile *outfptr, int incol, int outcol, 
           int create_col, int *status);
int ffcprw(fitsfile *infptr, fitsfile *outfptr, LONGLONG firstrow, 
           LONGLONG nrows, int *status);

/*--------------------- WCS Utilities ------------------*/
int ffgics(fitsfile *fptr, double *xrval, double *yrval, double *xrpix,
           double *yrpix, double *xinc, double *yinc, double *rot,
           char *type, int *status);
int ffgicsa(fitsfile *fptr, char version, double *xrval, double *yrval, double *xrpix,
           double *yrpix, double *xinc, double *yinc, double *rot,
           char *type, int *status);
int ffgtcs(fitsfile *fptr, int xcol, int ycol, double *xrval,
           double *yrval, double *xrpix, double *yrpix, double *xinc,
           double *yinc, double *rot, char *type, int *status);
int ffwldp(double xpix, double ypix, double xref, double yref,
           double xrefpix, double yrefpix, double xinc, double yinc,
           double rot, char *type, double *xpos, double *ypos, int *status);
int ffxypx(double xpos, double ypos, double xref, double yref, 
           double xrefpix, double yrefpix, double xinc, double yinc,
           double rot, char *type, double *xpix, double *ypix, int *status);

/*   WCS support routines (provide interface to Doug Mink's WCS library */
int ffgiwcs(fitsfile *fptr,  char **header, int *status); 
int ffgtwcs(fitsfile *fptr, int xcol, int ycol, char **header, int *status);

/*--------------------- lexical parsing routines ------------------*/
int fftexp( fitsfile *fptr, char *expr, int maxdim,
	    int *datatype, long *nelem, int *naxis,
	    long *naxes, int *status );

int fffrow( fitsfile *infptr, char *expr,
	    long firstrow, long nrows,
            long *n_good_rows, char *row_status, int *status);

int ffffrw( fitsfile *fptr, char *expr, long *rownum, int *status);

int fffrwc( fitsfile *fptr, char *expr, char *timeCol,    
            char *parCol, char *valCol, long ntimes,      
            double *times, char *time_status, int  *status );

int ffsrow( fitsfile *infptr, fitsfile *outfptr, char *expr, 
            int *status);

int ffcrow( fitsfile *fptr, int datatype, char *expr,
	    long firstrow, long nelements, void *nulval,
	    void *array, int *anynul, int *status );

int ffcalc_rng( fitsfile *infptr, char *expr, fitsfile *outfptr,
               char *parName, char *parInfo, int nRngs,
                 long *start, long *end, int *status );

int ffcalc( fitsfile *infptr, char *expr, fitsfile *outfptr,
            char *parName, char *parInfo, int *status );

  /* ffhist is not really intended as a user-callable routine */
  /* but it may be useful for some specialized applications   */
  /* ffhist2 is a newer version which is strongly recommended instead of ffhist */

int ffhist(fitsfile **fptr, char *outfile, int imagetype, int naxis,
           char colname[4][FLEN_VALUE],
           double *minin, double *maxin, double *binsizein,
           char minname[4][FLEN_VALUE], char maxname[4][FLEN_VALUE],
           char binname[4][FLEN_VALUE], 
           double weightin, char wtcol[FLEN_VALUE],
           int recip, char *rowselect, int *status);
int ffhist2(fitsfile **fptr, char *outfile, int imagetype, int naxis,
           char colname[4][FLEN_VALUE],
           double *minin, double *maxin, double *binsizein,
           char minname[4][FLEN_VALUE], char maxname[4][FLEN_VALUE],
           char binname[4][FLEN_VALUE], 
           double weightin, char wtcol[FLEN_VALUE],
           int recip, char *rowselect, int *status);

int fits_select_image_section(fitsfile **fptr, char *outfile,
           char *imagesection, int *status);
int fits_copy_image_section(fitsfile *infptr, fitsfile *outfile,
           char *imagesection, int *status);

int fits_calc_binning(fitsfile *fptr, int naxis, char colname[4][FLEN_VALUE], 
    double *minin, double *maxin,  double *binsizein,
    char minname[4][FLEN_VALUE],  char maxname[4][FLEN_VALUE], 
    char binname[4][FLEN_VALUE],  int *colnum,  long *haxes,  float *amin, 
    float *amax, float *binsize,  int *status);

int fits_write_keys_histo(fitsfile *fptr,  fitsfile *histptr, 
      int naxis, int *colnum, int *status);  
int fits_rebin_wcs( fitsfile *fptr, int naxis, float *amin,  float *binsize, 
      int *status);      
int fits_make_hist(fitsfile *fptr, fitsfile *histptr, int bitpix,int naxis,
     long *naxes,  int *colnum,  float *amin,  float *amax, float *binsize,
     float weight, int wtcolnum, int recip, char *selectrow, int *status);

typedef struct
{
	/* input(s) */
	int count;
	char ** path;
	char ** tag;
	fitsfile ** ifptr;

	char * expression;

	/* output control */
	int bitpix;
	long blank;
	fitsfile * ofptr;
	char keyword[FLEN_KEYWORD];
	char comment[FLEN_COMMENT];
} PixelFilter;


int fits_pixel_filter (PixelFilter * filter, int * status);


/*--------------------- grouping routines ------------------*/

int ffgtcr(fitsfile *fptr, char *grpname, int grouptype, int *status);
int ffgtis(fitsfile *fptr, char *grpname, int grouptype, int *status);
int ffgtch(fitsfile *gfptr, int grouptype, int *status);
int ffgtrm(fitsfile *gfptr, int rmopt, int *status);
int ffgtcp(fitsfile *infptr, fitsfile *outfptr, int cpopt, int *status);
int ffgtmg(fitsfile *infptr, fitsfile *outfptr, int mgopt, int *status);
int ffgtcm(fitsfile *gfptr, int cmopt, int *status);
int ffgtvf(fitsfile *gfptr, long *firstfailed, int *status);
int ffgtop(fitsfile *mfptr,int group,fitsfile **gfptr,int *status);
int ffgtam(fitsfile *gfptr, fitsfile *mfptr, int hdupos, int *status);
int ffgtnm(fitsfile *gfptr, long *nmembers, int *status);
int ffgmng(fitsfile *mfptr, long *nmembers, int *status);
int ffgmop(fitsfile *gfptr, long member, fitsfile **mfptr, int *status);
int ffgmcp(fitsfile *gfptr, fitsfile *mfptr, long member, int cpopt, 
	   int *status);
int ffgmtf(fitsfile *infptr, fitsfile *outfptr,	long member, int tfopt,	       
	   int *status);
int ffgmrm(fitsfile *fptr, long member, int rmopt, int *status);

/*--------------------- group template parser routines ------------------*/

int	fits_execute_template(fitsfile *ff, char *ngp_template, int *status);

int fits_img_stats_short(short *array,long nx, long ny, int nullcheck,   
    short nullvalue,long *ngoodpix, short *minvalue, short *maxvalue, double *mean,  
    double *sigma, double *noise1, double *noise2, double *noise3, double *noise5, int *status);
int fits_img_stats_int(int *array,long nx, long ny, int nullcheck,   
    int nullvalue,long *ngoodpix, int *minvalue, int *maxvalue, double *mean,  
    double *sigma, double *noise1, double *noise2, double *noise3, double *noise5, int *status);
int fits_img_stats_float(float *array, long nx, long ny, int nullcheck,   
    float nullvalue,long *ngoodpix, float *minvalue, float *maxvalue, double *mean,  
    double *sigma, double *noise1, double *noise2, double *noise3, double *noise5, int *status);

/*--------------------- image compression routines ------------------*/

int fits_set_compression_type(fitsfile *fptr, int ctype, int *status);
int fits_set_tile_dim(fitsfile *fptr, int ndim, long *dims, int *status);
int fits_set_noise_bits(fitsfile *fptr, int noisebits, int *status);
int fits_set_quantize_level(fitsfile *fptr, float qlevel, int *status);
int fits_set_hcomp_scale(fitsfile *fptr, float scale, int *status);
int fits_set_hcomp_smooth(fitsfile *fptr, int smooth, int *status);
int fits_set_quantize_dither(fitsfile *fptr, int dither, int *status);
int fits_set_dither_offset(fitsfile *fptr, int offset, int *status);
int fits_set_lossy_int(fitsfile *fptr, int lossy_int, int *status);

int fits_get_compression_type(fitsfile *fptr, int *ctype, int *status);
int fits_get_tile_dim(fitsfile *fptr, int ndim, long *dims, int *status);
int fits_get_quantize_level(fitsfile *fptr, float *qlevel, int *status);
int fits_get_noise_bits(fitsfile *fptr, int *noisebits, int *status);
int fits_get_hcomp_scale(fitsfile *fptr, float *scale, int *status);
int fits_get_hcomp_smooth(fitsfile *fptr, int *smooth, int *status);
int fits_get_dither_offset(fitsfile *fptr, int *offset, int *status);

int fits_img_compress(fitsfile *infptr, fitsfile *outfptr, int *status);
int fits_compress_img(fitsfile *infptr, fitsfile *outfptr, int compress_type,
         long *tilesize, int parm1, int parm2, int *status);
int fits_is_compressed_image(fitsfile *fptr, int *status);
int fits_is_reentrant(void);
int fits_decompress_img (fitsfile *infptr, fitsfile *outfptr, int *status);
int fits_img_decompress_header(fitsfile *infptr, fitsfile *outfptr, int *status);
int fits_img_decompress (fitsfile *infptr, fitsfile *outfptr, int *status);

/* H-compress routines */
int fits_hcompress(int *a, int nx, int ny, int scale, char *output, 
    long *nbytes, int *status);
int fits_hcompress64(LONGLONG *a, int nx, int ny, int scale, char *output, 
    long *nbytes, int *status);
int fits_hdecompress(unsigned char *input, int smooth, int *a, int *nx, 
       int *ny, int *scale, int *status);
int fits_hdecompress64(unsigned char *input, int smooth, LONGLONG *a, int *nx, 
       int *ny, int *scale, int *status);

int fits_transpose_table(fitsfile *infptr, fitsfile *outfptr, int *status);
int fits_compress_table_fast(fitsfile *infptr, fitsfile *outfptr, int *status);
int fits_compress_table_best(fitsfile *infptr, fitsfile *outfptr, int *status);
int fits_compress_table_rice(fitsfile *infptr, fitsfile *outfptr, int *status);
int fits_uncompress_table(fitsfile *infptr, fitsfile *outfptr, int *status);
int fits_gzip_datablocks(fitsfile *fptr, size_t *size, int *status);

/*  The following exclusion if __CINT__ is defined is needed for ROOT */
#ifndef __CINT__
#ifdef __cplusplus
}
#endif
#endif

#endif

