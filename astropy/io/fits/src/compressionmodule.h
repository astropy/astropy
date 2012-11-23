#ifndef _COMPRESSIONMODULE_H
#define _COMPRESSIONMODULE_H

/* Some defines for Python3 support--bytes objects should be used where */
/* strings were previously used                                         */
#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif

#ifdef IS_PY3K
#define PyString_FromString PyUnicode_FromString
#define PyInt_AsLong PyLong_AsLong
#endif


/* CFITSIO version-specific feature support */
#ifndef CFITSIO_MAJOR
    // Define a minimized version
    #define CFITSIO_MAJOR 0
    #ifdef _MSC_VER
        #pragma warning ( "CFITSIO_MAJOR not defined; your CFITSIO version may be too old; compile at your own risk" )
    #else
        #warning "CFITSIO_MAJOR not defined; your CFITSIO version may be too old; compile at your own risk"
    #endif
#endif

#ifndef CFITSIO_MINOR
    #define CFITSIO_MINOR 0
#endif


#if CFITSIO_MAJOR >= 3
    #if CFITSIO_MINOR >= 28
        #define CFITSIO_SUPPORTS_GZIPDATA
    #else
        #ifdef _MSC_VER
            #pragma warning ( "GZIP_COMPRESSED_DATA columns not supported" )
        #else
            #warning "GZIP_COMPRESSED_DATA columns not supported"
        #endif
    #endif
#endif


/* These defaults mirror the defaults in pyfits.hdu.compressed */
#define DEFAULT_COMPRESSION_TYPE "RICE_1"
#define DEFAULT_QUANTIZE_LEVEL 16.0
#define DEFAULT_HCOMP_SCALE 0
#define DEFAULT_HCOMP_SMOOTH 0
#define DEFAULT_BLOCK_SIZE 32
#define DEFAULT_BYTE_PIX 4

#endif
