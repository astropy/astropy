#ifndef _COMPRESSIONMODULE_H
#define _COMPRESSIONMODULE_H


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
    #if CFITSIO_MINOR >= 35
        #define CFITSIO_SUPPORTS_Q_FORMAT_COMPRESSION
        #define CFITSIO_SUPPORTS_SUBTRACTIVE_DITHER_2
    #else
        /* This constant isn't defined in older versions and has a different */
        /* value anyways. */
        #define NO_DITHER 0
    #endif
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


#define CFITSIO_LOSSLESS_COMP_SUPPORTED_VERS 3.22


/* These defaults mirror the defaults in io.fits.hdu.compressed */
#define DEFAULT_COMPRESSION_TYPE "RICE_1"
#define DEFAULT_QUANTIZE_LEVEL 16.0
#define DEFAULT_HCOMP_SCALE 0
#define DEFAULT_HCOMP_SMOOTH 0
#define DEFAULT_BLOCK_SIZE 32
#define DEFAULT_BYTE_PIX 4

/* This constant is defined by cfitsio in imcompress.c */
#define NO_QUANTIZE 9999

#endif
