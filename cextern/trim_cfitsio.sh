#!/bin/sh

set -euv

# This script should be run every time cfitsio is updated.
# This moves all the code needed for the actual library to lib
# and deletes everything else (except License.txt and doc/changes.txt)

# So, the standard update would be to execute, from this directory,
# rm -rf cfitsio
# tar xvf <PATH_TO_TAR>   # (e.g., cfitsio-4.2.0.tar.gz)
# mv cfitsio-?.?.? cfitsio  # (e.g., mv cfitsio-4.2.0 cfitsio)
# ./trim_cfitsio.sh

if [ ! -d cfitsio/lib ]; then
    mkdir cfitsio/lib
fi

mv cfitsio/fits_hcompress.c cfitsio/lib/
mv cfitsio/fits_hdecompress.c cfitsio/lib/
mv cfitsio/pliocomp.c cfitsio/lib/
mv cfitsio/quantize.c cfitsio/lib/
mv cfitsio/ricecomp.c cfitsio/lib/

rm -f cfitsio/README
rm -f cfitsio/configure
rm -f cfitsio/install-sh
rm -f cfitsio/docs/*.tex
rm -f cfitsio/docs/*.ps
rm -f cfitsio/docs/*.pdf
rm -f cfitsio/docs/*.doc
rm -f cfitsio/docs/*.toc
rm -rf cfitsio/[^L]*.*
rm -rf cfitsio/utilities

# We only use a very small subset of fitsio2.h, so here we generate that
# file. If there are compilation issues after updating, it may be that
# the definitions below need tweaking or that some definitions need to be
# removed or added.
cat <<EOF > cfitsio/lib/fitsio2.h
#ifndef LONGLONG_TYPE
typedef long long LONGLONG;
typedef unsigned long long ULONGLONG;
#define LONGLONG_TYPE
#endif

# define DATA_COMPRESSION_ERR 413
# define DATA_DECOMPRESSION_ERR 414

void ffpmsg(const char *err_message);

#define FFLOCK
#define FFUNLOCK

#define N_RANDOM 10000

#define MEMORY_ALLOCATION 113

#define NO_DITHER -1
#define SUBTRACTIVE_DITHER_1 1
#define SUBTRACTIVE_DITHER_2 2

int fits_init_randoms(void);
EOF


cat <<EOF >cfitsio/README.rst
This directory only contains the small subset of files from CFITSIO which are
required for the astropy.io.fits.hdu.compressed._tiled_compression package. All files are
copied verbatim from CFITSIO and can easily be updated if needed.
EOF
