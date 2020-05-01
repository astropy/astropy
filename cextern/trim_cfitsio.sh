#!/bin/sh

set -euv

# This script should be run every time cfitsio is updated.
# This moves all the code needed for the actual library to lib
# and deletes everything else (except License.txt and doc/changes.txt)

# So, the standard update would be to execute, from this directory,
# rm -rf cfitsio
# tar xvf <PATH_TO_TAR>   # (e.g., cfitsio3410.tar.gz)
# ./trim_cfitsio.sh


# This just gets CORE_SOURCES from Makefile.in, excluding anything beyond zlib
lib_files=`make -f cfitsio/Makefile.in cfitsioLibSrcs | sed 's/zlib\/.*//'`
# The include files cannot be directly inferred from Makefile.in
inc_files='fitsio.h fitsio2.h longnam.h drvrsmem.h eval_defs.h eval_tab.h region.h group.h simplerng.h grparser.h'

if [ ! -d cfitsio/lib ]; then
    mkdir cfitsio/lib
fi

for fil in $lib_files $inc_files; do
    if [ -f cfitsio/$fil ]; then
        mv cfitsio/$fil cfitsio/lib/
    fi
done

rm -f cfitsio/README
rm -f cfitsio/configure
rm -f cfitsio/install-sh
rm -f cfitsio/docs/*.tex
rm -f cfitsio/docs/*.ps
rm -f cfitsio/docs/*.pdf
rm -f cfitsio/docs/*.doc
rm -f cfitsio/docs/*.toc
rm -rf cfitsio/[^L]*.*

cat <<EOF >cfitsio/README.txt
Note: astropy only requires the CFITSIO library, and hence in this bundled version,
we removed all other files except the required license (License.txt) and changelog
(docs/changes.txt, which has the version number).
EOF
