#!/bin/sh

# This script should be run every time wcslib is updated.

# It patches it for the use of astropy.wcs and removes extra files
# that aren't needed.

cat patches/*.patch | patch -p0 -b -d src/wcslib
rm -rf src/wcslib/C/test
rm -rf src/wcslib/doxygen
rm -rf src/wcslib/Fortran
rm -rf src/wcslib/html
rm -rf src/wcslib/pgsbox
rm -rf src/wcslib/utils
rm src/wcslib/*.pdf