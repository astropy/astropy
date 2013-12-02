#!/bin/sh

# This script should be run every time wcslib is updated.

# This removes extra large files from wcslib that aren't needed.

rm -rf wcslib/C/test
rm -rf wcslib/doxygen
rm -rf wcslib/Fortran
rm -rf wcslib/html
rm -rf wcslib/pgsbox
rm -rf wcslib/utils
rm wcslib/*.pdf

patch -p1 < patches/cunit-fix-message.diff
patch -p1 < patches/unit-translate-information.diff
