#!/bin/sh

# This script should be run every time wcslib is updated.

# This removes extra large files from wcslib that aren't needed.

rm -rf wcslib/C/flexed/RCS
rm -rf wcslib/C/RCS
rm -rf wcslib/C/test
rm -rf wcslib/doxygen
rm -rf wcslib/Fortran
rm -rf wcslib/html
rm -rf wcslib/pgsbox
rm -rf wcslib/RCS/
rm -rf wcslib/utils
rm wcslib/*.pdf
rm wcslib/.clang-tidy
rm wcslib/clang-tidy

# Autotools/build infrastructure -- astropy builds wcslib via
# astropy/wcs/setup_package.py and writes its own wcsconfig.h, so none of
# the configure/make machinery is used.
rm -rf wcslib/config
rm wcslib/configure
rm wcslib/configure.ac
rm wcslib/GNUmakefile
rm wcslib/C/GNUmakefile
rm wcslib/makedefs.in
rm wcslib/flavours
rm wcslib/wcsconfig.h.in
rm wcslib/wcsconfig_f77.h.in
rm wcslib/wcsconfig_tests.h.in
rm wcslib/wcsconfig_utils.h.in
rm wcslib/wcslib.pc.in
rm wcslib/INSTALL
