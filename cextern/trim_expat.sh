#!/bin/sh

set -euv

# This script should be run every time expat is updated.

# So, the standard update would be to execute, from this directory,
# rm -rf expat
# tar xvf <PATH_TO_TAR>   # (e.g., expat-2.2.6.tar.bz2)
# ./trim_expat.sh

rm -rf expat/{conftools,doc,examples,m4,test,win32,CMake*,configure*,Makefile*,lib/*vcxproj*,tests,*.cmake,aclocal.m4,run.sh.in,test-driver-wrapper.sh,expat.*,xmlwf}

cat <<EOF >expat/README.txt
Note: astropy only requires the expat library, and hence in this bundled version,
we removed all other files except the required license and changelog.
EOF
