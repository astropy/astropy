#!/bin/sh

set -euv

# This script should be run every time expat is updated.

# So, the standard update would be to execute, from this directory:
#
# rm -rf expat
# tar xvf <PATH_TO_TAR>   # (e.g., expat-2.2.6.tar.bz2)
# cd expat
# ./buildconf.sh
# ./configure --without-getrandom --without-sys-getrandom
# cp expat_config.h ../../astropy/utils/xml/src/
# cd ..
# ./trim_expat.sh
# cd ..
# git add . -u
#
# And if the trim script missed anything, try to remove the extra files
# with a "git clean -xdf" command, do a local build, and run a full test suite
# with "pytest --remote-data" to make sure we didn't accidentally deleted an
# important file.

rm -rf expat/{conftools,doc,examples,m4,test,win32,CMake*,configure*,Makefile*,lib/*vcxproj*,tests,*.cmake,aclocal.m4,run.sh.in,test-driver-wrapper.sh,expat.*,xmlwf}

cat <<EOF >expat/README.txt
Note: astropy only requires the expat library, and hence in this bundled version,
we removed all other files except the required license and changelog.
EOF
