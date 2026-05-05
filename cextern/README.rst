External Packages/Libraries
===========================

This directory contains C libraries included with Astropy. Note that only C
libraries without python-specific code  should be included in this directory.
Cython or C code intended for use with Astropy or wrapper code should be in
the Astropy source tree.

CFITSIO
-------
The cfitsio directory only contains the small subset of files from CFITSIO which are
required for the astropy.io.fits.hdu.compressed._tiled_compression package. All files are
copied verbatim from CFITSIO and can easily be updated if needed.
