This directory only contains the small subset of files from CFITSIO which are
required for the astropy.io.fits._tiled_compression package. All files are
copied verbatim from CFITSIO and can easily be updated if needed, and the
``imcompress.c`` file only contains the functions needed for unquantization.

Each time a new version of CFITSIO is released, new versions of the files
present here can be copied over manually, and ``imcompress.c`` should be
edited accordingly.
