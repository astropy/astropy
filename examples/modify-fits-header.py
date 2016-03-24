# -*- coding: utf-8 -*-
"""
==================
Edit a FITS header
==================

This example describes how to read in, edit a FITS header, and then write it back out to disk.
For this example we're going to change the OBJECT keyword.

This example uses `~astropy.io.fits`, which was formerly released separately as pyfits. If you
have used  pyfits to manipulate FITS files then you may already be familiar with the features
and syntax of the package. We start by importing the subpackage into our local namespace, and
allows us to access the functions and classes as ``fits.name_of_function()``. For example, to access
the `~astropy.io.fits.getdata()` function, we don't have to do `~astropy.io.fits.getdata()` and can instead simple use
`~astropy.io.fits.getdata()`. You may run across old documentation or tutorials that use the name pyfits.
Such examples will begin with import pyfits and then the command `~astropy.io.fits.getdata() (for example)
would be written as pyfits.getdata().
"""

# Code source: Adrian Price-Whelan
# License: BSD

from astropy.io import fits

##############################################################################
# The following line is needed to download the example FITS files used here.

from astropy.utils.data import download_file
fits_file = download_file('http://data.astropy.org/tutorials/FITS-Header/input_file.fits',
                          cache=True)

##############################################################################
# `astropy.io.fits` provides a lot of flexibility for reading FITS files and
# headers, but most of the time the convenience functions are the easiest way
# to access the data. `~astropy.io.fits.getdata()` reads just the data from a FITS file,
# but with the ``header=True`` keyword argument will also read the header.

data, header = fits.getdata(fits_file, header=True)

##############################################################################
# There is also a dedicated function for reading just the header:

hdu_number = 0
hdr = fits.getheader(fits_file, hdu_number)
print(hdr['OBJECT'])

##############################################################################
# but `~astropy.io.fits.getdata()` can get both the data and the header, so it is a
# useful command to remember. Since the primary HDU of a FITS file must
# contain image data, the data is now stored in a numpy array. The header is
# stored in an object that acts like a standard Python dictionary.

print(type(data))
print(header['NAXIS'])

##############################################################################
# Now let's change the header to give it the correct object:

header['OBJECT'] = 'M31'

##############################################################################
# Finally, write out the FITS file. Again, the convenience function
# for this is the most useful command to remember:

fits.writeto('output_file.fits', data, header, clobber=True)

##############################################################################
# Two common more complicated cases are worth mentioning (but if your needs
# are much more complex, you should consult the full documentation in `astropy.io.fits`).
#
# The first complication is that a FITS 
# might have multiple HDU's (extensions), in which case the
# extension can be specified like this:

data,header = fits.getdata(fits_file, ext=1, header=True)

##############################################################################
# This will get the data and header associated with the index=1 extension
# in the FITS file. Without specifying a number, `~astropy.io.fits.getdata()` will get the
# 0th extension (equivalent to saying ``ext=0``).
#
# Use the clobber keyword argument to overwrite an existing FITS file.
# By default, `~astropy.io.fits.writeto()` will not allow this. 

fits.writeto('output_file.fits', data, header, clobber=True)
