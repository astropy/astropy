# -*- coding: utf-8 -*-
"""
==================
Edit a FITS header
==================

This tutorial describes how to read in, edit a FITS header, and then write it back out to disk.
For this example we're going to change the OBJECT keyword.

This tutorial uses `~astropy.io.fits`, which was formerly released separately as pyfits. If you
have used  pyfits to manipulate FITS files then you may already be familiar with the features
and syntax of the package. We start by importing the subpackage into our local namespace, and
allows us to access the functions and classes as ``fits.name_of_function()``. For example, to access
the getdata() function, we don't have to do `~astropy.io.fits.getdata()` and can instead simple use
`~astropy.io.fits.getdata()`. You may run across old documentation or tutorials that use the name pyfits.
Such examples will begin with import pyfits and then the command `astropy.io.fits.getdata() (for example)
would be written as pyfits.getdata().
"""

# Code source: Adrian Price-Whelan
# License: BSD

from astropy.io import fits

##############################################################################
# `astropy.io.fits` provides a lot of flexibility for reading FITS files and
# headers, but most of the time the convenience functions are the easiest way
# to access the data. `fits.getdata()` reads just the data from a FITS file,
# but with the ``header=True`` keyword argument will also read the header.

data, header = fits.getdata("input_file.fits", header=True)

##############################################################################
# There is also a dedicated function for reading just the header:

hdu_number = 0
fits.getheader('input_file.fits', hdu_number)

##############################################################################
# but `~fits.getdata()` can get both the data and the header, so it is a
# useful command to remember. Since the primary HDU of a FITS file must
# contain image data, the data is now stored in a numpy array. The header is
# stored in an object that acts like a standard Python dictionary.

print(type(data))
print(header["NAXIS"])

##############################################################################
# Now let's change the header to give it the correct object:

header['OBJECT'] = "M31"

##############################################################################
# Finally, we have to write out the FITS file. Again, the convenience function
# for this is the most useful command to remember:

fits.writeto('output_file.fits', data, header, clobber=True)

##############################################################################
# That's it; you're done!
#
# Two common more complicated cases are worth mentioning (but if your needs
# are much more complex, you should consult the full documentation).
#
# The first complication is that the FITS file you're examining and editing
# might have multiple HDU's (extensions), in which case you can specify the
# extension like this:

data,header = fits.getdata("input_file.fits", ext=1, header=True)

##############################################################################
# This will get you the data and header associated with the index=1 extension
# in the FITS file. Without specifying a number, `fits.getdata()` will get the
# 0th extension (equivalent to saying ``ext=0``).
#
# Another useful tip is if you want to overwrite an existing FITS file. By
# default, `~fits.writeto()` won't let you do this, so you need to explicitly
# give it permission using the clobber keyword argument:

fits.writeto('output_file.fits', data, header, clobber=True)

##############################################################################
# ## Exercise
#
# Read in the file you just wrote, and add three header keywords:
#
# * 'RA' for the Right Ascension of M31
# * 'DEC' for the Declination of M31
# * 'RADECSRC' with text indicating where you found the RA/Dec (web URL, textbook name, your photographic memory, etc.).

