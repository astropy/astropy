# -*- coding: utf-8 -*-
"""
==================
Edit a FITS header
==================

This example describes how to edit a value in a FITS header
using `astropy.io.fits`.

-------------------

*By: Adrian Price-Whelan*

*License: BSD*

-------------------

"""

from astropy.io import fits

##############################################################################
# Download a FITS file:

from astropy.utils.data import get_pkg_data_filename

fits_file = get_pkg_data_filename('tutorials/FITS-Header/input_file.fits')

##############################################################################
# Look at contents of the FITS file

fits.info(fits_file)

##############################################################################
# Look at the headers of the two extensions:

print("Before modifications:")
print()
print("Extension 0:")
print(repr(fits.getheader(fits_file, 0)))
print()
print("Extension 1:")
print(repr(fits.getheader(fits_file, 1)))

##############################################################################
# `astropy.io.fits` provides an object-oriented interface for reading and
# interacting with FITS files, but for small operations (like this example) it
# is often easier to use the
# `convenience functions <http://docs.astropy.org/en/latest/io/fits/index.html#convenience-functions>`_.
#
# To edit a single header value in the header for extension 0, use the
# `~astropy.io.fits.setval()` function. For example, set the OBJECT keyword
# to 'M31':

fits.setval(fits_file, 'OBJECT', value='M31')

##############################################################################
# With no extra arguments, this will modify the header for extension 0, but
# this can be changed using the ``ext`` keyword argument. For example, we can
# specify extension 1 instead:

fits.setval(fits_file, 'OBJECT', value='M31', ext=1)

##############################################################################
# This can also be used to create a new keyword-value pair ("card" in FITS
# lingo):

fits.setval(fits_file, 'ANEWKEY', value='some value')

##############################################################################
# Again, this is useful for one-off modifications, but can be inefficient
# for operations like editing multiple headers in the same file
# because `~astropy.io.fits.setval()` loads the whole file each time it
# is called. To make several modifications, it's better to load the file once:

with fits.open(fits_file, 'update') as f:
    for hdu in f:
        hdu.header['OBJECT'] = 'CAT'

print("After modifications:")
print()
print("Extension 0:")
print(repr(fits.getheader(fits_file, 0)))
print()
print("Extension 1:")
print(repr(fits.getheader(fits_file, 1)))
