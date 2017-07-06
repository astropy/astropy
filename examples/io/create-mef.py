# -*- coding: utf-8 -*-
"""
=====================================================
Create a multi-extension FITS (MEF) file from scratch
=====================================================

This example demonstrates how to create a multi-extension FITS (MEF)
file from scratch using `astropy.io.fits`.

-------------------

*By: Erik Bray*

*License: BSD*

-------------------

"""

import os

##############################################################################
# HDUList objects are used to hold all the HDUs in a FITS file. This
# ``HDUList`` class is a subclass of Python's builtin `list`. and can be
# created from scratch. For example, to create a FITS file with
# three extensions:

from astropy.io import fits
new_hdul = fits.HDUList()
new_hdul.append(fits.ImageHDU())
new_hdul.append(fits.ImageHDU())

##############################################################################
# Write out the new file to disk:

new_hdul.writeto('test.fits')

##############################################################################
# Alternatively, the HDU instances can be created first (or read from an
# existing FITS file).
#
# Create a multi-extension FITS file with two empty IMAGE extensions (a
# default PRIMARY HDU is prepended automatically if one is not specified;
# we use ``overwrite=True`` to overwrite the file if it already exists):

hdu1 = fits.PrimaryHDU()
hdu2 = fits.ImageHDU()
new_hdul = fits.HDUList([hdu1, hdu2])
new_hdul.writeto('test.fits', overwrite=True)

##############################################################################
# Finally, we'll remove the file we created:

os.remove('test.fits')
