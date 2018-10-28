# -*- coding: utf-8 -*-
"""
==========================================
Create a very large FITS file from scratch
==========================================

This example demonstrates how to create a large file (larger than will fit in
memory) from scratch using `astropy.io.fits`.

-------------------

*By: Erik Bray*

*License: BSD*

-------------------
"""

##############################################################################
#  Normally to create a single image FITS file one would do something like:

import os
import numpy as np
from astropy.io import fits
data = np.zeros((40000, 40000), dtype=np.float64)
hdu = fits.PrimaryHDU(data=data)

##############################################################################
# Then use the `astropy.io.fits.writeto()` method to write out the new
# file to disk

hdu.writeto('large.fits')

##############################################################################
# However, a 40000 x 40000 array of doubles is nearly twelve gigabytes! Most
# systems won't be able to create that in memory just to write out to disk. In
# order to create such a large file efficiently requires a little extra work,
# and a few assumptions.
#
# First, it is helpful to anticipate about how large (as in, how many keywords)
# the header will have in it. FITS headers must be written in 2880 byte
# blocks, large enough for 36 keywords per block (including the END keyword in
# the final block). Typical headers have somewhere between 1 and 4 blocks,
# though sometimes more.
#
# Since the first thing we write to a FITS file is the header, we want to write
# enough header blocks so that there is plenty of padding in which to add new
# keywords without having to resize the whole file. Say you want the header to
# use 4 blocks by default. Then, excluding the END card which Astropy will add
# automatically, create the header and pad it out to 36 * 4 cards.
#
# Create a stub array to initialize the HDU; its
# exact size is irrelevant, as long as it has the desired number of
# dimensions

data = np.zeros((100, 100), dtype=np.float64)
hdu = fits.PrimaryHDU(data=data)
header = hdu.header
while len(header) < (36 * 4 - 1):
    header.append()  # Adds a blank card to the end

##############################################################################
# Now adjust the NAXISn keywords to the desired size of the array, and write
# only the header out to a file. Using the ``hdu.writeto()`` method will cause
# astropy to "helpfully" reset the NAXISn keywords to match the size of the
# dummy array. That is because it works hard to ensure that only valid FITS
# files are written. Instead, we can write just the header to a file using the
# `astropy.io.fits.Header.tofile` method:

header['NAXIS1'] = 40000
header['NAXIS2'] = 40000
header.tofile('large.fits')

##############################################################################
# Finally, grow out the end of the file to match the length of the
# data (plus the length of the header). This can be done very efficiently on
# most systems by seeking past the end of the file and writing a single byte,
# like so:

with open('large.fits', 'rb+') as fobj:
    # Seek past the length of the header, plus the length of the
    # Data we want to write.
    # 8 is the number of bytes per value, i.e. abs(header['BITPIX'])/8
    # (this example is assuming a 64-bit float)
    # The -1 is to account for the final byte that we are about to
    # write:
    fobj.seek(len(header.tostring()) + (40000 * 40000 * 8) - 1)
    fobj.write(b'\0')

##############################################################################
# More generally, this can be written:

shape = tuple(header['NAXIS{0}'.format(ii)] for ii in range(1, header['NAXIS']+1))
with open('large.fits', 'rb+') as fobj:
    fobj.seek(len(header.tostring()) + (np.product(shape) * np.abs(header['BITPIX']//8)) - 1)
    fobj.write(b'\0')

##############################################################################
# On modern operating systems this will cause the file (past the header) to be
# filled with zeros out to the ~12GB needed to hold a 40000 x 40000 image. On
# filesystems that support sparse file creation (most Linux filesystems, but not
# the HFS+ filesystem used by most Macs) this is a very fast, efficient
# operation. On other systems your mileage may vary.
#
# This isn't the only way to build up a large file, but probably one of the
# safest. This method can also be used to create large multi-extension FITS
# files, with a little care.

##############################################################################
# Finally, we'll remove the file we created:

os.remove('large.fits')
