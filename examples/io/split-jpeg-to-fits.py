# -*- coding: utf-8 -*-
"""
=====================================================
Convert a 3-color image (JPG) to separate FITS images
=====================================================

This example opens an RGB JPEG image and writes out each channel as a separate
FITS (image) file.

This example uses `pillow <http://python-pillow.org>`_ to read the image,
`matplotlib.pyplot` to display the image, and `astropy.io.fits` to save FITS files.

-------------------

*By: Erik Bray, Adrian Price-Whelan*

*License: BSD*

-------------------

"""

import numpy as np
from PIL import Image
from astropy.io import fits

##############################################################################
# Set up matplotlib and use a nicer set of plot parameters

import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

##############################################################################
# Load and display the original 3-color jpeg image:

image = Image.open('Hs-2009-14-a-web.jpg')
xsize, ysize = image.size
print("Image size: {} x {}".format(xsize, ysize))
plt.imshow(image)

##############################################################################
# Split the three channels (RGB) and get the data as Numpy arrays. The arrays
# are flattened, so they are 1-dimensional:

r, g, b = image.split()
r_data = np.array(r.getdata()) # data is now an array of length ysize*xsize
g_data = np.array(g.getdata())
b_data = np.array(b.getdata())
print(r_data.shape)

##############################################################################
# Reshape the image arrays to be 2-dimensional:

r_data = r_data.reshape(ysize, xsize)
g_data = g_data.reshape(ysize, xsize)
b_data = b_data.reshape(ysize, xsize)

##############################################################################
# Write out the channels as separate FITS images

red = fits.PrimaryHDU(data=r_data)
red.header['LATOBS'] = "32:11:56" # add spurious header info
red.header['LONGOBS'] = "110:56"
red.writeto('red.fits')

green = fits.PrimaryHDU(data=g_data)
green.header['LATOBS'] = "32:11:56"
green.header['LONGOBS'] = "110:56"
green.writeto('green.fits')

blue = fits.PrimaryHDU(data=b_data)
blue.header['LATOBS'] = "32:11:56"
blue.header['LONGOBS'] = "110:56"
blue.writeto('blue.fits')

##############################################################################
# Delete the files created
import os
os.remove('red.fits')
os.remove('green.fits')
os.remove('blue.fits')
