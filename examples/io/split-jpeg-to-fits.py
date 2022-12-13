"""
=====================================================
Convert a 3-color image (JPG) to separate FITS images
=====================================================

This example opens an RGB JPEG image and writes out each channel as a separate
FITS (image) file.

This example uses `pillow <https://python-pillow.org>`_ to read the image,
`matplotlib.pyplot` to display the image, and `astropy.io.fits` to save FITS files.


*By: Erik Bray, Adrian Price-Whelan*

*License: BSD*


"""

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

from astropy.io import fits
from astropy.visualization import astropy_mpl_style

##############################################################################
# Set up matplotlib and use a nicer set of plot parameters

plt.style.use(astropy_mpl_style)

##############################################################################
# Load and display the original 3-color jpeg image:

image = Image.open('Hs-2009-14-a-web.jpg')
xsize, ysize = image.size
print(f"Image size: {ysize} x {xsize}")
print(f"Image bands: {image.getbands()}")
ax = plt.imshow(image)

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

r_data = r_data.reshape(ysize, xsize) # data is now a matrix (ysize, xsize)
g_data = g_data.reshape(ysize, xsize)
b_data = b_data.reshape(ysize, xsize)
print(r_data.shape)

##############################################################################
# Write out the channels as separate FITS images.
# Add and visualize header info

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

from pprint import pprint

pprint(red.header)

##############################################################################
# Delete the files created

import os

os.remove('red.fits')
os.remove('green.fits')
os.remove('blue.fits')
