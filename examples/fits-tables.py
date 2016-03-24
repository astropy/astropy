# -*- coding: utf-8 -*-
"""
========================
Fits Table example
========================

Demonstrates astropy.utils.data to download the file, astropy.io.fits to open
and view the file, matplotlib for making

"""

# Code source:  Lia R. Corrales

import numpy as np

##############################################################################
# Set up matplotlib and use a nicer set of plot parameters
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

##############################################################################
# The following line is needed to download the example FITS files used here.

from astropy.utils.data import download_file
from astropy.io import fits

image_file = download_file('http://data.astropy.org/tutorials/FITS-images/HorseHead.fits',
                           cache=True)

##############################################################################
# Open the FITS file and find out what it contains.

hdu_list = fits.open(image_file)
hdu_list.info()

##############################################################################
# Generally the image information is located in the PRIMARY block.
# The blocks are numbered and can be accessed by indexing hdu_list.

image_data = hdu_list[0].data

##############################################################################
# The data is now stored as a 2-D numpy array.
# Look at the shape of the array to see the dimensions of the image

print(type(image_data))
print(image_data.shape)

##############################################################################
# At this point, the FITS file can be closed since everything needed is now
# stored in variables.

hdu_list.close()

##############################################################################
# To simply load the image and not examine the header, use fits.getdata to
# bypass the previous steps.

image_data = fits.getdata(image_file)
print(type(image_data))
print(image_data.shape)

##############################################################################
# Viewing the image data and getting basic statistics

plt.imshow(image_data, cmap='gray')
plt.colorbar()

##############################################################################
# Get some basic statistics about the image

print('Min:', np.min(image_data))
print('Max:', np.max(image_data))
print('Mean:', np.mean(image_data))
print('Stdev:', np.std(image_data))
