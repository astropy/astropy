import numpy as np
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization.wcsaxes.frame import EllipticalFrame
from astropy.convolution import Gaussian2DKernel
from scipy.signal import convolve as scipy_convolve
from astropy.convolution import convolve,convolve_fft
import matplotlib.pyplot as plt

# Load the data from data.astropy.org
filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')
hdu = fits.open(filename)[0]

# Scale the file to have reasonable numbers
# (this is mostly so that colorbars don't have too many digits)
# Also, we crop it so you can see individual pixels
img = hdu.data[50:90,60:100] * 1e5

# This example is intended to demonstrate how astropy.convolve and
# scipy.convolve handle missing data, so we start by setting the brightest
# pixels to NaN to simulate a "saturated" data set
img[img > 2e1] = np.nan

# We also create a copy of the data and set those NaNs to zero.  We'll use this
# for the scipy convolution
img_zerod = img.copy()
img_zerod[np.isnan(img)] = 0

# We smooth with a Gaussian kernel with stddev=1
# It is a 9x9 array
kernel = Gaussian2DKernel(stddev=1)

# Convolution: scipy's direct convolution mode spreads out NaNs (see panel 2 below)
scipy_conv = scipy_convolve(img, kernel, mode='same', method='direct')

# scipy's direct convolution mode run on the 'zero'd' image will not have NaNs,
# but will have some very low value zones where the NaNs were
# (see panel 3 below)
scipy_conv_zerod = scipy_convolve(img_zerod, kernel, mode='same', method='direct')

# astropy's convolution replaces the NaN pixels with a kernel-weighted
# interpolation from their neighbors
astropy_conv = convolve(img, kernel)
