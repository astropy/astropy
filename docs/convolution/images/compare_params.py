import numpy as np
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization.wcsaxes.frame import EllipticalFrame
from astropy.convolution import Gaussian2DKernel
from scipy.signal import convolve as scipy_convolve
from astropy.convolution import convolve,convolve_fft
import matplotlib.pyplot as plt

filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')
hdu = fits.open(filename)[0]
img = hdu.data[50:90,60:100] * 1e5
img[img > 2e1] = np.nan
img_zerod = img.copy()
img_zerod[np.isnan(img)] = 0

kernel = Gaussian2DKernel(stddev=1)
scipy_conv = scipy_convolve(img, kernel, mode='same', method='direct')
scipy_conv_zerod = scipy_convolve(img_zerod, kernel, mode='same', method='direct')
astropy_conv = convolve(img, kernel)

plt.figure(1, figsize=(12,12)).clf()
ax1 = plt.subplot(2,2,1)
im = ax1.imshow(img, vmin=-2., vmax=2.e1, origin='lower',
                interpolation='nearest', cmap='viridis')
ax1.set_title("Original")
ax1.set_xticklabels([])
ax1.set_yticklabels([])

ax2 = plt.subplot(2,2,2)
im = ax2.imshow(scipy_conv, vmin=-2., vmax=2.e1, origin='lower',
                interpolation='nearest', cmap='viridis')
ax2.set_title("Scipy")
ax2.set_xticklabels([])
ax2.set_yticklabels([])

ax3 = plt.subplot(2,2,3)
im = ax3.imshow(scipy_conv_zerod, vmin=-2., vmax=2.e1, origin='lower',
                interpolation='nearest', cmap='viridis')
ax3.set_title("Scipy nan->zero")
ax3.set_xticklabels([])
ax3.set_yticklabels([])


ax4 = plt.subplot(2,2,4)
im = ax4.imshow(astropy_conv, vmin=-2., vmax=2.e1, origin='lower',
                    interpolation='nearest', cmap='viridis')
ax4.set_title("Default astropy")
ax4.set_xticklabels([])
ax4.set_yticklabels([])

plt.figure(2).clf()
plt.plot(img[:,25], label='input', drawstyle='steps-mid', linewidth=2, alpha=0.5)
plt.plot(scipy_conv[:,25], label='scipy', drawstyle='steps-mid',
         linewidth=2, alpha=0.5, marker='s')
plt.plot(scipy_conv_zerod[:,25], label='scipy nan->zero', drawstyle='steps-mid',
         linewidth=2, alpha=0.5, marker='s')
plt.plot(astropy_conv[:,25], label='astropy', drawstyle='steps-mid', linewidth=2, alpha=0.5)
plt.legend(loc='best')
plt.show()
