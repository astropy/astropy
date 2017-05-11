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
scipy_conv = scipy_convolve(img, kernel, mode='same')
scipy_conv_zerod = scipy_convolve(img_zerod, kernel, mode='same')
astropy_conv = convolve(img, kernel)
astropy_conv_intnan = convolve(img, kernel, nan_treatment='interpolate', normalize_kernel=False)
astropy_conv_intnan_norm = convolve(img, kernel, nan_treatment='interpolate',
                                    normalize_kernel=True)
astropy_conv_norm = convolve(img, kernel, nan_treatment='interpolate'e,
                                    normalize_kernel=True)
astropy_conv_fft = convolve_fft(img, kernel)
astropy_conv_intnan_fft = convolve_fft(img, kernel, nan_treatment='interpolate', normalize_kernel=False)
astropy_conv_intnan_fft_norm = convolve_fft(img, kernel,
                   nan_treatment='interpolate', normalize_kernel=True)
astropy_conv_fft_norm = convolve_fft(img, kernel,
                   nan_treatment='fill', normalize_kernel=True)

plt.figure(1, figsize=(12,12)).clf()
ax1 = plt.subplot(3,3,1)
im = ax1.imshow(img, vmin=-2., vmax=2.e1, origin='lower',
                interpolation='nearest', cmap='viridis')
ax1.set_title("Original")
ax1.set_xticklabels([])
ax1.set_yticklabels([])

ax2 = plt.subplot(3,3,2)
im = ax2.imshow(scipy_conv, vmin=-2., vmax=2.e1, origin='lower',
                interpolation='nearest', cmap='viridis')
ax2.set_title("Scipy")
ax2.set_xticklabels([])
ax2.set_yticklabels([])

ax3 = plt.subplot(3,3,3)
im = ax3.imshow(scipy_conv_zerod, vmin=-2., vmax=2.e1, origin='lower',
                interpolation='nearest', cmap='viridis')
ax3.set_title("Scipy nan->zero")
ax3.set_xticklabels([])
ax3.set_yticklabels([])

ax4 = plt.subplot(3,3,4)
im = ax4.imshow(astropy_conv_intnan, vmin=-2., vmax=2.e1,
                origin='lower', interpolation='nearest', cmap='viridis')
ax4.set_title("astropy interpolate_nan")
ax4.set_xticklabels([])
ax4.set_yticklabels([])

ax5 = plt.subplot(3,3,5)
im = ax5.imshow(astropy_conv_intnan_norm, vmin=-2., vmax=2.e1,
                origin='lower', interpolation='nearest', cmap='viridis')
ax5.set_title("astropy interpolate_nan norm")
ax5.set_xticklabels([])
ax5.set_yticklabels([])

ax6 = plt.subplot(3,3,6)
im = ax6.imshow(astropy_conv, vmin=-2., vmax=2.e1, origin='lower',
                    interpolation='nearest', cmap='viridis')
ax6.set_title("Default astropy")
ax6.set_xticklabels([])
ax6.set_yticklabels([])

ax7 = plt.subplot(3,3,7)
im = ax7.imshow(astropy_conv_intnan_fft, vmin=-2., vmax=2.e1,
                origin='lower', interpolation='nearest', cmap='viridis')
ax7.set_title("astropyfft interpolate_nan")
ax7.set_xticklabels([])
ax7.set_yticklabels([])

ax8 = plt.subplot(3,3,8)
im = ax8.imshow(astropy_conv_intnan_fft_norm, vmin=-2., vmax=2.e1,
                origin='lower', interpolation='nearest', cmap='viridis')
ax8.set_title("astropyfft interpolate_nan norm")
ax8.set_xticklabels([])
ax8.set_yticklabels([])

ax9 = plt.subplot(3,3,9)
im = ax9.imshow(astropy_conv_fft, vmin=-2., vmax=2.e1, origin='lower',
                    interpolation='nearest', cmap='viridis')
ax9.set_title("Default fft astropy")
ax9.set_xticklabels([])
ax9.set_yticklabels([])

pl.figure(2).clf()
pl.plot(img[:,25], label='input', drawstyle='steps-mid', linewidth=2, alpha=0.5)
pl.plot(scipy_conv[:,25], label='scipy', drawstyle='steps-mid',
        linewidth=2, alpha=0.5, marker='s')
pl.plot(scipy_conv_zerod[:,25], label='scipy nan->zero', drawstyle='steps-mid',
        linewidth=2, alpha=0.5, marker='s')
pl.plot(astropy_conv[:,25], label='astropy', drawstyle='steps-mid', linewidth=2, alpha=0.5)
pl.plot(astropy_conv_intnan_norm[:,25], label='astropy intnan norm', drawstyle='steps-mid', linewidth=2, alpha=0.5)
pl.plot(astropy_conv_intnan[:,25], label='astropy intnan', drawstyle='steps-mid', linewidth=2, alpha=0.5)
pl.legend(loc='best')
