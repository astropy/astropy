*************************************
Convolving with un-normalized kernels
*************************************

There are some tasks, such as source finding, where you want to apply a filter
with a kernel that is not normalized.

For data that are well-behaved (contain no missing or infinite values), this is
easy and can be done in one step::

    convolve(image, kernel)

For example, we can try to run a commonly-used peak enhancing kernel:

.. plot::
   :context: reset
   :include-source:
   :align: center

    import numpy as np
    import matplotlib.pyplot as plt

    from astropy.io import fits
    from astropy.utils.data import get_pkg_data_filename
    from astropy.convolution import CustomKernel
    from scipy.signal import convolve as scipy_convolve
    from astropy.convolution import convolve, convolve_fft


    # Load the data from data.astropy.org
    filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')
    hdu = fits.open(filename)[0]

    # Scale the file to have reasonable numbers
    # (this is mostly so that colorbars don't have too many digits)
    # Also, we crop it so you can see individual pixels
    img = hdu.data[50:90, 60:100] * 1e5

    kernel = CustomKernel([[-1,-1,-1], [-1, 8, -1], [-1,-1,-1]])

    astropy_conv = convolve(img, kernel, normalize_kernel=False, nan_treatment='fill')
    #astropy_conv_fft = convolve_fft(img, kernel, normalize_kernel=False, nan_treatment='fill')

    plt.figure(1, figsize=(12, 12)).clf()
    ax1 = plt.subplot(1, 2, 1)
    im = ax1.imshow(img, vmin=-6., vmax=5.e1, origin='lower',
                    interpolation='nearest', cmap='viridis')

    ax2 = plt.subplot(1, 2, 2)
    im = ax2.imshow(astropy_conv, vmin=-6., vmax=5.e1, origin='lower',
                    interpolation='nearest', cmap='viridis')
