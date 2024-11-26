.. _astropy_convolve:

*************************************************
Convolution and Filtering (`astropy.convolution`)
*************************************************

Introduction
============

`astropy.convolution` provides convolution functions and kernels that offer
improvements compared to the SciPy `scipy.ndimage` convolution routines,
including:

* Proper treatment of NaN values (ignoring them during convolution and
  replacing NaN pixels with interpolated values)

* Both direct and Fast Fourier Transform (FFT) versions

* Built-in kernels that are commonly used in Astronomy

The following thumbnails show the difference between SciPy and
Astropy's convolve functions on an astronomical image that contains NaN
values. Scipy's function returns NaN for all pixels that are within a
kernel-sized region of any NaN value, which is often not the desired
result.

.. plot::
    :show-source-link:

    import matplotlib.pyplot as plt
    import numpy as np
    from astropy.convolution import Gaussian2DKernel, convolve
    from astropy.io import fits
    from astropy.utils.data import get_pkg_data_filename
    from scipy.ndimage import convolve as scipy_convolve

    # Load the data from data.astropy.org
    filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')
    hdu = fits.open(filename)[0]

    # Scale the file to have reasonable numbers
    # (this is mostly so that colorbars do not have too many digits)
    # Also, we crop it so you can see individual pixels
    img = hdu.data[50:90, 60:100] * 1e5

    # This example is intended to demonstrate how astropy.convolve and
    # scipy.convolve handle missing data, so we start by setting the
    # brightest pixels to NaN to simulate a "saturated" data set
    img[img > 20] = np.nan

    # We also create a copy of the data and set those NaNs to zero.  We will
    # use this for the scipy convolution
    img_zerod = img.copy()
    img_zerod[np.isnan(img)] = 0

    # We smooth with a Gaussian kernel with x_stddev=1 (and y_stddev=1)
    # It is a 9x9 array
    kernel = Gaussian2DKernel(x_stddev=1)

    # Convolution: scipy's direct convolution mode spreads out NaNs (see
    # panel 2 below)
    scipy_conv = scipy_convolve(img, kernel)

    # scipy's direct convolution mode run on the 'zero'd' image will not
    # have NaNs, but will have some very low value zones where the NaNs were
    # (see panel 3 below)
    scipy_conv_zerod = scipy_convolve(img_zerod, kernel)

    # astropy's convolution replaces the NaN pixels with a kernel-weighted
    # interpolation from their neighbors
    astropy_conv = convolve(img, kernel)

    # Now we do a bunch of plots.  In the first two plots, the originally masked
    # values are marked with red X's
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(8, 8), layout="tight")
    fig.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95,
                        wspace=0.3, hspace=0.3)
    ax = ax.flatten()
    for axis in ax:
        axis.axis('off')

    im = ax[0].imshow(img, vmin=-2.0, vmax=20.0, origin='lower',
                      interpolation='nearest', cmap='viridis')
    y, x = np.where(np.isnan(img))
    ax[0].plot(x, y, 'rx', markersize=4)
    ax[0].set_title("Input Data")
    ax[0].set_xticklabels([])
    ax[0].set_yticklabels([])

    im = ax[1].imshow(scipy_conv, vmin=-2.0, vmax=20.0, origin='lower',
                      interpolation='nearest', cmap='viridis')
    ax[1].plot(x, y, 'rx', markersize=4)
    ax[1].set_title("Scipy convolved")
    ax[1].set_xticklabels([])
    ax[1].set_yticklabels([])

    im = ax[2].imshow(scipy_conv_zerod, vmin=-2.0, vmax=20.0, origin='lower',
                      interpolation='nearest', cmap='viridis')
    ax[2].set_title("Scipy convolved (NaN to zero)")
    ax[2].set_xticklabels([])
    ax[2].set_yticklabels([])

    im = ax[3].imshow(astropy_conv, vmin=-2.0, vmax=20.0, origin='lower',
                      interpolation='nearest', cmap='viridis')
    ax[3].set_title("Astropy convolved")
    ax[3].set_xticklabels([])
    ax[3].set_yticklabels([])

    # we make a second plot of the amplitudes vs offset position to more
    # clearly illustrate the value differences
    fig2, ax2 = plt.subplots(figsize=(8, 6))
    ax2.plot(img[:, 25], label='Input data', drawstyle='steps-mid', linewidth=2,
             alpha=0.5)
    ax2.plot(scipy_conv[:, 25], label='SciPy convolved', drawstyle='steps-mid',
             linewidth=2, alpha=0.5, marker='s')
    ax2.plot(scipy_conv_zerod[:, 25], label='SciPy convolved (NaN to zero)',
             drawstyle='steps-mid', linewidth=2, alpha=0.5, marker='s')
    ax2.plot(astropy_conv[:, 25], label='Astropy convolved', drawstyle='steps-mid',
             linewidth=2, alpha=0.5)
    ax2.set(xlabel="Pixel", ylabel="Amplitude")
    ax2.legend(loc='best')
    plt.show()


Getting Started
===============

Two convolution functions are provided. They are imported as::

    from astropy.convolution import convolve, convolve_fft

and are both used as::

    result = convolve(image, kernel)
    result = convolve_fft(image, kernel)

:func:`~astropy.convolution.convolve` is implemented as a direct convolution
algorithm, while :func:`~astropy.convolution.convolve_fft` uses a Fast Fourier
Transform (FFT). Thus, the former is better for small kernels, while the latter
is much more efficient for larger kernels.

Example
-------

..
  EXAMPLE START
  Convolution for User-Specified Kernels

To convolve a 1D dataset with a user-specified kernel, you can do::

    >>> from astropy.convolution import convolve
    >>> convolve([1, 4, 5, 6, 5, 7, 8], [0.2, 0.6, 0.2])  # doctest: +FLOAT_CMP
    array([1.4, 3.6, 5. , 5.6, 5.6, 6.8, 6.2])

The ``boundary`` keyword determines how the input array is extended
beyond its boundaries. The default value is ``'fill'``, meaning values
outside of the array boundary are set to the input ``fill_value``
(default is 0.0). Setting ``boundary='extend'`` causes values near the
edges to be extended using a constant extrapolation beyond the boundary.
The values at the end are computed assuming that any value below the
first point is ``1``, and any value above the last point is ``8``::

    >>> from astropy.convolution import convolve
    >>> convolve([1, 4, 5, 6, 5, 7, 8], [0.2, 0.6, 0.2], boundary='extend')  # doctest: +FLOAT_CMP
    array([1.6, 3.6, 5. , 5.6, 5.6, 6.8, 7.8])

For a more detailed discussion of boundary treatment, see :doc:`using`.

..
  EXAMPLE END

Example
-------

..
  EXAMPLE START
  Convolution for Built-In Kernels

The convolution module also includes built-in kernels that can be imported as,
for example::

    >>> from astropy.convolution import Gaussian1DKernel

To use a kernel, first create a specific instance of the kernel::

    >>> gauss = Gaussian1DKernel(stddev=2)

``gauss`` is not an array, but a kernel object. The underlying array can be
retrieved with::

    >>> gauss.array  # doctest: +FLOAT_CMP
    array([6.69162896e-05, 4.36349021e-04, 2.21596317e-03, 8.76430436e-03,
           2.69959580e-02, 6.47599366e-02, 1.20987490e-01, 1.76035759e-01,
           1.99474648e-01, 1.76035759e-01, 1.20987490e-01, 6.47599366e-02,
           2.69959580e-02, 8.76430436e-03, 2.21596317e-03, 4.36349021e-04,
           6.69162896e-05])

The kernel can then be used directly when calling
:func:`~astropy.convolution.convolve`:

.. plot::
    :show-source-link:

    import numpy as np
    import matplotlib.pyplot as plt

    from astropy.convolution import Gaussian1DKernel, convolve

    fig, ax = plt.subplots()

    # Generate fake data
    rng = np.random.default_rng(963)
    x = np.arange(1000).astype(float)
    y = np.sin(x / 100.) + rng.normal(0., 1., x.shape)
    y[::3] = np.nan

    # Create kernel
    g = Gaussian1DKernel(stddev=50)

    # Convolve data
    z = convolve(y, g)

    # Plot data before and after convolution
    ax.plot(x, y, label='Data')
    ax.plot(x, z, label='Convolved Data', linewidth=2)
    ax.legend(loc='best')
    plt.show()

..
  EXAMPLE END

Using ``astropy``'s Convolution to Replace Bad Data
---------------------------------------------------

``astropy``'s convolution methods can be used to replace bad data with values
interpolated from their neighbors. Kernel-based interpolation is useful for
handling images with a few bad pixels or for interpolating sparsely sampled
images.

The interpolation tool is implemented and used as::

    from astropy.convolution import interpolate_replace_nans
    result = interpolate_replace_nans(image, kernel)

Some contexts in which you might want to use kernel-based interpolation include:

 * Images with saturated pixels. Generally, these are the highest-intensity
   regions in the imaged area, and the interpolated values are not reliable,
   but this can be useful for display purposes.
 * Images with flagged pixels (e.g., a few small regions affected by cosmic
   rays or other spurious signals that require those pixels to be flagged out).
   If the affected region is small enough, the resulting interpolation will have
   a small effect on source statistics and may allow for robust source-finding
   algorithms to be run on the resulting data.
 * Sparsely sampled images such as those constructed with single-pixel
   detectors. Such images will only have a few discrete points sampled across
   the imaged area, but an approximation of the extended sky emission can still
   be constructed.

.. note::
    Care must be taken to ensure that the kernel is large enough to completely
    cover potential contiguous regions of NaN values.
    An ``AstropyUserWarning`` is raised if NaN values are detected post-
    convolution, in which case the kernel size should be increased.

Example
^^^^^^^

..
  EXAMPLE START
  Kernel Interpolation to Fill in Flagged-Out Pixels

The script below shows an example of kernel interpolation to fill in
flagged-out pixels:

.. plot::
   :context:
   :show-source-link:

   import numpy as np
   import matplotlib.pyplot as plt

   from astropy.io import fits
   from astropy.utils.data import get_pkg_data_filename
   from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans

   # Load the data from data.astropy.org
   filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')

   hdu = fits.open(filename)[0]
   img = hdu.data[50:90, 60:100] * 1e5

   # This example is intended to demonstrate how astropy.convolve and
   # scipy.convolve handle missing data, so we start by setting the brightest
   # pixels to NaN to simulate a "saturated" data set
   img[img > 2e1] = np.nan

   # We smooth with a Gaussian kernel with x_stddev=1 (and y_stddev=1)
   # It is a 9x9 array
   kernel = Gaussian2DKernel(x_stddev=1)

   # create a "fixed" image with NaNs replaced by interpolated values
   fixed_image = interpolate_replace_nans(img, kernel)

   # Now we do a bunch of plots.  In the first two plots, the originally masked
   # values are marked with red X's
   fig, axs = plt.subplots(ncols=2, figsize=(12, 6))

   ax1 = axs[0]
   im = ax1.imshow(img, vmin=-2., vmax=2.e1, origin='lower',
                   interpolation='nearest', cmap='viridis')
   y, x = np.where(np.isnan(img))
   ax1.plot(x, y, 'rx', markersize=4)
   ax1.set(
       autoscale_on=False,
       title="Original",
       xticklabels=[],
       yticklabels=[],
   )

   ax2 = axs[1]
   im = ax2.imshow(fixed_image, vmin=-2., vmax=2.e1, origin='lower',
                   interpolation='nearest', cmap='viridis')
   ax2.set(
       title="Fixed",
       xticklabels=[],
       yticklabels=[],
   )

..
  EXAMPLE END

Example
^^^^^^^

..
  EXAMPLE START
  Kernel Interpolation to Reconstruct Images from Sparse Sampling.

This script shows the power of this technique for reconstructing images from
sparse sampling. Note that the image is not perfect: the pointlike sources
are sometimes missed, but the extended structure is very well recovered by
eye.

.. plot::
   :context:
   :show-source-link:

   import numpy as np
   import matplotlib.pyplot as plt

   from astropy.io import fits
   from astropy.utils.data import get_pkg_data_filename
   from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans

   # Load the data from data.astropy.org
   filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')

   hdu = fits.open(filename)[0]
   img = hdu.data[50:90, 60:100] * 1e5

   rng = np.random.default_rng(1379)
   indices = rng.integers(low=0, high=img.size, size=300)

   sampled_data = img.flat[indices]

   # Build a new, sparsely sampled version of the original image
   new_img = np.tile(np.nan, img.shape)
   new_img.flat[indices] = sampled_data

   # We smooth with a Gaussian kernel with x_stddev=1 (and y_stddev=1)
   # It is a 9x9 array
   kernel = Gaussian2DKernel(x_stddev=1)

   # create a "reconstructed" image with NaNs replaced by interpolated values
   reconstructed_image = interpolate_replace_nans(new_img, kernel)

   # Now we do a bunch of plots.  In the first two plots, the originally masked
   # values are marked with red X's
   fig, axs = plt.subplots(ncols=3, figsize=(12, 6))

   ax1 = axs[0]
   im = ax1.imshow(img, vmin=-2., vmax=2.e1, origin='lower',
                   interpolation='nearest', cmap='viridis')
   y, x = np.where(np.isnan(img))
   ax1.set(
       autoscale_on=False,
       title="Original",
       xticklabels=[],
       yticklabels=[],
   )

   ax2 = axs[1]
   im = ax2.imshow(new_img, vmin=-2., vmax=2.e1, origin='lower',
                   interpolation='nearest', cmap='viridis')
   ax2.set(
       title="Sparsely Sampled",
       xticklabels=[],
       yticklabels=[],
   )

   ax3 = axs[2]
   im = ax3.imshow(reconstructed_image, vmin=-2., vmax=2.e1, origin='lower',
                   interpolation='nearest', cmap='viridis')
   ax3.set(
       title="Reconstructed",
       xticklabels=[],
       yticklabels=[],
   )

..
  EXAMPLE END

Using `astropy.convolution`
===========================

.. toctree::
   :maxdepth: 2

   using.rst
   kernels.rst
   non_normalized_kernels.rst

.. note that if this section gets too long, it should be moved to a separate
   doc page - see the top of performance.inc.rst for the instructions on how to
   do that
.. include:: performance.inc.rst

Reference/API
=============

.. toctree::
   :maxdepth: 2

   ref_api
