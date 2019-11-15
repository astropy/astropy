Convolution Kernels
*******************

Introduction and Concept
========================

The convolution module provides several built-in kernels to cover the most
common applications in astronomy. It is also possible to define custom kernels
from arrays or combine existing kernels to match specific applications.

Every filter kernel is characterized by its response function. For time series
we speak of an "impulse response function" or for images we call it "point
spread function." This response function is given for every kernel by a
`~astropy.modeling.FittableModel`, which is evaluated on a grid with
:func:`~astropy.convolution.discretize_model` to obtain a kernel
array, which can be used for discrete convolution with the binned data.


Examples
========

1D Kernels
----------

One application of filtering is to smooth noisy data. In this case we
consider a noisy Lorentz curve:

>>> import numpy as np
>>> from astropy.modeling.models import Lorentz1D
>>> from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
>>> lorentz = Lorentz1D(1, 0, 1)
>>> x = np.linspace(-5, 5, 100)
>>> data_1D = lorentz(x) + 0.1 * (np.random.rand(100) - 0.5)

Smoothing the noisy data with a `~astropy.convolution.Gaussian1DKernel`
with a standard deviation of 2 pixels:

>>> gauss_kernel = Gaussian1DKernel(2)
>>> smoothed_data_gauss = convolve(data_1D, gauss_kernel)

Smoothing the same data with a `~astropy.convolution.Box1DKernel` of width 5
pixels:

>>> box_kernel = Box1DKernel(5)
>>> smoothed_data_box = convolve(data_1D, box_kernel)

The following plot illustrates the results:

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.modeling.models import Lorentz1D
    from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel

    # Fake Lorentz data including noise
    lorentz = Lorentz1D(1, 0, 1)
    x = np.linspace(-5, 5, 100)
    data_1D = lorentz(x) + 0.1 * (np.random.rand(100) - 0.5)

    # Smooth data
    gauss_kernel = Gaussian1DKernel(2)
    smoothed_data_gauss = convolve(data_1D, gauss_kernel)
    box_kernel = Box1DKernel(5)
    smoothed_data_box = convolve(data_1D, box_kernel)

    # Plot data and smoothed data
    plt.plot(x, data_1D, label='Original')
    plt.plot(x, smoothed_data_gauss, label='Smoothed with Gaussian1DKernel')
    plt.plot(x, smoothed_data_box, label='Smoothed with Box1DKernel')
    plt.xlabel('x [a.u.]')
    plt.ylabel('amplitude [a.u.]')
    plt.xlim(-5, 5)
    plt.ylim(-0.1, 1.5)
    plt.legend(prop={'size':12})
    plt.show()


Beside the ``astropy`` convolution functions
`~astropy.convolution.convolve` and
`~astropy.convolution.convolve_fft`, it is also possible to use
the kernels with ``numpy`` or ``scipy`` convolution by passing the ``array``
attribute. This will be faster in most cases than the ``astropy`` convolution,
but will not work properly if NaN values are present in the data.

>>> smoothed = np.convolve(data_1D, box_kernel.array)

2D Kernels
----------

As all 2D kernels are symmetric, it is sufficient to specify the width in one
direction. Therefore the use of 2D kernels is basically the same as for 1D
kernels. Here we consider a small Gaussian-shaped source of amplitude 1 in the
middle of the image and add 10% noise:

>>> import numpy as np
>>> from astropy.convolution import convolve, Gaussian2DKernel, Tophat2DKernel
>>> from astropy.modeling.models import Gaussian2D
>>> gauss = Gaussian2D(1, 0, 0, 3, 3)
>>> # Fake image data including noise
>>> x = np.arange(-100, 101)
>>> y = np.arange(-100, 101)
>>> x, y = np.meshgrid(x, y)
>>> data_2D = gauss(x, y) + 0.1 * (np.random.rand(201, 201) - 0.5)

Smoothing the noisy data with a
:class:`~astropy.convolution.Gaussian2DKernel` with a standard
deviation of 2 pixels:

>>> gauss_kernel = Gaussian2DKernel(2)
>>> smoothed_data_gauss = convolve(data_2D, gauss_kernel)

Smoothing the noisy data with a
:class:`~astropy.convolution.Tophat2DKernel` of width 5 pixels:

>>> tophat_kernel = Tophat2DKernel(5)
>>> smoothed_data_tophat = convolve(data_2D, tophat_kernel)

This is what the original image looks like:

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.modeling.models import Gaussian2D
    gauss = Gaussian2D(1, 0, 0, 2, 2)
    # Fake image data including noise
    x = np.arange(-100, 101)
    y = np.arange(-100, 101)
    x, y = np.meshgrid(x, y)
    data_2D = gauss(x, y) + 0.1 * (np.random.rand(201, 201) - 0.5)
    plt.imshow(data_2D, origin='lower')
    plt.xlabel('x [pixels]')
    plt.ylabel('y [pixels]')
    plt.colorbar()
    plt.show()

The following plot illustrates the differences between several 2D kernels
applied to the simulated data. Note that it has a slightly different color
scale compared to the original image.

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt

    from astropy.convolution import *
    from astropy.modeling.models import Gaussian2D

    # Small Gaussian source in the middle of the image
    gauss = Gaussian2D(1, 0, 0, 2, 2)
    # Fake data including noise
    x = np.arange(-100, 101)
    y = np.arange(-100, 101)
    x, y = np.meshgrid(x, y)
    data_2D = gauss(x, y) + 0.1 * (np.random.rand(201, 201) - 0.5)

    # Setup kernels, including unity kernel for original image
    # Choose normalization for linear scale space for RickerWavelet

    kernels = [TrapezoidDisk2DKernel(11, slope=0.2),
               Tophat2DKernel(11),
               Gaussian2DKernel(11),
               Box2DKernel(11),
               11 ** 2 * RickerWavelet2DKernel(11),
               AiryDisk2DKernel(11)]

    fig, axes = plt.subplots(nrows=2, ncols=3)

    # Plot kernels
    for kernel, ax in zip(kernels, axes.flat):
        smoothed = convolve(data_2D, kernel, normalize_kernel=False)
        im = ax.imshow(smoothed, vmin=-0.01, vmax=0.08, origin='lower',
                       interpolation='None')
        title = kernel.__class__.__name__
        ax.set_title(title, fontsize=12)
        ax.set_yticklabels([])
        ax.set_xticklabels([])

    cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    fig.colorbar(im, cax=cax)
    plt.subplots_adjust(left=0.05, right=0.85, top=0.95, bottom=0.05)
    plt.show()


The Gaussian kernel has better smoothing properties compared to the Box and the
Top Hat. The Box filter is not isotropic and can produce artifacts (the source
appears rectangular). The Ricker Wavelet filter removes noise and slowly varying
structures (i.e., background), but produces a negative ring around the source.
The best choice for the filter strongly depends on the application.


Available Kernels
=================

.. currentmodule:: astropy.convolution

.. autosummary::

   AiryDisk2DKernel
   Box1DKernel
   Box2DKernel
   CustomKernel
   Gaussian1DKernel
   Gaussian2DKernel
   RickerWavelet1DKernel
   RickerWavelet2DKernel
   Model1DKernel
   Model2DKernel
   Ring2DKernel
   Tophat2DKernel
   Trapezoid1DKernel
   TrapezoidDisk2DKernel

Kernel Arithmetics
==================

Addition and Subtraction
------------------------

As convolution is a linear operation, kernels can be added or subtracted from
each other. They can also be multiplied with some number. One basic example
would be the definition of a Difference of Gaussian filter:

>>> from astropy.convolution import Gaussian1DKernel
>>> gauss_1 = Gaussian1DKernel(10)
>>> gauss_2 = Gaussian1DKernel(16)
>>> DoG = gauss_2 - gauss_1

Another application is to convolve faked data with an instrument response
function model. For example, if the response function can be described by
the weighted sum of two Gaussians:

>>> gauss_1 = Gaussian1DKernel(10)
>>> gauss_2 = Gaussian1DKernel(16)
>>> SoG = 4 * gauss_1 + gauss_2

Most times it will be necessary to normalize the resulting kernel by calling
explicitly:

>>> SoG.normalize()


Convolution
-----------

Furthermore, two kernels can be convolved with each other, which is useful when
data is filtered with two different kinds of kernels or to create a new,
special kernel:

>>> from astropy.convolution import Gaussian1DKernel, convolve
>>> gauss_1 = Gaussian1DKernel(10)
>>> gauss_2 = Gaussian1DKernel(16)
>>> broad_gaussian = convolve(gauss_2,  gauss_1)  # doctest: +IGNORE_WARNINGS

Or in case of multistage smoothing:

>>> import numpy as np
>>> from astropy.modeling.models import Lorentz1D
>>> from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
>>> lorentz = Lorentz1D(1, 0, 1)
>>> x = np.linspace(-5, 5, 100)
>>> data_1D = lorentz(x) + 0.1 * (np.random.rand(100) - 0.5)

>>> gauss = Gaussian1DKernel(3)
>>> box = Box1DKernel(5)
>>> smoothed_gauss = convolve(data_1D, gauss)
>>> smoothed_gauss_box = convolve(smoothed_gauss, box)

You would rather do the following:

>>> gauss = Gaussian1DKernel(3)
>>> box = Box1DKernel(5)
>>> smoothed_gauss_box = convolve(data_1D, convolve(box, gauss))  # doctest: +IGNORE_WARNINGS

Which, in most cases, will also be faster than the first method because only
one convolution with the often times larger data array will be necessary.

Discretization
==============

To obtain the kernel array for discrete convolution, the kernel's response
function is evaluated on a grid with
:func:`~astropy.convolution.discretize_model`. For the
discretization step the following modes are available:

* Mode ``'center'`` (default) evaluates the response function on the grid by
  taking the value at the center of the bin.

   >>> from astropy.convolution import Gaussian1DKernel
   >>> gauss_center = Gaussian1DKernel(3, mode='center')

* Mode ``'linear_interp'`` takes the values at the corners of the bin and
  linearly interpolates the value at the center:

  >>> gauss_interp = Gaussian1DKernel(3, mode='linear_interp')

* Mode ``'oversample'`` evaluates the response function by taking the mean on an
  oversampled grid. The oversample factor can be specified with the ``factor``
  argument. If the oversample factor is too large, the evaluation becomes slow.

 >>> gauss_oversample = Gaussian1DKernel(3, mode='oversample', factor=10)

* Mode ``'integrate'`` integrates the function over the pixel using
  ``scipy.integrate.quad`` and ``scipy.integrate.dblquad``. This mode is very
  slow and is only recommended when the highest accuracy is required.

.. doctest-requires:: scipy

    >>> gauss_integrate = Gaussian1DKernel(3, mode='integrate')

Especially in the range where the kernel width is in order of only a few pixels,
it can be advantageous to use the mode ``oversample`` or ``integrate`` to
conserve the integral on a subpixel scale.


Normalization
=============

The kernel models are normalized per default (i.e.,
:math:`\int_{-\infty}^{\infty} f(x) dx = 1`). But because of the limited kernel
array size, the normalization for kernels with an infinite response can differ
from one. The value of this deviation is stored in the kernel's ``truncation``
attribute.

The normalization can also differ from one, especially for small kernels, due
to the discretization step. This can be partly controlled by the ``mode``
argument, when initializing the kernel. (See also
:func:`~astropy.convolution.discretize_model`.) Setting the
``mode`` to ``'oversample'`` allows us to conserve the normalization even on the
subpixel scale.

The kernel arrays can be renormalized explicitly by calling either the
``normalize()`` method or by setting the ``normalize_kernel`` argument in the
:func:`~astropy.convolution.convolve` and
:func:`~astropy.convolution.convolve_fft` functions. The latter
method leaves the kernel itself unchanged but works with an internal normalized
version of the kernel.

Note that for :class:`~astropy.convolution.RickerWavelet1DKernel`
and :class:`~astropy.convolution.RickerWavelet2DKernel` there is
:math:`\int_{-\infty}^{\infty} f(x) dx = 0`. To define a proper normalization,
both filters are derived from a normalized Gaussian function.
