Convolution Kernels
===================

Introduction and Concept
------------------------
The convolution module provides several built-in kernels to cover the most common applications in astronomy.
It is also possible to define custom kernels from arrays or combine existing kernels to match specific 
applications. 

Every filter kernel is characterized by its response function. For time series we speak of an 
"impulse response function" or for images we call it "point spread function". 
This response function is given for every kernel by a `~astropy.modeling.core.ParametricModel`, 
which is evaluated on a grid with :func:`~astropy.nddata.convolution.utils.discretize_model` to obtain a 
kernel array, which can be used for discrete convolution with the binned data. 


Examples
--------

1D Kernels
^^^^^^^^^^

One application of filtering is to smooth noisy data. In this case we consider a noisy 
Lorentz curve: 

>>> import numpy as np
>>> from astropy.modeling.models import Lorentz1DModel
>>> from astropy.nddata.convolution import convolve, Gaussian1DKernel, Box1DKernel
>>> lorentz = Lorentz1DModel(1, 0, 1)
>>> x = np.linspace(-5, 5, 100)
>>> data = lorentz(x) + 0.1 * (np.random.rand(100) - 0.5)

Smoothing the noisy data with a `~astropy.nddata.convolution.kernels.Gaussian1DKernel` of width 2 pixels:

>>> gauss_kernel = Gaussian1DKernel(2)
>>> smoothed_data_gauss = convolve(data, gauss_kernel)

Smoothing the same data with a `~astropy.nddata.convolution.kernels.Box1DKernel` of width 5 pixels:

>>> box_kernel = Box1DKernel(5)
>>> smoothed_data_box = convolve(data, box_kernel)

The following plot illustrates the results:

.. plot::

	import numpy as np
	import matplotlib.pyplot as plt
	from astropy.modeling.models import Lorentz1DModel
	from astropy.nddata.convolution import convolve, Gaussian1DKernel, Box1DKernel
	
	# Fake Lorentz data including noise
	lorentz = Lorentz1DModel(1, 0, 1)
	x = np.linspace(-5, 5, 100)
	data = lorentz(x) + 0.1 * (np.random.rand(100) - 0.5)
	
	# Smooth data
	gauss_kernel = Gaussian1DKernel(2)
	smoothed_data_gauss = convolve(data, gauss_kernel)
	box_kernel = Box1DKernel(5)
	smoothed_data_box = convolve(data, box_kernel)

	# Plot data and smoothed data
	plt.plot(x, data, label='Original')
	plt.plot(x, smoothed_data_gauss, label='Smoothed with Gaussian1DKernel')
	plt.plot(x, smoothed_data_box, label='Smoothed with Box1DKernel')
	plt.xlabel('x [a.u.]')
	plt.ylabel('amplitude [a.u.]')
	plt.xlim(-5, 5)
	plt.ylim(-0.1, 1.5)
	plt.legend(prop={'size':12})
	plt.show()


Beside the astropy convolution functions  `~astropy.nddata.convolution.convolve.convolve` and 
`~astropy.nddata.convolution.convolve.convolve_fft`, it is also possible to use the kernels 
with Numpy or Scipy convolution by passing the ``array`` attribute. This will be faster in most
cases than the astropy convolution, but will not work properly if ``NaN`` values are present in the data.

>>> smoothed = np.convolve(data, box_kernel.array)

2D Kernels
^^^^^^^^^^
As all 2D kernels are symmetric it is sufficient to specify the width in one direction.
Therefore the use of 2D kernels is basically the same as for 1D kernels. We consider a 
small Gaussian shaped source of amplitude one in the middle of the image and add 10% noise: 

>>> import numpy as np
>>> from astropy.nddata.convolution import convolve, Gaussian2DKernel, TopHat2DKernel
>>> from astropy.modeling.models import Gaussian2DModel
>>> gauss = Gaussian2DModel(1, 0, 0, 3, 3)
>>> # Fake image data including noise
>>> x = np.arange(-100, 101)
>>> y = np.arange(-100, 101)
>>> x, y = np.meshgrid(x, y)
>>> data = gauss(x, y) + 0.1 * (np.random.rand(201, 201) - 0.5)

Smoothing the noisy data with a `~astropy.nddata.convolution.kernels.Gaussian2DKernel` of width 2 pixels:

>>> gauss_kernel = Gaussian2DKernel(2)
>>> smoothed_data_gauss = convolve(data, gauss_kernel)

Smoothing the noisy data with a `~astropy.nddata.convolution.kernels.Tophat2DKernel` of width 5 pixels:

>>> tophat_kernel = TopHat2DKernel(5)
>>> smoothed_data_tophat = convolve(data, tophat_kernel)

This is what the original image looks like:

.. plot::

	import numpy as np
	import matplotlib.pyplot as plt
	from astropy.modeling.models import Gaussian2DModel
	gauss = Gaussian2DModel(1, 0, 0, 3, 3)
	# Fake image data including noise
	x = np.arange(-100, 101)
	y = np.arange(-100, 101)
	x, y = np.meshgrid(x, y)
	data = gauss(x, y) + 0.1 * (np.random.rand(201, 201) - 0.5)
	plt.imshow(data, origin='lower')
	plt.xlabel('x [pixels]')
	plt.ylabel('y [pixels]')
	plt.colorbar()
	plt.show()

The following plot illustrates the differences between several 2D kernels applied to the simulated data. 
Note that it has a slightly different color scale compared to the original image.  

.. plot:: 
		
	import numpy as np
	import matplotlib.pyplot as plt
	
	from astropy.nddata.convolution import *
	from astropy.modeling.models import Gaussian2DModel

	# Small Gaussian source in the middle of the image
	gauss = Gaussian2DModel(1, 0, 0, 3, 3)
	# Fake data including noise
	x = np.arange(-100, 101)
	y = np.arange(-100, 101)
	x, y = np.meshgrid(x, y)
	data = gauss(x, y) + 0.1 * (np.random.rand(201, 201) - 0.5)
	
	# Setup kernels, including unity kernel for original image
	# Choose normalization for linear scale space for MexicanHat
	
	kernels = [TrapezoidDisk2DKernel(11, slope=0.2),
			   Tophat2DKernel(11),
			   Gaussian2DKernel(11),
			   Box2DKernel(22),
			   - 11 ** 2 * MexicanHat2DKernel(11),
			   AiryDisk2DKernel(22)]

	fig, axes = plt.subplots(nrows=2, ncols=3)
		
	# Plot kernels
	for kernel, ax in zip(kernels, axes.flat):
		smoothed = convolve(data, kernel)
		im = ax.imshow(smoothed, vmin=-0.01, vmax=0.15, origin='lower')
		title = kernel.__class__.__name__
		ax.set_title(title)
		ax.set_yticklabels([])
		ax.set_xticklabels([])
			
	cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
	fig.colorbar(im, cax=cax)
	plt.subplots_adjust(left=0.05, right=0.85, top=0.95, bottom=0.05)
	plt.show()
	

The Gaussian kernel has better smoothing properties compared to the Box and the Tophat. The Box filter is not isotropic
and can produce artifact (the source appears rectangular). The Mexican-Hat filter is almost noise free, but produces a negative
ring around the source. The best choice for the filter strongly depends on the application. 


Available Kernels
-----------------

.. automodsumm:: astropy.nddata.convolution.kernels
	:classes-only:


Kernel Arithmetics
------------------

As convolution is a linear operation, kernels can be added or subtracted from each other. They can also be multiplied with some
number. One basic example would be the definition of a Difference of Gaussian filter:

>>> gauss_1 = Gaussian1DKernel(10)
>>> gauss_2 = Gaussian2Dkernel(16)
>>> DoG = gauss_2 - gauss_1

Another application is to convolve faked data with an instrument response function model. 
E.g. if the response function can be be described by the weighted sum of two Gaussians:

>>> gauss_1 = Gaussian1DKernel(10)
>>> gauss_2 = Gaussian2Dkernel(16)
>>> SoG = 4 * gauss_1 + gauss_2

Most times it will be necessary to normalize the resulting kernel by calling explicitly:

>>> SoG.normalize()

Normalization
-------------

The kernel models are normalized per default, i.e. :math:`\int_{-\infty}^{\infty} f(x) dx = 1`. But because of the limited 
kernel array size the normalization for kernels with an infinite response can differ from one. 
The value of this deviation is stored in the kernel's ``truncation`` attribute.

The normalization can also differ from one, especially for small kernels, due to the discretization step.
This can be partly controlled by the ``mode`` argument, when initializing the kernel (See also 
:func:`~astropy.nddata.convolution.utils.discretize_model`). Setting the ``mode`` to ``'oversample'`` allows
to conserve the normalization even on the subpixel scale.
 
The kernel arrays can be renormalized explicitly by calling either the ``normalize()`` method or by setting
the ``normalize_kernel`` argument in the `~astropy.nddata.convolution.convolve.convolve` and 
`~astropy.nddata.convolution.convolve.convolve_fft` functions. The latter method leaves the kernel itself unchanged
but works with an internal normalized version of the kernel.   

Note that for `~astropy.nddata.convolution.kernels.MexicanHat1DKernel` 
and `~astropy.nddata.convolution.kernels.MexicanHat2DKernel` there is :math:`\int_{-\infty}^{\infty} f(x) dx = 0`. 
To define a proper normalization both filters are derived from a normalized Gaussian function. 

 
	 