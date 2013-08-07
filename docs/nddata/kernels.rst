Convolution Kernels
===================

There are a several built-in convolution kernels which should cover the most common applications.


Introduction and Concept
------------------------
A filter `~astropy.nddata.convolution.core.Kernel` is characterized by its response function. In the 1D
case we speak of "impulse response function", in the 2D case we call it "point spread function". 
This response function is given for every kernel by a `~astropy.modeling.core.ParametricModel`, 
which is evaluated on a grid with `~astropy.nddata.convolution.utils.discretize_model` to obtain a 
filter array, which can then be applied to the binned data.

**Note:** Currently only **symmetric** 2D kernels are supported.


Examples
--------



Kernel Arithmetics
------------------

As convolution is a linear operation, kernels can be added or subtracted from each other. They can aslo be multiplied with some
number. One basic example would be the definition of a Difference of Gaussian filter:

>>> gauss_1 = Gaussian1DKernel(10)
>>> gauss_2 = Gaussian2Dkernel(16)
>>> DoG = gauss_2 - gauss_1

Another application is to convolve faked data with an instrument response function model. 
E.g. if the response function can be be described by the weighted sum of two Gaussians:

>>> gauss_1 = Gaussian1DKernel(10)
>>> gauss_2 = Gaussian2Dkernel(16)
>>> SoG = 4 * gauss_1 + gauss_2

Most times it will be necessary to normalize the resulting kernel by calling:

>>> SoG.normalize()


Available Kernels
-----------------

.. automodsumm:: astropy.nddata.convolution.kernels
	:classes-only:
	

Normalization
-------------

The kernel models are normalized per default, i.e. :math:`\int_{-\infty}^{\infty} f(x) dx = 1`. But because of the limited 
kernel array size the normalization for kernels with an infinite response can differ from one. 
The value of this deviation is stored in the kernel's ``truncation`` attribute.

The normalization can also differ from one, especially for small kernels, due to the discretization process.
This can be controlled by the ``mode`` argument, when initializing the kernel (See also 
`~astropy.nddata.convolution.utils.discretize_model`).
 
The kernel arrays can be renormalized explicitly by calling either the ``normalize()`` method or by setting
the ``normalize_kernel`` argument in the `~astropy.nddata.convolution.convolve.convolve` and 
`~astropy.nddata.convolution.convolve.convolve_fft` functions. 
	 