Convolution Kernels
===================

There are a few different convolution kernels available, including standard Gaussian and Box filters,
but also Mexican Hat or Top Hat filters. A complete list can be found on 


Introduction and Concept
------------------------
A filter kernel is characterized by its response function. In the 1D
case we speak of "impulse response function", in the 2D case we call it
"point spread function". This response function is given for every kernel
by an astropy ParametricModel, which is evaluated on a grid to obtain a filter
array, which can then be applied to the binned data.

Currently only ``symmetric´´ 2D kernels are supported.


Normalization
-------------

The kernel models are normalized per default, but because of the discretization process, 
the normalization can differ from one, especially for small kernels. The value of
this deviation is stored in the truncation attribute. 
The kernel arrays can be renormalized explicitly (Note that this is also controlled by the 
normalize_kernel flag in the convolution functions).  


Examples
--------






Kernel Arithmetics
------------------

As convolution is a linear operation, kernels can be added subtracted or multiplied with some
number. One basic example would be the definition of a Difference of Gaussian filter:

>>> gauss_1 = Gaussian1DKernel(10)
>>> gauss_2 = Gaussian2Dkernel(16)
>>> DoG = gauss_2 - gauss_1

Another use case is to fake data and convolve it with an instrument response function, 
which can be modeled by the sum of two Gaussians of different weight:

>>> gauss_1 = Gaussian1DKernel(10)
>>> gauss_2 = Gaussian2Dkernel(16)
>>> SoG = 4 * gauss_1 + gauss_2

Most times it will be necessary to normalize the resulting kernel by calling:

>>> SoG.normalize()




 

 