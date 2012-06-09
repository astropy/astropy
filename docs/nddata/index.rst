.. _astropy_nddata:

=========================================
N-dimensional datasets (`astropy.nddata`)
=========================================

Introduction
------------

`astropy.nddata` provides the `~astropy.nddata.nddata.NDData`
class and related tools to manage n-dimensional array-based data (e.g.
CCD images, IFU Data, grid-based simulation data, ...). This is more than
just `numpy.ndarray` objects, because it provides metadata that cannot
be easily provided by a single array.

Getting started
---------------

At the moment, the `~astropy.nddata.nddata.NDData` class is still under development. You can however already make use of the convolution routines included in `astropy.nddata`::

    from astropy.nddata import convolve, convolve_fft, make_kernel

These routines differ from the Scipy convolution routines in that they offer a proper treatment of NaN values. Both functions are used as::

    result = convolve(image, kernel)

For example, you can smooth a 1D array with a custom kernel using::

    >>> convolve([1, 4, 5, 6, 5, 7, 8], [0.2, 0.6, 0.2])
    array([ 0. ,  3.4,  5. ,  5.6,  5.6,  5.2,  0. ])

You can also use `make_kernel` to generate n-dimensional kernels::

    >>> make_kernel([3,3], 1, 'boxcar')
    array([[ 0.  0.  0.]
           [ 0.  1.  0.]
           [ 0.  0.  0.]])

See the documentation below for more information.

Using `nddata`
--------------

.. toctree::
   :maxdepth: 2

   convolution.rst

Reference/API
-------------

.. automodapi:: astropy.nddata
    :no-inheritance-diagram:
