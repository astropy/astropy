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

This subpackage also provides new convolution routines that differ from
Scipy in that they offer a proper treatment of NaN values.


Getting started
---------------

The `~astropy.nddata.nddata.NDData` class is still under development, and
many of it's more advanced features are not yet implemented.  It already
functions as an array container with metadata, however::

    >>> from astropy.nddata import NDData
    >>> ndd = NDData(mydataarray, error=myerrorarray)
    >>> ndd['Exposure time(s)'] = 5

You can, also make use of the new convolution routines. For example, if your
data is 1D, you might smooth it by a gaussian kernel by doing::

    >>> from astropy.nddata import convolve, make_kernel
    >>> kernel = make_kernel((9,), 1.5, 'gaussian')
    >>> ndd.data = convolve(ndd.data, kernel)

The convolution routines can also be used on bare arrays.

`~astropy.nddata.nddata.NDData` objects can also be easily converted to
numpy arrays::

    >>> import numpy as np
    >>> arr = np.array(ndd)
    >>> np.all(arr == mydataarray)
    True

If a `mask` is defined, this will result in a `~numpy.ma.core.MaskedArray`, so
in all cases a useable `numpy.ndarray` or subclass will result. This allows
straightforward plotting of `~astropy.nddata.nddata.NDData` objects with
`matplotlib`::

    >>> from matplotlib import pyplot as plt
    >>> plt.plot(ndd)

This works because the `matplotlib` plotting functions automatically convert
their inputs using `numpy.array`.


Using `nddata`
--------------

.. toctree::
   :maxdepth: 2

   convolution.rst

Reference/API
-------------

.. automodapi:: astropy.nddata
    :no-inheritance-diagram:
