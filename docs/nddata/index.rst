.. _astropy_nddata:

*****************************************
N-dimensional datasets (`astropy.nddata`)
*****************************************

Introduction
============

`astropy.nddata` provides the `~astropy.nddata.nddata.NDData`
class and related tools to manage n-dimensional array-based data (e.g.
CCD images, IFU Data, grid-based simulation data, ...). This is more than
just `numpy.ndarray` objects, because it provides metadata that cannot
be easily provided by a single array.

This subpackage also provides new convolution routines that differ from
Scipy in that they offer a proper treatment of NaN values.


Getting started
===============

The `~astropy.nddata.nddata.NDData` class is still under development, and
support for WCS and units is not yet implemented.

Initializing
------------

An `~astropy.nddata.nddata.NDData` object can be instantiated by passing it an
n-dimensional Numpy array::

    >>> from astropy.nddata import NDData
    >>> array = np.random.random((12, 12, 12))  # a random 3-dimensional array
    >>> ndd = NDData(array)

This object has a few attributes in common with Numpy:

    >>> ndd.ndim
    3
    >>> ndd.shape
    (12, 12, 12)
    >>> ndd.dtype
    dtype('float64')

The underlying Numpy array can be accessed via the `data` attribute::

    >>> ndd.data
    array([[[ 0.05621944,  0.85569765,  0.71609697, ...,  0.76049288,
    ...

Mask
----

Values can be masked using the `mask` attribute, which should be a boolean
Numpy array with the same dimensions as the data, e.g.::

     >>> ndd.mask = ndd.data > 0.9

A mask value of `True` indicates a value that should be ignored, while a mask
value of `False` indicates a valid value.

Flags
-----

Values can be assigned one or more flags. The `flags` attribute is used to
store either a single Numpy array (of any type) with dimensions matching that
of the data, or a `~astropy.nddata.flag_collection.FlagCollection`, which is
essentially a dictionary of Numpy arrays (of any type) with the same shape as
the data. The following example demonstrates setting a single set of integer
flags::

    >>> ndd.flags = np.zeros(ndd.shape)
    >>> ndd.flags[ndd.data < 0.1] = 1
    >>> ndd.flags[ndd.data < 0.01] = 2

but one can also have multiple flag layers with different types::

    >>> from astropy.nddata import FlagCollection
    >>> ndd.flags = FlagCollection(shape=(12, 12, 12))
    >>> ndd.flags['photometry'] = np.zeros(ndd.shape, dtype=str)
    >>> ndd.flags['photometry'][ndd.data > 0.9] = 's'
    >>> ndd.flags['cosmic_rays'] = np.zeros(ndd.shape, dtype=int)
    >>> ndd.flags['cosmic_rays'][ndd.data > 0.99] = 99

and flags can easily be used to set the mask::

    >>> ndd.mask = ndd.flags['cosmic_rays'] == 99

Errors
------

`~astropy.nddata.nddata.NDData` objects have an `error` attribute that can be
used to set the uncertainty on the data values. This is done by using classes
to represent the errors of a given type. For example, to set standard
deviation errors on the pixel values, you can do::

    >>> from astropy.nddata import StandardDeviationError
    >>> ndd.error = StandardDeviationError(np.ones((12, 12, 12)) * 0.1)

Arithmetic
----------

Provided that the world coordinate system (WCS), units, and shape match, two
:class:`~astropy.nddata.nddata.NDData` instances can be added or subtracted from each
other, with error propagation, creating a new
:class:`~astropy.nddata.nddata.NDData` object::

    ndd3 = ndd1.add(ndd2)
    ndd4 = ndd1.subtract(ndd2)

The purpose of the :meth:`~astropy.nddata.nddata.NDData.add` and
:meth:`~astropy.nddata.nddata.NDData.subtract` methods is to allow the
combination of two data objects that have common WCS, units, and shape, with
consistent behavior for masks and flags, and with a framework to propagate
uncertainties. These methods are intended for use by sub-classes and functions
that deal with more complex combinations.

.. warning:: Error propagation is still experimental, and does not take into
             account correlated errors.

Convolution
-----------

`astropy.nddata` also contains convenience functions for n-dimensional
convolution that will later be integrated into the
`~astropy.nddata.nddata.NDData` class. For example, if your
data is 1D, you might smooth it by a gaussian kernel by doing::

    >>> from astropy.nddata import convolve, make_kernel
    >>> kernel = make_kernel((9,), 1.5, 'gaussian')
    >>> ndd.data = convolve(ndd.data, kernel)

The convolution routines can also be used on bare arrays.

Converting to Numpy arrays
--------------------------

`~astropy.nddata.nddata.NDData` objects can also be easily converted to
numpy arrays::

    >>> import numpy as np
    >>> arr = np.array(ndd)
    >>> np.all(arr == mydataarray)
    True

If a `mask` is defined, this will result in a `~numpy.ma.MaskedArray`, so
in all cases a useable `numpy.ndarray` or subclass will result. This allows
straightforward plotting of `~astropy.nddata.nddata.NDData` objects with 1-
and 2-dimensional datasets using `matplotlib`::

    >>> from matplotlib import pyplot as plt
    >>> plt.plot(ndd)

This works because the `matplotlib` plotting functions automatically convert
their inputs using `numpy.array`.


Using `nddata`
==============

.. toctree::
   :maxdepth: 2

   convolution.rst

Reference/API
=============

.. automodapi:: astropy.nddata
    :no-inheritance-diagram:
