.. _astropy_nddata:

*****************************************
N-dimensional datasets (`astropy.nddata`)
*****************************************

Introduction
============

`astropy.nddata` provides the `~astropy.nddata.nddata.NDData`
class and related tools to manage n-dimensional array-based data (e.g.
CCD images, IFU data, grid-based simulation data, ...). This is more than
just `numpy.ndarray` objects, because it provides metadata that cannot
be easily provided by a single array.

This subpackage also provides new convolution routines that differ from
Scipy in that they offer a proper treatment of NaN values.

.. note:: The `~astropy.nddata.nddata.NDData` class is still under
          development, and support for WCS and units is not yet implemented.

Getting started
===============

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

Values can be masked using the `mask` attribute, which should be a boolean
Numpy array with the same dimensions as the data, e.g.::

     >>> ndd.mask = ndd.data > 0.9

A mask value of `True` indicates a value that should be ignored, while a mask
value of `False` indicates a valid value.

Similarly, attributes are available to store generic meta-data, flags, and
uncertainties, and the `~astropy.nddata.nddata.NDData` class includes methods to
combine datasets with arithmetic operations (which include uncertainties propagation).
These are described in :doc:`nddata`.

Using `nddata`
==============

.. toctree::
   :maxdepth: 2

   nddata.rst
   convolution.rst
   subclassing.rst

Reference/API
=============

.. automodapi:: astropy.nddata
    :no-inheritance-diagram:
