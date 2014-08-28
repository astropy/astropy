.. _astropy_nddata:

*****************************************
N-dimensional datasets (`astropy.nddata`)
*****************************************

Introduction
============

`astropy.nddata` provides the `~astropy.nddata.NDData`
class and related tools to manage n-dimensional array-based data (e.g.
CCD images, IFU data, grid-based simulation data, ...). This is more than
just `numpy.ndarray` objects, because it provides metadata that cannot
be easily provided by a single array.

.. warning::

  The `~astropy.nddata.NDData` class has changed significantly in astropy 1.0. See the section :ref:`nddata_transition` for more information.

Getting started
===============

An `~astropy.nddata.NDData` object can be instantiated by passing it an
n-dimensional Numpy array::

    >>> import numpy as np
    >>> from astropy.nddata import NDData
    >>> array = np.zeros((12, 12, 12))  # a 3-dimensional array with all zeros
    >>> ndd = NDData(array)

This object has a few attributes in common with Numpy:

    >>> ndd.ndim
    3
    >>> ndd.shape
    (12, 12, 12)
    >>> ndd.dtype
    dtype('float64')

Although the underlying Numpy array can be accessed via the ``data`` attribute::

    >>> ndd.data
    array([[[ 0., 0., 0., ...
    ...

the preferred way to access that data (because it will incorporate mask
information) is either by using ``numpy.asarray`` or using the ``NDData`` object as if it were a numpy array::

    >>> np.asarray(ndd)
    array([[[ 0., 0., 0., ...
    ...

Values can be masked using the ``mask`` attribute, which should be a boolean
Numpy array with the same dimensions as the data, e.g.::

     >>> ndd.mask = ndd.data > 0.9

A mask value of `True` indicates a value that should be ignored, while a mask
value of `False` indicates a valid value.

Similarly, attributes are available to store generic meta-data, flags, and
uncertainties. The `~astropy.nddata.NDDataArithmetic` class includes methods to
combine datasets with arithmetic operations (which include limited support for
propagation of uncertainties). These are described in :doc:`nddata`.

.. _nddata_transition:

Transition from astropy pre-1.0 to 1.0
======================================

The most important change is that arithmetic is no longer included in the base
`~astropy.nddata.NDData` class. Code that only uses the metadata features of
`~astropy.nddata.NDData` should not need to be modified. Code that uses the
arithemtic methods that used to be included in `~astropy.nddata.NDData` should
instead subclass `~astropy.nddata.NDDataArithmetic`; that class is equivalent
to the original `~astropy.nddata.NDData` class.

Using ``nddata``
================

.. toctree::
   :maxdepth: 2

   nddata.rst
   subclassing.rst
   decorator.rst

Reference/API
=============

.. automodapi:: astropy.nddata
    :no-inheritance-diagram:

.. automodapi:: astropy.nddata.utils
    :no-inheritance-diagram:
