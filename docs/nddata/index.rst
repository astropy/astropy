.. _astropy_nddata:

*****************************************
N-dimensional datasets (`astropy.nddata`)
*****************************************

Introduction
============

The ``nddata`` package provides a uniform interface to N-dimensional datasets
(for tabulated data please have a look at `astropy.table`) in astropy through:

+ `~astropy.nddata.NDData`: a basic container for `numpy.ndarray`-like data.
+ `~astropy.nddata.NDDataRef`: like `~astropy.nddata.NDData` but with
  additional functionality like an reading/writing, simple arithmetic
  operations and slicing.
+ `~astropy.nddata.StdDevUncertainty` a class that can store and propagate
  uncertainties for a NDData object.
+ General utility functions (:ref:`nddata_utils`) for operations on these
  classes, including a decorator to facilitate writing functions for such
  classes.
+ A framework including different mixins and metaclasses to facilitate
  customizing subclasses.

Getting started
===============

NDData
------

The primary purpose of `~astropy.nddata.NDData` is to act as a *container* for
data, metadata, and other related information like a mask.

An `~astropy.nddata.NDData` object can be instantiated by passing it an
n-dimensional `numpy` array::

    >>> import numpy as np
    >>> from astropy.nddata import NDData
    >>> array = np.zeros((12, 12, 12))  # a 3-dimensional array with all zeros
    >>> ndd1 = NDData(array)

or something that can be converted to an `numpy.ndarray`::

    >>> ndd2 = NDData([1, 2, 3, 4])
    >>> ndd2
    NDData([1, 2, 3, 4])

and can be accessed again via the ``data`` attribute::

    >>> ndd2.data
    array([1, 2, 3, 4])

It also supports additional properties like a ``unit`` or ``mask`` for the
data, a ``wcs`` (world coordinate system) and ``uncertainty`` of the data and
additional ``meta`` attributes:

    >>> data = np.array([1,2,3,4])
    >>> mask = data > 2
    >>> unit = 'erg / s'
    >>> from astropy.nddata import StdDevUncertainty
    >>> uncertainty = StdDevUncertainty(np.sqrt(data)) # representing standard deviation
    >>> meta = {'object': 'fictional data.'}
    >>> from astropy.coordinates import SkyCoord
    >>> wcs = SkyCoord('00h42m44.3s', '+41d16m09s')
    >>> ndd = NDData(data, mask=mask, unit=unit, uncertainty=uncertainty,
    ...              meta=meta, wcs=wcs)
    >>> ndd
    NDData([1, 2, 3, 4])

The representation only displays the ``data``; the other attributes need to be
accessed directly, for example ``ndd.mask`` to access the mask.


NDDataRef
---------

Building upon this pure container `~astropy.nddata.NDDataRef` implements:

+ a ``read`` and ``write`` method to access astropy's unified file io interface.
+ simple arithmetics like addition, subtraction, division and multiplication.
+ slicing.

Instances are created in the same way::

    >>> from astropy.nddata import NDDataRef
    >>> ndd = NDDataRef(ndd)
    >>> ndd
    NDDataRef([1, 2, 3, 4])

But also support arithmetic (:ref:`nddata_arithmetic`) like addition::

    >>> import astropy.units as u
    >>> ndd2 = ndd.add([4, 3, 2, 1] * u.erg / u.s)
    >>> ndd2
    NDDataRef([ 5.,  5.,  5.,  5.])

Because these operations have a wide range of options these are not available
using arithmetic operators like ``+``.

Slicing or indexing (:ref:`nddata_slicing`) is possible (issuing warnings if
some attribute cannot be sliced)::

    >>> ndd2[2:]  # discard the first two elements
    INFO: wcs cannot be sliced. [astropy.nddata.mixins.ndslicing]
    NDDataRef([ 5.,  5.])
    >>> ndd2[1]   # get the second element
    INFO: wcs cannot be sliced. [astropy.nddata.mixins.ndslicing]
    NDDataRef(5.0)


StdDevUncertainty
-----------------

`~astropy.nddata.StdDevUncertainty` implements uncertainty based on standard
deviation and can propagate these using the arithmetic methods of
`~astropy.nddata.NDDataRef`::

    >>> from astropy.nddata import NDDataRef, StdDevUncertainty
    >>> import numpy as np

    >>> uncertainty = StdDevUncertainty(np.arange(5))
    >>> ndd = NDDataRef([5,5,5,5,5], uncertainty=uncertainty)

    >>> doubled_ndd = ndd.multiply(2)  # multiply by 2
    >>> doubled_ndd.uncertainty
    StdDevUncertainty([0, 2, 4, 6, 8])

    >>> ndd2 = ndd.add(doubled_ndd)    # add the doubled to the original
    >>> ndd2.uncertainty
    StdDevUncertainty([ 0.        ,  2.23606798,  4.47213595,  6.70820393,
                        8.94427191])

    >>> # or taking into account that both of these uncertainties are correlated
    >>> ndd3 = ndd.add(doubled_ndd, uncertainty_correlation=1)
    >>> ndd3.uncertainty
    StdDevUncertainty([  0.,   3.,   6.,   9.,  12.])

.. note::
    The "amount" of correlation must be given, so ``1`` means correlated, ``-1``
    anti-correlated and ``0`` (default) uncorrelated. See also
    :ref:`nddata_arithmetic` for more information about correlation handling.

Using ``nddata``
================

.. toctree::
   :maxdepth: 2

   nddata.rst
   decorator.rst
   mixins/index.rst
   subclassing.rst
   utils.rst

Reference/API
=============

.. automodapi:: astropy.nddata
    :no-inheritance-diagram:

.. automodapi:: astropy.nddata.utils
    :no-inheritance-diagram:

.. _APE 7: https://github.com/astropy/astropy-APEs/blob/master/APE7.rst
