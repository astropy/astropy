NDData overview
===============

Initializing
------------

An `~astropy.nddata.NDData` object can be instantiated by passing it an
n-dimensional Numpy array::

    >>> import numpy as np
    >>> from astropy.nddata import NDData
    >>> array = np.zeros((12, 12, 12))  # a 3-dimensional array with all zeros
    >>> ndd = NDData(array)

Note that the data in ``ndd`` is a reference to the original ``array``, so
changing the data in ``ndd`` will change the corresponding data in ``array``
in most circumstances.

An `~astropy.nddata.NDData` object can also be instantiated by passing it an
`~astropy.nddata.NDData` object:

    >>> ndd1 = NDData(array)
    >>> ndd2 = NDData(ndd1)

As above, the data in``ndd2`` is a reference to the data in ``ndd1``, so
changes to one will affect the other.

This object has a few attributes in common with Numpy:

    >>> ndd.ndim
    3
    >>> ndd.shape
    (12, 12, 12)
    >>> ndd.dtype
    dtype('float64')

The underlying Numpy array can be accessed via the ``data`` attribute::

    >>> ndd.data
    array([[[ 0.,  0.,  0., ...

Mask
----

Values can be masked using the ``mask`` attribute, which should be a boolean
Numpy array with the same dimensions as the data, e.g.::

     >>> ndd.mask = ndd.data > 0.9

A mask value of `True` indicates a value that should be ignored, while a mask
value of `False` indicates a valid value.

Flags
-----

Values can be assigned one or more flags. The ``flags`` attribute is used to
store either a single Numpy array (of any type) with dimensions matching that
of the data, or a `~astropy.nddata.FlagCollection`, which is
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

Uncertainties
-------------

`~astropy.nddata.NDData` objects have an ``uncertainty`` attribute that can be
used to set the uncertainty on the data values. This is done by using classes
to represent the uncertainties of a given type. For example, to set standard
deviation uncertainties on the pixel values, you can do::

    >>> from astropy.nddata import StdDevUncertainty
    >>> ndd.uncertainty = StdDevUncertainty(np.ones((12, 12, 12)) * 0.1)

.. note:: For information on creating your own uncertainty classes,
          see :doc:`subclassing`.

Arithmetic
----------

Provided that the world coordinate system (WCS) and shape match, and that the
units are consistent, two :class:`~astropy.nddata.NDData` instances can be
added, subtracted, multiplied or divided from each other, with uncertainty
propagation, creating a new :class:`~astropy.nddata.NDData` object::

    ndd3 = ndd1.add(ndd2)
    ndd4 = ndd1.subtract(ndd2)
    ndd5 = ndd1.multiply(ndd2)
    ndd6 = ndd1.divide(ndd2)

The purpose of the :meth:`~astropy.nddata.nddata.NDData.add`,
:meth:`~astropy.nddata.nddata.NDData.subtract`,
:meth:`~astropy.nddata.nddata.NDData.multiply` and
:meth:`~astropy.nddata.nddata.NDData.divide` methods is to allow the
combination of two data objects that have common WCS and shape and units
consistent with the operation performed, with consistent behavior for masks,
and with a framework to propagate uncertainties. Currently any flags on the
operands are dropped so that the result of the operation always has no flags.
These methods are intended for use by sub-classes and functions that deal with
more complex combinations.

Entries that are masked in either of the operands are also masked in the
result.

.. warning:: Uncertainty propagation is still experimental, and does not take
             into account correlated uncertainties.

Meta-data
---------

The :class:`~astropy.nddata.NDData` class includes a ``meta`` attribute
that defaults to an empty dictionary, and can be used to set overall meta-data
for the dataset::

    ndd.meta['exposure_time'] = 340.
    ndd.meta['filter'] = 'J'

Elements of the meta-data dictionary can be set to any valid Python object::

    ndd.meta['history'] = ['calibrated', 'aligned', 'flat-fielded']

Converting to Numpy arrays
--------------------------

`~astropy.nddata.NDData` objects can also be easily converted to
numpy arrays::

    >>> import numpy as np
    >>> arr = np.array(ndd)
    >>> np.all(arr == mydataarray)  # doctest: +SKIP
    True

If a ``mask`` is defined, this will result in a `~numpy.ma.MaskedArray`, so
in all cases a useable `numpy.ndarray` or subclass will result. This allows
straightforward plotting of `~astropy.nddata.NDData` objects with 1-
and 2-dimensional datasets using Matplotlib::

    >>> from matplotlib import pyplot as plt  # doctest: +SKIP
    >>> plt.plot(ndd)  # doctest: +SKIP

This works because the Matplotlib plotting functions automatically convert
their inputs using `numpy.array`.
