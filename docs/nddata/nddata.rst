.. _nddata_details:

NDData
======

Overview
--------

The `~astropy.nddata.NDData` class is a container for gridded N-dimensional
data. It has a ``data`` attribute, which can be any object that presents a
`~numpy.ndarray`-like interface, and optional attributes:

+ ``meta``, for metadata.
+ ``unit``, for the ``data`` unit.
+ ``uncertainty``, for the uncertainty in the data (which could be standard
  deviation, variance or something else).
+ ``mask`` for the ``data``.
+ ``wcs``, representing the relationship  between ``data`` and world
  coordinates.

Of these, only ``mask`` and  ``uncertainty`` may be changed after the NDData
object is created. Altough it *is* possible to alter attributes if they
are mutable by altering the attributes return.

Initializing
------------

An `~astropy.nddata.NDData` object can be instantiated by passing it an
n-dimensional Numpy array::

    >>> import numpy as np
    >>> from astropy.nddata import NDData
    >>> array = np.zeros((12, 12, 12))  # a 3-dimensional array with all zeros
    >>> ndd1 = NDData(array)

Note that the data in ``ndd`` is a reference to the original ``array``, so
changing the data in ``ndd`` will change the corresponding data in ``array``
in most circumstances. If you like to save a copy instead of the reference you
can specify the ``copy`` parameter to be `True`::

    >>> ndd2 = NDData(array, copy=True)

An `~astropy.nddata.NDData` object can also be instantiated by passing it an
`~astropy.nddata.NDData` object:

    >>> ndd1 = NDData(array)
    >>> ndd3 = NDData(ndd1)

As above, the data in``ndd2`` is a reference to the data in ``ndd1``, so
changes to one will affect the other if copy was not `True`.

It can also be instantiated by passing in an object that can be converted to a
numpy numerical array::

    >>> ndd4 = NDData([1, 2, 3, 4])

This way the data is always copied since converting a `list` to an
`~numpy.ndarray` enforces a copy.

Also certain subclasses of `~numpy.ndarray` can be used as input like
`~numpy.ma.MaskedArray` or `~astropy.units.Quantity` for these cases check
the following subsections.

The final way to instantiate an `~astropy.nddata.NDData` object is with a data
object that presents a numpy array-like interface. If the object passed to the
intializer has all four of the attributes ``shape``, ``__getitem__`` (so it
is indexable), ``dtype`` and ``__array__`` (so that it can act like a numpy
array in expression) then the ``data`` attribute will be set to that object.

The purpose of this mechanism is to allow considerable flexibility in the
objects used to store the data while providing a useful to default (numpy
array).

Mask
----

Values can be masked using the ``mask`` attribute.  One straightforward way to
provide a mask is to use a boolean numpy array::

    >>> ndd_masked = NDData(ndd4, mask = ndd4.data > 1.5)

Another is to simply initialize an `~astropy.nddata.NDData` object  with a
masked numpy array (`~numpy.ma.MaskedArray`)::

    >>> masked_array = np.ma.array([1, 2, 3, 4], mask=[1, 0, 0, 1])
    >>> ndd_masked = NDData(masked_array)
    >>> ndd_masked.mask
    array([ True, False, False,  True], dtype=bool)

A mask value of `True` indicates a value that should be ignored, while a mask
value of `False` indicates a valid value.

There is no requirement that the mask actually be a numpy array; for example,
a function which evaluates a mask value as needed is acceptable as long as it
follows the convention that `True` indicates a value that should be ignored.

Unit
----

The unit of the data can be set by either explicitly providing an astropy unit
(as `~astropy.units.Unit` or as `str`) when creating the ``NDData`` object::

    >>> import astropy.units as u
    >>> ndd_unit = NDData([1, 2, 3, 4], unit=u.m)
    >>> ndd_unit.unit
    Unit("m")
    >>> ndd_unit2 = NDData([1, 2, 3, 4], unit="meter")
    >>> ndd_unit2.unit
    Unit("m")

or by initializing with data that is an astropy `~astropy.units.Quantity`::

    >>> q = [1, 2, 3, 4] * u.meter
    >>> ndd_unit3 = NDData(q)
    >>> ndd_unit3.unit
    Unit("m")

Uncertainties
-------------

`~astropy.nddata.NDData` objects have an ``uncertainty`` attribute that can be
used to set the uncertainty on the data values. The ``uncertainty`` should have
an attribute ``uncertainty_type`` which is a string.

While not a requirement, the following ``uncertainty_type`` strings
are strongly recommended for common ways of specifying normal
distributions:

+ ``"std"``: if ``uncertainty`` stores the standard deviation/sigma
  (either a single value or on a per-pixel basis).
+ ``"var"``: if ``uncertainty`` stores the variance (either a single
  value or on a per-pixel basis).
+ ``"ivar"``: if ``uncertainty`` stores the inverse variance (either a
  single value or on a per-pixel basis).

.. note:: For information on creating your own uncertainty classes,
          see :doc:`subclassing`.

Meta-data
---------

The :class:`~astropy.nddata.NDData` class includes a ``meta`` attribute that
defaults to an empty ordered dictionary, and can be used to set overall meta-
data for the dataset::

    ndd.meta['exposure_time'] = 340.
    ndd.meta['filter'] = 'J'

Elements of the meta-data dictionary can be set to any valid Python object::

    ndd.meta['history'] = ['calibrated', 'aligned', 'flat-fielded']

The metadata can be any python object that presents a dict-like interface. For
example, a FITS header can be used as the metadata::

    >>> from astropy.io import fits
    >>> header = fits.Header()
    >>> header['observer'] = 'Edwin Hubble'
    >>> ndd = NDData(np.zeros([10, 10]), meta=header)
    >>> ndd.meta['observer']
    'Edwin Hubble'

WCS
---

At the moment the ``wcs`` attribute can be set to any object, though in the
future it may be restricted to an `~astropy.wcs.WCS` object once a generalized
WCS object is developed.

Alteration of Attributs
-----------------------

Most attributes (except ``mask`` and ``uncertainty``) do not allow for
overriding attributes as a whole but since most attributes (except ``unit``)
can be stored as mutable types you may alter what these attributes contain.

To change the data you may::

    >>> ndd1 = NDData([1, 2, 3, 4])
    >>> ndd1.data = [4, 3, 2, 1]        # doctest: +SKIP
    AttributeError ...
    ...
    AttributeError: can't set attribute
    >>> ndd1.data[:] = [4, 3, 2, 1]
    >>> ndd1.data
    array([4, 3, 2, 1])

similar operations work also on the ``meta`` (just use `dict` operations there)
and for the ``wcs`` depending on the type of object you stored there.

Converting to Numpy arrays
--------------------------

Data should be accessed through the ``data`` attribute::

    >>> array = np.asarray(ndd.data)

Though using ``np.asarray`` is not required it will ensure that an additional
copy of the data is not made if the data is a numpy array.

Note that if the data is masked you must explicitly construct a numpy masked
array like this::

    >>> input_array = np.ma.array([1, 2, 3, 4], mask=[1, 0, 0, 1])
    >>> ndd_masked = NDData(input_array)
    >>> masked_array = np.ma.array(ndd_masked.data, mask=ndd_masked.mask)

The same applies to extracting the data as Quantity::

    >>> ndd_quantity = NDData([1, 2, 3, 4], unit="meter")
    >>> quantity = u.Quantity(ndd_quantity.data, unit=ndd_quantity.unit)
