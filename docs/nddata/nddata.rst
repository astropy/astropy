.. _nddata_details:

NDData
======

Overview
--------

The `~astropy.nddata.NDData` class is a container for gridded N-dimensional
data. It has a ``data`` attribute, which can be any object that presents an
array-like interface, and optional attributes:

+  ``meta``, for metadata
+ ``unit`` for the ``data`` unit
+ ``uncertainty`` for the uncertainty of the data (which could be standard
  deviation, variance or something else),
+ ``mask`` for the ``data``
+ ``wcs``, representing the relationship  between ``data`` and world
  coordinates.

Of these, only ``mask`` and  ``uncertainty`` may be changed after the NDData
object is created.

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

It can also be instantiated by passing in an object that can be converted to a
numpy numerical array::

    >>> ndd3 = NDData([1, 2, 3, 4])

The final way to instantiate an `~astropy.nddata.NDData` object is with a data
object that presents a numpy array-like interface. If the object passed to the
intializer has all three of the attributes ``shape``, ``__getitem__`` (so it
is indexable) and ``__array__`` (so that it can act like a numpy array in
expression) then the ``data`` attribute will be set to that object.

The purpose of this mechanism is to allow considerable flexibility in the
objects used to store the data while providing a useful to default (numpy
array).

Mask
----

Values can be masked using the ``mask`` attribute.  One straightforward way to
provide a mask is to use a boolean numpy array::

     >>> ndd_masked = NDData(ndd, mask = ndd.data > 0.9)

Another is to simply initialize an `~astropy.nddata.NDData` object  with a
masked numpy array::

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
when creating the ``NDData`` object::

    >>> import astropy.units as u
    >>> ndd_unit = NDData([1, 2, 3, 4], unit="meter")
    >>> ndd_unit.unit
    Unit("m")

or by initializing with data that is an astropy `~astropy.units.Quantity`::

    >>> q = [1, 2, 3, 4] * u.meter
    >>> ndd_unit2 = NDData(q)
    >>> ndd_unit2.unit
    Unit("m")

Uncertainties
-------------

`~astropy.nddata.NDData` objects have an ``uncertainty`` attribute that can be
used to set the uncertainty on the data values. The ``uncertainty`` must have
an attribute ``uncertainty_type`` otherwise it is wrapped inside an
`~astropy.nddata.UnknownUncertainty`.

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

Initialization with copy
------------------------

The default way to create an `~astropy.nddata.NDData` instance is to try saving
the parameters as reference rather than as copy. Sometimes this is not possible
because the internal mechanics doesn't allow for this. For example if the
``data`` is a `list` then during initialization this is copied while converting
to a `~numpy.ndarray`. But it is also possible to enforce copies
during initialization by setting the ``copy`` parameter to ``True``::

    >>> array = np.array([1, 2, 3, 4])
    >>> ndd = NDData(array)
    >>> ndd.data[2] = 10
    >>> array[2]  # Original array is changed
    10
    >>> ndd2 = NDData(array, copy=True)
    >>> ndd2.data[2] = 3
    >>> array[2]  # Original array is not changed.
    10

Conflicting parameters
----------------------

It is possible to initialize `~astropy.nddata.NDData` with an explicit
parameter in the initialization call and an implicit one in the ``data``
parameter. If such a combination is detected a warning is
issued and the explicit parameter is kept (without any attempt of conversion).

For example with quantities::

    >>> import astropy.units as u
    >>> quantity = np.array([1, 2, 3, 4]) * u.m
    >>> ndd = NDData(quantity, unit="cm")  # Explicit unit is cm and implicit is m, keep explicit.
    INFO: Overwriting Quantity's current unit with specified unit [astropy.nddata.nddata]
    >>> ndd.unit
    Unit("cm")

But this can affect the mask if the ``data`` parameter is a `~numpy.ma.MaskedArray`
or even all attributes if the ``data`` is another `~astropy.nddata.NDData`
instance.

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
