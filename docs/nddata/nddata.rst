.. _nddata_details:

NDData
******

Overview
========

:class:`~astropy.nddata.NDData` is based on `numpy.ndarray`-like ``data`` with
additional meta attributes:

+  ``meta`` for general metadata
+ ``unit`` represents the physical unit of the data
+ ``uncertainty`` for the uncertainty of the data
+ ``mask`` indicates invalid points in the data
+ ``wcs`` represents the relationship between the data grid and world
  coordinates

Each of these attributes can be set during initialization or directly on the
instance. Only the ``data`` cannot be directly set after creating the instance.

Data
====

The data is the base of `~astropy.nddata.NDData` and is required to be
`numpy.ndarray`-like. It is the only property that is required to create an
instance and it cannot be directly set on the instance.

Example
-------

..
  EXAMPLE START
  Creating Instances with NumPy NDarray-like Data

To create an instance::

    >>> import numpy as np
    >>> from astropy.nddata import NDData
    >>> array = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
    >>> ndd = NDData(array)
    >>> ndd
    NDData([[0, 1, 0],
            [1, 0, 1],
            [0, 1, 0]])

And access by the ``data`` attribute::

    >>> ndd.data
    array([[0, 1, 0],
           [1, 0, 1],
           [0, 1, 0]])

As already mentioned, it is not possible to set the data directly. So
``ndd.data = np.arange(9)`` will raise an exception. But the data can be
modified in place::

    >>> ndd.data[1,1] = 100
    >>> ndd.data
    array([[  0,   1,   0],
           [  1, 100,   1],
           [  0,   1,   0]])

..
  EXAMPLE END

Data During Initialization
--------------------------

During initialization it is possible to provide data that is not a
`numpy.ndarray` but convertible to one.

Examples
^^^^^^^^

..
  EXAMPLE START
  Data Convertible to a NumPy NDarray During Initialization

To provide data that is convertible to a `numpy.ndarray`, you can pass a `list`
containing numerical values::

    >>> alist = [1, 2, 3, 4]
    >>> ndd = NDData(alist)
    >>> ndd.data  # data will be a numpy-array:
    array([1, 2, 3, 4])

A nested `list` or `tuple` is possible, but if these contain non-numerical
values the conversion might fail.

Besides input that is convertible to such an array, you can also use the
``data`` parameter to pass implicit additional information. For example, if the
data is another `~astropy.nddata.NDData` object it implicitly uses its
properties::

    >>> ndd = NDData(ndd, unit = 'm')
    >>> ndd2 = NDData(ndd)
    >>> ndd2.data  # It has the same data as ndd
    array([1, 2, 3, 4])
    >>> ndd2.unit  # but it also has the same unit as ndd
    Unit("m")

Another possibility is to use a `~astropy.units.Quantity` as a ``data``
parameter::

    >>> import astropy.units as u
    >>> quantity = np.ones(3) * u.cm  # this will create a Quantity
    >>> ndd3 = NDData(quantity)
    >>> ndd3.data  # doctest: +FLOAT_CMP
    array([1., 1., 1.])
    >>> ndd3.unit
    Unit("cm")

Or a `numpy.ma.MaskedArray`::

    >>> masked_array = np.ma.array([5,10,15], mask=[False, True, False])
    >>> ndd4 = NDData(masked_array)
    >>> ndd4.data
    array([ 5, 10, 15])
    >>> ndd4.mask
    array([False,  True, False]...)

If such an implicitly passed property conflicts with an explicit parameter, the
explicit parameter will be used and an info message will be issued::

    >>> quantity = np.ones(3) * u.cm
    >>> ndd6 = NDData(quantity, unit='m')
    INFO: overwriting Quantity's current unit with specified unit. [astropy.nddata.nddata]
    >>> ndd6.data  # doctest: +FLOAT_CMP
    array([1., 1., 1.])
    >>> ndd6.unit
    Unit("m")

The unit of the `~astropy.units.Quantity` is being ignored and the unit is set
to the explicitly passed one.

It might be possible to pass other classes as a ``data`` parameter as long as
they have the properties ``shape``, ``dtype``, ``__getitem__``, and
``__array__``.

The purpose of this mechanism is to allow considerable flexibility in the
objects used to store the data while providing a useful default (``numpy``
array).

..
  EXAMPLE END

Mask
====

The ``mask`` is being used to indicate if data points are valid or invalid.
`~astropy.nddata.NDData` does not restrict this mask in any way but it is
expected to follow the `numpy.ma.MaskedArray` convention in that the mask:

+ Returns ``True`` for data points that are considered **invalid**.
+ Returns ``False`` for those points that are **valid**.

Examples
--------

..
  EXAMPLE START
  Masks Used to Indicate Valid or Invalid Data Points in NDData

One possibility is to create a mask by using ``numpy``'s comparison operators::

    >>> array = np.array([0, 1, 4, 0, 2])

    >>> mask = array == 0  # Mask points containing 0
    >>> mask
    array([ True, False, False,  True, False]...)

    >>> other_mask = array > 1  # Mask points with a value greater than 1
    >>> other_mask
    array([False, False,  True, False,  True]...)

And initialize the `~astropy.nddata.NDData` instance using the ``mask``
parameter::

    >>> ndd = NDData(array, mask=mask)
    >>> ndd.mask
    array([ True, False, False,  True, False]...)

Or by replacing the mask::

    >>> ndd.mask = other_mask
    >>> ndd.mask
    array([False, False,  True, False,  True]...)

There is no requirement that the mask actually be a ``numpy`` array; for
example, a function which evaluates a mask value as needed is acceptable as
long as it follows the convention that ``True`` indicates a value that should
be ignored.

..
  EXAMPLE END

Unit
====

The ``unit`` represents the unit of the data values. It is required to be
`~astropy.units.Unit`-like or a string that can be converted to such a
`~astropy.units.Unit`::

    >>> import astropy.units as u
    >>> ndd = NDData([1, 2, 3, 4], unit="meter")  # using a string
    >>> ndd.unit
    Unit("m")

..note::
    Setting the ``unit`` on an instance is not possible.

Uncertainties
=============

The ``uncertainty`` represents an arbitrary representation of the error of the
data values. To indicate which kind of uncertainty representation is used, the
``uncertainty`` should have an ``uncertainty_type`` property. If no such
property is found it will be wrapped inside a
`~astropy.nddata.UnknownUncertainty`.

The ``uncertainty_type`` should follow the `~astropy.nddata.StdDevUncertainty`
convention in that it returns a short string like ``"std"`` for an uncertainty
given in standard deviation. Other examples are
`~astropy.nddata.VarianceUncertainty` and `~astropy.nddata.InverseVariance`.

Examples
--------

..
  EXAMPLE START
  Setting Uncertainties During Initialization in NDData

Like the other properties the ``uncertainty`` can be set during
initialization::

    >>> from astropy.nddata import StdDevUncertainty
    >>> array = np.array([10, 7, 12, 22])
    >>> uncert = StdDevUncertainty(np.sqrt(array))
    >>> ndd = NDData(array, uncertainty=uncert)
    >>> ndd.uncertainty  # doctest: +FLOAT_CMP
    StdDevUncertainty([3.16227766, 2.64575131, 3.46410162, 4.69041576])

Or on the instance directly::

    >>> other_uncert = StdDevUncertainty([2,2,2,2])
    >>> ndd.uncertainty = other_uncert
    >>> ndd.uncertainty
    StdDevUncertainty([2, 2, 2, 2])

But it will print an info message if there is no ``uncertainty_type``::

    >>> ndd.uncertainty = np.array([5, 1, 2, 10])
    INFO: uncertainty should have attribute uncertainty_type. [astropy.nddata.nddata]
    >>> ndd.uncertainty
    UnknownUncertainty([ 5,  1,  2, 10])

..
  EXAMPLE END

WCS
---

The ``wcs`` should contain a mapping from the gridded data to world
coordinates. There are no restrictions placed on the property currently but it
may be restricted to an `~astropy.wcs.WCS` object or a more generalized WCS
object in the future.

.. note::
    Like the unit the ``wcs`` cannot be set on an instance.

Metadata
=========

The ``meta`` property contains all further meta information that does not fit
any other property.

Examples
--------

..
  EXAMPLE START
  Metadata in NDData

If the ``meta`` property is given it must be `dict`-like::

    >>> ndd = NDData([1,2,3], meta={'observer': 'myself'})
    >>> ndd.meta
    {'observer': 'myself'}

`dict`-like means it must be a mapping from some keys to some values. This
also includes `~astropy.io.fits.Header` objects::

    >>> from astropy.io import fits
    >>> header = fits.Header()
    >>> header['observer'] = 'Edwin Hubble'
    >>> ndd = NDData(np.zeros([10, 10]), meta=header)
    >>> ndd.meta['observer']
    'Edwin Hubble'

If the ``meta`` property is not provided or explicitly set to ``None``, it will
default to an empty `collections.OrderedDict`::

    >>> ndd.meta = None
    >>> ndd.meta
    OrderedDict()

    >>> ndd = NDData([1,2,3])
    >>> ndd.meta
    OrderedDict()

The ``meta`` object therefore supports adding or updating these values::

    >>> ndd.meta['exposure_time'] = 340.
    >>> ndd.meta['filter'] = 'J'

Elements of the metadata dictionary can be set to any valid Python object::

    >>> ndd.meta['history'] = ['calibrated', 'aligned', 'flat-fielded']

..
  EXAMPLE END

Initialization with Copy
========================

The default way to create an `~astropy.nddata.NDData` instance is to try saving
the parameters as references to the original rather than as copy. Sometimes
this is not possible because the internal mechanics do not allow for this.

Examples
--------

..
  EXAMPLE START
  Creating an NDData Instance with Copy

If the ``data`` is a `list` then during initialization this is copied
while converting to a `~numpy.ndarray`. But it is also possible to enforce
copies during initialization by setting the ``copy`` parameter to ``True``::

    >>> array = np.array([1, 2, 3, 4])
    >>> ndd = NDData(array)
    >>> ndd.data[2] = 10
    >>> array[2]  # Original array has changed
    10

    >>> ndd2 = NDData(array, copy=True)
    >>> ndd2.data[2] = 3
    >>> array[2]  # Original array hasn't changed.
    10

.. note::
    In some cases setting ``copy=True`` will copy the ``data`` twice. Known
    cases are if the ``data`` is a `list` or `tuple`.

..
  EXAMPLE END

Converting NDData to Other Classes
==================================

There is limited support to convert a `~astropy.nddata.NDData` instance to
other classes. In the process some properties might be lost.

    >>> data = np.array([1, 2, 3, 4])
    >>> mask = np.array([True, False, False, True])
    >>> unit = 'm'
    >>> ndd = NDData(data, mask=mask, unit=unit)

`numpy.ndarray`
---------------

Converting the ``data`` to an array::

    >>> array = np.asarray(ndd.data)
    >>> array
    array([1, 2, 3, 4])

Though using ``np.asarray`` is not required, in most cases it will ensure that
the result is always a `numpy.ndarray`

`numpy.ma.MaskedArray`
----------------------

Converting the ``data`` and ``mask`` to a MaskedArray::


    >>> masked_array = np.ma.array(ndd.data, mask=ndd.mask)
    >>> masked_array
    masked_array(data=[--, 2, 3, --],
                 mask=[ True, False, False,  True],
           fill_value=999999)

`~astropy.units.Quantity`
-------------------------

Converting the ``data`` and ``unit`` to a Quantity::

    >>> quantity = u.Quantity(ndd.data, unit=ndd.unit)
    >>> quantity  # doctest: +FLOAT_CMP
    <Quantity [1., 2., 3., 4.] m>

.. note::
    Ideally, you would construct masked quantities, but these are not properly
    supported: many operations on them fail.
