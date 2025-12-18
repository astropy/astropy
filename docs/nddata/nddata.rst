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
+ ``psf`` holds an image representation of the point spread function (PSF)

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
    array([0.01, 0.01, 0.01])
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

    >>> from astropy.nddata import StdDevUncertainty, InverseVariance
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

It is also possible to convert between uncertainty types::

    >>> uncert.represent_as(InverseVariance)
    InverseVariance([0.1       , 0.14285714, 0.08333333, 0.04545455])

..
  EXAMPLE END

Covariance
----------

A `~astropy.nddata.Covariance` uncertainty type is also implemented; however,
its functionality is generally limited to construction and storage of sparse
covariance matrices.  Additional functionality will be implemented as requested.
See :ref:`nddata-covariance` for more description and example usage.

WCS
===

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
    np.int64(10)

    >>> ndd2 = NDData(array, copy=True)
    >>> ndd2.data[2] = 3
    >>> array[2]  # Original array hasn't changed.
    np.int64(10)

.. note::
    In some cases setting ``copy=True`` will copy the ``data`` twice. Known
    cases are if the ``data`` is a `list` or `tuple`.

..
  EXAMPLE END


Collapsing an NDData object along one or more axes
==================================================

..
  EXAMPLE START
  Collapsing an NDData object along one or more axes

A common operation on an `~numpy.ndarray` is to take the sum, mean,
maximum, or minimum along one or more axes, reducing the dimensions
of the output. These four operations are implemented on
`~astropy.nddata.NDData` with appropriate propagation of uncertainties,
masks, and units.

For example, let's work on the following ``data`` with a mask, unit, and
(uniform) uncertainty::

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from astropy.nddata import NDDataArray, StdDevUncertainty
    >>>
    >>> data = [
    ...     [1, 2, 3],
    ...     [2, 3, 4]
    ... ]
    >>> mask = [
    ...     [True, False, False],
    ...     [False, False, False]
    ... ]
    >>> uncertainty = StdDevUncertainty(np.ones_like(data))
    >>> nddata = NDDataArray(data=data, uncertainty=uncertainty, mask=mask, unit='m')

The sum along axis ``1`` gives one result per row::

    >>> sum_axis_1 = nddata.sum(axis=1)  # this is a new NDDataArray
    >>> print(np.asanyarray(sum_axis_1))  # this converts data to a numpy masked array. doctest: +FLOAT_CMP
    [-- 9.0]
    >>> print(sum_axis_1.uncertainty)  # doctest: +FLOAT_CMP
    StdDevUncertainty([1.41421356, 1.73205081])

The result has one masked value derived from the logical OR of the original mask
along ``axis=1``. The uncertainties are the square-root of the sum of the squares
of the input uncertainties. Since the original uncertainties were all unity, the
result is the square root of the number of unmasked data entries,
:math:`[\sqrt{2},\,\sqrt{3}]`.

We can similarly take the mean along ``axis=1``::

    >>> mean_axis_1 = nddata.mean(axis=1)
    >>> print(np.asanyarray(mean_axis_1))  # doctest: +FLOAT_CMP
    [2.5 3.0]
    >>> print(mean_axis_1.uncertainty)  # doctest: +FLOAT_CMP
    StdDevUncertainty([0.70710678, 0.57735027])

The result is the mean of the values where ``mask==False``, and in this example,
the result would only have ``mask==True`` if an entire row was masked. Since the
uncertainties were given as `~astropy.nddata.StdDevUncertainty`, the propagated
uncertainties decrease proportional to the number of unmasked measurements in each
row, following :math:`[2^{-1/2},\,3^{-1/2}]`.

There's no single, correct way of defining the uncertainties associated
with the ``min`` or ``max`` of a set of measurements, so
`~astropy.nddata.NDData` resists the temptation to guess, and returns
the minimum data value along the axis/axes, and the propagated mask, but
no uncertainties::

    >>> min_axis_1 = nddata.min(axis=1)
    >>> print(np.asanyarray(min_axis_1))  # doctest: +FLOAT_CMP
    [2.0 2.0]
    >>> print(min_axis_1.uncertainty)
    None

For some use cases, it may be helpful to return the uncertainty
at the same index as the minimum/maximum ``data`` value, so that
the original ``data`` retains its uncertainty. You can get this
behavior with::

    >>> min_axis_1 = nddata.min(axis=1, propagate_uncertainties=True)

    >>> print(np.asanyarray(min_axis_1))  # doctest: +FLOAT_CMP
    [2.0 2.0]
    >>> print(min_axis_1.uncertainty)  # doctest: +FLOAT_CMP
    StdDevUncertainty([1, 1])

Finally, in some cases it may be useful to do perform a collapse
operation only on the unmasked values, and only return a masked
result when all of the input values are masked. If we refer back to
the first example in this section, we see that the underlying
``data`` attribute has been summed over all values, including
masked ones::

    >>> sum_axis_1  # doctest: +FLOAT_CMP
    NDDataArray([——, 9.], unit='m')

where the first data element is masked. We can instead get the sum
for only unmasked values with the ``operation_ignores_mask`` option::

    >>> nddata.sum(axis=1, operation_ignores_mask=True)
    NDDataArray([5, 9], unit='m')

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

MaskedQuantity
--------------

Converting the ``data``, ``unit``, and ``mask`` to a ``MaskedQuantity``::

    >>> from astropy.utils.masked import Masked
    >>> Masked(u.Quantity(ndd.data, ndd.unit), ndd.mask)  # doctest: +FLOAT_CMP
    <MaskedQuantity [——, 2., 3., ——] m>
