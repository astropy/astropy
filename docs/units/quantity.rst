.. _quantity:

Quantity
********

The |Quantity| object is meant to represent a value that has some unit
associated with the number.

Creating Quantity Instances
===========================

|Quantity| objects are normally created through multiplication with
:class:`~astropy.units.Unit` objects.

Examples
--------

.. EXAMPLE START: Creating Quantity Instances Through Multiplication

To create a |Quantity| to represent 15 m/s:

    >>> import astropy.units as u
    >>> 15 * u.m / u.s  # doctest: +FLOAT_CMP
    <Quantity 15. m / s>

This extends as expected to division by a unit, or using ``numpy`` arrays or
`Python sequences <https://docs.python.org/3/library/stdtypes.html#typesseq>`_:

    >>> 1.25 / u.s
    <Quantity 1.25 1 / s>
    >>> [1, 2, 3] * u.m  # doctest: +FLOAT_CMP
    <Quantity [1., 2., 3.] m>
    >>> import numpy as np
    >>> np.array([1, 2, 3]) * u.m  # doctest: +FLOAT_CMP
    <Quantity [1., 2., 3.] m>

.. EXAMPLE END

.. EXAMPLE START: Creating Quantity Instances Using the Quantity Constructor

You can also create instances using the |Quantity| constructor directly, by
specifying a value and unit:

    >>> u.Quantity(15, u.m / u.s)  # doctest: +FLOAT_CMP
    <Quantity 15. m / s>

The constructor gives a few more options. In particular, it allows you to
merge sequences of |Quantity| objects (as long as all of their units are
equivalent), and to parse simple strings (which may help, for example, to parse
configuration files, etc.):

    >>> qlst = [60 * u.s, 1 * u.min]
    >>> u.Quantity(qlst, u.minute)  # doctest: +FLOAT_CMP
    <Quantity [1.,  1.] min>
    >>> u.Quantity('15 m/s')  # doctest: +FLOAT_CMP
    <Quantity 15. m / s>

The current unit and value can be accessed via the
`~astropy.units.quantity.Quantity.unit` and
`~astropy.units.quantity.Quantity.value` attributes:

    >>> q = 2.5 * u.m / u.s
    >>> q.unit
    Unit("m / s")
    >>> q.value
    2.5

.. note:: |Quantity| objects are converted to float by default. Furthermore, any
          data passed in are copied, which for large arrays may not be optimal.
          As discussed :ref:`further below <astropy-units-quantity-no-copy>`,
          you can instead obtain a `view
          <https://numpy.org/doc/stable/glossary.html#term-view>`_ by passing
          ``copy=False`` to |Quantity| or by using the ``<<`` operator.

.. EXAMPLE END

.. _quantity_unit_conversion:

Converting to Different Units
=============================

|Quantity| objects can be converted to different units using the
:meth:`~astropy.units.quantity.Quantity.to` method.

Examples
--------

.. EXAMPLE START: Converting Quantity Objects to Different Units

To convert |Quantity| objects to different units:

    >>> q = 2.3 * u.m / u.s
    >>> q.to(u.km / u.h)  # doctest: +FLOAT_CMP
    <Quantity 8.28 km / h>

For convenience, the :attr:`~astropy.units.quantity.Quantity.si` and
:attr:`~astropy.units.quantity.Quantity.cgs` attributes can be used to convert
the |Quantity| to base `SI
<https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-EN.pdf>`_ or `CGS
<https://en.wikipedia.org/wiki/Centimetre-gram-second_system_of_units>`_ units:

    >>> q = 2.4 * u.m / u.s
    >>> q.si  # doctest: +FLOAT_CMP
    <Quantity 2.4 m / s>
    >>> q.cgs  # doctest: +FLOAT_CMP
    <Quantity 240. cm / s>

If you want the value of the quantity in a different unit, you can use
:meth:`~astropy.units.Quantity.to_value` as a shortcut:

    >>> q = 2.5 * u.m
    >>> q.to_value(u.cm)
    250.0

.. note:: You could get the value in ``cm`` also by using ``q.to(u.cm).value``.
          The difference is that :meth:`~astropy.units.Quantity.to_value` does
          no copying if the unit is already the correct one, instead
          returning a `view
          <https://numpy.org/doc/stable/glossary.html#term-view>`_  of the data
          (just as if you had done ``q.value``). In contrast,
          :meth:`~astropy.units.Quantity.to` always returns a copy (which also
          means it is slower for the case where no conversion is necessary).
          As discussed :ref:`further below <astropy-units-quantity-no-copy>`,
          you can avoid the copying by using the ``<<`` operator.

Comparing Quantities
====================

The equality of |Quantity| objects is best tested using the
:func:`~astropy.units.allclose` and :func:`~astropy.units.isclose` functions,
which are unit-aware analogues of the ``numpy`` functions with the same name::

    >>> u.allclose([1, 2] * u.m, [100, 200] * u.cm)
    True
    >>> u.isclose([1, 2] * u.m, [100, 20] * u.cm)
    array([ True, False])

The use of `Python comparison operators
<https://docs.python.org/3/reference/expressions.html#comparisons>`_ is also
supported::

    >>> 1*u.m < 50*u.cm
    False

Plotting Quantities
===================

|Quantity| objects can be conveniently plotted using |Matplotlib| — see
:ref:`plotting-quantities` for more details.

.. _quantity_arithmetic:

Arithmetic
==========

Addition and Subtraction
------------------------

Addition or subtraction between |Quantity| objects is supported when their
units are equivalent.

Examples
^^^^^^^^

.. EXAMPLE START: Addition and Subtraction Between Quantity Objects

When the units are equal, the resulting object has the same unit:

    >>> 11 * u.s + 30 * u.s  # doctest: +FLOAT_CMP
    <Quantity 41. s>
    >>> 30 * u.s - 11 * u.s  # doctest: +FLOAT_CMP
    <Quantity 19. s>

If the units are equivalent, but not equal (e.g., kilometer and meter), the
resulting object **has units of the object on the left**:

    >>> 1100.1 * u.m + 13.5 * u.km
    <Quantity 14600.1 m>
    >>> 13.5 * u.km + 1100.1 * u.m  # doctest: +FLOAT_CMP
    <Quantity 14.6001 km>
    >>> 1100.1 * u.m - 13.5 * u.km
    <Quantity -12399.9 m>
    >>> 13.5 * u.km - 1100.1 * u.m  # doctest: +FLOAT_CMP
    <Quantity 12.3999 km>

Addition and subtraction are not supported between |Quantity| objects and basic
numeric types, except for dimensionless quantities (see `Dimensionless
Quantities`_) or special values like zero and infinity::

    >>> 13.5 * u.km + 19.412  # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
      ...
    UnitConversionError: Can only apply 'add' function to dimensionless
    quantities when other argument is not a quantity (unless the
    latter is all zero/infinity/nan)

.. EXAMPLE END

Multiplication and Division
---------------------------

Multiplication and division are supported between |Quantity| objects with any
units, and with numeric types. For these operations between objects with
equivalent units, the **resulting object has composite units**.

Examples
^^^^^^^^

.. EXAMPLE START: Multiplication and Division Between Quantity Objects

To perform these operations on |Quantity| objects:

    >>> 1.1 * u.m * 140.3 * u.cm  # doctest: +FLOAT_CMP
    <Quantity 154.33 cm m>
    >>> 140.3 * u.cm * 1.1 * u.m  # doctest: +FLOAT_CMP
    <Quantity 154.33 cm m>
    >>> 1. * u.m / (20. * u.cm)  # doctest: +FLOAT_CMP
    <Quantity 0.05 m / cm>
    >>> 20. * u.cm / (1. * u.m)  # doctest: +FLOAT_CMP
    <Quantity 20. cm / m>

For multiplication, you can change how to represent the resulting object by
using the :meth:`~astropy.units.quantity.Quantity.to` method:

    >>> (1.1 * u.m * 140.3 * u.cm).to(u.m**2)  # doctest: +FLOAT_CMP
    <Quantity 1.5433 m2>
    >>> (1.1 * u.m * 140.3 * u.cm).to(u.cm**2)  # doctest: +FLOAT_CMP
    <Quantity 15433. cm2>

For division, if the units are equivalent, you may want to make the resulting
object dimensionless by reducing the units. To do this, use the
:meth:`~astropy.units.quantity.Quantity.decompose()` method:

    >>> (20. * u.cm / (1. * u.m)).decompose()  # doctest: +FLOAT_CMP
    <Quantity 0.2>

This method is also useful for more complicated arithmetic:

    >>> 15. * u.kg * 32. * u.cm * 15 * u.m / (11. * u.s * 1914.15 * u.ms)  # doctest: +FLOAT_CMP
    <Quantity 0.34195097 cm kg m / (ms s)>
    >>> (15. * u.kg * 32. * u.cm * 15 * u.m / (11. * u.s * 1914.15 * u.ms)).decompose()  # doctest: +FLOAT_CMP
    <Quantity 3.41950973 m2 kg / s2>

.. EXAMPLE END

.. _quantity_and_numpy:

NumPy Functions
===============

|Quantity| objects are actually full ``numpy`` arrays (the |Quantity| class
inherits from and extends :class:`numpy.ndarray`), and we have tried to ensure
that ``numpy`` functions behave properly with quantities:

    >>> q = np.array([1., 2., 3., 4.]) * u.m / u.s
    >>> np.mean(q)
    <Quantity 2.5 m / s>
    >>> np.std(q)  # doctest: +FLOAT_CMP
    <Quantity 1.11803399 m / s>

This includes functions that only accept specific units such as angles:

    >>> q = 30. * u.deg
    >>> np.sin(q)  # doctest: +FLOAT_CMP
    <Quantity 0.5>

Or `Dimensionless Quantities`_::

    >>> from astropy.constants import h, k_B
    >>> nu = 3 * u.GHz
    >>> T = 30 * u.K
    >>> np.exp(-h * nu / (k_B * T))  # doctest: +FLOAT_CMP
    <Quantity 0.99521225>

.. note:: Support for functions from other packages, such as |SciPy|, is more
          incomplete (contributions to improve this are welcomed!).

Dimensionless Quantities
========================

Dimensionless quantities have the characteristic that if they are
added to or subtracted from a Python scalar or unitless `~numpy.ndarray`,
or if they are passed to a ``numpy`` function that takes dimensionless
quantities, the units are simplified so that the quantity is
dimensionless and scale-free. For example:

    >>> 1. + 1. * u.m / u.km  # doctest: +FLOAT_CMP
    <Quantity 1.001>

Which is different from:

    >>> 1. + (1. * u.m / u.km).value
    2.0

In the latter case, the result is ``2.0`` because the unit of ``(1. * u.m /
u.km)`` is not scale-free by default:

    >>> q = (1. * u.m / u.km)
    >>> q.unit
    Unit("m / km")
    >>> q.unit.decompose()
    Unit(dimensionless with a scale of 0.001)

However, when combining with an object that is not a |Quantity|, the unit is
automatically decomposed to be scale-free, giving the expected result.

This also occurs when passing dimensionless quantities to functions that take
dimensionless quantities:

    >>> nu = 3 * u.GHz
    >>> T = 30 * u.K
    >>> np.exp(- h * nu / (k_B * T))  # doctest: +FLOAT_CMP
    <Quantity 0.99521225>

The result is independent from the units in which the different quantities were
specified:

    >>> nu = 3.e9 * u.Hz
    >>> T = 30 * u.K
    >>> np.exp(- h * nu / (k_B * T))  # doctest: +FLOAT_CMP
    <Quantity 0.99521225>

Converting to Plain Python Scalars
==================================

Converting |Quantity| objects does not work for non-dimensionless quantities:

    >>> float(3. * u.m)
    Traceback (most recent call last):
      ...
    TypeError: only dimensionless scalar quantities can be converted
    to Python scalars

Only dimensionless values can be converted to plain Python scalars:

    >>> float(3. * u.m / (4. * u.m))
    0.75
    >>> float(3. * u.km / (4. * u.m))
    750.0
    >>> int(6. * u.km / (2. * u.m))
    3000

Functions that Accept Quantities
================================

If a function accepts a |Quantity| as an argument then it can be a good idea to
check that the provided |Quantity| belongs to one of the expected
:ref:`physical_types`. This can be done with the `decorator
<https://docs.python.org/3/glossary.html#term-decorator>`_
:func:`~astropy.units.quantity_input`.

The decorator does not convert the input |Quantity| to the desired unit, say
arcseconds to degrees in the example below, it merely checks that such a
conversion is possible, thus verifying that the `~astropy.units.Quantity`
argument can be used in calculations.

Keyword arguments to :func:`~astropy.units.quantity_input` specify which
arguments should be validated and what unit they are expected to be compatible
with.

Examples
--------

.. EXAMPLE START: Functions that Accept Quantities

To verify if a |Quantity| argument can be used in calculations::

    >>> @u.quantity_input(myarg=u.deg)
    ... def myfunction(myarg):
    ...     return myarg.unit

    >>> myfunction(100*u.arcsec)
    Unit("arcsec")
    >>> myfunction(2*u.m)  # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ...
    UnitsError: Argument 'myarg' to function 'myfunction' must be in units
    convertible to 'deg'.

It is also possible to instead specify the :ref:`physical type
<physical_types>` of the desired unit::

    >>> @u.quantity_input(myarg='angle')
    ... def myfunction(myarg):
    ...     return myarg.unit

    >>> myfunction(100*u.arcsec)
    Unit("arcsec")

Optionally, `None` keyword arguments are also supported; for such cases, the
input is only checked when a value other than `None` is passed::

    >>> @u.quantity_input(a='length', b='angle')
    ... def myfunction(a, b=None):
    ...     return a, b

    >>> myfunction(1.*u.km)  # doctest: +FLOAT_CMP
    (<Quantity 1. km>, None)
    >>> myfunction(1.*u.km, 1*u.deg)  # doctest: +FLOAT_CMP
    (<Quantity 1. km>, <Quantity 1. deg>)

Alternatively, you can use the `annotations syntax
<https://docs.python.org/3/library/typing.html>`_ to provide the units.
While the raw unit or string can be used, the preferred method is with the
unit-aware Quantity-annotation syntax.

``Quantity[unit or "string", metadata, ...]``

    >>> @u.quantity_input
    ... def myfunction(myarg: u.Quantity[u.arcsec]):
    ...     return myarg.unit
    >>>
    >>> myfunction(100*u.arcsec)
    Unit("arcsec")

You can also annotate for different types in non-unit expecting arguments:

    >>> @u.quantity_input
    ... def myfunction(myarg: u.Quantity[u.arcsec], nice_string: str):
    ...     return myarg.unit, nice_string
    >>> myfunction(100*u.arcsec, "a nice string")
    (Unit("arcsec"), 'a nice string')

The output can be specified to have a desired unit with a function annotation,
for example

    >>> @u.quantity_input
    ... def myfunction(myarg: u.Quantity[u.arcsec]) -> u.deg:
    ...     return myarg*1000
    >>>
    >>> myfunction(100*u.arcsec)  # doctest: +FLOAT_CMP
    <Quantity 27.77777778 deg>

This both checks that the return value of your function is consistent with what
you expect and makes it much neater to display the results of the function.

.. EXAMPLE END

Specifying a list of valid equivalent units or :ref:`physical_types` is
supported for functions that should accept inputs with multiple valid units:

    >>> @u.quantity_input(a=['length', 'speed'])
    ... def myfunction(a):
    ...     return a.unit

    >>> myfunction(1.*u.km)
    Unit("km")
    >>> myfunction(1.*u.km/u.s)
    Unit("km / s")

Representing Vectors with Units
===============================

|Quantity| objects can, like ``numpy`` arrays, be used to represent vectors or
matrices by assigning specific dimensions to represent the coordinates or
matrix elements, but that implies tracking those dimensions carefully. For
vectors :ref:`astropy-coordinates-representations` can be more convenient as
doing so allows you to use representations other than Cartesian (such as
spherical or cylindrical), as well as simple vector arithmetic.

.. _astropy-units-quantity-no-copy:

Creating and Converting Quantities without Copies
=================================================

When creating a |Quantity| using multiplication with a unit, a copy of the
underlying data is made. This can be avoided by passing on ``copy=False`` in
the initializer.

Examples
--------

.. EXAMPLE START: Creating and Converting Quantities without Copies

To avoid duplication using ``copy=False``::

    >>> a = np.arange(5.)
    >>> q = u.Quantity(a, u.m, copy=False)
    >>> q  # doctest: +FLOAT_CMP
    <Quantity [0., 1., 2., 3., 4.] m>
    >>> np.may_share_memory(a, q)
    True
    >>> a[0] = -1.
    >>> q  # doctest: +FLOAT_CMP
    <Quantity [-1.,  1.,  2.,  3.,  4.] m>

This may be particularly useful in functions which do not change their input
while ensuring that if a user passes in a |Quantity| then it will be converted
to the desired unit.

.. EXAMPLE END

As a shortcut, you can "shift" to the requested unit using the ``<<``
operator::

    >>> q = a << u.m
    >>> np.may_share_memory(a, q)
    True
    >>> q  # doctest: +FLOAT_CMP
    <Quantity [-1.,  1.,  2.,  3.,  4.] m>

The operator works identically to the initialization with ``copy=False``
mentioned above::

    >>> q << u.cm  # doctest: +FLOAT_CMP
    <Quantity [-100.,  100.,  200.,  300.,  400.] cm>

It can also be used for in-place conversion::

    >>> q <<= u.cm
    >>> q  # doctest: +FLOAT_CMP
    <Quantity [-100.,  100.,  200.,  300.,  400.] cm>
    >>> a  # doctest: +FLOAT_CMP
    array([-100.,  100.,  200.,  300.,  400.])


The `numpy.dtype` of a Quantity
===============================

|Quantity| subclasses `numpy.ndarray` and similarly accepts a ``dtype``
argument.

    >>> q = u.Quantity(1.0, dtype=np.float32)
    >>> q.dtype
    dtype('float32')

Like for `numpy.ndarray`, ``dtype`` does not have to be specified, in which case
the data is inspected to find the best ``dtype``. For `numpy` this means
integers remain integers, while |Quantity| instead upcasts integers to floats.

    >>> v = np.array(1)
    >>> np.issubdtype(v.dtype, np.integer)
    True

    >>> q = u.Quantity(1)
    >>> np.issubdtype(q.dtype, np.integer)
    False

|Quantity| promotes integer to floating types because it has a different default
value for ``dtype`` than `numpy` -- `numpy.inexact` versus `None`. For |Quantity|
to use the same ``dtype`` inspection as `numpy`, use ``dtype=None``.

    >>> q = u.Quantity(1, dtype=None)
    >>> np.issubdtype(q.dtype, np.integer)
    True

Note that `numpy.inexact` is a deprecated ``dtype`` argument for
`numpy.ndarray`. |Quantity| changes `numpy.inexact` to `numpy.float64`, but does
not change data that are already floating point or complex.


QTable
======

It is possible to use |Quantity| objects as columns in :mod:`astropy.table`.
See :ref:`quantity_and_qtable` for more details.

Subclassing Quantity
====================

To subclass |Quantity|, you generally proceed as you would when subclassing
|ndarray| (i.e., you typically need to override ``__new__()``, rather than
``__init__()``, and use the ``numpy.ndarray.__array_finalize__()`` method to
update attributes). For details, see the `NumPy documentation on subclassing
<https://numpy.org/doc/stable/user/basics.subclassing.html>`_.  To get a sense
of what is involved, have a look at |Quantity| itself, where, for example, the
``astropy.units.Quantity.__array_finalize__()`` method is used to pass on the
``unit``, at :class:`~astropy.coordinates.Angle`, where strings are parsed as
angles in the ``astropy.coordinates.Angle.__new__()`` method and at
:class:`~astropy.coordinates.Longitude`, where the
``astropy.coordinates.Longitude.__array_finalize__()`` method is used to pass
on the angle at which longitudes wrap.

Another method that is meant to be overridden by subclasses, specific to
|Quantity|, is ``astropy.units.Quantity.__quantity_subclass__()``. This is
called to decide which type of subclass to return, based on the unit of the
|Quantity| that is to be created. It is used, for example, in
:class:`~astropy.coordinates.Angle` to return a |Quantity| if a calculation
returns a unit other than an angular one. The implementation of this is via
:class:`~astropy.units.SpecificTypeQuantity`, which more generally allows users
to construct |Quantity| subclasses that have methods that are useful only for a
specific physical type.
