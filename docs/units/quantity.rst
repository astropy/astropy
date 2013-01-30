Quantity
========

The :class:`~astropy.units.quantity.Quantity` object is meant to represent a
value that has some unit associated with the number.

Creating Quantity instances
---------------------------

:class:`~astropy.units.quantity.Quantity` objects are created through
multiplication or divison with :class:`~astropy.units.core.Unit` objects. For
example, to create a :class:`~astropy.units.quantity.Quantity` to represent 15
m/s:

    >>> import astropy.units as u
    >>> 15 * u.m / u.s
    <Quantity 15 m / (s)>

or 1.14/s:

    >>> 1.14 / u.s
    <Quantity 1.14 1 / (s)>

You can also create instances using the
:class:`~astropy.units.quantity.Quantity` constructor directly, by specifying
a value and unit:

    >>> u.Quantity(15, u.m / u.s)
    <Quantity 15 m / (s)>

:class:`~astropy.units.quantity.Quantity` objects can also be created
automatically from Numpy arrays:

    >>> import numpy as np
    >>> np.array([1,2,3]) * u.m
    <Quantity [1 2 3] m>

Finally, the current unit and value can be accessed via the ``unit`` and
``value`` attributes:

    >>> q = 2.3 * u.m / u.s
    >>> q.unit
    Unit("m / (s)")
    >>> q.value
    2.3

Converting to different units
-----------------------------

:class:`~astropy.units.quantity.Quantity` objects can be converted to
different units using the :meth:`~astropy.units.quantity.Quantity.to` method::

    >>> q = 2.3 * u.m / u.s
    >>> q.to(u.km / u.h)
    <Quantity 8.28 km / (h)>

For convenience, the `si` and `cgs` attributes can be used to convert the
:class:`~astropy.units.quantity.Quantity` to base S.I. or c.g.s units:

    >>> q = 2.4 * u.m / u.s
    >>> q.si
    <Quantity 2.4 m / (s)>
    >>> q.cgs
    <Quantity 240.0 cm / (s)>

Arithmetic
----------

Addition and Subtraction
~~~~~~~~~~~~~~~~~~~~~~~~

Addition or subtraction between :class:`~astropy.units.quantity.Quantity`
objects is supported when their units are equivalent. When the units are
equal, the resulting object has the same unit:

    >>> 11 * u.s + 30 * u.s
    <Quantity 41 s>
    >>> 30 * u.s - 11 * u.s
    <Quantity 19 s>

If the units are equivalent, but not equal (e.g. kilometer and meter), the
resulting object **has units of the object on the left**:

    >>> 1100.1 * u.m + 13.5 * u.km
    <Quantity 14600.1 m>
    >>> 13.5 * u.km + 1100.1 * u.m
    <Quantity 14.6001 km>
    >>> 1100.1 * u.m - 13.5 * u.km
    <Quantity -12399.9 m>
    >>> 13.5 * u.km - 1100.1 * u.m
    <Quantity 12.3999 km>

Addition and subtraction is not supported between
:class:`~astropy.units.quantity.Quantity` objects and basic numeric types:

    >>> 13.5 * u.km + 19.412
    TypeError: Object of type '<type 'float'>' cannot be added with a Quantity
    object. Addition is only supported between Quantity objects.

Multiplication and Division
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Multiplication and division is supported between
:class:`~astropy.units.quantity.Quantity` objects with any units, and with
numeric types. For these operations between objects with equivalent units, the
**resulting object has composite units**:

    >>> 1.1 * u.m * 140.3 * u.cm
    <Quantity 154.33 cm m>
    >>> 140.3 * u.cm * 1.1 * u.m
    <Quantity 154.33 cm m>
    >>> 1. * u.m / (20. * u.cm)
    <Quantity 0.05 m / (cm)>
    >>> 20. * u.cm / (1. * u.m)
    <Quantity 20.0 cm / (m)>

For multiplication, you can change how to represent the resulting object by
using the `~astropy.units.quantity.Quantity.to` method:

    >>> (1.1 * u.m * 140.3 * u.cm).to(u.m**2)
    <Quantity 1.5433 m2>
    >>> (1.1 * u.m * 140.3 * u.cm).to(u.cm**2)
    <Quantity 15433.0 cm2>

For division, if the units are equivalent, you may want to make the resulting
object dimensionless by reducing the units. To do this, use the
:meth:`~astropy.units.quantity.Quantity.decompose()` method:

    >>> (20. * u.cm / (1. * u.m)).decompose()
    <Quantity 0.2 >

This method is also useful for more complicated arithmetic:

    >>> 15. * u.kg * 32. * u.cm * 15 * u.m / (11. * u.s * 1914.15 * u.ms)
    <Quantity 0.341950972779 cm kg m / (ms s)>
    >>> (15. * u.kg * 32. * u.cm * 15 * u.m / (11. * u.s * 1914.15 * u.ms)).decompose()
    <Quantity 3.41950972779 kg m2 / (s2)>
