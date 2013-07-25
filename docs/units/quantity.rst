Quantity
========

.. |quantity| replace:: :class:`~astropy.units.quantity.Quantity`

The |quantity| object is meant to represent a value that has some unit
associated with the number.

Creating Quantity instances
---------------------------

|quantity| objects are created through multiplication or divison with
:class:`~astropy.units.core.Unit` objects. For example, to create a |quantity|
to represent 15 m/s:

    >>> import astropy.units as u
    >>> 15 * u.m / u.s
    <Quantity 15 m / s>

or 1.14/s:

    >>> 1.14 / u.s
    <Quantity 1.14 1 / s>

You can also create instances using the |quantity| constructor directly, by
specifying a value and unit:

    >>> u.Quantity(15, u.m / u.s)
    <Quantity 15 m / s>

|quantity| objects can also be created automatically from Numpy arrays:

    >>> import numpy as np
    >>> np.array([1,2,3]) * u.m
    <Quantity [1 2 3] m>

|quantity| objects can also be created from sequences of |quantity| objects,
and will automatically convert to Numpy arrays.

    >>> qlst = [60 * u.s, 120 * u.s]
    >>> u.Quantity(qlst, u.minute)
    <Quantity [ 1.  2.] min>

Finally, the current unit and value can be accessed via the ``unit`` and
``value`` attributes:

    >>> q = 2.3 * u.m / u.s
    >>> q.unit
    Unit("m / s")
    >>> q.value
    2.3

Converting to different units
-----------------------------

|quantity| objects can be converted to different units using the
:meth:`~astropy.units.quantity.Quantity.to` method::

    >>> q = 2.3 * u.m / u.s
    >>> q.to(u.km / u.h)
    <Quantity 8.28 km / h>

For convenience, the `si` and `cgs` attributes can be used to convert the
|quantity| to base S.I. or c.g.s units:

    >>> q = 2.4 * u.m / u.s
    >>> q.si
    <Quantity 2.4 m / s>
    >>> q.cgs
    <Quantity 240.0 cm / s>

The value of converting to another unit can also be obtained, by specifying
the desired unit as an attribute::

    >>> q = 180.0 * u.degree
    >>> q.radian
    3.141592653589793
    >>> # The above is shorthand for...
    >>> q.to(u.radian).value
    3.141592653589793

In code, it is considered best practice to use the `q.SOME_UNIT` rather than
`q.value`, because then the value will always be in the requested unit.

Arithmetic
----------

Addition and Subtraction
~~~~~~~~~~~~~~~~~~~~~~~~

Addition or subtraction between |quantity| objects is supported when their
units are equivalent. When the units are equal, the resulting object has the
same unit:

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

Addition and subtraction is not supported between |quantity| objects and basic
numeric types:

    >>> 13.5 * u.km + 19.412
    TypeError: Object of type '<type 'float'>' cannot be added with a Quantity
    object. Addition is only supported between Quantity objects.

except for dimensionless quantities (see `Dimensionless quantities`_).

Multiplication and Division
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Multiplication and division is supported between |quantity| objects with any
units, and with numeric types. For these operations between objects with
equivalent units, the **resulting object has composite units**:

    >>> 1.1 * u.m * 140.3 * u.cm
    <Quantity 154.33 cm m>
    >>> 140.3 * u.cm * 1.1 * u.m
    <Quantity 154.33 cm m>
    >>> 1. * u.m / (20. * u.cm)
    <Quantity 0.05 m / cm>
    >>> 20. * u.cm / (1. * u.m)
    <Quantity 20.0 cm / m>

For multiplication, you can change how to represent the resulting object by
using the :meth:`~astropy.units.quantity.Quantity.to` method:

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
    <Quantity 3.41950972779 kg m2 / s2>


Numpy functions
---------------

|quantity| objects are actually full Numpy arrays (the |quantity| object class
inherits from and extends the ``np.ndarray`` class), and we have tried to
ensure that most Numpy functions behave properly with units:

    >>> q = np.array([1., 2., 3., 4.]) * u.m / u.s
    >>> np.mean(q)
    <Quantity 2.5 m / s>
    >>> np.std(q)
    <Quantity 1.11803398875 m / s>

including functions that only accept specific units such as angles::

    >>> q = 30. * u.deg
    >>> np.sin(q)
    <Quantity 0.5 >

or dimensionless quantities::

    >>> nu = 3 * u.GHz
    >>> T = 30 * u.K
    >>> np.exp(- h * nu / (k_B * T))
    <Quantity 0.995212254619 >

(see `Dimensionless quantities`_ for more details).

Dimensionless quantities
------------------------

Dimensionless quantities have the characteristic that if they are added or
subtracted from a Python scalar or unitless ``ndarray``, or if they are passed
to a Numpy function that takes dimensionless quantities, the units are
simplified so that the quantity is dimensionless and scale-free. For example::

    >>> 1. + 1. * u.m / u.km
    <Quantity 1.001 >

which is different from::

    >>> 1. + (1. * u.m / u.km).value
    2.0

In the latter case, the result is ``2.0`` because the unit of ``(1. * u.m /
u.km)`` is not scale-free by default::

    >>> q = (1. * u.m / u.km)
    >>> q.unit
    Unit("m / km")
    >>> q.unit.decompose()
    Unit(dimensionless with a scale of 0.001)

However, when combining with a non-quantity object, the unit is automatically
decomposed to be scale-free, giving the expected result.

This also occurs when passing dimensionless quantities to functions that take
dimensionless quantities::

    >>> nu = 3 * u.GHz
    >>> T = 30 * u.K
    >>> np.exp(- h * nu / (k_B * T))
    <Quantity 0.995212254619 >

The result is independent from the units the different quantities were specified in::

    >>> nu = 3.e9 * u.Hz
    >>> T = 30 * u.K
    >>> np.exp(- h * nu / (k_B * T))
    <Quantity 0.995212254619 >

Converting to plain Python scalars or Numpy arrays
--------------------------------------------------

Converting |quantity| objects does not work for non-dimensionless quantities::

    >>> float(3. * u.m)
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "/Users/tom/Library/Python/2.7/lib/python/site-packages/astropy-0.3.dev4691-py2.7-macosx-10.6-x86_64.egg/astropy/units/quantity.py", line 627, in __float__
        raise TypeError('Only dimensionless scalar quantities can be '
    TypeError: Only dimensionless scalar quantities can be converted to Python scalars

Instead, only dimensionless values can be converted to plain Python scalars::

    >>> float(3. * u.m / (4. * u.m))
    0.75

Note that scaled dimensionless quantities such as ``m / km`` also do not work::

    >>> float(3. * u.m / (4. * u.km))
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "/Users/tom/Library/Python/2.7/lib/python/site-packages/astropy-0.3.dev4691-py2.7-macosx-10.6-x86_64.egg/astropy/units/quantity.py", line 627, in __float__
        raise TypeError('Only dimensionless scalar quantities can be '
    TypeError: Only dimensionless scalar quantities can be converted to Python scalars

If you want to simplify e.g. dimensionless quantities to their true
dimensionless value, then you can make use of the
:meth:`~astropy.units.quantity.Quantity.decompose` method:

    >>> q = 3. * u.m / (4. * u.km)
    >>> q
    <Quantity 0.75 m / km>
    >>> q.decompose()
    <Quantity 0.00075 >
    >>> float(q.decompose())
    0.00075

Similarly, ``int`` and ``long`` work, but only for dimensionless unscaled quantities::

    >>> int(6. * u.m / (2. * u.m))
    3

Since |quantity| objects are actual Numpy arrays, we are not able to ensure
that only dimensionless quantities are converted to Numpy arrays::

    >>> np.array([1, 2, 3] * u.m)
    array([1, 2, 3])

