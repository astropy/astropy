Quantity
========

.. |quantity| replace:: :class:`~astropy.units.Quantity`

The |quantity| object is meant to represent a value that has some unit
associated with the number.

Creating Quantity instances
---------------------------

|quantity| objects are created through multiplication or division with
:class:`~astropy.units.Unit` objects. For example, to create a |quantity|
to represent 15 m/s:

    >>> import astropy.units as u
    >>> 15 * u.m / u.s
    <Quantity 15.0 m / s>

.. note:: |quantity| objects are converted to float by default.

As another example:

    >>> 1.25 / u.s
    <Quantity 1.25 1 / s>

You can also create instances using the |quantity| constructor directly, by
specifying a value and unit:

    >>> u.Quantity(15, u.m / u.s)
    <Quantity 15.0 m / s>

|quantity| objects can also be created automatically from Numpy arrays
or Python sequences:

    >>> [1, 2, 3] * u.m
    <Quantity [ 1., 2., 3.] m>
    >>> import numpy as np
    >>> np.array([1, 2, 3]) * u.m
    <Quantity [ 1., 2., 3.] m>

|quantity| objects can also be created from sequences of |quantity|
objects, as long as all of their units are equivalent, and will
automatically convert to Numpy arrays:

    >>> qlst = [60 * u.s, 1 * u.min]
    >>> u.Quantity(qlst, u.minute)
    <Quantity [ 1.,  1.] min>

Finally, the current unit and value can be accessed via the
`~astropy.units.quantity.Quantity.unit` and
`~astropy.units.quantity.Quantity.value` attributes:

    >>> q = 2.5 * u.m / u.s
    >>> q.unit
    Unit("m / s")
    >>> q.value
    2.5

Converting to different units
-----------------------------

|quantity| objects can be converted to different units using the
:meth:`~astropy.units.quantity.Quantity.to` method:

    >>> q = 2.3 * u.m / u.s
    >>> q.to(u.km / u.h)
    <Quantity 8.2... km / h>

For convenience, the `~astropy.units.quantity.Quantity.si` and
`~astropy.units.quantity.Quantity.cgs` attributes can be used to
convert the |quantity| to base S.I. or c.g.s units:

    >>> q = 2.4 * u.m / u.s
    >>> q.si
    <Quantity 2... m / s>
    >>> q.cgs
    <Quantity 240.0 cm / s>

Arithmetic
----------

Addition and Subtraction
~~~~~~~~~~~~~~~~~~~~~~~~

Addition or subtraction between |quantity| objects is supported when their
units are equivalent. When the units are equal, the resulting object has the
same unit:

    >>> 11 * u.s + 30 * u.s
    <Quantity 41.0 s>
    >>> 30 * u.s - 11 * u.s
    <Quantity 19.0 s>

If the units are equivalent, but not equal (e.g. kilometer and meter), the
resulting object **has units of the object on the left**:

    >>> 1100.1 * u.m + 13.5 * u.km
    <Quantity 14600.1 m>
    >>> 13.5 * u.km + 1100.1 * u.m
    <Quantity 14.600... km>
    >>> 1100.1 * u.m - 13.5 * u.km
    <Quantity -12399.9 m>
    >>> 13.5 * u.km - 1100.1 * u.m
    <Quantity 12.399... km>

Addition and subtraction is not supported between |quantity| objects and basic
numeric types:

    >>> 13.5 * u.km + 19.412
    Traceback (most recent call last):
      ...
    UnitsError: Can only apply 'add' function to dimensionless
    quantities when other argument is not a quantity (unless the
    latter is all zero/infinity/nan)

except for dimensionless quantities (see `Dimensionless quantities`_).

Multiplication and Division
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Multiplication and division is supported between |quantity| objects with any
units, and with numeric types. For these operations between objects with
equivalent units, the **resulting object has composite units**:

    >>> 1.1 * u.m * 140.3 * u.cm
    <Quantity 154.33... cm m>
    >>> 140.3 * u.cm * 1.1 * u.m
    <Quantity 154.33... cm m>
    >>> 1. * u.m / (20. * u.cm)
    <Quantity 0.05... m / cm>
    >>> 20. * u.cm / (1. * u.m)
    <Quantity 20.0 cm / m>

For multiplication, you can change how to represent the resulting object by
using the :meth:`~astropy.units.quantity.Quantity.to` method:

    >>> (1.1 * u.m * 140.3 * u.cm).to(u.m**2)
    <Quantity 1.5433... m2>
    >>> (1.1 * u.m * 140.3 * u.cm).to(u.cm**2)
    <Quantity 15433.0... cm2>

For division, if the units are equivalent, you may want to make the resulting
object dimensionless by reducing the units. To do this, use the
:meth:`~astropy.units.quantity.Quantity.decompose()` method:

    >>> (20. * u.cm / (1. * u.m)).decompose()
    <Quantity 0.2...>

This method is also useful for more complicated arithmetic:

    >>> 15. * u.kg * 32. * u.cm * 15 * u.m / (11. * u.s * 1914.15 * u.ms)
    <Quantity 0.341950972... cm kg m / (ms s)>
    >>> (15. * u.kg * 32. * u.cm * 15 * u.m / (11. * u.s * 1914.15 * u.ms)).decompose()
    <Quantity 3.41950972... kg m2 / s2>


Numpy functions
---------------

|quantity| objects are actually full Numpy arrays (the |quantity|
object class inherits from and extends the ``numpy.ndarray`` class), and
we have tried to ensure that most Numpy functions behave properly with
quantities:

    >>> q = np.array([1., 2., 3., 4.]) * u.m / u.s
    >>> np.mean(q)
    <Quantity 2.5 m / s>
    >>> np.std(q)
    <Quantity 1.118033... m / s>

including functions that only accept specific units such as angles:

    >>> q = 30. * u.deg
    >>> np.sin(q)
    <Quantity 0.4999999...>

or dimensionless quantities:

    >>> from astropy.constants import h, k_B
    >>> nu = 3 * u.GHz
    >>> T = 30 * u.K
    >>> np.exp(-h * nu / (k_B * T))
    <Quantity 0.99521225...>

(see `Dimensionless quantities`_ for more details).

Dimensionless quantities
------------------------

Dimensionless quantities have the characteristic that if they are
added or subtracted from a Python scalar or unitless `~numpy.ndarray`,
or if they are passed to a Numpy function that takes dimensionless
quantities, the units are simplified so that the quantity is
dimensionless and scale-free. For example:

    >>> 1. + 1. * u.m / u.km
    <Quantity 1.00...>

which is different from:

    >>> 1. + (1. * u.m / u.km).value
    2.0

In the latter case, the result is ``2.0`` because the unit of ``(1. * u.m /
u.km)`` is not scale-free by default:

    >>> q = (1. * u.m / u.km)
    >>> q.unit
    Unit("m / km")
    >>> q.unit.decompose()
    Unit(dimensionless with a scale of 0.001)

However, when combining with a non-quantity object, the unit is automatically
decomposed to be scale-free, giving the expected result.

This also occurs when passing dimensionless quantities to functions that take
dimensionless quantities:

    >>> nu = 3 * u.GHz
    >>> T = 30 * u.K
    >>> np.exp(- h * nu / (k_B * T))
    <Quantity 0.99521225...>

The result is independent from the units the different quantities were specified in:

    >>> nu = 3.e9 * u.Hz
    >>> T = 30 * u.K
    >>> np.exp(- h * nu / (k_B * T))
    <Quantity 0.99521225...>

Converting to plain Python scalars
----------------------------------

Converting |quantity| objects does not work for non-dimensionless quantities:

    >>> float(3. * u.m)
    Traceback (most recent call last):
      ...
    TypeError: Only dimensionless scalar quantities can be converted
    to Python scalars

Instead, only dimensionless values can be converted to plain Python scalars:

    >>> float(3. * u.m / (4. * u.m))
    0.75
    >>> float(3. * u.km / (4. * u.m))
    750.0
    >>> int(6. * u.km / (2. * u.m))
    3000

Functions Accepting Quantities
------------------------------

Validation of quantity arguments to functions can lead to many repetitons 
of the same checking code. A decorator is provided which verifies that certain
arguments to a function are `~astropy.units.Quantity` objects and that the units
are compatible with a desired unit.

The decorator does not convert the unit to the desired unit, say arcseconds
to degrees, it merely checks that such a conversion is possible, thus verifying 
that the `~astropy.units.Quantity` argument can be used in calculations.

The decorator `~astropy.units.quantity_input` accepts keyword arguments to 
spcifiy which arguments should be validated and what unit they are expected to 
be compatible with:

    >>> @u.quantity_input(myarg=u.deg)
    ... def myfunction(myarg):
    ...     return myarg.unit

    >>> myfunction(100*u.arcsec)
    Unit("arcsec")

Under Python 3 you can use the annotations syntax to provide the units:

    >>> @u.quantity_input  # doctest: +SKIP
    ... def myfunction(myarg: u.arcsec):
    ...     return myarg.unit

    >>> myfunction(100*u.arcsec)  # doctest: +SKIP
    Unit("arcsec")

Known issues with conversion to numpy arrays
--------------------------------------------

Since |quantity| objects are Numpy arrays, we are not able to ensure
that only dimensionless quantities are converted to Numpy arrays:

    >>> np.array([1, 2, 3] * u.m)
    array([ 1., 2., 3.])

Similarly, while most numpy functions work properly, a few have :ref:`known
issues <quantity_issues>`, either ignoring the unit (e.g., ``np.dot``) or
not reinitializing it properly (e.g., ``np.hstack``).  This propagates to
more complex functions such as ``np.linalg.norm`` and
``scipy.integrate.odeint``.

Subclassing Quantity
--------------------

To subclass |quantity|, one generally proceeds as one would when subclassing
:class:`~numpy.ndarray`, i.e., one typically needs to override ``__new__``
(rather than ``__init__``) and uses the ``numpy.ndarray.__array_finalize__``
method to update attributes.  For details, see the `numpy documentation on
subclassing
<http://docs.scipy.org/doc/numpy/user/basics.subclassing.html>`__.  For
examples, one can look at |quantity| itself, where, e.g., the
``astropy.units.Quantity.__array_finalize__`` method is used to pass on the
``unit``, at :class:`~astropy.coordinates.Angle`, where strings are parsed
as angles in the ``astropy.coordinates.Angle.__new__`` method and at
:class:`~astropy.coordinates.Longitude`, where the
``astropy.coordinates.Longitude.__array_finalize__`` method is used to pass
on the angle at which longitudes wrap.

Another method that is meant to be overridden by subclasses, one specific to
|quantity|, is ``astropy.units.Quantity.__quantity_subclass__``.  This is
called to decide which type of subclass to return, based on the unit of the
quantity that is to be created.  It is used, e.g., in
:class:`~astropy.coordinates.Angle` to return a |quantity| if a calculation
returns a unit other than an angular one.
