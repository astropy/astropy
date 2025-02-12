Combining and Defining Units
****************************

Basic example
=============

.. EXAMPLE START: Combining Units and Quantities

Units and quantities can be combined together using the regular Python
numeric operators::

  >>> from astropy import units as u
  >>> fluxunit = u.erg / (u.cm ** 2 * u.s)
  >>> fluxunit
  Unit("erg / (s cm2)")
  >>> 52.0 * fluxunit  # doctest: +FLOAT_CMP
  <Quantity  52. erg / (s cm2)>
  >>> 52.0 * fluxunit / u.s  # doctest: +FLOAT_CMP
  <Quantity  52. erg / (cm2 s2)>

.. EXAMPLE END

Fractional powers
=================

.. EXAMPLE START: Using Fractional Powers with Units

Units support fractional powers, which retain their precision through
complex operations. To do this, it is recommended to use
:class:`fractions.Fraction` objects::

  >>> from fractions import Fraction
  >>> Franklin = u.g ** Fraction(1, 2) * u.cm ** Fraction(3, 2) * u.s ** -1

.. note::

    Floating-point powers that are effectively the same as fractions
    with a denominator less than 10 are implicitly converted to
    `~fractions.Fraction` objects under the hood. Therefore, the
    following are equivalent::

        >>> x = u.m ** Fraction(1, 3)
        >>> x.powers
        [Fraction(1, 3)]
        >>> x = u.m ** (1. / 3.)
        >>> x.powers
        [Fraction(1, 3)]

.. EXAMPLE END

Defining units
==============

.. EXAMPLE START: Defining New Units

Users are free to define new units, either fundamental or compound,
using the :func:`~astropy.units.def_unit` function::

  >>> bakers_fortnight = u.def_unit('bakers_fortnight', 13 * u.day)

The addition of a string gives the new unit a name that will show up
when the unit is printed::

  >>> 10. * bakers_fortnight  # doctest: +FLOAT_CMP
  <Quantity  10. bakers_fortnight>

Creating a new fundamental unit is also possible::

  >>> titter = u.def_unit('titter')
  >>> chuckle = u.def_unit('chuckle', 5 * titter)
  >>> laugh = u.def_unit('laugh', 4 * chuckle)
  >>> guffaw = u.def_unit('guffaw', 3 * laugh)
  >>> rofl = u.def_unit('rofl', 4 * guffaw)
  >>> death_by_laughing = u.def_unit('death_by_laughing', 10 * rofl)
  >>> (1. * rofl).to(titter)  # doctest: +FLOAT_CMP
  <Quantity  240. titter>

Users can see the definition of a unit and its :ref:`decomposition
<decomposing>` via::

  >>> rofl.represents
  Unit("4 guffaw")
  >>> rofl.decompose()
  Unit("240 titter")

By default, custom units are not searched by methods such as
:meth:`~astropy.units.core.UnitBase.find_equivalent_units`. However, they
can be enabled by calling :func:`~astropy.units.add_enabled_units`::

  >>> kmph = u.def_unit('kmph', u.km / u.h)
  >>> (u.m / u.s).find_equivalent_units()
  There are no equivalent units
  >>> u.add_enabled_units([kmph])
  <astropy.units.core._UnitContext object at ...>
  >>> (u.m / u.s).find_equivalent_units()
    Primary name | Unit definition | Aliases
  [
    kmph         | 0.277778 m / s  |         ,
  ]

.. testcleanup::

    >>> u.core._unit_registries.pop()  # doctest: +IGNORE_OUTPUT

If new units are defined with prefixes enabled, the prefixed units must be
explicitly enabled as well, e.g., by using the ``namespace`` argument::

  >>> new_units = dict()
  >>> foo = u.def_unit(['Fo', 'foo'], prefixes=True, namespace=new_units)
  >>> u.add_enabled_units(new_units)
  <astropy.units.core._UnitContext object at ...>

Now, the prefixed units can be parsed etc::

  >>> print(u.Unit("megafoo").find_equivalent_units())
  Primary name | Unit definition | Aliases
  [
    Fo           | irreducible     | foo     ,
  ]
  >>> print(u.Unit("megafoo").to(u.Unit("kFo")))
  1000.0

.. testcleanup::

    >>> u.core._unit_registries.pop()  # doctest: +IGNORE_OUTPUT

.. EXAMPLE END
