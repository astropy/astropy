Composing and defining units
============================

Units can be composed together using the regular Python numeric
operators.  For example::

  >>> fluxunit = u.erg / (u.cm ** 2 * u.s)
  >>> fluxunit
  Unit("erg / (s cm2)")

Users are free to define new units, either fundamental or compound
using the `def_unit` function.  For example::

  >>> bakers_fortnight = u.def_unit('bakers_fortnight', 13 * u.day)

The addition of a string gives the new unit a name that will show up
when the unit is printed.

Creating a new fundamental unit is simple::

  >>> titter = u.def_unit('titter')
  >>> chuckle = u.def_unit('chuckle', 5 * titter)
  >>> laugh = u.def_unit('laugh', 4 * chuckle)
  >>> guffaw = u.def_unit('guffaw', 3 * laugh)
  >>> rofl = u.def_unit('rofl', 4 * guffaw)
  >>> death_by_laughing = u.def_unit('death_by_laughing', 10 * rofl)
  >>> rofl.to(titter, 1)
  240

Reducing a unit to its irreducible parts
----------------------------------------

A unit can be decomposed into its irreducible parts using the `decompose`
method::

  >>> u.Ry
  Unit("Ry")
  >>> u.Ry.decompose()
  Unit("2.18e-18 m2 kg / (s2)")

.. TODO: Add function and description to factor units into high-level pieces
