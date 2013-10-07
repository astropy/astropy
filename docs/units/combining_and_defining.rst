Combining and defining units
============================

Units can be combined together using the regular Python numeric
operators.  For example::

  >>> from astropy import units as u
  >>> fluxunit = u.erg / (u.cm ** 2 * u.s)
  >>> fluxunit
  Unit("erg / (cm2 s)")

Users are free to define new units, either fundamental or compound
using the `~astropy.units.core.def_unit` function.  For example::

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

By default, custom units are not searched by methods such as
`UnitBase.find_equivalent_units`.  However, they can be enabled by
calling `~astropy.units.add_enabled_units`::

  >>> kmph = u.def_unit('kmph', u.km / u.h)
  >>> (u.m / u.s).find_equivalent_units()
  # ... nothing ...
  >>> u.add_enabled_units([kmph])
  >>> (u.m / u.s).find_equivalent_units()
    Primary name | Unit definition | Aliases
  [
    kmph         | 0.277778 m / s  |         ,
  ]
