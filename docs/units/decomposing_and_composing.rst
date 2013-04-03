Decomposing and composing units
===============================

Reducing a unit to its irreducible parts
----------------------------------------

A unit can be decomposed into its irreducible parts using the
`~astropy.units.core.UnitBase.decompose` method::

  >>> from astropy import units as u
  >>> u.Ry
  Unit("Ry")
  >>> u.Ry.decompose()
  Unit("2.17987e-18 m2 kg / s2")

You can limit the selection of units that you want to decompose to
using the `bases` keyword argument::

  >>> u.Ry.decompose(bases=[u.m, u.N])
  Unit("2.17987e-18 m N")

This is also useful to decompose to a particular system.  For example,
to decompose the Rydberg unit in terms of CGS units::

  >>> u.Ry.decompose(bases=u.cgs.bases)
  Unit("2.17987e-11 cm2 g / s2")

Automatically composing a unit into more complex units
------------------------------------------------------

Conversely, a unit may be recomposed back into more complex units
using the `~astropy.units.core.UnitBase.compose` method.  Since there
may be multiple equally good results, a list is always returned::

  >>> x = u.Ry.decompose()
  >>> x.compose()
  [Unit("1 Ry"),
   Unit("5.21002e-22 kcal"),
   Unit("2.06612e-21 BTU"),
   Unit("5.21002e-19 cal"),
   Unit("2.17987e-18 J"),
   Unit("2.17987e-11 erg"),
   Unit("13.6057 eV")]

Some other interesting examples::

   >>> (u.s ** -1).compose()
   [Unit("Hz"),
    Unit("1 / s"),
    Unit("60 / min"),
    Unit("3600 / h"),
    Unit("86164.1 / sday"),
    Unit("86400 / d"),
    Unit("604800 / wk"),
    Unit("1.2096e+06 / fortnight"),
    Unit("3.15576e+07 / yr"),
    Unit("3.15576e+07 / a")]

Composition can be combined with :ref:`unit_equivalencies`::

   >>> (u.s ** -1).compose(equivalencies=u.spectral())
   [Unit("m"),
    Unit("Hz"),
    Unit("J"),
    Unit("1 / s"),
    Unit("3.24078e-17 pc"),
    Unit("1.057e-16 lyr"),
    Unit("6.68459e-12 AU"),
    Unit("1.4378e-09 solRad"),
    Unit("0.000239006 kcal"),
    Unit("0.000621371 mi"),
    Unit("0.000947817 BTU"),
    Unit("0.239006 cal"),
    Unit("1.09361 yd"),
    Unit("3.28084 ft"),
    Unit("39.3701 inch"),
    Unit("100 cm"),
    Unit("1e+06 micron"),
    Unit("1e+07 erg"),
    Unit("1e+10 Angstrom"),
    Unit("4.58743e+17 Ry"),
    Unit("6.24151e+18 eV"),
    Unit("60 / min"),
    Unit("3600 / h"),
    Unit("86164.1 / sday"),
    Unit("86400 / d"),
    Unit("604800 / wk"),
    Unit("1.2096e+06 / fortnight"),
    Unit("3.15576e+07 / a"),
    Unit("3.15576e+07 / yr")]

Obviously a name doesn't exist for every arbitrary derived unit
imaginable.  In that case, the system will do its best to reduce the
unit to the fewest possible symbols::

   >>> (u.cd * u.sr * u.V * u.s).compose()
   [Unit("lm Wb")]

Converting between systems
--------------------------

Built on top of this functionality is a convenience method to convert
between unit systems.

   >>> u.Pa.to_system(u.cgs)
   [Unit("10 Ba")]

This is equivalent to decomposing into the new system and then
composing into the most complex units possible, though `to_system`
adds some extra logic to return the results sorted in the most useful
order::

   >>> u.Pa.decompose(bases=u.cgs.bases)
   Unit("10 g / (cm s2)")
   >>> _.compose(units=u.cgs)
   [Unit("10 Ba")]
