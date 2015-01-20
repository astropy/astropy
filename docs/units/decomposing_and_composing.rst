Decomposing and composing units
===============================

Reducing a unit to its irreducible parts
----------------------------------------

A unit or quantity can be decomposed into its irreducible parts using
the `Unit.decompose <astropy.units.core.UnitBase.decompose>` or
`Quantity.decompose <astropy.units.quantity.Quantity.decompose>`
methods::

  >>> from astropy import units as u
  >>> u.Ry
  Unit("Ry")
  >>> u.Ry.decompose()
  Unit("2.17987e-18 kg m2 / s2")

You can limit the selection of units that you want to decompose to
using the ``bases`` keyword argument::

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
  [Unit("Ry"),
   Unit("2.17987e-18 J"),
   Unit("2.17987e-11 erg"),
   Unit("13.6057 eV")]

Some other interesting examples::

   >>> (u.s ** -1).compose()  # doctest: +SKIP
   [Unit("Bq"), Unit("Hz"), Unit("3.7e+10 Ci")]

Composition can be combined with :ref:`unit_equivalencies`::

   >>> (u.s ** -1).compose(equivalencies=u.spectral())  # doctest: +SKIP
   [Unit("m"),
    Unit("Hz"),
    Unit("J"),
    Unit("Bq"),
    Unit("3.24078e-17 pc"),
    Unit("1.057e-16 lyr"),
    Unit("6.68459e-12 AU"),
    Unit("1.4378e-09 solRad"),
    Unit("0.01 k"),
    Unit("100 cm"),
    Unit("1e+06 micron"),
    Unit("1e+07 erg"),
    Unit("1e+10 Angstrom"),
    Unit("3.7e+10 Ci"),
    Unit("4.58743e+17 Ry"),
    Unit("6.24151e+18 eV")]

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
   [Unit("10 P / s"), Unit("10 Ba")]

There is also a shorthand for this which only returns the first of
many possible matches::

   >>> u.Pa.cgs
   Unit("10 P / s")

This is equivalent to decomposing into the new system and then
composing into the most complex units possible, though
`~astropy.units.core.UnitBase.to_system` adds some extra logic to
return the results sorted in the most useful order::

   >>> u.Pa.decompose(bases=u.cgs.bases)
   Unit("10 g / (cm s2)")
   >>> _.compose(units=u.cgs)
   [Unit("10 Ba"), Unit("10 P / s")]
