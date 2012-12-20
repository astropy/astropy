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
  Unit("2.18e-18 m2 kg / (s2)")

You can limit the selection of units that you want to decompose to
using the `bases` keyword argument::

  >>> u.Ry.decompose(bases=[u.m, u.N])
  Unit("2.179872e-18 m N")

This is also useful to decompose to a particular system.  For example,
to decompose the Rydberg unit in terms of CGS units::

  >>> u.Ry.decompose(bases=u.cgs.bases)
  Unit("2.179872e-11 cm2 g / (s2)")

Automatically composing a unit into more complex units
------------------------------------------------------

Conversely, a unit may be recomposed back into more complex units
using the `~astropy.units.core.UnitBase.compose` method.  Since there
may be multiple equally good results, a list is always returned::

  >>> x = u.Ry.decompose()
  >>> x.compose()
  [Unit("1.000000e+00 Ry"), Unit("2.179872e-18 J"), Unit("2.066120e-21 BTU"),
   Unit("5.210019e-22 kcal"), Unit("5.210019e-19 cal"),
   Unit("1.360569e+01 eV"), Unit("2.179872e-11 erg")]

Some other interesting examples::

   >>> (u.s ** -1).compose()
   [Unit("Hz")]

Composition can be combined with :ref:`unit-equivalencies`::

   >>> (u.s ** -1).compose(equivs=u.spectral())
   [Unit("Hz"), Unit("1.437798e-09 solRad"), Unit("1.057001e-16 lyr"),
    Unit("1.000000e+06 micron"), Unit("m"), Unit("3.240779e-17 pc"),
    Unit("1.000000e+10 Angstrom"), Unit("6.684587e-12 AU"),
    Unit("1.000000e+02 cm"), Unit("3.937008e+01 inch"),
    Unit("3.280840e+00 ft"), Unit("1.093613e+00 yd"), Unit("6.213712e-04 mi"),
    Unit("4.587425e+17 Ry"), Unit("J"), Unit("9.478171e-04 BTU"),
    Unit("2.390057e-04 kcal"), Unit("2.390057e-01 cal"),
    Unit("6.241509e+18 eV"), Unit("1.000000e+07 erg")]

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
   [Unit("1.000000e+01 Ba")]

This is equivalent to decomposing into the new system and then
composing into the most complex units possible::

   >>> u.Pa.decompose(bases=u.cgs.bases)
   Unit("1.000000e+01 g / (cm s2)")
   >>> _.compose(units=u.cgs)
   [Unit("1.000000e+01 Ba")]
