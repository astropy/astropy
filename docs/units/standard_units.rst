Standard units
==============

Standard units are defined in the `astropy.units` package as object
instances.

All units are defined in term of basic 'irreducible' units. The
irreducible units include:

  - Length (meter)
  - Time (second)
  - Mass (kilogram)
  - Current (ampere)
  - Temperature (Kelvin)
  - Angular distance (radian)
  - Solid angle (steradian)
  - Luminous intensity (candela)
  - Stellar magnitude (mag)
  - Amount of substance (mole)
  - Photon count (photon)

(There are also some more obscure base units required by the FITS
standard that are no longer recommended for use.)

Units that involve combinations of fundamental units are instances of
`~astropy.units.core.CompositeUnit`. In most cases, one does not need
to worry about the various kinds of unit classes unless one wants to
design a more complex case.

There are many units already predefined in the module. One may use the
following function to list all the existing predefined units of a
given type::

  >>> from astropy import units as u
  >>> u.g.find_equivalent_units()
    Primary name | Unit definition | Aliases
  [
    M_e          | 9.109383e-31 kg |            ,
    M_p          | 1.672622e-27 kg |            ,
    g            | 1.000000e-03 kg | gram       ,
    kg           | irreducible     | kilogram   ,
    lb           | 4.535924e-01 kg | pound      ,
    oz           | 2.834952e-02 kg | ounce      ,
    solMass      | 1.989100e+30 kg |            ,
    t            | 1.000000e+03 kg | tonne      ,
    ton          | 9.071847e+02 kg |            ,
    u            | 1.660539e-27 kg | Da, Dalton ,
  ]

The dimensionless unit
----------------------

In addition to these units, `astropy.units` includes the concept of
the dimensionless unit, used to indicate quantities that don't have a
physical dimension.  This is distinct in concept from a unit that is
equal to `None`: that indicates that no unit was specified in the data
or by the user.

For convenience, there is a unit that is both dimensionless and
unscaled: the `~astropy.units.dimensionless_unscaled` object::

   >>> from astropy import units as u
   >>> u.dimensionless_unscaled
   Unit(dimensionless)

Dimensionless quantities are often defined as products or ratios of
quantities that are not dimensionless, but whose dimensions cancel out
when their powers are multiplied.  For example::

   >>> u.m / u.m
   Unit(dimensionless)

For compatibility with the supported unit string formats, this is
equivalent to ``Unit('')`` and ``Unit(1)``, though using
``u.dimensionless_unscaled`` in Python code is preferred for
readability::

   >>> u.dimensionless_unscaled == u.Unit('')
   True
   >>> u.dimensionless_unscaled == u.Unit(1)
   True

Note that in many cases, a dimensionless unit may also have a scale.
For example::

   >>> (u.km / u.m).decompose()
   Unit(dimensionless with a scale of 1000.0)
   >>> (u.km / u.m).decompose() == u.dimensionless_unscaled
   False

To determine if a unit is dimensionless (but regardless of the scale),
use the `physical_type` property::

   >>> (u.km / u.m).physical_type
   u'dimensionless'
   # This also has a scale, so it is not the same as u.dimensionless_unscaled
   >>> (u.km / u.m) == u.dimensionless_unscaled
   False
   # However, (u.m / u.m) has a scale of 1.0, so it is the same
   >>> (u.m / u.m) == u.dimensionless_unscaled
   True
