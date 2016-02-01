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
`~astropy.units.CompositeUnit`. In most cases, one does not need
to worry about the various kinds of unit classes unless one wants to
design a more complex case.

There are many units already predefined in the module. One may use the
`~astropy.units.core.UnitBase.find_equivalent_units` method to list
all the existing predefined units of a given type::

  >>> from astropy import units as u
  >>> u.g.find_equivalent_units()
    Primary name | Unit definition | Aliases
  [
    M_e          | 9.10938e-31 kg  |                                  ,
    M_p          | 1.67262e-27 kg  |                                  ,
    earthMass    | 5.9742e+24 kg   | M_earth, Mearth                  ,
    g            | 0.001 kg        | gram                             ,
    jupiterMass  | 1.8987e+27 kg   | M_jup, Mjup, M_jupiter, Mjupiter ,
    kg           | irreducible     | kilogram                         ,
    solMass      | 1.9891e+30 kg   | M_sun, Msun                      ,
    t            | 1000 kg         | tonne                            ,
    u            | 1.66054e-27 kg  | Da, Dalton                       ,
  ]

The dimensionless unit
----------------------

In addition to these units, `astropy.units` includes the concept of
the dimensionless unit, used to indicate quantities that don't have a
physical dimension.  This is distinct in concept from a unit that is
equal to `None`: that indicates that no unit was specified in the data
or by the user.

For convenience, there is a unit that is both dimensionless and
unscaled: the ``dimensionless_unscaled`` object::

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
use the `~astropy.units.core.UnitBase.physical_type` property::

   >>> (u.km / u.m).physical_type
   u'dimensionless'
   >>> # This also has a scale, so it is not the same as u.dimensionless_unscaled
   >>> (u.km / u.m) == u.dimensionless_unscaled
   False
   >>> # However, (u.m / u.m) has a scale of 1.0, so it is the same
   >>> (u.m / u.m) == u.dimensionless_unscaled
   True

.. _enabling-other-units:

Enabling other units
--------------------

By default, only the "default" units are searched by
`~astropy.units.core.UnitBase.find_equivalent_units` and similar
methods that do searching.  This includes SI, CGS and astrophysical
units.  However, one may wish to enable the imperial or other
user-defined units.

For example, to enable Imperial units, simply do::

    >>> from astropy.units import imperial
    >>> imperial.enable()  # doctest: +SKIP
    >>> u.m.find_equivalent_units()  # doctest: +SKIP
      Primary name | Unit definition | Aliases
    [
      AU           | 1.49598e+11 m   | au, astronomical_unit ,
      Angstrom     | 1e-10 m         | AA, angstrom          ,
      cm           | 0.01 m          | centimeter            ,
      ft           | 0.3048 m        | foot                  ,
      fur          | 201.168 m       | furlong               ,
      inch         | 0.0254 m        |                       ,
      lyr          | 9.46073e+15 m   | lightyear             ,
      m            | irreducible     | meter                 ,
      mi           | 1609.34 m       | mile                  ,
      micron       | 1e-06 m         |                       ,
      mil          | 2.54e-05 m      | thou                  ,
      nmi          | 1852 m          | nauticalmile, NM      ,
      pc           | 3.08568e+16 m   | parsec                ,
      solRad       | 6.95508e+08 m   | R_sun, Rsun           ,
      yd           | 0.9144 m        | yard                  ,
    ]


This may also be used with the ``with`` statement, to temporarily
enable additional units::

    >>> from astropy import units as u
    >>> from astropy.units import imperial
    >>> with imperial.enable():
    ...     u.m.find_equivalent_units()  # doctest: +SKIP
    ...

To enable just specific units, use `~astropy.units.add_enabled_units`::

    >>> from astropy import units as u
    >>> from astropy.units import imperial
    >>> with u.add_enabled_units_context([imperial.knot]):
    ...     u.m.find_equivalent_units()  # doctest: +SKIP
    ...
