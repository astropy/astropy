.. _doc_standard_units:

Standard units
**************

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
    earthMass    | 5.97217e+24 kg  | M_earth, Mearth                  ,
    g            | 0.001 kg        | gram                             ,
    jupiterMass  | 1.89812e+27 kg  | M_jup, Mjup, M_jupiter, Mjupiter ,
    kg           | irreducible     | kilogram                         ,
    solMass      | 1.98841e+30 kg  | M_sun, Msun                      ,
    t            | 1000 kg         | tonne                            ,
    u            | 1.66054e-27 kg  | Da, Dalton                       ,
  ]


Prefixes
========

Most units can be used with prefixes, with both the standard SI prefixes and
the IEEE 1514 binary prefixes (for ``bit`` and ``byte``) supported:

+------------------------------+
|  Available decimal prefixes  |
+--------+-------------+-------+
| Symbol |    Prefix   | Value |
+========+=============+=======+
|    Y   |    yotta-   |  1e24 |
+--------+-------------+-------+
|    Z   |    zetta-   |  1e21 |
+--------+-------------+-------+
|    E   |     exa-    |  1e18 |
+--------+-------------+-------+
|    P   |    peta-    |  1e15 |
+--------+-------------+-------+
|    T   |    tera-    |  1e12 |
+--------+-------------+-------+
|    G   |    giga-    |  1e9  |
+--------+-------------+-------+
|    M   |    mega-    |  1e6  |
+--------+-------------+-------+
|    k   |    kilo-    |  1e3  |
+--------+-------------+-------+
|    h   |    hecto-   |  1e2  |
+--------+-------------+-------+
|   da   | deka-, deca |  1e1  |
+--------+-------------+-------+
|    d   |    deci-    |  1e-1 |
+--------+-------------+-------+
|    c   |    centi-   |  1e-2 |
+--------+-------------+-------+
|    m   |    milli-   |  1e-3 |
+--------+-------------+-------+
|    u   |    micro-   |  1e-6 |
+--------+-------------+-------+
|    n   |    nano-    |  1e-9 |
+--------+-------------+-------+
|    p   |    pico-    | 1e-12 |
+--------+-------------+-------+
|    f   |    femto-   | 1e-15 |
+--------+-------------+-------+
|    a   |    atto-    | 1e-18 |
+--------+-------------+-------+
|    z   |    zepto-   | 1e-21 |
+--------+-------------+-------+
|    y   |    yocto-   | 1e-24 |
+--------+-------------+-------+

+---------------------------+
| Available binary prefixes |
+--------+--------+---------+
| Symbol | Prefix |  Value  |
+========+========+=========+
|   Ki   |  kibi- | 2 ** 10 |
+--------+--------+---------+
|   Mi   |  mebi- | 2 ** 20 |
+--------+--------+---------+
|   Gi   |  gibi- | 2 ** 30 |
+--------+--------+---------+
|   Ti   |  tebi- | 2 ** 40 |
+--------+--------+---------+
|   Pi   |  pebi- | 2 ** 50 |
+--------+--------+---------+
|   Ei   |  exbi- | 2 ** 60 |
+--------+--------+---------+


.. _doc_dimensionless_unit:

The dimensionless unit
======================

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

As an example of why you might want to create a scaled dimensionless
quantity, say you will be doing many calculations with some big
unitless number, ``big_unitless_num = 20000000  # 20 million``,
but you want all of your answers to be in multiples of a million. This
can be done by simply dividing ``big_unitless_num`` by ``1e6``, but this
requires you to remember that this scaling factor has been applied,
which may be difficult to do after many calculations. Instead, create
a scaled dimensionless quantity by multiplying a value by ``Unit(scale)``
to keep track of the scaling factor, e.g.::

   >>> scale = 1e6
   >>> big_unitless_num = 20 * u.Unit(scale)  # 20 million

   >>> some_measurement = 5.0 * u.cm
   >>> some_measurement * big_unitless_num  # doctest: +FLOAT_CMP
   <Quantity 100. 1e+06 cm>

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
====================

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
      solRad       | 6.957e+08 m     | R_sun, Rsun           ,
      yd           | 0.9144 m        | yard                  ,
    ]


This may also be used with the ``with`` statement, to temporarily
enable additional units::

    >>> from astropy import units as u
    >>> from astropy.units import imperial
    >>> with imperial.enable():
    ...     print(u.m.find_equivalent_units())
          Primary name | Unit definition | Aliases
    ...

To enable just specific units, use `~astropy.units.add_enabled_units`::

    >>> from astropy import units as u
    >>> from astropy.units import imperial
    >>> with u.add_enabled_units([imperial.knot]):
    ...     print(u.m.find_equivalent_units())
          Primary name | Unit definition | Aliases
    ...
