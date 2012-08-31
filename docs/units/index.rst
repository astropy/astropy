***********************
Units (`astropy.units`)
***********************

.. currentmodule:: astropy.units

Introduction
============

``astropy.units`` is a Python package to handle defining and converting
between physical units

Getting Started
===============

  >>> from astropy import units as u
  >>> u.pc.to(u.m, 1)
  3.0856776e+16
  >>> speed_unit = u.cm / u.s
  >>> speed_unit.to(u.mile / u.hour, 1)
  0.02236936292054402
  >>> speed_converter = speed_unit.get_converter("mile hour^-1")
  >>> speed_converter([1., 1000., 5000.])
  array([  2.23693629e-02,   2.23693629e+01,   1.11846815e+02])

Using `astropy.units`
=====================

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
`CompositeUnit`. In most cases, one does not need to worry about the
various kinds of unit classes unless one wants to design a more
complex case.

There are many units already predefined in the module. One may use the
following function to list all the existing predefined units of a
given type::

  >>> u.print_equivalent_units(u.g)
  Primary name | Unit definition | Aliases
  g            | 1.00e-03 kg     | gram
  kg           | irreducible     | kilogram
  lb           | 4.54e-01 kg     | pound
  m_e          | 9.11e-31 kg     |
  m_p          | 1.67e-27 kg     |
  oz           | 2.83e-02 kg     | ounce
  solMass      | 1.99e+30 kg     |
  t            | 1.00e+03 kg     | tonne
  ton          | 9.07e+02 kg     |
  u            | 1.66e-27 kg     |

Composing units
---------------

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

Defining new units
------------------

Converting units to and from strings
------------------------------------

Units can be created from strings using the `Unit` factory function::

  >>> u.Unit("m")
  Unit("m")
  >>> u.Unit("erg / (s cm2)")
  Unit("erg / (s cm2)")

Units can be converted to strings using the `to_string` method::

  >>> fluxunit = u.erg / (u.cm ** 2 * u.s)
  >>> fluxunit.to_string()
  u'erg / (s cm2)'

By default, the string format used is referred to as the "generic"
format, which is based on the FITS standard's format for representing
units, but includes all of the units defined within the
`astropy.units` framework, including user-defined units.  The `Unit`
and `to_string` functions also take an optional `format` parameter to
select a different format.  This parameter may be either a string or a
`astropy.units.format.Base` instance.

`astropy.units` includes support for parsing and writing the following
formats:

  - `"fits"`: This is the format defined in the Units section of the
    `FITS Standard <http://fits.gsfc.nasa.gov/fits_standard.html>`_.
    Unlike the "generic" string format, this will only accept or
    generate units defined in the FITS standard.

  - `"vounit"`: The `proposed IVOA standard
    <http://www.ivoa.net/Documents/VOUnits/>`_ for representing units
    in the VO.  Again, based on the FITS syntax, but the collection of
    supported units is different.

.. These are to-be-implemented

  - OGIP Units: A standard for storing units in `OGIP FITS files
    <http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/>`_.

  - `Standards for astronomical catalogues
    <http://cds.u-strasbg.fr/doc/catstd-3.2.htx>`_: This is the
    standard used, for example, by VOTable versions 1.2 and earlier.

`astropy.units` also is able to write units in the following formats:

  - ``"latex"``: Writes units out using LaTeX math syntax using the
    `IAU Style Manual
    <http://www.iau.org/static/publications/stylemanual1989.pdf>`_
    recommendations for unit presentation.  This format is
    automatically used when printing a unit in the IPython notebook::

      >>> fluxunit

    .. math::

       \mathrm{\frac{erg}{s\ cm^{2}}}

  - ``"console"``: Writes a multi-line representation of the unit
    useful for display in a text console::

      >>> print fluxunit.to_string('console')
       erg
      ------
      s cm^2

  - ``"unicode"``: Same as ``"console"``, except uses Unicode
    characters::

      >>> print u.Ry.decompose().to_string('unicode')
                 m² kg
      2.18×10-¹⁸ ─────
                  s²

Unit Conversion
---------------

There are two ways of handling conversions between units.

Direct Conversion
^^^^^^^^^^^^^^^^^

In this case, one give a unit both the new unit to convert to,
and the value or values to be converted; the value(s) in the new
units is(are) returned.

  >>> u.pc.to(u.m, 3.26)
  1.0059317615e+17

This converts 3.26 parsecs to meters. The first argument is the new
unit desired, the second is the value to be converted.

Arrays are permitted as arguments.

  >>> u.h.to(u.s, [1, 2, 5, 10.1])
  array([  3600.,   7200.,  18000.,  36360.])

Obtaining a Conversion Function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, one may obtain a function that can be used to convert to the
new unit. Normally this may seem like overkill when all one needs to
do is multiply by a scale factor, but there are cases where it is not
quite as simple as multiplying by a scale factor, for example when
there are equivalencies in use.

Conversion to different units involves obtaining a conversion function
and then applying it to the value, or values to be converted.

  >>> speed_unit = u.cm / u.s
  >>> speed_converter = speed_unit.get_converter(u.mile / u.hour)
  >>> speed_converter(100.)
  2.2366936292054402
  >>> speed_converter([1000, 2000])
  array([ 22.36936292,  44.73872584])

Incompatible Conversions
^^^^^^^^^^^^^^^^^^^^^^^^

If you attempt to convert to a incompatible unit, an exception will result:

  >>> speed_unit.to(u.mile)
  ...
  UnitsException: 'cm / (s)' and 'mi' are not convertible

You can check whether a particular conversion is possible using the
`is_equivalent` method::

  >>> u.m.is_equivalent(u.foot)
  True
  >>> u.m.is_equivalent("second")
  False
  >>> (u.m**2).is_equivalent(u.acre)
  True

.. _equivalencies:

Equivalencies
-------------

The unit module has machinery for supporting equivalences between
different units in certain contexts. Namely when equations can
uniquely relate a value in one unit to a different unit. A good
example is the equivalence between wavelength, frequency and energy
for specifying a wavelength of radiation. Normally these units are not
convertable, but when understood as representing light, they are
convertable.  This will describe how to use two of the equivalencies
include in `astropy.units` and then describe how to define new
equivalencies.

Equivalencies are used by passing a list of equivalency pairs to the
`to` or `get_converter` methods.

Spectral Units
^^^^^^^^^^^^^^

The `spectral` (`sp`) is a function that returns an equivalency list
to handle conversions between wavelength, frequency and energy.

Length and frequency are not normally convertible, so
`to` raises an exception::

  >>> u.nm.to(u.Hz, [1000, 2000])
  UnitsException: 'nm' and 'Hz' are not convertible

However, when passing the result of `sp` as the third argument to the
`to`, function, wavelength, frequency and energy can be
converted.

  >>> u.nm.to(u.Hz, [1000, 2000], u.sp())
  array([  2.99792458e+14,   1.49896229e+14])
  >>> u.nm.to(u.eV, [1000, 2000], u.sp())
  array([ 1.23984201,  0.61992101])

Spectral Flux Density Units
^^^^^^^^^^^^^^^^^^^^^^^^^^^

There is also support for spectral flux density units. Their use is
more complex, since it is necessary to also supply the location in the
spectrum for which the conversions will be done, and the units of
those spectral locations. The function that handles these unit
conversions is `spectral_density` (aliased to `sd`).  This function
takes as its arguments the unit and value for the spectral location.
For example::

  >>> u.flam.to(u.fnu, 1, u.sd(u.AA, 3500))
  4.086160166177361e-12
  >>> u.flam.get_converter(u.Jy, u.sd(u.eV, 2.2))(0.0001)
  105941625.20578358

Writing new equivalencies
^^^^^^^^^^^^^^^^^^^^^^^^^

An equivalence list is just a list of tuples, where each tuple has 4
elements::

  (from_unit, to_unit, forward, backward)

`from_unit` and `to_unit` are the equivalent units.  `forward` and
`backward` are functions that convert values between those units.

For example, in an old definition of the metric liter, it was defined
as the volume of 1000 grams of water at 4∘C, a standard pressure and a
bunch of other assumptions. Volumes and masses are not normally
directly convertible, but if we hold those assumptions as true, we
could build an presumably want to build an equivalency for them::

  >>> liters_water = [
         (u.l, u.g, lambda x: 1000.0 * x, lambda x: x / 1000.0)
      ]
  >>> u.l.to(u.kg, 1, liters_water)
  1.0

Note that the equivalency can be used with any other compatible units::

  >>> u.gallon.to(u.pound, 1, liters_water)
  8.345404463333525

Displaying available equivalencies
----------------------------------

The `print_equivalent_units` function also understands equivalencies.
For example, without passing equivalencies, there are no compatible
units for `Hz` in the standard set::

  >>> u.print_equivalent_units(u.Hz)
  Primary name | Unit definition | Aliases
  Hz           | 1 / (s)         | Hertz, hertz

However, when passing the spectral equivalency, you can see there are
all kinds of things that `Hz` can be converted to::

  >>>u.print_equivalent_units(u.Hz, u.sp())
  Primary name | Unit definition       | Aliases
  AA           | 1.00e-10 m            | Angstrom, angstrom
  AU           | 1.50e+11 m            |
  BTU          | 1.06e+03 m2 kg / (s2) | btu
  Hz           | 1 / (s)               | Hertz, hertz
  J            | m2 kg / (s2)          | Joule, joule
  Ry           | 2.18e-18 m2 kg / (s2) | rydberg
  cal          | 4.18e+00 m2 kg / (s2) | calorie
  eV           | 1.60e-19 m2 kg / (s2) | electron_volt
  erg          | 1.00e-07 m2 kg / (s2) |
  ft           | 3.05e-01 m            | foot
  inch         | 2.54e-02 m            |
  kcal         | 4.18e+03 m2 kg / (s2) | Cal, Calorie, kilocal, kilocalorie
  lyr          | 9.46e+15 m            |
  m            | irreducible           | meter
  mi           | 1.61e+03 m            | mile
  micron       | 1.00e-06 m            |
  pc           | 3.09e+16 m            | parsec
  solRad       | 6.96e+08 m            |
  yd           | 9.14e-01 m            | yard

See Also
========

- `FITS Standard <http://fits.gsfc.nasa.gov/fits_standard.html>`_ for
  units in FITS.

- The `proposed IVOA standard
  <http://www.ivoa.net/Documents/VOUnits/>`_ for representing units in
  the VO.

- OGIP Units: A standard for storing units in `OGIP FITS files
  <http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/>`_.

- `Standards for astronomical catalogues units
  <http://cds.u-strasbg.fr/doc/catstd-3.2.htx>`_.

- `IAU Style Manual
  <http://www.iau.org/static/publications/stylemanual1989.pdf>`_.

- `A table of astronomical unit equivalencies
  <http://astro.wku.edu/strolger/UNITS.txt>`_

Reference/API
=============

.. automodapi:: astropy.units.core

.. automodapi:: astropy.units.format

.. currentmodule:: astropy.units.standard_units

.. autosummary::
    m
    kg
    C

.. TODO: Autogenerate the above list

Acknowledgments
===============

astropy.units was adopted from the pynbody units module (with a
multitude of changes; so do not expect it to behave in the same way or
use the same names).
