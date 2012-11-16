.. _unit_equivalencies:

Equivalencies
=============

The unit module has machinery for supporting equivalences between
different units in certain contexts. Namely when equations can
uniquely relate a value in one unit to a different unit. A good
example is the equivalence between wavelength, frequency and energy
for specifying a wavelength of radiation. Normally these units are not
convertible, but when understood as representing light, they are
convertible.  This will describe how to use two of the equivalencies
include in `astropy.units` and then describe how to define new
equivalencies.

Equivalencies are used by passing a list of equivalency pairs to the
`equivs` keyword argument of `~astropy.units.core.UnitBase.to` or
`~astropy.units.core.UnitBase.get_converter` methods.

Built-in equivalencies
----------------------

Spectral Units
^^^^^^^^^^^^^^

`~astropy.units.equivalencies.spectral` is a function that returns an
equivalency list to handle conversions between wavelength, frequency
and energy.

Length and frequency are not normally convertible, so
`~astropy.units.core.UnitBase.to` raises an exception::

  >>> from astropy import units as u
  >>> u.nm.to(u.Hz, [1000, 2000])
  UnitsException: 'nm' (length) and 'Hz' (frequency) are not convertible

However, when passing the result of `~astropy.units.equivalencies.spectral`
as the third argument to the `~astropy.units.core.UnitBase.to` method,
wavelength, frequency and energy can be converted.

  >>> u.nm.to(u.Hz, [1000, 2000], equivs=u.spectral())
  array([  2.99792458e+14,   1.49896229e+14])
  >>> u.nm.to(u.eV, [1000, 2000], equivs=u.spectral())
  array([ 1.23984201,  0.61992101])

These equivalencies even work with non-base units::

  >>> # Inches to calories
  >>> u.inch.to(u.Cal, 1, equivs=u.spectral())
  1.869180759162485e-27

Spectral Flux Density Units
^^^^^^^^^^^^^^^^^^^^^^^^^^^

There is also support for spectral flux density units. Their use is more
complex, since it is necessary to also supply the location in the spectrum for
which the conversions will be done, and the units of those spectral locations.
The function that handles these unit conversions is
`~astropy.units.equivalencies.spectral_density`. This function takes as its
arguments the unit and value for the spectral location. For example::

  >>> u.Jy.to(u.erg / u.cm**2 / u.s / u.Hz, 1., equivs=u.spectral_density(u.AA, 3500))
  1.0000000000000001e-23

  >>> u.Jy.to(u.erg / u.cm**2 / u.s / u.micron, 1., equivs=u.spectral_density(u.AA, 3500))
  2.4472853714285712e-08

Writing new equivalencies
-------------------------

An equivalence list is just a list of tuples, where each tuple has 4
elements::

  (from_unit, to_unit, forward, backward)

`from_unit` and `to_unit` are the equivalent units.  `forward` and
`backward` are functions that convert values between those units.

For example, until 1964 the metric liter was defined as the volume of
1kg of water at 4°C at 760mm mercury pressure.  Volumes and masses are
not normally directly convertible, but if we hold the constants in the
1964 definition of the liter as true, we could build an equivalency
for them::

  >>> liters_water = [
         (u.l, u.g, lambda x: 1000.0 * x, lambda x: x / 1000.0)
      ]
  >>> u.l.to(u.kg, 1, equivs=liters_water)
  1.0

Note that the equivalency can be used with any other compatible units::

  >>> u.gallon.to(u.pound, 1, equivs=liters_water)
  8.345404463333525

And it also works in the other direction::

  >>> u.lb.to(u.pint, 1, equivs=liters_water)
  0.9586114172355458

Displaying available equivalencies
----------------------------------

The `find_equivalent_units` function also understands equivalencies.
For example, without passing equivalencies, there are no compatible
units for `Hz` in the standard set::

  >>> u.Hz.find_equivalent_units()
    Primary name | Unit definition | Aliases
  [
    Hz           | 1 / (s)         | Hertz, hertz ,
  ]

However, when passing the spectral equivalency, you can see there are
all kinds of things that `Hz` can be converted to::

  >>> u.Hz.find_equivalent_units(equivs=u.spectral())
    Primary name | Unit definition           | Aliases
  [
    AU           | 1.495979e+11 m            | au                                 ,
    Angstrom     | 1.000000e-10 m            | AA, angstrom                       ,
    BTU          | 1.055056e+03 kg m2 / (s2) | btu                                ,
    Hz           | 1 / (s)                   | Hertz, hertz                       ,
    J            | kg m2 / (s2)              | Joule, joule                       ,
    Ry           | 2.179872e-18 kg m2 / (s2) | rydberg                            ,
    cal          | 4.184000e+00 kg m2 / (s2) | calorie                            ,
    eV           | 1.602177e-19 kg m2 / (s2) | electronvolt                       ,
    erg          | 1.000000e-07 kg m2 / (s2) |                                    ,
    ft           | 3.048000e-01 m            | foot                               ,
    inch         | 2.540000e-02 m            |                                    ,
    kcal         | 4.184000e+03 kg m2 / (s2) | Cal, Calorie, kilocal, kilocalorie ,
    lyr          | 9.460730e+15 m            |                                    ,
    m            | irreducible               | meter                              ,
    mi           | 1.609344e+03 m            | mile                               ,
    micron       | 1.000000e-06 m            |                                    ,
    pc           | 3.085678e+16 m            | parsec                             ,
    solRad       | 6.955080e+08 m            |                                    ,
    yd           | 9.144000e-01 m            | yard                               ,
  ]
