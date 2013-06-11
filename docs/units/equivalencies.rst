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
`equivalencies` keyword argument of `~astropy.units.core.UnitBase.to`
or `~astropy.units.core.UnitBase.get_converter` methods.

Built-in equivalencies
----------------------

Parallax Units
^^^^^^^^^^^^^^
`~astropy.units.equivalencies.parallax` is a function that returns an
equivalency list to handle conversions between angles and length.

Length and angles are not normally convertible, so
`~astropy.units.core.UnitBase.to` raises an exception::

  >>> from astropy import units as u
  >>> u.arcsec.to(u.parsec, 8)
  UnitsException: 'arcsec' (angle) and 'pc' (length) are not convertible

However, when passing the result of `~astropy.units.equivalencies.parallax`
as the third argument to the `~astropy.units.core.UnitBase.to` method,
angles can be converted into units of length (and vice versa).

    >>> u.arcsec.to(u.parsec, 8, equivalencies=u.parallax())
    0.125
    >>> u.AU.to(u.arcminute, 1, equivalencies=u.parallax())
    3437.7467707580054

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

As mentioned above with parallax units, we simply pass the proper conversion
function (in this case `~astropy.units.equivalencies.spectral`) as the third
argument to the `~astropy.units.core.UnitBase.to` method and wavelength,
frequency and energy can be converted.

  >>> u.nm.to(u.Hz, [1000, 2000], equivalencies=u.spectral())
  array([  2.99792458e+14,   1.49896229e+14])
  >>> u.nm.to(u.eV, [1000, 2000], equivalencies=u.spectral())
  array([ 1.23984201,  0.61992101])

These equivalencies even work with non-base units::

  >>> # Inches to calories
  >>> u.inch.to(u.Cal, 1, equivalencies=u.spectral())
  1.869180759162485e-27

Spectral Flux Density Units
^^^^^^^^^^^^^^^^^^^^^^^^^^^

There is also support for spectral flux density units. Their use is more
complex, since it is necessary to also supply the location in the spectrum for
which the conversions will be done, and the units of those spectral locations.
The function that handles these unit conversions is
`~astropy.units.equivalencies.spectral_density`. This function takes as its
arguments the unit and value for the spectral location. For example::

  >>> u.Jy.to(u.erg / u.cm**2 / u.s / u.Hz, 1.,
              equivalencies=u.spectral_density(u.AA, 3500))
  1.0000000000000001e-23

  >>> u.Jy.to(u.erg / u.cm**2 / u.s / u.micron, 1.,
              equivalencies=u.spectral_density(u.AA, 3500))
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
  >>> u.l.to(u.kg, 1, equivalencies=liters_water)
  1.0

Note that the equivalency can be used with any other compatible units::

  >>> u.gallon.to(u.pound, 1, equivalencies=liters_water)
  8.345404463333525

And it also works in the other direction::

  >>> u.lb.to(u.pint, 1, equivalencies=liters_water)
  0.9586114172355458

Writing Spectral (velocity) equivalencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Spectral equivalencies allow you to convert between wavelength, frequency, and
energy, but not to velocity, which is frequently the quantity of interest.

It is fairly straightforward to define the equivalency, but note that there are
different `conventions <http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html>`__.  
In these conventions :math:`f_0` is the rest frequency, :math:`f` is the observed frequency,
:math:`V` is the velocity, and :math:`c` is the speed of light:
        
    * Radio         :math:`V = c \frac{f_0 - f}{f_0}  ;  f(V) = f_0 ( 1 - V/c )`
    * Optical       :math:`V = c \frac{f_0 - f}{f  }  ;   f(V) = f_0 ( 1 + V/c )^{-1}`
    * Redshift      :math:`z =   \frac{f_0 - f}{f  }  ;  f(V) = f_0 ( 1 + z )^{-1}`
    * Relativistic  :math:`V = c (f_0^2 - f^2)/(f_0^2 + f^2) ;  f(V) = f_0 \frac{\left(1 - (V/c)^2\right)^{1/2}}{(1+V/c)}`

To define an equivalency using the radio convention for CO 1-0::

    >>> restfreq = 115.27120  # rest frequency of 12 CO 1-0 in GHz
    >>> ghz_kms = [(u.GHz, u.km/u.s, 
        lambda x: (x-restfreq) / restfreq * c.c.to('km/s').value,
        lambda x: (x/c.c.to('km/s').value) * restfreq + restfreq)]
    >>> u.Hz.to(u.km/u.s,116e9,equivalencies=ghz_kms)
    1895.43219287

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

  >>> u.Hz.find_equivalent_units(equivalencies=u.spectral())
  Primary name | Unit definition        | Aliases
  [
    AU           | 1.49598e+11 m          | au                                 ,
    Angstrom     | 1e-10 m                | AA, angstrom                       ,
    BTU          | 1055.06 kg m2 / s2     | btu                                ,
    Hz           | 1 / s                  | Hertz, hertz                       ,
    J            | kg m2 / s2             | Joule, joule                       ,
    Ry           | 2.17987e-18 kg m2 / s2 | rydberg                            ,
    a            | 3.15576e+07 s          | annum                              ,
    cal          | 4.184 kg m2 / s2       | calorie                            ,
    cm           | 0.01 m                 | centimeter                         ,
    d            | 86400 s                | day                                ,
    eV           | 1.60218e-19 kg m2 / s2 | electronvolt                       ,
    erg          | 1e-07 kg m2 / s2       |                                    ,
    fortnight    | 1.2096e+06 s           |                                    ,
    ft           | 0.3048 m               | foot                               ,
    h            | 3600 s                 | hour, hr                           ,
    inch         | 0.0254 m               |                                    ,
    kcal         | 4184 kg m2 / s2        | Cal, Calorie, kilocal, kilocalorie ,
    lyr          | 9.46073e+15 m          | lightyear                          ,
    m            | irreducible            | meter                              ,
    mi           | 1609.34 m              | mile                               ,
    micron       | 1e-06 m                |                                    ,
    min          | 60 s                   | minute                             ,
    pc           | 3.08568e+16 m          | parsec                             ,
    s            | irreducible            | second                             ,
    sday         | 86164.1 s              |                                    ,
    solRad       | 6.95508e+08 m          | R_sun                              ,
    wk           | 604800 s               | week                               ,
    yd           | 0.9144 m               | yard                               ,
    yr           | 3.15576e+07 s          | year                               ,
  ]
