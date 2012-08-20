******************************
A description of astropy.units
******************************

Introduction
============

This is a proposed API for a physical units package. Most, but not
all, of what is proposed here is implemented.

``astropy.units`` is a Python package to handle defining and converting
between physical units

Getting Started
===============

  >>> from astropy import units as u
  >>> speed_unit = u.cm / u.s
  >>> speed_unit.convert_to(u.mile/u.hour, 1)
  0.02236936292054402
  >>> speed_converter = speed_unit.converter_to("mile hour^-1")
  >>> speed_converter([1.,1000.,5000.])
  array([  2.23693629e-02,   2.23693629e+01,   1.11846815e+02])

Using ``astropy.units``
=======================

Standard units are defined in the module as object instances. 
Of course, any redefinition of any the unit instances will break 
unit behavior.

All units are defined in term of basic 'irreducible' units. The following 
irreducible units are currently implemented:

Length (meter)
Time (second)
Mass (kilogram)
Charge (coulomb)

Units that involve combinations of fundamental units are instances of 
CompositeUnits. In most cases, one does not need to worry about the 
various kinds of unit classes unless one wants to design a more complex
case (such as spectral units).

There are many units already predefined in the module. One may use the 
following function to list all the existing predefined units of a given 
type:

  >>> u.list_like(u.g)
  u.list_like(u.g)
  Primary name | Unit definition | Aliases
  Msol            1.99e+30 kg     []
  m_p             1.67e-27 kg     []
  kg              irreducible     ['kilogram']
  g               1.00e-03 kg     ['gram']
  m_e             9.11e-31 kg     []
  
Examples of defining new units
------------------------------

  >>> from astropy import units as u
  >>> speed_unit = u.cm / u.s
  >>> speed_unit = u.unit("cm s^-1")
  >>> fluxunit = u.unit("erg cm^-2 s^-1")
  >>> fluxunit = u.Unit("bozos", erg/cm**2/s)
  
Note the last example give a name 'bozos' to the unit that can be referred
to by other machinery (currently doesn't work outside the module package).

Unit Conversion
---------------

There are three ways of handling conversions between units

Direct Conversion
^^^^^^^^^^^^^^^^^

In this case, one give a unit both the new unit to convert to, 
and the value or values to be converted; the value(s) in the new
units is(are) returned.

  >>> u.pc.convert_to(u.m, 3.26)
  1.0059317615e+17
  
This converts 3.26 parsecs to meters. The first argument is the new unit
desired, the second is the value to be converted.

Arrays are permitted as arguments.

  >>> u.hr.convert_to(u.s, [1,2,5,10.1])
  array([  3600.,   7200.,  18000.,  36360.])

Obtaining The Conversion Scale Factor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One can obtain the needed scale factor for converting between two units
instead.

 >>> u.hr.scale_to(u.s)
 3600
 
Obtaining a Conversion Function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, one may obtain a function that can be used to convert to the 
new unit. Normally this may seem like overkill when all one needs to 
do is multiply by a scale factor, but there are cases where it is not
quite as simple as multiplying by a scale factor. (e.g., see 
the section on spectral units.)

Conversion to different units involves obtaining a conversion function
and then applying it to the value, or values to be converted.

  >>> speed_unit = u.cm/u.s
  >>> speed_converter = speed_unit.converter_to(u.mile/u.hour)
  >>> speed_converter(100.)
  2.2366936292054402
  >>> speed_converter([1000,2000])
  array([ 22.36936292,  44.73872584])

Incompatible Conversions
^^^^^^^^^^^^^^^^^^^^^^^^

If you attempt to convert to a incompatible unit, an exception will result:

  >>> speed_unit.scale_to(u.miles)
  ...
  UnitsException: Not convertible
  

Users are free to define new units, either fundamental or compound

e.g.:

  >>> bakers_fortnight = u.Unit('bakers_fortnight',13 * u.day) 

The addition of a string gives the new unit a name that will show up when
the unit is printed.


Creating a new fundamental unit is simple
 
  >>> titter = u.IrreducibleUnit('titter')
  >>> chuckle = u.Unit('chuckle',5 * titter)
  >>> laugh = u.Unit('laugh',4 * chuckle)
  >>> guffaw = u.Unit('guffaw',3 * laugh)
  >>> rofl = u.Unit('rofl',4 * guffaw)
  >>> death_by_laughing = u.Unit('death_by_laughing',10 * rofl)
  >>> rofl.scale_to(titter)
  240

Using strings to define units
-----------------------------

Units may be specified or combined by their string representations.
For example:

  >>> mile = u.Unit("mile") # same as mile = u.mile
  >>> speed = u.Unit("mile hour**-1") # same as speed = u.mile/u.hour

Checking for unit consistency (not implemented yet)
---------------------------------------------------

  >>> u.m.consistent_with(u.foot)
  True
  >>> u.m.consistent_with("second")
  False
  >>> (u.m**2).consistent_with(u.area_unit)
  True

Equivalence Units
-----------------

The unit module has machinery for supporting equivalences between 
different units in certain contexts. Namely when equations can 
uniquely relate a value in one unit to a different unit. A good
example is the equivalence between wavelength, frequency and energy
for specifying a wavelength of radiation. Normally these units are
not convertable, but when understood as representing light, they
are convertable. This won't describe the means of adding new kinds
of such units, but will describe using two cases already implemented.

Spectral Units
^^^^^^^^^^^^^^

There is a special unit class called SpectralUnit that handles unit
equivalences between wavelength, frequency, and energy. The unit module
defines special spectral units in these terms, all of which have 'sp_'
prepended, e.g., sp_nm for spectral nanometers. These units can be 
converted to other forms. Examples are the easiest way to show how it
works.

  >>> u.sp_nm.convert_to(u.sp_Hz, [1000, 2000])
  array([  2.99792458e+14,   1.49896229e+14])
  >>> u.sp_nm.convert_to(u.sp_eV, [1000,2000])
  array([ 1.23984201,  0.61992101])
  

Note that one can convert to an ordinary unit

  >>> u.sp_nm.converter_to(u.Hz)(1000)
  299792457999999.94
  
Unlike ordinary units, one cannot form composite units with these
(other than applying simple scaling factors). For example:

  >>> u.sp_nm * u.m
  TypeError: can only multiply by scalars
  
Although these units cannot be combined with other units, they are intended
to simplify the user interface for any module that involves spectral units.
It allow the user to use any one of the equivalent representations without
the module author being burdened with the bookkeeping of checking which unit
the user has has supplied and doing an explicit conversion for each 
alternative, and similarly, handing output units in the form the user wishes.


Spectral Flux Density Units
^^^^^^^^^^^^^^^^^^^^^^^^^^^

There is also support for Spectral Flux Density Units. Their use is
more complex, since it is necessary to also supply the location in the 
spectrum for which the conversions will be done, and the units of those
spectral locations. The class that handles this unit is SpectralDensityUnit
and all the predefined units of this type are prefixed with 'sd_'

  >>> u.sd_flam.convert_to(u.fnu,1,u.sp_A,3500)
  4.086160166177361e-12
  >>> u.sd_flam.converter_to(u.Jy)(0.0001,u.sp_eV,2.2)
  105941625.20578358
  

Acknowledgments
===============

astropy.units was adopted from the pynbody units module (with a number of
changes; so do not expect it to behave in the same way or use
the same names for everything)