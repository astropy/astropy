#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

from astropy.coordinates import Angle
from astropy.units import Units as u
'''
The emphasis of this package is that common tasks should
be performed in as simple and intuitive manner as possible.
Complex or unusual tasks are still possible, but should not
interfere with the simplicity of common tasks.

One design principle is that objects should silently accept parameters
of any data type that can unambiguously be determined. Where data types
are accepted, they should be specified as constants, not as strings.
If an unknown unit is specified as a string, the code certainly
cannot handle it.

A note on units: Units (as with everything else in Python) are objects.
Any instance of a "unit" keyword will only accept "Unit" objects, not
bare strings. All commonly used units will be predefined in astropy.units
(which is still under development). The exact syntax might change,
but the ideas will not. Units may be specified as part of a string in
an initializer (subscribing to the "accept and parse values that are
unambiguous to a human reader" philosophy), but will not be accepted as
a bare string anywhere else.

On arrays: The interface outlined here only accepts scalar angles (and hence
coordinates) - we will likely later expand it to allow arrays to be stored in
these objects, but only after the scalar interfaces is complete and working.
'''

# ============
# Base Objects
# ============
'''
The "angle" is a fundamental object. The internal representation
is stored in radians, but this is transparent to the user. Units
*must* be specified rather than a default value be assumed. This is
as much for self-documenting code as anything else.

Angle objects simply represent a single angular coordinate. More specific
angular coordinates (e.g. RA, Dec) are subclasses of Angle.

'''
# Creating Angles
# ---------------
angle = Angle(54.12412, unit=u.DEGREE)
angle = Angle("54.12412", unit=u.DEGREE)
angle = Angle("54:07:26.832", unit=u.DEGREE)
angle = Angle("54.12412 deg")
angle = Angle("54.12412 degrees")
angle = Angle("54.12412°") # because we like Unicode

angle = Angle(3.60827466667, unit=u.HOUR)
angle = Angle("3:36:29.7888000120", unit=u.HOUR)

angle = Angle(0.944644098745, unit=u.RADIAN)

angle = Angle(54.12412)
#raises an exception because this is ambiguous

# Angle operations
'''
Angles can be added and subtracted. Multiplication and division by
a scalar is also permitted. A negative operator is also valid.
All of these operate in a single dimension. Attempting to
multiply or divide two Angle objects will raise an exception.
'''
a1 = Angle(3.60827466667, unit=u.HOUR)
a2 = Angle("54:07:26.832", unit=u.DEGREE)
a3 = a1 + a2 # creates new Angle object
a4 = a1 - a2
a5 = -a1

a6 = a1 / a2 # raises an exception
a6 = a1 * a2 # raises an exception

a7 = Angle(a1) # makes a *copy* of the object, but identical content as a1


# Bounds
# ------
'''
By default the Angle object can accept any value, but will return
values in [-360,360] (retaining the sign that was specified).
	
One can also set artificial bounds for custom applications and range checking.
As Angle objects are intended to be immutable (the angle value anyway), this provides a bound check
upon creation. The units given in the "bounds" keyword must match the units
specified in the scalar.
'''
a1 = Angle(13343, unit=u.DEGREE)
print(a1.degrees)
# 23

a2 = Angle(-50, unit=u.DEGREE)
print(a2.degrees)
# -50

a3 = Angle(-361, unit=u.DEGREE)
print (a3.degrees)
# -1

# custom bounds

a4 = Angle(66, unit=u.DEGREE, bounds=(-45,45))
# RangeError

a5 = Angle(390, unit=u.DEGREE, bounds=(-75,75))
print(a5.degrees)
# 30, no RangeError because while 390>75, 30 is within the bounds

a6 = Angle(390, unit=u.DEGREE, bounds=(-720, 720))
print(a6.degrees)
# 390

a7 = Angle(1020, unit=u.DEGREE, bounds=None)
print(a7.degrees)
# 1020

# bounds and operations

a8 = a5 + a6
# ValueError - the bounds don't match

a9 = a5 + a5
print(a9.bounds) 
# (-75, 75) - if the bounds match, there is no error and the bound is kept
print(a9.degrees)
# 60

#To get the default bounds back, you need to create a new object with the equivalent angle
a9 = Angle(a5.degrees + a5.degrees, unit=u.DEGREE) 

a10 = Angle(a5.degrees + a6.degrees, unit=u.DEGREE, bounds=[-180,180])
print(a10.degrees)
# 60 - if they don't match and you want to combine, just re-assign the bounds yourself

a11 = a7 - a7
print(a11.degrees)
# 0 - bounds of None can also be operated on without complaint


# Converting units
# ----------------
angle = Angle("54.12412", unit=u.DEGREE)

print("Angle in hours: {0}.".format(angle.hours))
# Angle in hours: 3.60827466667.

print("Angle in radians: {0}.".format(angle.radians))
# Angle in radians: 0.944644098745.

print("Angle in degrees: {0}.".format(angle.degrees))
# Angle in degrees: 54.12412.

print("Angle in HMS: {0}".format(angle.hms)) # returns a tuple, e.g. (12, 21, 2.343)
# Angle in HMS: (3, 36, 29.78879999999947)

print("Angle in DMS: {0}".format(angle.dms)) # returns a tuple, e.g. (12, 21, 2.343)
# Angle in DMS: (54, 7, 26.831999999992036)


# String formatting
# -----------------
'''
The string method of Angle has this signature:
def string(self, unit=DEGREE, decimal=False, sep=" ", precision=5, pad=False):

The "decimal" parameter defaults to False since if you need to print the
Angle as a decimal, there's no need to use the "to_string" method (see above).
'''
print("Angle as HMS: {0}".format(angle.to_string(unit=u.HOUR)))
# Angle as HMS: 3 36 29.78880

print("Angle as HMS: {0}".format(angle.to_string(unit=u.HOUR, sep=":")))
# Angle as HMS: 3:36:29.78880

print("Angle as HMS: {0}".format(angle.to_string(unit=u.HOUR, sep=":", precision=2)))
# Angle as HMS: 3:36:29.79

# Note that you can provide one, two, or three separators passed as a tuple or list
#
print("Angle as HMS: {0}".format(angle.string(unit=u.HOUR, sep=("h","m","s"), precision=4)))
# Angle as HMS: 3h36m29.7888s

print("Angle as HMS: {0}".format(angle.string(unit=u.HOUR, sep=["-","|"], precision=4)))
# Angle as HMS: 3-36|29.7888

print("Angle as HMS: {0}".format(angle.string(unit=u.HOUR, sep="-", precision=4)))
# Angle as HMS: 3-36-29.7888

# The "pad" parameter will add leading zeros to numbers less than 10.
#
print("Angle as HMS: {0}".format(angle.string(unit=u.HOUR, precision=4, pad=True)))
# Angle as HMS: 03 36 29.7888

# RA/Dec Objects
# --------------
from astropy.coordinates import RA, Dec
'''
RA and Dec are objects that are subclassed from Angle. As with Angle, RA and Dec can
parse any unambiguous format (tuples, formatted strings, etc.).

The intention is not to create an Angle subclass for every possible coordinate object
(e.g. galactic l, galactic b). However, equatorial RA/dec are so prevalent in astronomy
that it's worth creating ones for these units. They will be noted as "special" in the
docs and use of the just the Angle class is to be used for other coordinate systems.
'''
ra = RA("4:08:15.162342") # error - hours or degrees?
ra = RA("26:34:65.345634") # unambiguous

# Units can be specified
ra = RA("4:08:15.162342", unit=u.HOUR)

# Where RA values are commonly found in hours or degrees, declination is nearly always
# specified in degrees, so this is the default.
dec = Dec("-41:08:15.162342")
dec = Dec("-41:08:15.162342", unit=u.DEGREE) # same as above

# The Dec object will have bounds hard-coded at [-90,90] degrees.

# Coordinates
# -----------
'''
A coordinate marks a position on the sky. This is an object that contains two Angle
objects. There are a wide array of coordinate systems that should be implemented, and
it will be easy to subclass to create custom user-made coordinates with conversions to
standard coordinates.
'''
from astropy.coordinates import ICRSCoordinate, GalacticCoordinate, HorizontalCoordinate, Coordinate
from astropy.coordinates import GALACTIC, EQUATORIAL #constants, not classes

# A coordinate in the ICRS standard frame (~= J2000)
c = ICRSCoordinate(ra, dec) #ra and dec are RA and Dec objects, or Angle objects
c = ICRSCoordinate("54.12412 deg", "-41:08:15.162342") #both strings are unambiguous

dec = c.dec
# dec is a Dec object
print(dec.degrees)
# -41.137545095

# It would be convenient to accept both (e.g. ra, dec) coordinates as a single
# string in the initializer. This can lead to ambiguities particularly when both
# are different units. The solution is to accept an array for the "unit" keyword
# that one would expect to sent to Angle:
c = ICRSCoordinate('4 23 43.43  +23 45 12.324', unit=[u.HOUR])
# The first element in 'units' refers to the first coordinate.
# DEGREE is assumed for the second coordinate

c = ICRSCoordinate('4 23 43.43  +23 45 12.324', unit=[u.HOUR, u.DEGREE])
# Both can be specified and should be when there is ambiguity.


# Other types of coordinate systems have their own classes
c = GalacticCoordinates(l, b) #this only accepts Angles, *not* RA and Dec objects

c.l
# this is an Angle object, *not* RA or Dec

#some coordinates require an epoch - the argument will be an astropy.time.Time object
#and the syntax for that is still being ironed out
c = HorizontalCoordinates(alt, az, epoch=timeobj)


# Coordinate Factory
# ------------------
'''
To simplify usage, syntax will be provided to figure out the type of coordinates the user
wants without requiring them to know exactly which frame they want.  The coordinate
system is determined from the keywords used when the object is initialized or can
be explicitly specified using "a1" and "a2" as angle parameters.

Question for the community: Should this be a *class* that overrides __new__ to provide
a new object, or a generator *function*. Erik prefers a function because he thinks it's
less magical to create an object from a function. Demitri prefers an abstract class
(Coordinate) that will return a new object that is subclassed from Coordinate, as
it's a more object-oriented solution and Demitri doesn't like functions in Python.
'''
#Note: `Coordinate` as used below would be `coordinate` if a function were used
c = Coordinate(ra=ra, dec=dec) # equatorial coordinate system - returns ICRSCoordinate object
c = Coordinate(l=158.558650, b=-43.350066, unit=u.DEGREE) # galactic coordinates - returns GalacticCoordinate object

c = Coordinate(a1=139.686111, a2=4.875278, unit=u.DEGREE, system=GALACTIC)
c = Coordinate(a1=139.686111, a2=4.875278, unit=u.DEGREE, system=EQUATORIAL)

# Any acceptable input for RA() is accepted in Coordinate, etc.
c = Coordinate(ra="24:08:15.162342", dec=-41.432345)

# Mismatched keywords produce an error
c = Coordinate(ra="24:08:15.162342", b=-43.350066) # error

# Angle objects also accepted
ra = RA("4:08:15.162342", unit=u.HOUR)
dec = Dec("-41:08:15.162342")
c = Coordinate(ra=ra, dec=dec)



# Coordinate Conversions
# ----------------------
'''
Coordinate conversion occurs on-demand internally
'''
# using the c from above, which is RA,Dec = 4:08:15, -41:08:15.2
print(c.galactic)
# GalacticCoordinate(245.28098,-47.554501)

print(c.galactic.l, c.galactic.b) #the `galactic` result will be cached to speed this up
# Angle(158.558650) Angle(-43.350066)

# can also explicitly specify a coordinate class to convert to
gal = c.convert_to(GalacticCoordinate)

# can still convert back to equatorial using the shorthand
print(gal.equatorial.ra.to_string(unit=u.HOUR, sep=":", precision=2))
# 4:08:15.16

horiz = c.convert_to(HorizontalCoordinate)
# ConvertError - there's no way to convert to alt/az without a specified location


# users can specify their own coordinates and conversions
class CustomCoordinate(BaseCoordinate):
    coordsysname = 'my_coord'
    ... specify conversions somewhere in the class ...

# allows both ways of converting
mycoord1 = c.my_coord
mycoord2 = c.convert_to(CustomCoordinate)


# Separations
# -----------
'''
Angular separations between two points on a sphere are supported via the 
`separation` method.
'''
c1 = Coordinate(ra=0, dec=0, unit=u.DEGREE)
c2 = Coordinate(ra=0, dec=1, unit=u.DEGREE)

sep = c2.separation(c1)
#returns an AngularSeparation object (a subclass of Angle)

print(sep.degrees)
# 1.0
print(sep.arcmin)
# 60.0

c1 + c2
c1 - c2
# TypeError - these operations have ambiguous interpretations for points on a sphere

c3 = Coordinate(l=0, b=0, unit=u.DEGREE) #Galactic Coordinates

# if there is a defined conversion between the relevant coordinate systems, it will
# be automatically performed to get the right angular separation
print c3.separation(c1).degrees
# 62.8716627659 - distance from the north galactic pole to celestial pole

c4 = CustomCoordinate(0, 0, unit=u.DEGREE)
c4.separation(c1) 
# raises an error if no conversion from the custom to equatorial 
# coordinates is available


# Distances
# ---------
'''
Distances can also be specified, and allow for a full 3D definition of the coordinate.
'''
from astropy.coordinates import Distance

d = Distance(12, u.PARSECS)
# need to provide a unit

#standard units are pre-defined
print(distance.light_years)
# 39.12

print(distance.km)
# 3.7e14

print(distance.z) # redshift, assuming "current" cosmology
# (very small for 12 pc)
print(distance(<cosmology>).z) # custom cosmology possible

#Coordinate objects can be assigned a distance object, giving them a full 3D position:
c.distance = Distance(12, u.PARSECS)

# Coordinate objects can be initialized with a distance using special syntax
c1 = Coordinate(l=158.558650, b=-43.350066, unit=u.DEGREE, distance=12 * u.KPC)

# Coordinate objects can be instantiated with cartesian coordinates
# Internally they will immediately be converted to two angles + a distance
c2 = ICRSCoordinates(x=2, y=4, z=8, unit=u.PARSEC)

c1.separation3d(c2, cosmology=<cosmology>) #if cosmology isn't given, use current
# returns a *3d* distance between the c1 and c2 coordinates with units


# Cartesian Coordinates
# ----------------------
'''
All spherical coordinate systems with distances can be converted to
cartesian coordinates.
'''

(x, y, z) = (c.x, c.y, c.z)
 #this only computes the CartesianPoint *once*, and then caches it

c.cartesian
# returns CartesianPoint object, raise exception if no distance was specified

# CartesianPoint objects can be added and subtracted, which are vector/elementwise
# they can also be given as arguments to a coordinate system
csum = ICRSCoordinates(c1.cartesian + c2.cartesian)