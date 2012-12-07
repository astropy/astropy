# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import print_function

from numpy import testing as npt
from ...tests.helper import pytest
raises = pytest.raises

from ... import units as u
from ..errors import *


try:
    import scipy
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True

# notes from the original api document:
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


def test_create_angles():
    """
    Tests creating and accessing Angle objects
    """
    from .. import Angle

    ''' The "angle" is a fundamental object. The internal
    representation is stored in radians, but this is transparent to the user.
    Units *must* be specified rather than a default value be assumed. This is
    as much for self-documenting code as anything else.

    Angle objects simply represent a single angular coordinate. More specific
    angular coordinates (e.g. RA, Dec) are subclasses of Angle.'''

    a1 = Angle(54.12412, unit=u.degree)
    a2 = Angle("54.12412", unit=u.degree)
    a3 = Angle("54:07:26.832", unit=u.degree)
    a4 = Angle("54.12412 deg")
    a5 = Angle("54.12412 degrees")
    a6 = Angle("54.12412°") # because we like Unicode
    a7 = Angle((54, 7, 26.832), unit=u.degree)
    # (deg,min,sec) *tuples* are acceptable, but lists/arrays are *not*
    # because of the need to eventually support arrays of coordinates
    with raises(NotImplementedError):
        Angle([54, 7, 26.832], unit=u.degree)

    a10 = Angle(3.60827466667, unit=u.hour)
    a11 = Angle("3:36:29.7888000120", unit=u.hour)
    a12 = Angle((3, 36, 29.7888000120), unit=u.hour)  # *must* be a tuple

    Angle(0.944644098745, unit=u.radian)

    with raises(UnitsError):
        Angle(54.12412)
        #raises an exception because this is ambiguous

    with raises(ValueError):
    	a13 = Angle(12.34, unit="not a unit")

    a14 = Angle("12h43m32") # no trailing 's', but unambiguous

    a15 = Angle("5h4m3s") # single digits, no decimal

    #ensure the above angles that should match do
    a1 == a2 == a3 == a4 == a5 == a6 == a7
    npt.assert_almost_equal(a1.radians, a2.radians)
    npt.assert_almost_equal(a2.degrees, a3.degrees)
    npt.assert_almost_equal(a3.radians, a4.radians)
    npt.assert_almost_equal(a4.radians, a5.radians)
    npt.assert_almost_equal(a5.radians, a6.radians)
    npt.assert_almost_equal(a6.radians, a7.radians)
    #assert a10 == a11 == a12

    # check for illegal ranges / values
    with raises(IllegalSecondError):
        a = Angle("12 32 99", unit=u.degree)

    with raises(IllegalMinuteError):
        a = Angle("12 99 23", unit=u.degree)

    with raises(IllegalSecondError):
        a = Angle("12 32 99", unit=u.hour)

    with raises(IllegalMinuteError):
        a = Angle("12 99 23", unit=u.hour)

    with raises(IllegalHourError):
        a = Angle("99 25 51.0", unit=u.hour)

    with raises(ValueError):
        a = Angle("12 25 51.0xxx", unit=u.hour)

    with raises(ValueError):
        a = Angle("12h34321m32.2s")

    assert a1 is not None

def test_angle_ops():
    """
    Tests operations on Angle objects
    """
    from .. import Angle
    import numpy.testing as npt

    '''
    Angles can be added and subtracted. Multiplication and division by
    a scalar is also permitted. A negative operator is also valid.
    All of these operate in a single dimension. Attempting to
    multiply or divide two Angle objects will raise an exception.
    '''
    a1 = Angle(3.60827466667, unit=u.hour)
    a2 = Angle("54:07:26.832", unit=u.degree)
    a1 + a2  # creates new Angle object
    a1 - a2
    -a1

    # division and multiplication have no unambiguous meaning here
    with raises(NotImplementedError):
        a1 / a2

    with raises(NotImplementedError):
        a1 * a2

    (a1 * 2).hours == 2 * 3.60827466667
    (a1 / 3.123456).hours == 3.60827466667 / 3.123456

    # commutativity
    (2 * a1).hours == (a1 * 2).hours

    a3 = Angle(a1)  # makes a *copy* of the object, but identical content as a1
    npt.assert_almost_equal(a1.radians, a3.radians)
    assert a1 is not a3

    a4 = abs(-a1)
    assert a4.radians == a1.radians

    a5 = Angle(5.0, unit=u.hour)
    assert a5 > a1
    assert a5 >= a1
    assert a1 < a5
    assert a1 <= a5

def test_angle_bounds():
    """
    Tests setting and obeying of bounds for Angle objects, as well as
    how operations interact with bounds
    """
    from .. import Angle, RangeError, BoundsError
    import numpy.testing as npt

    '''
    By default the Angle object can accept any value, but will return
    values in [-360,360] (retaining the sign that was specified).

    One can also set artificial bounds for custom applications and range
    checking. As Angle objects are intended to be immutable (the angle value
    anyway), this provides a bound check upon creation. The units given in the
    `bounds` keyword must match the units specified in the scalar.
    '''

    a1 = Angle(13343, unit=u.degree)
    npt.assert_almost_equal(a1.degrees, 23)

    a2 = Angle(-50, unit=u.degree)
    assert a2.degrees == -50

    a3 = Angle(-361, unit=u.degree)
    npt.assert_almost_equal(a3.degrees, -1)

    # custom bounds

    with raises(BoundsError):
        Angle(66, unit=u.degree, bounds=(-45, 45))

    a4 = Angle(390, unit=u.degree, bounds=(-75, 75))
    npt.assert_almost_equal(a4.degrees, 30)
    # no BoundsError because while 390>75, 30 is within the bounds

    a5 = Angle(390, unit=u.degree, bounds=(-720, 720))
    a5.degrees == 390

    a6 = Angle(1020, unit=u.degree, bounds=None)
    a6.degrees == 1020

    # bounds and operations

    with raises(ValueError):
        a4 + a5
        # ValueError - the bounds don't match

    a7 = a4 + a4
    assert a7.bounds == (-75, 75)
    # if the bounds match, there is no error and the bound is kept
    npt.assert_almost_equal(a7.degrees, 60)

    a8 = a4 - a4
    assert a8.bounds == (-75, 75)
    # To get the default bounds back, you need to create a new object with the
    # equivalent angle
    Angle(a4.degrees + a4.degrees, unit=u.degree)

    a9 = Angle(a4.degrees + a5.degrees, unit=u.degree, bounds=[-180, 180])
    npt.assert_almost_equal(a9.degrees, 60)
    # if they don't match and you want to combine, just re-assign the bounds
    # yourself

    # bounds of None can also be operated on without complaint
    a10 = a6 - a6
    a10.degrees == 0

    with raises(AttributeError):
	    a10.bounds = (0,34)

def test_angle_convert():
    """
    Test unit conversion of Angle objects
    """

    from .. import Angle
    import numpy.testing as npt

    angle = Angle("54.12412", unit=u.degree)

    npt.assert_almost_equal(angle.hours, 3.60827466667)
    npt.assert_almost_equal(angle.radians, 0.944644098745)
    npt.assert_almost_equal(angle.degrees, 54.12412)

    assert isinstance(angle.hms, tuple)
    assert angle.hms[0] == 3
    assert angle.hms[1] == 36
    npt.assert_almost_equal(angle.hms[2], 29.78879999999947)

    assert isinstance(angle.dms, tuple)
    assert angle.dms[0] == 54
    assert angle.dms[1] == 7
    npt.assert_almost_equal(angle.dms[2], 26.831999999992036)

    assert isinstance(angle.dms[0], float)
    assert isinstance(angle.hms[0], float)


def test_angle_formatting():
    """
    Tests string formatting for Angle objects
    """
    from .. import Angle

    '''
    The string method of Angle has this signature:
    def string(self, unit=DEGREE, decimal=False, sep=" ", precision=5,
               pad=False):

    The "decimal" parameter defaults to False since if you need to print the
    Angle as a decimal, there's no need to use the "format" method (see
    above).
    '''

    angle = Angle("54.12412", unit=u.degree)

    #__str__ is the default `format`
    assert str(angle) == angle.format()

    res = 'Angle as HMS: 3h36m29.78880s'
    assert "Angle as HMS: {0}".format(angle.format(unit=u.hour)) == res

    res = 'Angle as HMS: 3:36:29.78880'
    assert "Angle as HMS: {0}".format(angle.format(unit=u.hour, sep=":")) == res

    res = 'Angle as HMS: 3:36:29.79'
    assert "Angle as HMS: {0}".format(angle.format(unit=u.hour, sep=":",
                                      precision=2)) == res

    # Note that you can provide one, two, or three separators passed as a
    # tuple or list

    res = 'Angle as HMS: 3h36m29.7888s'
    assert "Angle as HMS: {0}".format(angle.format(unit=u.hour,
                                                   sep=("h", "m", "s"),
                                                   precision=4)) == res

    res = 'Angle as HMS: 3-36|29.7888'
    assert "Angle as HMS: {0}".format(angle.format(unit=u.hour, sep=["-", "|"],
                                                   precision=4)) == res

    res = 'Angle as HMS: 3-36-29.7888'
    assert "Angle as HMS: {0}".format(angle.format(unit=u.hour, sep="-",
                                                    precision=4)) == res

    res = 'Angle as HMS: 03h36m29.7888s'
    assert "Angle as HMS: {0}".format(angle.format(unit=u.hour, precision=4,
                                                  pad=True)) == res

    # Same as above, in degrees

    angle = Angle("3 36 29.78880", unit=u.degree)

    res = 'Angle as DMS: 3d36m29.78880s'
    assert "Angle as DMS: {0}".format(angle.format(unit=u.degree)) == res

    res = 'Angle as DMS: 3:36:29.78880'
    assert "Angle as DMS: {0}".format(angle.format(unit=u.degree, sep=":")) == res

    res = 'Angle as DMS: 3:36:29.79'
    assert "Angle as DMS: {0}".format(angle.format(unit=u.degree, sep=":",
                                      precision=2)) == res

    # Note that you can provide one, two, or three separators passed as a
    # tuple or list

    res = 'Angle as DMS: 3d36m29.7888s'
    assert "Angle as DMS: {0}".format(angle.format(unit=u.degree,
                                                   sep=("d", "m", "s"),
                                                   precision=4)) == res

    res = 'Angle as DMS: 3-36|29.7888'
    assert "Angle as DMS: {0}".format(angle.format(unit=u.degree, sep=["-", "|"],
                                                   precision=4)) == res

    res = 'Angle as DMS: 3-36-29.7888'
    assert "Angle as DMS: {0}".format(angle.format(unit=u.degree, sep="-",
                                                    precision=4)) == res

    res = 'Angle as DMS: 03d36m29.7888s'
    assert "Angle as DMS: {0}".format(angle.format(unit=u.degree, precision=4,
                                                  pad=True)) == res


    # check negative angles

    angle = Angle(-1.23456789, unit=u.degree)

    assert angle.format() == '-1d14m04.44440s'
    assert angle.format(unit=u.hour) == '-0h04m56.29629s'
    assert angle.format(unit=u.radian, decimal=True) == '-0.021547'

def test_angle_format_roundtripping():
    """
    Ensures that the string represtation of an angle can be used to create a
    new valid Angle.
    """
    from .. import Angle, RA, Dec

    a1 = Angle(0, unit=u.radian)
    a2 = Angle(10, unit=u.degree)
    a3 = Angle(0.543, unit=u.degree)
    a4 = Angle('1d2m3.4s')

    assert Angle(str(a1)).degrees == a1.degrees
    assert Angle(str(a2)).degrees == a2.degrees
    assert Angle(str(a3)).degrees == a3.degrees
    assert Angle(str(a4)).degrees == a4.degrees

    #also check RA/Dec
    ra = RA('1h2m3.4s')
    dec = Dec('1d2m3.4s')

    assert Angle(str(ra)).degrees == ra.degrees
    assert Angle(str(dec)).degrees == dec.degrees


def test_radec():
    """
    Tests creation/operations of RA and Dec objects
    """
    from .. import RA, Dec, Angle
    from ...time import Time

    '''
    RA and Dec are objects that are subclassed from Angle. As with Angle, RA
    and Dec can parse any unambiguous format (tuples, formatted strings, etc.).

    The intention is not to create an Angle subclass for every possible
    coordinate object (e.g. galactic l, galactic b). However, equatorial RA/dec
    are so prevalent in astronomy that it's worth creating ones for these
    units. They will be noted as "special" in the docs and use of the just the
    Angle class is to be used for other coordinate systems.
    '''

    with raises(UnitsError):
        ra = RA("4:08:15.162342")  # error - hours or degrees?
    with raises(RangeError):
        ra = RA("-4:08:15.162342")  # same, should check sign

    ra = RA("26:34:15.345634")  # unambiguous b/c hours don't go past 24
    npt.assert_almost_equal(ra.degrees, 26.570929342)

    with raises(ValueError):
	    ra = RA("garbage containing a d and no units")

    ra = RA(68)

    with raises(UnitsError):
        ra = RA(12)

    ra = RA("12h43m23s")
    npt.assert_almost_equal(ra.hours, 12.7230555556)

    ra = RA((56,14,52.52))		# can accept tuples
    with raises(ValueError):
	    ra = RA((12,14,52)) # ambiguous w/o units
    ra = RA((12,14,52), unit=u.hour)

    with raises(ValueError):
        ra = RA([56,64,52.2])	# ...but not arrays (yet)

    # Units can be specified
    ra = RA("4:08:15.162342", unit=u.hour)

    # Where RA values are commonly found in hours or degrees, declination is
    # nearly always specified in degrees, so this is the default.
    dec = Dec("-41:08:15.162342")
    dec = Dec("-41:08:15.162342", unit=u.degree)  # same as above

    # The RA and Dec objects have bounds hard-coded at (0,360) and (-90,90)
    # degrees, respectively.
    assert ra.bounds == (0, 360)
    with raises(AttributeError):
        ra.bounds = (-45, 45)
    assert dec.bounds == (-90, 90)
    with raises(AttributeError):
        dec.bounds = (-45, 45)


    #RA objects can also compute hour angle and local siderial times
    ra = RA("1:00:00", unit=u.hour)
    ha1 = ra.hour_angle(Angle(1.5, u.hour))
    assert isinstance(ha1, Angle)
    npt.assert_almost_equal(ha1.hours, .5)
    ha2 = ra.hour_angle(Time('2012-1-1 3:00:00', scale='utc'))
    npt.assert_almost_equal(ha2.hours, 23.125)

    lst = ra.lst(Angle(1.5, u.hour))
    assert isinstance(lst, Angle)
    npt.assert_almost_equal(lst.hours, 2.5)

def test_create_coordinate():
    """
    Tests creation and basic attributes of coordinates
    """

    '''
    A coordinate marks a position on the sky. This is an object that contains
    two Angle objects. There are a wide array of coordinate systems that should
    be implemented, and it will be easy to subclass to create custom user-made
    coordinates with conversions to standard coordinates.
    '''

    from .. import Angle, RA, Dec, ICRSCoordinates, GalacticCoordinates
    from .. import HorizontalCoordinates
    import numpy.testing as npt

    ra = RA("4:08:15.162342", unit=u.hour)
    dec = Dec("-41:08:15.162342", unit=u.degree)

    # ra and dec are RA and Dec objects, or Angle objects
    c = ICRSCoordinates(ra, dec)
    assert isinstance(c, ICRSCoordinates)

    c = ICRSCoordinates("54.12412 deg", "-41:08:15.162342 deg")
    assert isinstance(c.dec, Dec) # dec is a Dec object

    npt.assert_almost_equal(dec.degrees, -41.137545095)

    # We should be really robust in what we accept.
    with raises(ValueError):
        c = ICRSCoordinates("12 34 56  -56 23 21") # ambiguous

    with raises(TypeError):
	    c = ICRSCoordinates() # not allowed

    c = ICRSCoordinates(ra="12 43 12", dec=dec, unit=(u.hour, u.hour))

    with raises(TypeError):
        c = ICRSCoordinates(ra="12 43 12", unit=(u.hour,))

    with raises(TypeError):
        c = ICRSCoordinates(ra="12h43m32", b="12:32:43")

    with raises(TypeError):
        c = ICRSCoordinates(ra="12h43m32")

    with raises(TypeError):
        c = ICRSCoordinates(dec="12 32 54")

    # It would be convenient to accept both (e.g. ra, dec) coordinates as a
    # single string in the initializer. This can lead to ambiguities
    # particularly when both are different units. The solution is to accept a
    # sequence for the "unit" keyword  that one would expect to sent to Angle.
    # This must be a 2-sequence with the two elements' units.
    ICRSCoordinates('4 23 43.43  +23 45 12.324', unit=(u.hour, u.degree))

    # If one of them is None, try to guess for unambiguous cases
    ICRSCoordinates('12h43m32 +23:45:12.324', unit=(None, u.degree))

    # unit=None is the same as unit=(None, None)
    ICRSCoordinates('12h43m32 +23d45m12.324s', unit=None)


    # Other types of coordinate systems have their own classes
    l = Angle(123.4, unit=u.degree)
    b = Angle(76.5, unit=u.degree)
    c = GalacticCoordinates(l, b)  # only accepts Angle objects *not RA/Dec
    with raises(TypeError):
        GalacticCoordinates(ra, dec)

    assert isinstance(c.l, Angle)  # *not* RA or Dec
    assert isinstance(c.b, Angle)  # *not* RA or Dec

    #some coordinates require an equinox - this is given as an astropy.time.Time
    from ...time import Time

    alt = Angle(20.5, unit=u.degree)
    az = Angle(45, unit=u.degree)
    timeobj = Time('J2000', scale='utc')
    HorizontalCoordinates(alt, az, equinox=timeobj)

    #some also have an option for an observation time
    ICRSCoordinates(12, 13, unit=(u.hour, u.degree), obstime=timeobj)

    #passing in a non-time object give a TypeError
    with raises(TypeError):
        ICRSCoordinates(12, 13, unit=(u.hour, u.degree), obstime=2000.)



def test_convert_api():
    """
    Tests the basic coordinate conversion functionality.
    """

    from .. import Angle, RA, Dec, ICRSCoordinates, GalacticCoordinates, HorizontalCoordinates
    from ..transformations import coordinate_alias, transform_function, master_transform_graph

    '''
    Coordinate conversion occurs on-demand internally
    '''

    ra = RA("4:08:15.162342", unit=u.hour)
    dec = Dec("-41:08:15.162342")
    c = ICRSCoordinates(ra=ra, dec=dec)

    npt.assert_almost_equal(c.galactic.l.degrees, -114.71902, 5)
    assert isinstance(c.galactic.b, Angle)
    npt.assert_almost_equal(c.galactic.b.degrees, -47.554501, 5)

    # can also explicitly specify a coordinate class to convert to
    gal = c.transform_to(GalacticCoordinates)

    # can still convert back to equatorial using the shorthand
    assert gal.icrs.ra.format(unit=u.hour, sep=":",
                              precision=2) == '4:08:15.16'

    with raises(ConvertError):
        # there's no way to convert to alt/az without a specified location
        c.transform_to(HorizontalCoordinates)

    # users can specify their own coordinates and conversions
    @coordinate_alias('my_coord')
    class CustomCoordinates(ICRSCoordinates):
        pass

    @transform_function(ICRSCoordinates, CustomCoordinates)
    def icrs_to_custom(icrs_coo):
        return CustomCoordinates(icrs_coo.ra.degrees,
                                 icrs_coo.dec.degrees + 2.5,
                                 unit=(u.degree, u.degree))

    try:
        # allows both ways of converting
        mycoord1 = c.my_coord
        mycoord2 = c.transform_to(CustomCoordinates)

        assert isinstance(mycoord1, CustomCoordinates)
        assert isinstance(mycoord2, CustomCoordinates)
    finally:
        #be sure to remove the registered transform
        master_transform_graph._graph[ICRSCoordinates][CustomCoordinates].unregister()


def test_proj_separations():
    """
    Test angular separation functionality
    """

    from .. import ICRSCoordinates, GalacticCoordinates, AngularSeparation, coordinate_alias

    '''
    Angular separations between two points on a sphere are supported via the
    `separation` method.
    '''

    c1 = ICRSCoordinates(ra=0, dec=0, unit=(u.degree, u.degree))
    c2 = ICRSCoordinates(ra=0, dec=1, unit=(u.degree, u.degree))

    sep = c2.separation(c1)
    #returns an AngularSeparation object (a subclass of Angle)
    assert isinstance(sep, AngularSeparation)

    assert sep.degrees == 1
    assert sep.arcmins == 60.

    # these operations have ambiguous interpretations for points on a sphere
    with raises(TypeError):
        c1 + c2
    with raises(TypeError):
        c1 - c2

    ngp = GalacticCoordinates(l=0, b=90, unit=(u.degree, u.degree))
    ncp = ICRSCoordinates(ra=0, dec=90, unit=(u.degree, u.degree))

    # if there is a defined conversion between the relevant coordinate systems,
    # it will be automatically performed to get the right angular separation
    npt.assert_almost_equal(ncp.separation(ngp.icrs).degrees, ncp.separation(ngp).degrees)

    # distance from the north galactic pole to celestial pole
    npt.assert_almost_equal(ncp.separation(ngp.icrs).degrees, 62.8716627659)


    @coordinate_alias('my_coord2')
    class CustomCoordinate(ICRSCoordinates):
        pass  # does not specify a coordinate transform

    c4 = CustomCoordinate(0, 0, unit=(u.degree, u.degree))
    with raises(ConvertError):
        # raises an error if no conversion from the custom to equatorial
        # coordinates is available
        c4.separation(c1)


def test_distances():
    """
    Tests functionality for Coordinate class distances and cartesian
    transformations.
    """
    from .. import Distance, ICRSCoordinates, GalacticCoordinates, CartesianPoints
    from ...cosmology import WMAP5, WMAP7

    '''
    Distances can also be specified, and allow for a full 3D definition of a
    coordinate.
    '''

    #try all the different ways to initialize a Distance
    distance = Distance(12, u.parsec)
    d2 = Distance(40, unit=u.au)
    d3 = Distance(value=5, unit=u.kpc)

    # need to provide a unit
    with raises(UnitsError):
        Distance(12)

    # standard units are pre-defined
    npt.assert_almost_equal(distance.lightyear, 39.13876728075561)
    npt.assert_almost_equal(distance.km, 370281309776063.0)

    # Coordinate objects can be assigned a distance object, giving them a full
    # 3D position
    c = GalacticCoordinates(l=158.558650, b=-43.350066, unit=(u.degree, u.degree))
    c.distance = Distance(12, u.parsec)

    #can also set distances using tuple-format
    c.distance = (12, u.parsec)
    c.distance.parsec = 12

    #or initialize distances via redshifts - this is actually tested in the
    #function below that checks for scipy. This is kept here as an example
    #c.distance = Distance(z=0.2)  # uses current cosmology
    #with whatever your preferred cosmology may be
    #c.distance = Distance(z=0.2, cosmology=WMAP5)


    # Coordinate objects can be initialized with a distance using special
    # syntax
    #TODO: use this commented line once quantity is in
    #c1 = GalacticCoordinates(l=158.558650, b=-43.350066, unit=u.degree, distance=12 * u.kpc)
    c1 = GalacticCoordinates(l=158.558650, b=-43.350066,
        unit=(u.degree, u.degree), distance=Distance(12, u.kpc))

    # Coordinate objects can be instantiated with cartesian coordinates
    # Internally they will immediately be converted to two angles + a distance
    c2 = GalacticCoordinates(x=2, y=4, z=8, unit=u.parsec)

    sep12 = c1.separation_3d(c2)
    # returns a *3d* distance between the c1 and c2 coordinates
    # not that this does *not*
    assert isinstance(sep12, Distance)
    npt.assert_almost_equal(sep12.pc, 12005.784163916317, 10)

    '''
    All spherical coordinate systems with distances can be converted to
    cartesian coordinates.
    '''

    (x, y, z) = (c2.x, c2.y, c2.z)
    #this only computes the CartesianPoints *once*, and then caches it
    npt.assert_almost_equal(x, 2)
    npt.assert_almost_equal(y, 4)
    npt.assert_almost_equal(z, 8)

    cpt = c2.cartesian
    assert isinstance(cpt, CartesianPoints)
    npt.assert_almost_equal(cpt.x, 2)
    npt.assert_almost_equal(cpt.y, 4)
    npt.assert_almost_equal(cpt.z, 8)

    # with no distance, the unit sphere is assumed when converting to cartesian
    c3 = GalacticCoordinates(l=158.558650, b=-43.350066, unit=(u.degree, u.degree), distance=None)
    unitcart = c3.cartesian
    npt.assert_almost_equal((unitcart.x**2 + unitcart.y**2 + unitcart.z**2)**0.5, 1.0)

    # CartesianPoints objects can be added and subtracted, which are
    # vector/elementwise they can also be given as arguments to a coordinate
    # system
    csum = ICRSCoordinates(c1.cartesian + c2.cartesian)

    npt.assert_almost_equal(csum.x, -8.12016610185)
    npt.assert_almost_equal(csum.y, 3.19380597435)
    npt.assert_almost_equal(csum.z, -8.2294483707)
    npt.assert_almost_equal(csum.ra.degrees, 158.529401774)
    npt.assert_almost_equal(csum.dec.degrees, -43.3235825777)
    npt.assert_almost_equal(csum.distance.kpc, 11.9942200501)


@pytest.mark.skipif('not HAS_SCIPY')
def test_distances_scipy():
    """
    The distance-related tests that require scipy due to the cosmology
    module needing scipy integration routines
    """
    from .. import Distance
    from ...cosmology import WMAP5

    #try all the different ways to initialize a Distance
    d4 = Distance(z=0.23)  # uses default cosmology - as of writing, WMAP7
    d5 = Distance(z=0.23, cosmology=WMAP5)

    assert abs(d4.z - 0.23) < 1e-8  # redshift, assuming "current" cosmology
    assert abs(d5.compute_z(WMAP5) - 0.23) < 1e-8 # specifying a cosmology possible