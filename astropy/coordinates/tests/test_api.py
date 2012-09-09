# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import print_function

from pytest import raises

# the commented unts import should work once units are merged, but for now use
# the temporary module in coordinates
#from ...units import Units as u
from .. import fakeunits as u
from .. import Angle

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
    with raises(TypeError):
        Angle([54, 7, 26.832], unit=u.degree)

    a10 = Angle(3.60827466667, unit=u.hour)
    a11 = Angle("3:36:29.7888000120", unit=u.hour)
    a12 = Angle((3, 36, 29.7888000120), unit=u.hour)  # *must* be a tuple

    Angle(0.944644098745, unit=u.radian)

    with raises(ValueError):
        Angle(54.12412)
        #raises an exception because this is ambiguous

    #ensure the above angles that should match do
    assert a1 == a2 == a3 == a4 == a5 == a6 == a7
    assert a10 == a11 == a12


def test_angle_ops():
    """
    Tests operations on Angle objects
    """

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
    with raises(TypeError):
        a1 / a2

    with raises(TypeError):
        a1 * a2

    a3 = Angle(a1)  # makes a *copy* of the object, but identical content as a1
    assert a1 == a3
    assert a1 is not a3


def test_angle_bounds():
    """
    Tests setting and obeying of bounds for Angle objects, as well as
    how operations interact with bounds
    """
    from .. import RangeError

    '''
    By default the Angle object can accept any value, but will return
    values in [-360,360] (retaining the sign that was specified).

    One can also set artificial bounds for custom applications and range
    checking. As Angle objects are intended to be immutable (the angle value
    anyway), this provides a bound check upon creation. The units given in the
    `bounds` keyword must match the units specified in the scalar.
    '''

    a1 = Angle(13343, unit=u.degree)
    assert a1.degrees == 23

    a2 = Angle(-50, unit=u.degree)
    assert a2.degrees == -50

    a3 = Angle(-361, unit=u.degree)
    assert a3.degrees == -1

    # custom bounds

    with raises(RangeError):
        Angle(66, unit=u.degree, bounds=(-45, 45))

    a4 = Angle(390, unit=u.degree, bounds=(-75, 75))
    assert a4.degrees == 30
    # no RangeError because while 390>75, 30 is within the bounds

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
    assert a7.degrees == 60

    # To get the default bounds back, you need to create a new object with the
    # equivalent angle
    Angle(a4.degrees + a4.degrees, unit=u.degree)

    a9 = Angle(a4.degrees + a5.degrees, unit=u.degree, bounds=[-180, 180])
    assert a9.degrees == 60
    # if they don't match and you want to combine, just re-assign the bounds
    # yourself

    # bounds of None can also be operated on without complaint
    a10 = a6 - a6
    a10.degrees == 0


def test_angle_convert():
    """
    Test unit conversion of Angle objects
    """
    from math import abs

    angle = Angle("54.12412", unit=u.degree)

    assert abs(angle.hours - 3.60827466667) < 1e-11
    assert abs(angle.radians - 0.944644098745) < 1e-12
    assert angle.degrees == 54.12412

    assert isinstance(angle.hms, tuple)
    assert angle.hms[0] == 3
    assert angle.hms[1] == 36
    assert (angle.hms[2] - 29.78879999999947) < 1e-13

    assert isinstance(angle.dms, tuple)
    assert angle.dms[0] == 54
    assert angle.dms[1] == 7
    assert (angle.dms[2] - 26.831999999992036) < 1e-13


def test_angle_formatting():
    """
    Tests string formatting for Angle objects
    """

    '''
    The string method of Angle has this signature:
    def string(self, unit=DEGREE, decimal=False, sep=" ", precision=5,
               pad=False):

    The "decimal" parameter defaults to False since if you need to print the
    Angle as a decimal, there's no need to use the "to_string" method (see
    above).
    '''

    angle = Angle("54.12412", unit=u.degree)

    res = 'Angle as HMS: 3 36 29.78880'
    assert "Angle as HMS: {0}".format(angle.to_string(unit=u.hour)) == res

    res = 'Angle as HMS: 3:36:29.78880'
    print("Angle as HMS: {0}".format(angle.to_string(unit=u.hour, sep=":")))

    res = 'Angle as HMS: 3:36:29.79'
    assert "Angle as HMS: {0}".format(angle.to_string(unit=u.hour, sep=":",
                                      precision=2)) == res

    # Note that you can provide one, two, or three separators passed as a
    # tuple or list

    res = 'Angle as HMS: 3h36m29.7888s'
    assert "Angle as HMS: {0}".format(angle.string(unit=u.hour,
                                                   sep=("h", "m", "s"),
                                                   precision=4)) == res

    res = 'Angle as HMS: 3-36|29.7888'
    assert "Angle as HMS: {0}".format(angle.string(unit=u.hour, sep=["-", "|"],
                                                   precision=4)) == res

    res = 'Angle as HMS: 3-36-29.7888'
    assert "Angle as HMS: {0}".format(angle.string(unit=u.hour, sep="-",
                                                    precision=4)) == res

    res = 'Angle as HMS: 03 36 29.7888'
    assert "Angle as HMS: {0}".format(angle.string(unit=u.hour, precision=4,
                                                  pad=True)) == res


def test_radec():
    """
    Tests creation/operations of RA and Dec objects
    """
    from .. import RA, Dec
    '''
    RA and Dec are objects that are subclassed from Angle. As with Angle, RA
    and Dec can parse any unambiguous format (tuples, formatted strings, etc.).

    The intention is not to create an Angle subclass for every possible
    coordinate object (e.g. galactic l, galactic b). However, equatorial RA/dec
    are so prevalent in astronomy that it's worth creating ones for these
    units. They will be noted as "special" in the docs and use of the just the
    Angle class is to be used for other coordinate systems.
    '''

    with raises(ValueError):
        ra = RA("4:08:15.162342")  # error - hours or degrees?
    ra = RA("26:34:65.345634")  # unambiguous b/c hours don't go past 24

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


def test_create_coordinate():
    """
    Tests creation and basic attributes of Coordinate and subclasses
    """

    '''
    A coordinate marks a position on the sky. This is an object that contains
    two Angle objects. There are a wide array of coordinate systems that should
    be implemented, and it will be easy to subclass to create custom user-made
    coordinates with conversions to standard coordinates.
    '''
    from math import abs

    from .. import RA, Dec, ICRSCoordinate, GalacticCoordinate
    from .. import HorizontalCoordinate
    from ...time import Time

    ra = RA("4:08:15.162342", unit=u.hour)
    dec = Dec("-41:08:15.162342", unit=u.degree)

    # ra and dec are RA and Dec objects, or Angle objects
    c = ICRSCoordinate(ra, dec)
    #both strings below are unambiguous
    c = ICRSCoordinate("54.12412 deg", "-41:08:15.162342")

    assert isinstance(c.dec, Dec)
    # dec is a Dec object
    assert abs(dec.degrees - -41.137545095) < 1e-8

    # It would be convenient to accept both (e.g. ra, dec) coordinates as a
    # single string in the initializer. This can lead to ambiguities
    # particularly when both are different units. The solution is to accept a
    # sequence for the "unit" keyword  that one would expect to sent to Angle.
    # The first element in 'units' refers to the first coordinate.
    # DEGREE is assumed for the second coordinate, unless specified
    c1 = ICRSCoordinate('4 23 43.43  +23 45 12.324', unit=[u.hour])

    # Both can be specified and should be when there is ambiguity.
    c2 = ICRSCoordinate('4 23 43.43  +23 45 12.324', unit=[u.hour, u.degree])

    assert c1 == c2

    # Other types of coordinate systems have their own classes
    l = Angle(123.4)
    b = Angle(76.5)
    c = GalacticCoordinate(l, b)  # only accepts Angle objects *not RA/Dec
    with raises(TypeError):
        GalacticCoordinate(ra, dec)

    assert isinstance(c.l, Angle)  # *not* RA or Dec
    assert isinstance(c.b, Angle)  # *not* RA or Dec

    #some coordinates require an epoch - this is given as an astropy.time.Time
    alt = Angle(20.5)
    az = Angle(45)
    timeobj = Time('J2000')
    HorizontalCoordinate(alt, az, epoch=timeobj)


def test_coord_factory():
    """
    Tests the coordinate factory class.
    """
    from .. import Coordinate, RA, Dec, ICRSCoordinate, GalacticCoordinate

    '''
    To simplify usage, syntax will be provided to figure out the type of
    coordinates the user wants without requiring them to know exactly which
    frame they want.  The coordinate system is determined from the keywords
    used when the object is initialized or can be explicitly specified using
    `a1` and `a2` as angle parameters.
    '''

    c1 = Coordinate(ra="12:43:53", dec=-23, angle1_unit=u.hour,
                   angle2_unit=u.degree)
    # The ra and dec keywords imply equatorial coordinates, which will default
    # to ICRS hence this returns an ICRSCoordinate
    assert isinstance(c1, ICRSCoordinate)

    # l and b are for galactic coordinates, so this returns a
    # GalacticCoordinate object
    c2 = Coordinate(l=158.558650, b=-43.350066, angle1_unit=u.degree,
                    angle2_unit=u.degree)
    assert isinstance(c2, GalacticCoordinate)

    # Any acceptable input for RA() is accepted in Coordinate, etc.
    Coordinate(ra="24:08:15.162342", dec=-41.432345, angle1_unit=u.hour,
               angle2_unit=u.degree)

    # Mismatched keywords produce an error
    with raises(ValueError):
        Coordinate(ra="24:08:15.162342", b=-43.350066, angle1_unit=u.hour,
                   angle2_unit=u.degree)  # error

    # Angle objects also accepted, and thus do not require units
    ra = RA("4:08:15.162342", unit=u.hour)
    dec = Dec("-41:08:15.162342")
    Coordinate(ra=ra, dec=dec)


def test_convert_api():
    """
    Tests the basic coordinate conversion functionality.
    """
    from math import abs

    from .. import Angle, RA, Dec, Coordinate, GalacticCoordinate
    from .. import HorizontalCoordinate, ConvertError, BaseCoordinate

    '''
    Coordinate conversion occurs on-demand internally
    '''

    ra = RA("4:08:15.162342", unit=u.hour)
    dec = Dec("-41:08:15.162342")
    c = Coordinate(ra=ra, dec=dec)

    assert c.galactic == GalacticCoordinate(245.28098, -47.554501)

    #the `galactic` result will be cached to speed this up
    assert isinstance(c.galactic.l, Angle)
    assert abs(c.galactic.l.degree - 158.558650) < 1e-5
    assert isinstance(c.galactic.b, Angle)
    assert abs(c.galactic.b.degree - -43.350066) < 1e-5

    # can also explicitly specify a coordinate class to convert to
    gal = c.convert_to(GalacticCoordinate)

    # can still convert back to equatorial using the shorthand
    assert gal.equatorial.ra.to_string(unit=u.hour, sep=":",
                                       precision=2) == '4:08:15.16'

    with raises(ConvertError):
        # there's no way to convert to alt/az without a specified location
        c.convert_to(HorizontalCoordinate)

    # users can specify their own coordinates and conversions
    class CustomCoordinate(BaseCoordinate):
        coordsysname = 'my_coord'
        #TODO: specify conversion rules

    # allows both ways of converting
    mycoord1 = c.my_coord
    mycoord2 = c.convert_to(CustomCoordinate)

    assert isinstance(mycoord1, CustomCoordinate)
    assert isinstance(mycoord2, CustomCoordinate)


def test_separations():
    """
    Test angular separation functionality
    """

    from .. import Coordinate, BaseCoordinate, ConvertError, AngularSeparation

    '''
    Angular separations between two points on a sphere are supported via the
    `separation` method.
    '''

    c1 = Coordinate(ra=0, dec=0, unit=u.degree)
    c2 = Coordinate(ra=0, dec=1, unit=u.degree)

    sep = c2.separation(c1)
    #returns an AngularSeparation object (a subclass of Angle)
    assert isinstance(sep, AngularSeparation)

    assert sep.degrees == 1
    assert sep.arcmin == 60.

    # these operations have ambiguous interpretations for points on a sphere
    with raises(TypeError):
        c1 + c2
    with raises(TypeError):
        c1 - c2

    c3 = Coordinate(l=0, b=0, unit=u.degree)  # Galactic Coordinates

    # if there is a defined conversion between the relevant coordinate systems,
    # it will be automatically performed to get the right angular separation

    # distance from the north galactic pole to celestial pole
    assert c3.separation(c1).degrees == 62.8716627659

    class CustomCoordinate(BaseCoordinate):
        coordsysname = 'my_coord2'
        # does not specify a coordinate transform

    c4 = CustomCoordinate(0, 0, unit=u.degree)
    with raises(ConvertError):
        # raises an error if no conversion from the custom to equatorial
        # coordinates is available
        c4.separation(c1)


def test_distances():
    """
    Tests functionality for Coordinate class distances and cartesian
    transformations.
    """
    from math import abs
    from .. import Distance, Coordinate, ICRSCoordinates, CartesianPoint
    from ...comology import WMAP5

    '''
    Distances can also be specified, and allow for a full 3D definition of a
    coordinate.
    '''

    distance = Distance(12, u.parsec)
    # need to provide a unit
    with raises(ValueError):
        Distance(12)

    # standard units are pre-defined
    assert distance.light_years == 39.12
    assert abs(distance.km - 3.7e14) < 1e13

    distance.z  # redshift, assuming "current" cosmology
    distance.get_z(WMAP5())  # specifying a cosmology possible

    # Coordinate objects can be assigned a distance object, giving them a full
    # 3D position
    c = Coordinate(l=158.558650, b=-43.350066, unit=u.degree)
    c.distance = Distance(12, u.parsec)

    # Coordinate objects can be initialized with a distance using special
    # syntax
    c1 = Coordinate(l=158.558650, b=-43.350066, unit=u.degree,
                    distance=12 * u.kpc)

    # Coordinate objects can be instantiated with cartesian coordinates
    # Internally they will immediately be converted to two angles + a distance
    c2 = ICRSCoordinates(x=2, y=4, z=8, unit=u.parsec)

    sep12 = c1.separation3d(c2)
    # returns a *3d* distance between the c1 and c2 coordinates, assuming
    # current cosmology if needed
    assert isinstance(sep12, Distance)
    assert sep12.kpc == 0  # TODO: actually put the right number in here

    # can also specify a cosmology
    c1.separation3d(c2, cosmology=WMAP5())

    '''
    All spherical coordinate systems with distances can be converted to
    cartesian coordinates.
    '''

    (x, y, z) = (c2.x, c2.y, c2.z)
    #this only computes the CartesianPoint *once*, and then caches it
    assert x == 2
    assert y == 4
    assert z == 8

    cpt = c2.cartesian
    assert isinstance(cpt, CartesianPoint)
    assert cpt.x == 2
    assert cpt.y == 4
    assert cpt.z == 8

    # returns CartesianPoint object, raise exception if no distance
    with raises(ValueError):
        c3 = Coordinate(l=158.558650, b=-43.350066, unit=u.degree)
        c3.cartesian

    # CartesianPoint objects can be added and subtracted, which are
    # vector/elementwise they can also be given as arguments to a coordinate
    # system
    csum = ICRSCoordinates(c1.cartesian + c2.cartesian)

    assert csum.x == 0  # TODO: fill in correct values
    assert csum.y == 0
    assert csum.z == 0
    assert csum.ra.d == 0
    assert csum.dec.d == 0
    assert csum.distance.pc == 0
