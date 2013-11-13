# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from numpy import testing as npt
from ...tests.helper import pytest
raises = pytest.raises

from ...extern import six

from ... import units as u
from ..errors import ConvertError, IllegalSecondError, IllegalMinuteError, IllegalHourError


try:
    import scipy  # pylint: disable=W0611
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
    angular coordinates (e.g. Longitude, Latitude) are subclasses of Angle.'''

    a1 = Angle(54.12412, unit=u.degree)
    a2 = Angle("54.12412", unit=u.degree)
    a3 = Angle("54:07:26.832", unit=u.degree)
    a4 = Angle("54.12412 deg")
    a5 = Angle("54.12412 degrees")
    a6 = Angle("54.12412°") # because we like Unicode
    a7 = Angle((54, 7, 26.832), unit=u.degree)
    a8 = Angle("54°07'26.832\"")
    # (deg,min,sec) *tuples* are acceptable, but lists/arrays are *not*
    # because of the need to eventually support arrays of coordinates
    a9 = Angle([54, 7, 26.832], unit=u.degree)
    npt.assert_allclose(a9.value, [54, 7, 26.832])
    assert a9.unit is u.degree

    a10 = Angle(3.60827466667, unit=u.hour)
    a11 = Angle("3:36:29.7888000120", unit=u.hour)
    a12 = Angle((3, 36, 29.7888000120), unit=u.hour)  # *must* be a tuple

    Angle(0.944644098745, unit=u.radian)

    with raises(u.UnitsError):
        Angle(54.12412)
        #raises an exception because this is ambiguous

    with raises(ValueError):
        a13 = Angle(12.34, unit="not a unit")

    a14 = Angle("12h43m32") # no trailing 's', but unambiguous

    a15 = Angle("5h4m3s") # single digits, no decimal

    a16 = Angle("1 d")
    a17 = Angle("1 degree")
    assert a16.degree == 1
    assert a17.degree == 1

    #ensure the above angles that should match do
    assert a1 == a2 == a3 == a4 == a5 == a6 == a7
    npt.assert_allclose(a1.radian, a2.radian)
    npt.assert_allclose(a2.degree, a3.degree)
    npt.assert_allclose(a3.radian, a4.radian)
    npt.assert_allclose(a4.radian, a5.radian)
    npt.assert_allclose(a5.radian, a6.radian)
    npt.assert_allclose(a6.radian, a7.radian)
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

    # Angles can be added and subtracted. Multiplication and division by a
    # scalar is also permitted. A negative operator is also valid.  All of
    # these operate in a single dimension. Attempting to multiply or divide two
    # Angle objects will raise an exception.

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

    npt.assert_allclose((a1 * 2).hour, 2 * 3.6082746666700003)
    assert abs((a1 / 3.123456).hour - 3.60827466667 / 3.123456) < 1e-10

    # commutativity
    assert (2 * a1).hour == (a1 * 2).hour

    a3 = Angle(a1)  # makes a *copy* of the object, but identical content as a1
    npt.assert_allclose(a1.radian, a3.radian)
    assert a1 is not a3

    a4 = abs(-a1)
    assert a4.radian == a1.radian

    a5 = Angle(5.0, unit=u.hour)
    assert a5 > a1
    assert a5 >= a1
    assert a1 < a5
    assert a1 <= a5


def test_angle_convert():
    """
    Test unit conversion of Angle objects
    """

    from .. import Angle

    angle = Angle("54.12412", unit=u.degree)

    npt.assert_allclose(angle.hour, 3.60827466667)
    npt.assert_allclose(angle.radian, 0.944644098745)
    npt.assert_allclose(angle.degree, 54.12412)

    assert len(angle.hms) == 3
    assert isinstance(angle.hms, tuple)
    assert angle.hms[0] == 3
    assert angle.hms[1] == 36
    npt.assert_allclose(angle.hms[2], 29.78879999999947)

    assert len(angle.dms) == 3
    assert isinstance(angle.dms, tuple)
    assert angle.dms[0] == 54
    assert angle.dms[1] == 7
    npt.assert_allclose(angle.dms[2], 26.831999999992036)

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
    assert str(angle) == angle.to_string()

    res = 'Angle as HMS: 3h36m29.7888s'
    assert "Angle as HMS: {0}".format(angle.to_string(unit=u.hour)) == res

    res = 'Angle as HMS: 3:36:29.7888'
    assert "Angle as HMS: {0}".format(angle.to_string(unit=u.hour, sep=":")) == res

    res = 'Angle as HMS: 3:36:29.79'
    assert "Angle as HMS: {0}".format(angle.to_string(unit=u.hour, sep=":",
                                      precision=2)) == res

    # Note that you can provide one, two, or three separators passed as a
    # tuple or list

    res = 'Angle as HMS: 3h36m29.7888s'
    assert "Angle as HMS: {0}".format(angle.to_string(unit=u.hour,
                                                   sep=("h", "m", "s"),
                                                   precision=4)) == res

    res = 'Angle as HMS: 3-36|29.7888'
    assert "Angle as HMS: {0}".format(angle.to_string(unit=u.hour, sep=["-", "|"],
                                                   precision=4)) == res

    res = 'Angle as HMS: 3-36-29.7888'
    assert "Angle as HMS: {0}".format(angle.to_string(unit=u.hour, sep="-",
                                                    precision=4)) == res

    res = 'Angle as HMS: 03h36m29.7888s'
    assert "Angle as HMS: {0}".format(angle.to_string(unit=u.hour, precision=4,
                                                  pad=True)) == res

    # Same as above, in degrees

    angle = Angle("3 36 29.78880", unit=u.degree)

    res = 'Angle as DMS: 3d36m29.7888s'
    assert "Angle as DMS: {0}".format(angle.to_string(unit=u.degree)) == res

    res = 'Angle as DMS: 3:36:29.7888'
    assert "Angle as DMS: {0}".format(angle.to_string(unit=u.degree, sep=":")) == res

    res = 'Angle as DMS: 3:36:29.79'
    assert "Angle as DMS: {0}".format(angle.to_string(unit=u.degree, sep=":",
                                      precision=2)) == res

    # Note that you can provide one, two, or three separators passed as a
    # tuple or list

    res = 'Angle as DMS: 3d36m29.7888s'
    assert "Angle as DMS: {0}".format(angle.to_string(unit=u.degree,
                                                   sep=("d", "m", "s"),
                                                   precision=4)) == res

    res = 'Angle as DMS: 3-36|29.7888'
    assert "Angle as DMS: {0}".format(angle.to_string(unit=u.degree, sep=["-", "|"],
                                                   precision=4)) == res

    res = 'Angle as DMS: 3-36-29.7888'
    assert "Angle as DMS: {0}".format(angle.to_string(unit=u.degree, sep="-",
                                                    precision=4)) == res

    res = 'Angle as DMS: 03d36m29.7888s'
    assert "Angle as DMS: {0}".format(angle.to_string(unit=u.degree, precision=4,
                                                  pad=True)) == res

    res = 'Angle as rad: 0.0629763rad'
    assert "Angle as rad: {0}".format(angle.to_string(unit=u.radian)) == res

    res = 'Angle as rad decimal: 0.0629763'
    assert "Angle as rad decimal: {0}".format(angle.to_string(unit=u.radian, decimal=True)) == res


    # check negative angles

    angle = Angle(-1.23456789, unit=u.degree)
    angle2 = Angle(-1.23456789, unit=u.hour)

    assert angle.to_string() == '-1d14m04.4444s'
    assert angle.to_string(pad=True) == '-01d14m04.4444s'
    assert angle.to_string(unit=u.hour) == '-0h04m56.2963s'
    assert angle2.to_string(unit=u.hour, pad=True) == '-01h14m04.4444s'
    assert angle.to_string(unit=u.radian, decimal=True) == '-0.0215473'

def test_angle_format_roundtripping():
    """
    Ensures that the string represtation of an angle can be used to create a
    new valid Angle.
    """
    from .. import Angle, Longitude, Latitude

    a1 = Angle(0, unit=u.radian)
    a2 = Angle(10, unit=u.degree)
    a3 = Angle(0.543, unit=u.degree)
    a4 = Angle('1d2m3.4s')

    assert Angle(str(a1)).degree == a1.degree
    assert Angle(str(a2)).degree == a2.degree
    assert Angle(str(a3)).degree == a3.degree
    assert Angle(str(a4)).degree == a4.degree

    #also check Longitude/Latitude
    ra = Longitude('1h2m3.4s')
    dec = Latitude('1d2m3.4s')

    npt.assert_allclose(Angle(str(ra)).degree, ra.degree)
    npt.assert_allclose(Angle(str(dec)).degree, dec.degree)


def test_radec():
    """
    Tests creation/operations of Longitude and Latitude objects
    """
    from .. import Longitude, Latitude, Angle
    from ...time import Time

    '''
    Longitude and Latitude are objects that are subclassed from Angle. As with Angle, Longitude
    and Latitude can parse any unambiguous format (tuples, formatted strings, etc.).

    The intention is not to create an Angle subclass for every possible
    coordinate object (e.g. galactic l, galactic b). However, equatorial Longitude/Latitude
    are so prevalent in astronomy that it's worth creating ones for these
    units. They will be noted as "special" in the docs and use of the just the
    Angle class is to be used for other coordinate systems.
    '''

    with raises(u.UnitsError):
        ra = Longitude("4:08:15.162342")  # error - hours or degrees?
    with raises(u.UnitsError):
        ra = Longitude("-4:08:15.162342")

    # the "smart" initializer allows >24 to automatically do degrees, but the
    #Angle-based one does not
    #TODO: adjust in 0.3 for whatever behavior is decided on

    #ra = Longitude("26:34:15.345634")  # unambiguous b/c hours don't go past 24
    #npt.assert_allclose(ra.degree, 26.570929342)
    with raises(u.UnitsError):
        ra = Longitude("26:34:15.345634")

    #ra = Longitude(68)
    with raises(u.UnitsError):
        ra = Longitude(68)

    with raises(u.UnitsError):
        ra = Longitude(12)

    with raises(ValueError):
        ra = Longitude("garbage containing a d and no units")

    ra = Longitude("12h43m23s")
    npt.assert_allclose(ra.hour, 12.7230555556)

    ra = Longitude((56, 14, 52.52), unit=u.degree)      # can accept tuples
    #TODO: again, fix based on >24 behavior
    #ra = Longitude((56,14,52.52))
    with raises(u.UnitsError):
        ra = Longitude((56, 14, 52.52))
    with raises(u.UnitsError):
        ra = Longitude((12, 14, 52))  # ambiguous w/o units
    ra = Longitude((12, 14, 52), unit=u.hour)

    ra = Longitude([56, 64, 52.2], unit=u.degree)  # ...but not arrays (yet)

    # Units can be specified
    ra = Longitude("4:08:15.162342", unit=u.hour)

    #TODO: this was the "smart" initializer behavior - adjust in 0.3 appropriately
    ## Where Longitude values are commonly found in hours or degrees, declination is
    ## nearly always specified in degrees, so this is the default.
    #dec = Latitude("-41:08:15.162342")
    with raises(u.UnitsError):
        dec = Latitude("-41:08:15.162342")
    dec = Latitude("-41:08:15.162342", unit=u.degree)  # same as above


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

    from .. import Angle, Longitude, Latitude, ICRS, Galactic
    from .. import AltAz

    ra = Longitude("4:08:15.162342", unit=u.hour)
    dec = Latitude("-41:08:15.162342", unit=u.degree)

    # ra and dec are Longitude and Latitude objects, or Angle objects
    c = ICRS(ra, dec)
    c = ICRS(Angle(4.137545095, u.hour), Angle(-41.137545095, u.degree))

    c = ICRS("54.12412 deg", "-41:08:15.162342 deg")
    assert isinstance(c.dec, Latitude)  # dec is a Latitude object

    #make sure ICRS has a working repr
    repr(c)

    npt.assert_allclose(dec.degree, -41.137545095)

    # We should be really robust in what we accept.
    with raises(u.UnitsError):
        c = ICRS("12 34 56  -56 23 21") # ambiguous

    with raises(TypeError):
        c = ICRS() # not allowed

    c = ICRS(ra="12 43 12", dec=dec, unit=(u.hour, u.hour))

    with raises(TypeError):
        c = ICRS(ra="12 43 12", unit=(u.hour,))

    with raises(TypeError):
        c = ICRS(ra="12h43m32", b="12:32:43")

    with raises(TypeError):
        c = ICRS(ra="12h43m32")

    with raises(TypeError):
        c = ICRS(dec="12 32 54")

    # It would be convenient to accept both (e.g. ra, dec) coordinates as a
    # single string in the initializer. This can lead to ambiguities
    # particularly when both are different units. The solution is to accept a
    # sequence for the "unit" keyword  that one would expect to sent to Angle.
    # This must be a 2-sequence with the two elements' units.
    ICRS('4 23 43.43  +23 45 12.324', unit=(u.hour, u.degree))

    # If one of them is None, try to guess for unambiguous cases
    ICRS('12h43m32 +23:45:12.324', unit=(None, u.degree))

    # unit=None is the same as unit=(None, None)
    ICRS('12h43m32 +23d45m12.324s', unit=None)


    # Other types of coordinate systems have their own classes
    l = Angle(123.4, unit=u.degree)
    b = Angle(76.5, unit=u.degree)
    c = Galactic(l, b)  # accepts Angle objects and Longitude/Latitude
    d = Galactic(ra, dec)

    #make sure Galactic has a working repr
    repr(c)

    assert isinstance(c.l, Angle)  # *not* Longitude or Latitude
    assert isinstance(c.b, Angle)  # *not* Longitude or Latitude

    #some coordinates require an equinox - this is given as an astropy.time.Time
    from ...time import Time

    alt = Angle(20.5, unit=u.degree)
    az = Angle(45, unit=u.degree)
    timeobj = Time('J2000', scale='utc')
    AltAz(alt, az, equinox=timeobj)

    #some also have an option for an observation time
    ICRS(12, 13, unit=(u.hour, u.degree), obstime=timeobj)

    #passing in a non-time object give a TypeError
    with raises(TypeError):
        ICRS(12, 13, unit=(u.hour, u.degree), obstime=2000.)

    #should also be able to initialize coordinates from other coord objects
    ICRS(ICRS(ra, dec))
    newi = ICRS(c) # including from other types like Galactic
    # and it will auto-convert
    assert isinstance(newi, ICRS)
    assert newi.ra != c.l
    assert newi.dec != c.b


def test_convert_api():
    """
    Tests the basic coordinate conversion functionality.
    """

    from .. import Angle, Longitude, Latitude, ICRS, Galactic, AltAz
    from ..transformations import coordinate_alias, transform_function, master_transform_graph

    '''
    Coordinate conversion occurs on-demand internally
    '''

    ra = Longitude("4:08:15.162342", unit=u.hour)
    dec = Latitude("-41:08:15.162342", unit=u.degree)
    c = ICRS(ra=ra, dec=dec)

    npt.assert_allclose(c.galactic.l.degree, 245.28098, 5)
    assert isinstance(c.galactic.b, Angle)
    npt.assert_allclose(c.galactic.b.degree, -47.554501, 5)

    # can also explicitly specify a coordinate class to convert to
    gal = c.transform_to(Galactic)

    # can still convert back to equatorial using the shorthand
    assert gal.icrs.ra.to_string(unit=u.hour, sep=":",
                                 precision=2) == '4:08:15.16'

    with raises(ConvertError):
        # there's no way to convert to alt/az without a specified location
        c.transform_to(AltAz)

    # users can specify their own coordinates and conversions
    @coordinate_alias('my_coord')
    class CustomCoordinates(ICRS):
        pass

    @transform_function(ICRS, CustomCoordinates)
    def icrs_to_custom(icrs_coo):
        return CustomCoordinates(icrs_coo.ra.degree,
                                 icrs_coo.dec.degree + 2.5,
                                 unit=(u.degree, u.degree))

    try:
        # allows both ways of converting
        mycoord1 = c.my_coord
        mycoord2 = c.transform_to(CustomCoordinates)

        assert isinstance(mycoord1, CustomCoordinates)
        assert isinstance(mycoord2, CustomCoordinates)
    finally:
        #be sure to remove the registered transform
        master_transform_graph._graph[ICRS][CustomCoordinates].unregister()


def test_proj_separations():
    """
    Test angular separation functionality
    """

    from .. import ICRS, Galactic, coordinate_alias, Angle

    '''
    Angular separations between two points on a sphere are supported via the
    `separation` method.
    '''

    c1 = ICRS(ra=0, dec=0, unit=(u.degree, u.degree))
    c2 = ICRS(ra=0, dec=1, unit=(u.degree, u.degree))

    sep = c2.separation(c1)
    #returns an Angle object
    assert isinstance(sep, Angle)

    assert sep.degree == 1
    npt.assert_allclose(sep.arcminute, 60.)

    # these operations have ambiguous interpretations for points on a sphere
    with raises(TypeError):
        c1 + c2
    with raises(TypeError):
        c1 - c2

    ngp = Galactic(l=0, b=90, unit=(u.degree, u.degree))
    ncp = ICRS(ra=0, dec=90, unit=(u.degree, u.degree))

    # if there is a defined conversion between the relevant coordinate systems,
    # it will be automatically performed to get the right angular separation
    npt.assert_allclose(ncp.separation(ngp.icrs).degree, ncp.separation(ngp).degree)

    # distance from the north galactic pole to celestial pole
    npt.assert_allclose(ncp.separation(ngp.icrs).degree, 62.8716627659)


    @coordinate_alias('my_coord2')
    class CustomCoordinate(ICRS):
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
    from .. import Distance, ICRS, Galactic, CartesianPoints
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
    with raises(u.UnitsError):
        Distance(12)

    # standard units are pre-defined
    npt.assert_allclose(distance.lyr, 39.138765325702551)
    npt.assert_allclose(distance.km, 370281309776063.0)

    # Coordinate objects can be assigned a distance object, giving them a full
    # 3D position
    c = Galactic(l=158.558650, b=-43.350066, unit=(u.degree, u.degree))
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
    #c1 = Galactic(l=158.558650, b=-43.350066, unit=u.degree, distance=12 * u.kpc)
    c1 = Galactic(l=158.558650, b=-43.350066,
        unit=(u.degree, u.degree), distance=Distance(12, u.kpc))

    # Coordinate objects can be instantiated with cartesian coordinates
    # Internally they will immediately be converted to two angles + a distance
    c2 = Galactic(x=2, y=4, z=8, unit=u.parsec)

    sep12 = c1.separation_3d(c2)
    # returns a *3d* distance between the c1 and c2 coordinates
    # not that this does *not*
    assert isinstance(sep12, Distance)
    npt.assert_allclose(sep12.pc, 12005.784163916317, 10)

    '''
    All spherical coordinate systems with distances can be converted to
    cartesian coordinates.
    '''

    (x, y, z) = (c2.x, c2.y, c2.z)
    #this only computes the CartesianPoints *once*, and then caches it
    #note that the x/y/z are Quantity objects
    assert isinstance(x, u.Quantity)
    npt.assert_allclose(x.value, 2)
    npt.assert_allclose(y.value, 4)
    npt.assert_allclose(z.value, 8)

    cpt = c2.cartesian
    assert isinstance(cpt, CartesianPoints)
    npt.assert_allclose(cpt.x.value, 2)
    npt.assert_allclose(cpt.y.value, 4)
    npt.assert_allclose(cpt.z.value, 8)

    # with no distance, the unit sphere is assumed when converting to cartesian
    c3 = Galactic(l=158.558650, b=-43.350066, unit=(u.degree, u.degree), distance=None)
    unitcart = c3.cartesian
    npt.assert_allclose(((unitcart.x**2 + unitcart.y**2 +
                          unitcart.z**2)**0.5).value, 1.0)

    # CartesianPoints objects can be added and subtracted, which are
    # vector/elementwise they can also be given as arguments to a coordinate
    # system
    csum = ICRS(c1.cartesian + c2.cartesian)

    npt.assert_allclose(csum.x.value, -8.12016610185)
    npt.assert_allclose(csum.y.value, 3.19380597435)
    npt.assert_allclose(csum.z.value, -8.2294483707)
    npt.assert_allclose(csum.ra.degree, 158.529401774)
    npt.assert_allclose(csum.dec.degree, -43.3235825777)
    npt.assert_allclose(csum.distance.kpc, 11.9942200501)


@pytest.mark.skipif(str('not HAS_SCIPY'))
def test_distances_scipy():
    """
    The distance-related tests that require scipy due to the cosmology
    module needing scipy integration routines
    """
    from .. import Distance
    from ...cosmology import WMAP5

    #try different ways to initialize a Distance
    d4 = Distance(z=0.23)  # uses default cosmology - as of writing, WMAP7
    npt.assert_allclose(d4.z, 0.23, rtol=1e-8)

    d5 = Distance(z=0.23, cosmology=WMAP5)
    npt.assert_allclose(d5.compute_z(WMAP5), 0.23, rtol=1e-8)

    d6 = Distance(z=0.23, cosmology=WMAP5, unit=u.km)
    npt.assert_allclose(d6.value, 3.5417046898762366e+22)


def test_unicode():
    """
    This test could only possibly fail when `sys.getdefaultencoding()`
    is not `utf-8` -- but given a recent fix, is expected to pass on
    all platforms.
    """
    from .. import ICRS
    u = six.text_type("12h46m11.086s -00d30m11.99s")
    c = ICRS(u)
