*******************************************************
Astronomical Coordinate Systems (`astropy.coordinates`)
*******************************************************

Introduction
============

The `~astropy.coordinates` package provides classes for representing celestial
coordinates, as well as tools for converting between standard systems in a
uniform way.

.. note::
    The current `~astropy.coordinates` framework only accepts scalar
    coordinates, i.e. one coordinate per object.  In the next release it will
    be expanded to accept arrays of coordinates.

.. warning::
    `~astropy.coordinates` is currently a work-in-progress, and thus it is
    possible there will be significant API changes in later versions of
    Astropy.


Getting Started
===============

Coordinate objects are instantiated with a flexible and natural approach that
supports both numeric angle values and (limited) string parsing::

    >>> from astropy import coordinates as coord
    >>> from astropy import units as u
    >>> coord.ICRSCoordinates(ra=10.68458, dec=41.26917, unit=(u.degree, u.degree))
    <ICRSCoordinates RA=10.68458 deg, Dec=41.26917 deg>
    >>> coord.ICRSCoordinates('00h42m44.3s +41d16m9s')
    <ICRSCoordinates RA=10.68458 deg, Dec=41.26917 deg>

The individual components of a coordinate are `~astropy.coordinates.angles.Angle`
objects, and their values are accessed using special attributes::

    >>> c = coord.ICRSCoordinates(ra=10.68458, dec=41.26917, unit=(u.degree, u.degree))
    >>> c.ra
    <RA 10.68458 deg>
    >>> c.ra.hours
    0.7123053333333333
    >>> c.ra.hms
    (0.0, 42, 44.2992000000001)
    >>> c.dec
    <Dec 41.26917 deg>
    >>> c.dec.radians
    0.7202828960652683

To convert to some other coordinate system, the easiest method is to use
attribute-style access with short names for the built-in systems, but explicit
transformations via the `transform_to` method are also available::

    >>> c.galactic
    <GalacticCoordinates l=121.17422 deg, b=-21.57283 deg>
    >>> c.transform_to(coord.GalacticCoordinates)
    <GalacticCoordinates l=121.17422 deg, b=-21.57283 deg>

Distances from the origin (which is system-dependent, but often the Earth
center) can also be assigned to a coordinate. This specifies a unique point
in 3D space, which also allows conversion to cartesian coordinates::

    >>> c = coord.ICRSCoordinates(ra=10.68458, dec=41.26917, unit=(u.degree, u.degree), distance=coord.Distance(770, u.kpc))
    >>> c.x
    568.7128654235232
    >>> c.y
    107.3008974042025
    >>> c.z
    507.88994291875713


Using `astropy.coordinates`
===========================

An exhaustive resource for the capabilities of this package is the
`astropy.coordinates.tests.test_api` testing file. It showcases most of the
major capabilities of the package, and hence is a useful supplement to this
document.  You can see it by either looking at it directly if you downloaded a
copy of the astropy source code, or typing the following in an IPython session::

    In [1]: from astropy.coordinates.tests import test_api
    In [2]: test_api??


Creating Coordinate Objects
---------------------------

Creating new coordinate objects is of course crucial to using
`~astropy.coordinates`.  The typical way to create a new coordinate object
is to directly initialize your preferred coordinate system using standard
python class creation, using the name of the class representing that
system.  For example::

    >>> from astropy.coordinates import ICRSCoordinates, FK4Coordinates, GalacticCoordinates
    >>> ICRSCoordinates(187.70592, 12.39112, unit=(u.degree, u.degree))
    <ICRSCoordinates RA=187.70592 deg, Dec=12.39112 deg>
    >>> FK4Coordinates(187.07317, 12.66715, unit=(u.degree, u.degree))
    <FK4Coordinates RA=187.07317 deg, Dec=12.66715 deg>
    >>> GalacticCoordinates(-76.22237, 74.49108, unit=(u.degree, u.degree))
    <GalacticCoordinates l=-76.22237 deg, b=74.49108 deg>

Note that if you do not provide units explicitly, this will fail::

    >>> ICRSCoordinates(23, 1)
    UnitsError: No unit was specified in Angle initializer; the unit parameter should be an object from the  astropy.units module (e.g. 'from astropy import units as u', then use 'u.degree').

While the above example uses python numerical types, you can also provide strings to create coordinates.
If the `unit` parameter is ``(None, None)`` (the default), strings will be interpreted using the `Angle` 
class' parsing scheme, and has a guiding principal of being able to interpret any *unambiguous* string 
specifying an angle. For the exact rules for how each string is parsed, see the 
`~astropy.coordinates.angles.Angle` documentation.  Some examples::

    >>> ICRSCoordinates("3h36m29.7888s -41d08m15.162342s", unit=(None, None))
    <ICRSCoordinates RA=54.12412 deg, Dec=-41.13755 deg>
    >>> ICRSCoordinates("3h36m29.7888s -41d08m15.162342s")
    <ICRSCoordinates RA=54.12412 deg, Dec=-41.13755 deg>
    >>> ICRSCoordinates("14.12412 hours", "-41:08:15.162342 degrees")
    <ICRSCoordinates RA=211.86180 deg, Dec=-41.13755 deg>
    >>> ICRSCoordinates("14.12412 -41:08:15.162342")
    UnitsError: Could not infer Angle units from provided string 14.12412

You can also directly specify the units for both to resolve ambiguities::

    >>> ICRSCoordinates("14.12412 -41:08:15.162342", unit=(u.hour, u.degree))
    <ICRSCoordinates RA=211.86180 deg, Dec=-41.13755 deg>
    >>> ICRSCoordinates("54:7:26.832 -41:08:15.162342", unit=(u.degree, u.degree))
    <ICRSCoordinates RA=54.12412 deg, Dec=-41.13755 deg>
    >>> ICRSCoordinates('3 4 5 +6 7 8', unit=(u.hour, u.degree))
    <ICRSCoordinates RA=46.02083 deg, Dec=6.11889 deg>
    >>> ICRSCoordinates('3h4m5s +6d7m8s', unit=(u.hour, u.degree))
    <ICRSCoordinates RA=46.02083 deg, Dec=6.11889 deg>

This will also give you an error if you give a string with units that conflict with your desired units::

    >>> ICRSCoordinates('3d4m5s +6h7m8s', unit=(u.hour, u.degree))
    ValueError: parse_hours: Invalid input string, can't parse to HMS. (3d4m5s)



Working with Angles
-------------------

The angular components of a coordinate are represented by objects of the
`~astropy.coordinates.angles.Angle` class. These objects can be instantiated on
their own anywhere a representation of an angle is needed, and support a variety
of ways of representing the value of the angle::

    >>> from astropy.coordinates import Angle
    >>> a = Angle(1, u.radian)
    >>> a
    <astropy.coordinates.angles.Angle 57.29578 deg>
    >>> a.radians
    1
    >>> a.degrees
    57.29577951308232
    >>> a.hours
    3.819718634205488
    >>> a.hms
    (3.0, 49, 10.987083139757061)
    >>> a.dms
    (57.0, 17, 44.80624709636231)
    >>> a.format()
    '57d17m44.80625s'
    >>> a.format(sep=':')
    '57:17:44.80625'
    >>> a.format(sep=('deg','m','s'))
    '57deg17m44.80625s'
    >>> a.format(u.hour)
    '3h49m10.98708s'

`~astropy.corodinates.angles.Angle` objects can also have bounds.  These specify
either a limited range in which the angle is valid (if it's <360 degrees), or
the limit at which an angle is wrapped back around to 0::

    >>> Angle(90, unit=u.degree, bounds=(0,180))
    <Angle 90.00000 deg>
    >>> Angle(-270, unit=u.degree, bounds=(0,180))
    <Angle 90.00000 deg>
    >>> Angle(181, unit=u.degree, bounds=(0,180))
    BoundsError: The angle given falls outside of the specified bounds.
    >>> Angle(361, unit=u.degree, bounds=(0,360))
    <Angle 1.00000 deg>

Angles will also behave correctly for appropriate arithmetic operations::

    >>> a = Angle(1, u.radian)
    >>> a + a
    <Angle 114.59156 deg>
    >>> a - a
    <Angle 0.00000 deg>
    >>> a == a
    True
    >>> a == (a + a)
    False

Angle objects can also be used for creating coordinate objects::

    >>> ICRSCoordinates(Angle(1, u.radian), Angle(2, u.radian))
    <ICRSCoordinates RA=57.29578 deg, Dec=114.59156 deg>
    >>> ICRSCoordinates(Angle(1, u.radian), Angle(2, u.radian))
    ValueError: Two angles were provided ('<astropy.coordinates.angles.Angle 57.29578 deg>', '<astropy.coordinates.angles.Angle 114.59156 deg>'), but the coordinate system was not provided. Specify the system via keywords or use the corresponding class (e.g. GalacticCoordinate).
    >>> ICRSCoordinates(RA(1, u.radian), Dec(2, u.radian))
    <ICRSCoordinates RA=57.29578 deg, Dec=114.59156 deg>


Separations
-----------

The on-sky separation is easily computed with the `separation` method::

    >>> c1 = ICRSCoordinates('5h23m34.5s -69d45m22s')
    >>> c2 = ICRSCoordinates('0h52m44.8s -72d49m43s')
    >>> sep = c1.separation(c2)
    >>> sep
    <AngularSeparation 20.74612 deg>

The `~astropy.coordinates.angles.AngularSeparation` object is a subclass of
`~astropy.coordinates.angles.Angle`, so it can be accessed in the same ways,
along with a few additions::

    >>> sep.radians
    0.36208807374669766
    >>> sep.hours
    1.383074562513832
    >>> sep.arcmins
    1244.7671062624488
    >>> sep.arcsecs
    74686.02637574692


Distances
---------

Coordinates can also have line-of-sight distances.  If these are provided, a
coordinate object becomes a full-fledged point in three-dimensional space.  If
not (i.e., the `distance` attribute of the coordinate object is `None`), the
point is interpreted as lying on the (dimensionless) unit sphere.

The `~astropy.coordinates.distances.Distance` class is provided to represent a
line-of-sight distance for a coordinate.  It must include a length unit to be
valid.::

    >>> from astropy.coordinates import Distance
    >>> d = Distance(770)
    UnitsError: A unit must be provided for distance.
    >>> d = Distance(770, u.kpc)
    >>> c = ICRSCoordinates('00h42m44.3s +41d16m9s', distance=d)
    >>> c
    <ICRSCoordinates RA=10.68458 deg, Dec=41.26917 deg, Distance=7.7e+02 kpc>

If a distance is available, the coordinate can be converted into cartesian
coordinates using the `x`/`y`/`z` attributes::

    >>> c.x
    568.7128882165681
    >>> c.y
    107.3009359688103
    >>> c.z
    507.8899092486349

.. note::

    The location of the origin is different for different coordinate
    systems, but for common celestial coordinate systems it is often
    the Earth center (or for precision work, the Earth/Moon barycenter).

The cartesian coordinates can also be accessed via the
`~astropy.coordinates.distances.CartesianCoordinates` object, which has
additional capabilities like arithmetic operations::

    >>> cp = c.cartesian
    >>> cp
    <CartesianPoints (568.712888217, 107.300935969, 507.889909249) kpc>
    >>> cp.x
    568.7128882165681
    >>> cp.y
    107.3009359688103
    >>> cp.z
    507.8899092486349
    >>> cp.unit
    Unit("kpc")
    >>> cp + cp
    <CartesianPoints (1137.42577643, 214.601871938, 1015.7798185) kpc>
    >>> cp - cp
    <CartesianPoints (0.0, 0.0, 0.0) kpc>

This cartesian representation can also be used to create a new coordinate
object::

    >>> ICRSCoordinates(x=568.7129, y=107.3009, z=507.8899, unit=u.kpc)
    <ICRSCoordinates RA=10.68458 deg, Dec=41.26917 deg, Distance=7.7e+02 kpc>

Finally, two coordinates with distances can be used to derive a real-space
distance (i.e., non-projected separation)::

    >>> c1 = ICRSCoordinates('5h23m34.5s -69d45m22s', distance=Distance(49, u.kpc))
    >>> c2 = ICRSCoordinates('0h52m44.8s -72d49m43s', distance=Distance(61, u.kpc))
    >>> sep3d = c1.separation_3d(c2)
    >>> sep3d
    <Distance 23.05685 kpc>
    >>> sep3d.kpc
    23.05684814695706
    >>> sep3d.Mpc
    0.02305684814695706
    >>> sep3d.au
    4755816315.663559


Transforming Between Systems
----------------------------

`astropy.coordinates` supports a rich system for transforming coordinates from
one system to another.  The key concept is that a registry of all the
transformations is used to determine which coordinates can convert to others.
When you ask for a transformation, the registry (or "transformation graph") is
searched for the shortest path from your starting coordinate to your target, and
it applies all of the transformations in that path in series.   This allows only
the simplest transformations to be defined, and the package will automatically
determine how to combine those transformations to get from one system to
another.

As described above, there are two ways of transforming coordinates.  Coordinates
that have an alias (created with
`~astropy.coordinates.transformations.coordinate_alias`) can be converted by
simply using attribute style access to any other coordinate system::

    >>> gc = GalacticCoordinates(l=0, b=45, unit=u.degree)
    >>> gc.fk5
    <FK5Coordinates RA=229.27250 deg, Dec=-1.12842 deg>
    >>> ic = ICRSCoordinates(ra=0, dec=45, unit=u.degree)
    >>> ic.fk5
    <FK5Coordinates RA=0.00001 deg, Dec=45.00000 deg>

While this appears to be simple attribute-style access, it is actually just
syntactic sugar for the `transform_to` method::

    >>> gc.transform_to(FK5Coordinates)
    <FK5Coordinates RA=229.27250 deg, Dec=-1.12842 deg>
    >>> ic.transform_to(FK5Coordinates)
    <FK5Coordinates RA=0.00001 deg, Dec=45.00000 deg>

The full list of supported coordinate systems and transformations is in the
`astropy.coordinates` API documentation below.

Additionally, some coordinate systems support precessing the coordinate to
produce a new coordinate in the same system but at a different equinox::

    >>> fk5c = FK5Coordinates('02h31m49.09s +89d15m50.8s', equinox=Time('J2000', scale='utc'))
    >>> fk5c
    <FK5Coordinates RA=37.95454 deg, Dec=89.26411 deg>
    >>> fk5c.precess_to(Time(2100, format='jyear', scale='utc'))
    <FK5Coordinates RA=88.32396 deg, Dec=89.54057 deg>

Coordinate systems do not necessarily all support an equinox nor precession, as it is a
meaningless action for coordinate systems that do not depend on a particular equinox.

Furthermore, coordinates typically have an `obstime` attribute, intended to record the
time of the observation.  Some systems (especially FK4) require this information due
to being non-inertial frames (i.e., they rotate over time due to motions of the
defining stars).


Designing Coordinate Systems
----------------------------

New coordinate systems can easily be added by users by simply subclassing the
`~astropy.coordinates.coordsystems.SphericalCoordinatesBase` object.
Detailed instructions for subclassing are in the docstrings for that class.  If
defining a latitude/longitude style of coordinate system, the
`_initialize_latlong` method and `_init_docstring_param_templ` attribute are
helpful for automated processing of the inputs.

To define transformations to and from this coordinate, the easiest method is to
define a function that accepts an object in one coordinate system and returns
the other.  Decorate this function with
`~astropy.coordinates.transformations.transform_function` function decorator,
supplying the information to determine which coordinates the function transforms
to or from.  This will register the transformation, allowing any other
coordinate object to use this converter.  You can also use the
`~astropy.coordinates.transformations.static_transform_matrix` and
`~astropy.coordinates.transformations.dynamic_transform_matrix` decorators to
specify the transformation in terms of 3 x 3 cartesian coordinate transformation
matrices (typically rotations).


See Also
========

Some references particularly useful in understanding subtleties of the
coordinate systems implemented here include:

* `Standards Of Fundamental Astronomy <http://www.iausofa.org/>`_
    The definitive implementation of IAU-defined algorithms.  The "SOFA Tools
    for Earth Attitude" document is particularly valuable for understanding
    the latest IAU standards in detail.
* `USNO Circular 179 <http://www.usno.navy.mil/USNO/astronomical-applications/publications/circ-179>`_
    A useful guide to the IAU 2000/2003 work surrounding ICRS/IERS/CIRS and
    related problems in precision coordinate system work.
* Meeus, J. "Astronomical Algorithms"
    A valuable text describing details of a wide range of coordinate-related
    problems and concepts.



Reference/API
=============

.. automodapi:: astropy.coordinates
