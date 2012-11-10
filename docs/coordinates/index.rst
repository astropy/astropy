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


Getting Started
===============

Coordinate objects are intantiated with a flexible and natural approach::

    >>> from astropy import coordinates as apc
    >>> from astropy import units as u
    >>> apc.Coordinates(ra=10.68458, dec=41.26917, unit=u.degree)
    <ICRSCoordinates RA=10.68458 deg, Dec=41.26917 deg>
    >>> apc.ICRSCoordinates('00h42m44.3s +41d16m9s')
    <ICRSCoordinates RA=10.68458 deg, Dec=41.26917 deg>

The individual components of a coordinate are `~astropy.coordinates.angles.Angle`
objects, and their values are acceessed using special attributes::

    >>> c = apc.Coordinates(ra=10.68458, dec=41.26917, unit=u.degree)
    >>> c.ra.hours
    0.7123053333333333
    >>> c.dec.radians
    0.7202828960652683
    >>> c.ra.hms
    (0.0, 42, 44.2992000000001)

To convert to some other coordinate system, the easiest method is to use
attribute-style access with short names for the builtin systems, but explicit
transformations via the `transform_to` method are also available::

    >>> c.galactic
    <GalacticCoordinates l=121.17422 deg, b=-21.57283 deg>
    >>> c.transform_to(apc.GalacticCoordinates)
    <GalacticCoordinates l=121.17422 deg, b=-21.57283 deg>

Distances can also be assigned to a coordinate, defining a unique point in 3D
space, which also allows conversion to cartesian coordinates::

    >>> c = apc.Coordinates(ra=10.68458, dec=41.26917, unit=u.degree, distance=apc.Distance(770, u.kpc))
    >>> c.x
    568.7128654235232
    >>> c.y
    107.3008974042025
    >>> c.z
    507.88994291875713


Using `astropy.coordinates`
===========================

An exhaustive resource for the capabilties of this package is the
`astropy.coordinates.tests.test_api` testing file. It showcases most of the
major capabilities of the package, and hence is a useful supplement to this
document.  You can see it by either looking at it directly if you downloaded a
copy of the astropy source code, or typing the following in an IPython session::

    In [1]: from astropy.coordinates.tests import test_api
    In [2]: test_api??


Creating Coordinate Objects
---------------------------

There are two basic ways to create coordinates.  The simplest is to use the
`~astropy.coordsystems.Coordinates` factory class. You provide this with
appropriate keywords, and it will determine from those what sort of coordinate
system you are likely to want.  E.g.::

    >>> Coordinates(ra='12h30m49.42s', dec='+12d23m28.044s')
    <ICRSCoordinates RA=187.70592 deg, Dec=12.39112 deg>
    >>> Coordinates(l=-76.22237, b=74.49108, unit=u.degree)
    <GalacticCoordinates l=-76.22237 deg, b=74.49108 deg>
    >>> Coordinates(az=45.8, el=42.3, unit=u.degree)
    <HorizontalCoordinates el=42.30000 deg, az=45.80000 deg>

.. warning::
    `~astropy.coordinates.coordsystems.Coordinates` is *not* an actual class
    that ever gets instantiated.  Do not do ``isinstance(coord, Coordinates)``
    in your code to check if something is a coordinate.  Instead, use
    `~astropy.coordinates.coordsystems.SphericalCoordinatesBase`, which is the
    true base class of all astropy coordinate classes.

The second method is to directly initialize your preferred coordinate system by
the name of the class representing that system.  E.g.,

    >>> ICRSCoordinates(187.70592, 12.39112, unit=u.degree)
    <ICRSCoordinates RA=187.70592 deg, Dec=12.39112 deg>
    >>> FK4Coordinates(187.07317, 12.66715, unit=u.degree)
    <FK4Coordinates RA=187.07317 deg, Dec=12.66715 deg>
    >>> GalacticCoordinates(-76.22237, 74.49108, unit=u.degree)
    <GalacticCoordinates l=-76.22237 deg, b=74.49108 deg>

Note that for these you don't need to specify `ra`/`dec` or `l`/`b` keywords,
because those are already the appropriate longitude/latitude angles for the
requested coordinate systems.

With this method, the parsing of coordinates is quite flexible, designed to
interpret any *unambiguous* string a human could interpret.  For example::

    >>> ICRSCoordinates(54.12412, "-41:08:15.162342", unit=u.degree)
    <ICRSCoordinates RA=54.12412 deg, Dec=-41.13755 deg>
    >>> ICRSCoordinates("3:36:29.7888 -41:08:15.162342", unit=u.hour)
    <ICRSCoordinates RA=54.12412 deg, Dec=-41.13755 deg>
    >>> ICRSCoordinates("54.12412 -41:08:15.162342")
    <ICRSCoordinates RA=54.12412 deg, Dec=-41.13755 deg>
    >>>apc.ICRSCoordinates("14.12412 -41:08:15.162342")
    UnitsError: No units were specified, and the angle value was ambiguous between hours and degrees.

Working with Angles
-------------------

The angular components of a coordinate are represented by objects of the
`~astropy.coordinates.angles.Angle` class. These objects can be instantiated on
their own anywhere a representation of an angle is needed, and support a variety
of ways of representing the value of the angle::

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
    >>> a.string()
    '57 17 44.80625'
    >>> a.string(sep=':')
    '57:17:44.80625'
    >>> a.string(sep='dms')
    '57d17m44.80625s'
    >>> a.string(u.hour, sep='hms')
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

Angle objects can also be used for creating coordinate objects::

    >>> ICRSCoordinates(Angle(1, u.radian), Angle(2, u.radian))
    <ICRSCoordinates RA=57.29578 deg, Dec=114.59156 deg>
    >>> Coordinates(Angle(1, u.radian), Angle(2, u.radian))
    ValueError: Two angles were provided ('<astropy.coordinates.angles.Angle 57.29578 deg>', '<astropy.coordinates.angles.Angle 114.59156 deg>'), but the coordinate system was not provided. Specify the system via keywords or use the corresponding class (e.g. GalacticCoordinate).
    >>> Coordinates(RA(1, u.radian), Dec(2, u.radian))
    <ICRSCoordinates RA=57.29578 deg, Dec=114.59156 deg>



Distances
---------

Content

Transforming Between Systems
----------------------------

Content

Custom Coordinate Classes
-------------------------

Content


See Also
========

Some references particularly useful in understanding subtleties of the
coordinate systems implemented here include:

* `Standards Of Funamental Astronomy <http://www.iausofa.org/>`_
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
