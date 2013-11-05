.. _astropy-coordinates:

*******************************************************
Astronomical Coordinate Systems (`astropy.coordinates`)
*******************************************************

Introduction
============

The `~astropy.coordinates` package provides classes for representing celestial
coordinates, as well as tools for converting between standard systems in a
uniform way.

.. warning::
    `~astropy.coordinates` is currently a work-in-progress, and thus it is
    possible there will be significant API changes in later versions of
    Astropy. If you have specific ideas for how it might be improved,
    feel free to let us know on the `astropy-dev mailing list`_ or at
    http://feedback.astropy.org


Getting Started
===============

Coordinate objects are instantiated with a flexible and natural approach that
supports both numeric angle values, (limited) string parsing, and can optionally
include lists of multiple coordinates in one object::

    >>> from astropy.coordinates import ICRS, Galactic
    >>> from astropy import units as u
    >>> ICRS(ra=10.68458, dec=41.26917, unit=(u.degree, u.degree))
    <ICRS RA=10.68458 deg, Dec=41.26917 deg>
    >>> ICRS('00h42m44.3s +41d16m9s')
    <ICRS RA=10.68458 deg, Dec=41.26917 deg>
    >>> ICRS('00h42m44.3s +41d16m9s')
    <ICRS RA=10.68458 deg, Dec=41.26917 deg>
    >>> ICRS(ra=[10.68458, 83.82208], dec=[41.26917, -5.39111], unit=(u.degree, u.degree))
    <ICRS RA=[ 10.68458  83.82208] deg, Dec=[ 41.26917  -5.39111] deg>

The individual components of a coordinate are `~astropy.coordinates.angles.Longitude`
or `~astropy.coordinates.angles.Latitude` objects, which are specialized versions
of the general `~astropy.coordinates.angles.Angle` class.  The component values are
accessed using aptly named attributes::

    >>> c = ICRS(ra=10.68458, dec=41.26917,
    ...          unit=(u.degree, u.degree))
    >>> c.ra
    <Longitude 10.684579999999983 deg>
    >>> c.ra.hour
    0.7123053333333323
    >>> c.ra.hms
    (0.0, 42.0, 44.299199999996262)
    >>> c.dec
    <Latitude 41.26917 deg>
    >>> c.dec.radian
    0.7202828960652683

Coordinates can easily be converted to strings using the ``to_string`` method::

    >>> c.to_string()
    '0h42m44.2992s 41d16m09.012s'
    >>> c.to_string(precision=1)
    '0h42m44.3s 41d16m09.0s'
    >>> c.to_string(precision=1, sep=' ')
    '0 42 44.3 41 16 09.0'

To convert to some other coordinate system, the easiest method is to
use attribute-style access with short names for the built-in systems,
but explicit transformations via the `transform_to` method are also
available::

    >>> c.galactic
    <Galactic l=121.17430 deg, b=-21.57280 deg>
    >>> c.transform_to(Galactic)
    <Galactic l=121.17430 deg, b=-21.57280 deg>

Distance from the origin (which is system-dependent, but often the
Earth center) can also be assigned to a coordinate. This specifies a
unique point in 3D space, which also allows conversion to Cartesian
coordinates::

    >>> from astropy.coordinates import Distance
    >>> c = ICRS(ra=10.68458, dec=41.26917,
    ...          unit=(u.degree, u.degree),
    ...          distance=Distance(770, u.kpc))
    >>> c.x
    <Quantity 568.7128654235232 kpc>
    >>> c.y
    <Quantity 107.30089740420232 kpc>
    >>> c.z
    <Quantity 507.88994291875713 kpc>

Coordinate objects can also store arrays of coordinates instead of a
single coordinate.  This has a major performance advantage over
transforming many individual coordinate objects separtely.  It also
allows coordinate objects to be used to find matches between two sets
of coordinates::

    >>> #assume ra1/dec1 and ra2/dec2 are arrays loaded from some file
    >>> c = ICRS(ra1, dec1, unit=(u.degree, u.degree))  # doctest: +SKIP
    >>> catalog = ICRS(ra2, dec2, unit=(u.degree, u.degree))  # doctest: +SKIP
    >>> idx, d2d, d3d = c1.match_to_catalog_sky(catalog)  # doctest: +SKIP

These array coordinates can also be indexed in the same way as numpy
arrays::

    >>> len(c[0].ra) # doctest: +SKIP
    TypeError: 'Longitude' object with a scalar value has no len()
    >>> len(c[1:5].ra) # doctest: +SKIP
    4
    >>> matches = catalog[idx]  # doctest: +SKIP
    >>> len(matches) == len(c)  # doctest: +SKIP
    True


The `astropy.coordinates` subpackage also provides a quick way to get
coordinates for named objects (if you have an active internet
connection). All subclasses of
`~astropy.coordinates.coordsystems.SphericalCoordinatesBase` have a
special class method, `from_name()`, that accepts a string and queries
`Sesame <http://cds.u-strasbg.fr/cgi-bin/Sesame>`_ to retrieve
coordinates for that object::

    >>> c = ICRS.from_name("M42")
    >>> c.ra, c.dec
    (<Longitude 83.82208... deg>, <Latitude -5.39111... deg>)

This works for any subclass of
`~astropy.coordinates.coordsystems.SphericalCoordinatesBase`::

    >>> c = Galactic.from_name("M42")
    >>> c.l, c.b
    (<Longitude 3.64797... rad>, <Latitude -0.33827... rad>)

.. note::

    This is intended to be a convenience, and is very simple. If you
    need precise coordinates for an object you should find the
    appropriate reference for that measurement and input the
    coordinates manually.


Using `astropy.coordinates`
===========================

More details of using `astropy.coordinates` are provided in the following sections:

.. toctree::
   :maxdepth: 1

   angles
   creating
   formatting
   separations
   distances
   transforming
   matching
   designing
   sgr-example


In addition, another resource for the capabilities of this package is the
`astropy.coordinates.tests.test_api` testing file. It showcases most of the
major capabilities of the package, and hence is a useful supplement to this
document.  You can see it by either looking at it directly if you downloaded a
copy of the astropy source code, or typing the following in an IPython session::

    In [1]: from astropy.coordinates.tests import test_api
    In [2]: test_api??


See Also
========

Some references particularly useful in understanding subtleties of the
coordinate systems implemented here include:

* `Standards Of Fundamental Astronomy <http://www.iausofa.org/>`_
    The definitive implementation of IAU-defined algorithms.  The "SOFA Tools
    for Earth Attitude" document is particularly valuable for understanding
    the latest IAU standards in detail.
* `USNO Circular 179 <http://aa.usno.navy.mil/publications/docs/Circular_179.php>`_
    A useful guide to the IAU 2000/2003 work surrounding ICRS/IERS/CIRS and
    related problems in precision coordinate system work.
* Meeus, J. "Astronomical Algorithms"
    A valuable text describing details of a wide range of coordinate-related
    problems and concepts.



Reference/API
=============

.. automodapi:: astropy.coordinates
  :skip: ICRSCoordinates
  :skip: FK5Coordinates
  :skip: FK4Coordinates
  :skip: FK4NoETermCoordinates
  :skip: GalacticCoordinates
  :skip: HorizontalCoordinates

.. the ":skip:"s above are to not document the v0.3 backwards-compatibility names.  They will be removed in the next version
