.. _astropy-coordinates:

*******************************************************
Astronomical Coordinate Systems (`astropy.coordinates`)
*******************************************************

Introduction
============

The `~astropy.coordinates` package provides classes for representing celestial
coordinates, as well as tools for converting between standard systems in a
uniform way.


Getting Started
===============

.. todo:: update to introduce SkyCoordinate first

Coordinate objects are instantiated with a flexible and natural approach that
supports both numeric angle values, (limited) string parsing, and can optionally
include lists of multiple coordinates in one object::

    >>> from astropy.coordinates import ICRS, Galactic
    >>> from astropy import units as u
    >>> ICRS(ra=10.68458*u.degree, dec=41.26917*u.degree)
    <ICRS Coordinate: ra=10.68458 deg, dec=41.26917 deg>
    >>> ICRS('00h42m44.3s', '+41d16m9s')
    <ICRS Coordinate: ra=10.6845833333 deg, dec=41.2691666667 deg>
    >>> ICRS(ra=[10.68458, 83.82208]*u.degree, dec=[41.26917, -5.39111]*u.degree)
    <ICRS Coordinate: (ra, dec) in deg
        [(10.684579999999983, 41.26917), (83.82208000000003, -5.39111)]>

The individual components of a coordinate are `~astropy.coordinates.Longitude`
or `~astropy.coordinates.Latitude` objects, which are specialized versions
of the general `~astropy.coordinates.Angle` class.  The component values are
accessed using aptly named attributes::

    >>> c = ICRS(ra=10.68458*u.degree, dec=41.26917*u.degree)
    >>> c.ra
    <Longitude 10.68457999999999 deg>
    >>> c.ra.hour
    0.712305333...
    >>> c.ra.hms
    hms_tuple(h=0.0, m=42.0, s=44.29919999999...)
    >>> c.dec
    <Latitude 41.26917... deg>
    >>> c.dec.radian
    0.7202828960652...

Coordinates can easily be converted to strings using the
:meth:`~astropy.coordinates.Angle.to_string` method::

.. todo:: update to show formatting with SkyCoordinate

To convert to some other coordinate system, the easiest method is to use
attribute-style access with short names for the built-in systems, but
explicit transformations via the
:meth:`~astropy.coordinates.SphericalCoordinatesBase.transform_to` method
are also available::

.. todo:: SkyCoordinate
.. >>> c.galactic
.. <Galactic l=121.17430 deg, b=-21.57280 deg>
.. >>> c.transform_to(Galactic)
.. <Galactic l=121.17430 deg, b=-21.57280 deg>

Distance from the origin (which is system-dependent, but often the
Earth center) can also be assigned to a coordinate. This specifies a
unique point in 3D space, which also allows conversion to Cartesian
coordinates::

    >>> from astropy.coordinates import Distance
    >>> c = ICRS(ra=10.68458*u.degree, dec=41.26917*u.degree,
    ...          distance=770*u.kpc)
    >>> c.cartesian.x
    <Quantity 568.712865423523... kpc>
    >>> c.cartesian.y
    <Quantity 107.300897404202... kpc>
    >>> c.cartesian.z
    <Quantity 507.889942918757... kpc>

Coordinate objects can also store arrays of coordinates instead of a
single coordinate.  This has a major performance advantage over
transforming many individual coordinate objects separtely.  It also
allows coordinate objects to be used to find matches between two sets
of coordinates::

.. todo:: update with SkyCoordinate
.. >>> #assume ra1/dec1 and ra2/dec2 are arrays loaded from some file
.. >>> c = ICRS(ra1*u.degree, dec1*u.degree)  # doctest: +SKIP
.. >>> catalog = ICRS(ra2*u.degree, dec2*u.degree)  # doctest: +SKIP
.. >>> idx, d2d, d3d = c.match_to_catalog_sky(catalog)  # doctest: +SKIP

These array coordinates can also be indexed in the same way as numpy
arrays::

    >>> len(c[0].ra) # doctest: +SKIP
    TypeError: 'Longitude' object with a scalar value has no len()
    >>> len(c[1:5].ra) # doctest: +SKIP
    4
    >>> matches = catalog[idx]  # doctest: +SKIP
    >>> len(matches) == len(c)  # doctest: +SKIP
    True


.. todo:: update for SkyCoordinate

The `astropy.coordinates` subpackage also provides a quick way to get
coordinates for named objects (if you have an active internet
connection). All subclasses of
`~astropy.coordinates.SphericalCoordinatesBase` have a special class method,
:meth:`~astropy.coordinates.SphericalCoordinatesBase.from_name`, that
accepts a string and queries `Sesame
<http://cds.u-strasbg.fr/cgi-bin/Sesame>`_ to retrieve coordinates for that
object::

.. >>> c = ICRS.from_name("M42")  # doctest: +REMOTE_DATA
.. >>> c.ra, c.dec  # doctest: +SKIP
.. (<Longitude 83.82208... deg>, <Latitude -5.39111... deg>)

This works for any subclass of
`~astropy.coordinates.SphericalCoordinatesBase`::

.. >>> c = Galactic.from_name("M42")  # doctest: +REMOTE_DATA
.. >>> c.l, c.b  # doctest: +SKIP
.. (<Longitude 3.64797... rad>, <Latitude -0.33827... rad>)

.. note::

    This is intended to be a convenience, and is very simple. If you
    need precise coordinates for an object you should find the
    appropriate reference for that measurement and input the
    coordinates manually.


Using `astropy.coordinates`
===========================

More details of using `astropy.coordinates` are provided in the following
sections:

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
``'astropy.coordinates.tests.test_api'`` testing file. It showcases most of
the major capabilities of the package, and hence is a useful supplement to
this document.  You can see it by either looking at it directly if you
downloaded a copy of the astropy source code, or typing the following in an
IPython session::

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
