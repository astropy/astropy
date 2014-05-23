.. include:: references.txt

.. _astropy-coordinates:

*******************************************************
Astronomical Coordinate Systems (`astropy.coordinates`)
*******************************************************

Introduction
============

The `~astropy.coordinates` package provides classes for representing a
variety of celestial/spatial  coordinates, as well as tools for
converting between common coordinate systems in a uniform way.


Getting Started
===============

The simplest way to use `~astropy.coordinates` is to use the |skycoord|
class. |skycoord| objects are instantiated with a flexible and natural
approach that supports angle values provided as |quantity| objects, as
well as string parsing::

    >>> from astropy import units as u
    >>> from astropy.coordinates import SkyCoord
    >>> SkyCoord(ra=10.68458*u.degree, dec=41.26917*u.degree)
    <SkyCoord (NoFrame): ra=10.68458... deg, dec=41.269... deg>
    >>> SkyCoord('00h42m44.3s', '+41d16m9s')
    <SkyCoord (NoFrame): ra=10.68458... deg, dec=41.269... deg>

The individual components of equatorial coordinates are
`~astropy.coordinates.Longitude` or `~astropy.coordinates.Latitude`
objects, which are specialized versions of the general
`~astropy.coordinates.Angle` class.  The component values are accessed
using aptly named attributes::

    >>> c = SkyCoord(ra=10.68458*u.degree, dec=41.26917*u.degree)
    >>> c.ra
    <Longitude 10.68458 deg>
    >>> c.ra.hour
    0.712305333...
    >>> c.ra.hms
    hms_tuple(h=0.0, m=42.0, s=44.299...)
    >>> c.dec
    <Latitude 41.2691... deg>
    >>> c.dec.degree
    41.2691...
    >>> c.dec.radian
    0.720282896065...

Coordinates can easily be converted to strings using the
:meth:`~astropy.coordinates.SkyCoord.to_string` method::

    >>> c = SkyCoord(ra=10.68458*u.degree, dec=41.26917*u.degree)
    >>> c.to_string('decimal')
    '10.6846 41.2692'
    >>> c.to_string('dms')
    '10d41m04.488s 41d16m09.012s'
    >>> c.to_string('hmsdms')
    '00h42m44.2992s +41d16m09.012s'

For more control over the string formatting, use the
`~astropy.coordinates.Angle.to_string` method of the individual
components::

    >>> c.ra.to_string(decimal=True)
    '10.6846'
    >>> c.dec.to_string(format='latex')
    '$41^\\circ16{}^\\prime09.012{}^{\\prime\\prime}$'
    >>> msg = 'My coordinates are: ra="{0}"" dec="{1}"'
    >>> msg.format(c.ra.to_string(sep=':'), c.dec.to_string(sep=':'))
    'My coordinates are: ra="10:41:04.488"" dec="41:16:09.012"'


The above examples used a "NoFrame" |skycoord|, because a coordinate
frame was not given.  However, to use the full power of
`~astropy.coordinates`,  you should specify the reference frame your
coordinates are defined in::

    >>> c_icrs = SkyCoord(ra=10.68458*u.degree, dec=41.26917*u.degree, frame='icrs')

Once you've defined the frame of your coordinates, you can transform from that
frame to another frame.  You can do this a few different ways: if you just want
the default version of that frame, you can use attribute-style access.  For
more control, you can use the `~astropy.coordinates.SkyCoord.transform_to` method,
which accepts  either the name of that frame, or an instance of the frame object::

    >>> from astropy.coordinates import FK5
    >>> c_icrs.galactic
    <SkyCoord (Galactic): l=121.17430... deg, b=-21.572... deg>
    >>> c_fk5 = c_icrs.transform_to('fk5')  # c_icrs.fk5 does the same thing
    >>> c_fk5
    <SkyCoord (FK5): equinox=J2000.000, ra=10.68459... deg, dec=41.269... deg>
    >>> c_fk5.transform_to(FK5(equinox='J1975'))  # precess to a different equinox
    <SkyCoord (FK5): equinox=J1975.000, ra=10.34209... deg, dec=41.132... deg>

|skycoord| (and all other `~astropy.coordinates`) objects also support
array coordinates.  These work the same as single-value coordinates, but
they store multiple coordinates in a single object.  When you're going
to apply the same operation to many different coordinates (say, from a
catalog), this is a better choice than a list of |skycoord| objects,
because it will be *much* faster than applying the operation to each
|skycoord| in a for loop.

    >>> SkyCoord(ra=[10, 11]*u.degree, dec=[41, -5]*u.degree)
    <SkyCoord (NoFrame): (ra, dec) in deg
        [(10.0, 41.0), (11.0, -5.0)]>


|skycoord| defines a number of convinience methods, as well, like on-sky
separation between two coordinates (and catalog matching, detailed in
:ref:`astropy-coordinates-matching`)::

    >>> c1 = SkyCoord(ra=10*u.degree, dec=9*u.degree, frame='icrs')
    >>> c2 = SkyCoord(ra=11*u.degree, dec=10*u.degree, frame='icrs')
    >>> c1.separation(c2)
    <Angle 1.4045397278... deg>


Distance from the origin (which is system-dependent, but often the Earth
center) can also be assigned to a |skycoord|. With two angles and a
distance, a unique point in 3D space is available, which also allows
conversion to the Cartesian representation of this location::

    >>> from astropy.coordinates import Distance
    >>> c = SkyCoord(ra=10.68458*u.degree, dec=41.26917*u.degree, distance=770*u.kpc)
    >>> c.cartesian.x
    <Quantity 568.712865423... kpc>
    >>> c.cartesian.y
    <Quantity 107.300897404... kpc>
    >>> c.cartesian.z
    <Quantity 507.889942918... kpc>

With distances assigned, |skycoord| convinience methods are more powerful, as
they can make use of the 3d information. For example::

    >>> c1 = SkyCoord(ra=10*u.degree, dec=9*u.degree, distance=10*u.pc, frame='icrs')
    >>> c2 = SkyCoord(ra=11*u.degree, dec=10*u.degree,distance=11.5*u.pc, frame='icrs')
    >>> c1.separation_3d(c2)
    <Distance 1.522860241... pc>


Finally, the `astropy.coordinates` subpackage also provides a quick way to get
coordinates for named objects (if you have an active internet
connection). The `~astropy.coordinates.SkyCoord.from_name` method of |skycoord|
uses  `Sesame
<http://cds.u-strasbg.fr/cgi-bin/Sesame>`_ to retrieve coordinates for a particular named object
object::

    >>> SkyCoord.from_name("M42")  # doctest: +REMOTE_DATA
    <SkyCoord (ICRS): ra=83.8220... deg, dec=-5.391... deg>

.. note::

    `~astropy.coordinates.SkyCoord.from_name` is intended to be a convenience,
    and is rather simple. If you need precise coordinates for an object you
    should find the appropriate reference for that measurement and input the
    coordinates manually.


.. _astropy-coordinates-overview:

Overview of `astropy.coordinates` concepts
==========================================

.. note ::
    The `~astropy.coordinates` package from v0.4 onward builds from
    previous versions of  the package, and more detailed information and
    justification of the design is available in `APE (Astropy Proposal for Enhancement) 5 <https://github.com/astropy/astropy- APEs/blob/master/APE5.rst>`_.

Here we provide an overview of the package and associated framework.
This background information is not necessary for simply using
`~astropy.coordinates`, particularly if you use the |skycoord| high-
level class, but it is helpful for more advanced usage, particularly
creating your own frame, transformations, or representations. Another
useful piece of background infromation are some
:ref:`astropy-coordinates-definitions` as they are used in
`~astropy.coordinates`.

`~astropy.coordinates` is built on a three-tired system of objects:
representations, frames, and a high-level class.  Representations
classes are a particular way of storing a three-dimensional data point
(or points), such as Cartesian coordinates or spherical polar
coordinates. Frames are particular reference frames like FK5 or ICRS,
which may store their data in different representations, but have well-
defined transformations between each other. These transformations are
all stored in the `astropy.coordinates.frame_transform_graph`, and new
transformations can be created by users. Finally, the high-level class
(|skycoord|) uses the frame classes, but provides a more accessible
interface to these objects as well as various convenience methods and
more string-parsing capabilities.

Separating these concepts makes it easier to extend the functionality of
`~astropy.coordinates`.  It allows representations, frames, and
transformations to be defined or extended separately, while still
preserving the high-level capabilities and simplicity of the |skycoord|
class.


Using `astropy.coordinates`
===========================

More detailed information on using the package is provided on separate pages,
listed below.

.. toctree::
   :maxdepth: 1

   angles
   skycoord
   transforming
   formatting
   matchsep
   representations
   frames
   sgr-example
   definitions


In addition, another resource for the capabilities of this package is the
``'astropy.coordinates.tests.test_api_ape5'`` testing file. It showcases most of
the major capabilities of the package, and hence is a useful supplement to
this document.  You can see it by either looking at it directly if you
downloaded a copy of the astropy source code, or typing the following in an
IPython session::

    In [1]: from astropy.coordinates.tests import test_api_ape5
    In [2]: test_api_ape5??


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
