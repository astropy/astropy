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
approach that supports inputs provided in a number of convenient
formats.  The following ways of initializing a coordinate are all
equivalent::

    >>> from astropy import units as u
    >>> from astropy.coordinates import SkyCoord

    >>> c = SkyCoord(ra=10.625*u.degree, dec=41.2*u.degree, frame='icrs')
    >>> c = SkyCoord(10.625, 41.2, frame='icrs', unit='deg')
    >>> c = SkyCoord('00h42m30s', '+41d12m00s', frame='icrs')
    >>> c = SkyCoord('00h42.5m', '+41d12m')
    >>> c = SkyCoord('00 42 30 +41 12 00', unit=(u.hourangle, u.deg))
    >>> c = SkyCoord('00:42.5 +41:12', unit=(u.hourangle, u.deg))
    >>> c
    <SkyCoord (ICRS): (ra, dec) in deg
        (10.625, 41.2)>

The examples above illustrate a few simple rules to follow when creating a coordinate
object:

- Coordinate values can be provided either as unnamed positional arguments or
  via keyword arguments like ``ra``, ``dec``, ``l``, or ``b`` (depending on the frame).
- Coordinate ``frame`` keyword is optional and defaults to ICRS.
- Angle units must be specified, either in the values themselves
  (e.g. ``10.5*u.degree`` or ``'+41d12m00s'``) or via the ``unit`` keyword.

|skycoord| and all other `~astropy.coordinates` objects also support
array coordinates.  These work the same as single-value coordinates, but
they store multiple coordinates in a single object.  When you're going
to apply the same operation to many different coordinates (say, from a
catalog), this is a better choice than a list of |skycoord| objects,
because it will be *much* faster than applying the operation to each
|skycoord| in a for loop.
::

    >>> c = SkyCoord(ra=[10, 11]*u.degree, dec=[41, -5]*u.degree)
    >>> c
    <SkyCoord (ICRS): (ra, dec) in deg
        [(10.0, 41.0), (11.0, -5.0)]>
    >>> c[1]
    <SkyCoord (ICRS): (ra, dec) in deg
        (11.0, -5.0)>

Coordinate access
-----------------

Once you have a coordinate object you can now access the components of that
coordinate (e.g. RA, Dec) and get a specific string representation of the full
coordinate.

The component values are accessed using aptly named attributes::

    >>> c = SkyCoord(ra=10.68458*u.degree, dec=41.26917*u.degree)
    >>> c.ra  # doctest: +FLOAT_CMP
    <Longitude 10.68458 deg>
    >>> c.ra.hour  # doctest: +FLOAT_CMP
    0.7123053333333335
    >>> c.ra.hms  # doctest: +FLOAT_CMP
    hms_tuple(h=0.0, m=42.0, s=44.299200000000525)
    >>> c.dec  # doctest: +FLOAT_CMP
    <Latitude 41.26917 deg>
    >>> c.dec.degree  # doctest: +FLOAT_CMP
    41.26917
    >>> c.dec.radian  # doctest: +FLOAT_CMP
    0.7202828960652683

Coordinates can easily be converted to strings using the
:meth:`~astropy.coordinates.SkyCoord.to_string` method::

    >>> c = SkyCoord(ra=10.68458*u.degree, dec=41.26917*u.degree)
    >>> c.to_string('decimal')
    '10.6846 41.2692'
    >>> c.to_string('dms')
    '10d41m04.488s 41d16m09.012s'
    >>> c.to_string('hmsdms')
    '00h42m44.2992s +41d16m09.012s'

For additional information see the section on :ref:`working_with_angles`.

Transformation
--------------

The simplest way to transform to a new coordinate frame is by accessing
the appropriately-named attribute.  For instance to get the coordinate in
the Galactic frame use::

    >>> c_icrs = SkyCoord(ra=10.68458*u.degree, dec=41.26917*u.degree, frame='icrs')
    >>> c_icrs.galactic  # doctest: +FLOAT_CMP
    <SkyCoord (Galactic): (l, b) in deg
        (121.174241811, -21.5728855724)>

For more control, you can use the `~astropy.coordinates.SkyCoord.transform_to`
method, which accepts a frame name, frame class, or frame instance::

    >>> c_fk5 = c_icrs.transform_to('fk5')  # c_icrs.fk5 does the same thing
    >>> c_fk5  # doctest: +FLOAT_CMP
    <SkyCoord (FK5: equinox=J2000.000): (ra, dec) in deg
        (10.6845915393, 41.2691714591)>

    >>> from astropy.coordinates import FK5
    >>> c_fk5.transform_to(FK5(equinox='J1975'))  # precess to a different equinox  # doctest: +FLOAT_CMP
    <SkyCoord (FK5: equinox=J1975.000): (ra, dec) in deg
        (10.3420913461, 41.1323211229)>

This form of `~astropy.coordinates.SkyCoord.transform_to` also makes it
straightforward to convert from celestial coordinates to
`~astropy.coordinates.AltAz` coordinates, allowing the use of |skycoord|
as a tool for planning observations.  For a more complete example of
this, see :doc:`observing-example`.

Representation
--------------

So far we have been using a spherical coordinate representation in the all the
examples, and this is the default for the built-in frames.  Frequently it is
convenient to initialize or work with a coordinate using a different
representation such as cartesian or cylindrical.  This can be done by setting
the ``representation`` for either |SkyCoord| objects or low-level frame
coordinate objects::

    >>> c = SkyCoord(x=1, y=2, z=3, unit='kpc', representation='cartesian')
    >>> c
    <SkyCoord (ICRS): (x, y, z) in kpc
        (1.0, 2.0, 3.0)>
    >>> c.x, c.y, c.z
    (<Quantity 1.0 kpc>, <Quantity 2.0 kpc>, <Quantity 3.0 kpc>)

    >>> c.representation = 'cylindrical'
    >>> c  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (rho, phi, z) in (kpc, deg, kpc)
        (2.2360679775, 63.4349488229, 3.0)>

For all the details see :ref:`astropy-skycoord-representations`.

Distance
--------

Distance from the origin (which is system-dependent, but often the Earth
center) can also be assigned to a |skycoord|. With two angles and a
distance, a unique point in 3D space is available, which also allows
conversion to the Cartesian representation of this location::

    >>> from astropy.coordinates import Distance
    >>> c = SkyCoord(ra=10.68458*u.degree, dec=41.26917*u.degree, distance=770*u.kpc)
    >>> c.cartesian.x  # doctest: +FLOAT_CMP
    <Quantity 568.7128654235232 kpc>
    >>> c.cartesian.y  # doctest: +FLOAT_CMP
    <Quantity 107.3008974042025 kpc>
    >>> c.cartesian.z  # doctest: +FLOAT_CMP
    <Quantity 507.88994291875713 kpc>

With distances assigned, |skycoord| convenience methods are more powerful, as
they can make use of the 3D information. For example::

    >>> c1 = SkyCoord(ra=10*u.degree, dec=9*u.degree, distance=10*u.pc, frame='icrs')
    >>> c2 = SkyCoord(ra=11*u.degree, dec=10*u.degree, distance=11.5*u.pc, frame='icrs')
    >>> c1.separation_3d(c2)  # doctest: +FLOAT_CMP
    <Distance 1.5228602415117989 pc>

Convenience methods
-------------------

|skycoord| defines a number of convenience methods as well, like on-sky
separation between two coordinates and catalog matching (detailed in
:ref:`astropy-coordinates-matching`)::

    >>> c1 = SkyCoord(ra=10*u.degree, dec=9*u.degree, frame='icrs')
    >>> c2 = SkyCoord(ra=11*u.degree, dec=10*u.degree, frame='fk5')
    >>> c1.separation(c2)  # Differing frames handled correctly  # doctest: +FLOAT_CMP
    <Angle 1.4045335865905868 deg>

The `astropy.coordinates` subpackage also provides a quick way to get
coordinates for named objects assuming you have an active internet
connection. The `~astropy.coordinates.SkyCoord.from_name` method of |skycoord|
uses `Sesame <http://cds.u-strasbg.fr/cgi-bin/Sesame>`_ to retrieve coordinates
for a particular named object::

    >>> SkyCoord.from_name("M42")  # doctest: +REMOTE_DATA +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec) in deg
        (83.82208, -5.39111)>

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
useful piece of background information are some
:ref:`astropy-coordinates-definitions` as they are used in
`~astropy.coordinates`.

`~astropy.coordinates` is built on a three-tiered system of objects:
representations, frames, and a high-level class.  Representations
classes are a particular way of storing a three-dimensional data point
(or points), such as Cartesian coordinates or spherical polar
coordinates. Frames are particular reference frames like FK5 or ICRS,
which may store their data in different representations, but have well-
defined transformations between each other. These transformations are
all stored in the ``astropy.coordinates.frame_transform_graph``, and new
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
   observing-example
   formatting
   matchsep
   representations
   frames
   sgr-example
   galactocentric
   definitions


In addition, another resource for the capabilities of this package is the
``astropy.coordinates.tests.test_api_ape5`` testing file. It showcases most of
the major capabilities of the package, and hence is a useful supplement to
this document.  You can see it by either looking at it directly if you
downloaded a copy of the astropy source code, or typing the following in an
IPython session::

    In [1]: from astropy.coordinates.tests import test_api_ape5
    In [2]: test_api_ape5??


.. _astropy-coordinates-seealso:

See Also
========

Some references particularly useful in understanding subtleties of the
coordinate systems implemented here include:

* `USNO Circular 179 <http://aa.usno.navy.mil/publications/docs/Circular_179.php>`_
    A useful guide to the IAU 2000/2003 work surrounding ICRS/IERS/CIRS and
    related problems in precision coordinate system work.
* `Standards Of Fundamental Astronomy <http://www.iausofa.org/>`_
    The definitive implementation of IAU-defined algorithms.  The "SOFA Tools
    for Earth Attitude" document is particularly valuable for understanding
    the latest IAU standards in detail.
* `IERS Conventions (2010) <http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html>`_
    An exhaustive reference covering the ITRS, the IAU2000 celestial coordinates
    framework, and other related details of modern coordinate conventions.
* Meeus, J. "Astronomical Algorithms"
    A valuable text describing details of a wide range of coordinate-related
    problems and concepts.


.. _astropy-coordinates-api:

Reference/API
=============

.. automodapi:: astropy.coordinates
