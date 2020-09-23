.. include:: references.txt

.. We call EarthLocation.of_site here first to force the downloading
.. of sites.json so that future doctest output isn't cluttered with
.. "Downloading ... [done]". This can be removed once we have a better
.. way of ignoring output lines based on pattern-matching, e.g.:
.. https://github.com/astropy/pytest-doctestplus/issues/11

.. testsetup::
    >>> from astropy.coordinates import EarthLocation
    >>> EarthLocation.of_site('greenwich') # doctest: +IGNORE_OUTPUT +IGNORE_WARNINGS

.. _astropy-coordinates:

*******************************************************
Astronomical Coordinate Systems (`astropy.coordinates`)
*******************************************************

Introduction
============

The `~astropy.coordinates` package provides classes for representing a variety
of celestial/spatial coordinates and their velocity components, as well as tools
for converting between common coordinate systems in a uniform way.

Getting Started
===============

The best way to start using `~astropy.coordinates` is to use the |skycoord|
class. |skycoord| objects are instantiated by passing in positions (and
optional velocities) with specified units and a coordinate frame. Sky positions
are commonly passed in as `~astropy.units.Quantity` objects and the frame is
specified with the string name. As an example of creating a |skycoord| to
represent an ICRS (Right ascension [RA], Declination [Dec]) sky position::

    >>> from astropy import units as u
    >>> from astropy.coordinates import SkyCoord
    >>> c = SkyCoord(ra=10.625*u.degree, dec=41.2*u.degree, frame='icrs')

The initializer for |skycoord| is very flexible and supports inputs provided in
a number of convenient formats. The following ways of initializing a coordinate
are all equivalent to the above::

    >>> c = SkyCoord(10.625, 41.2, frame='icrs', unit='deg')
    >>> c = SkyCoord('00h42m30s', '+41d12m00s', frame='icrs')
    >>> c = SkyCoord('00h42.5m', '+41d12m')
    >>> c = SkyCoord('00 42 30 +41 12 00', unit=(u.hourangle, u.deg))
    >>> c = SkyCoord('00:42.5 +41:12', unit=(u.hourangle, u.deg))
    >>> c  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec) in deg
        (10.625, 41.2)>

The examples above illustrate a few rules to follow when creating a
coordinate object:

- Coordinate values can be provided either as unnamed positional arguments or
  via keyword arguments like ``ra`` and ``dec``, or  ``l`` and ``b`` (depending
  on the frame).
- The coordinate ``frame`` keyword is optional because it defaults to
  `~astropy.coordinates.ICRS`.
- Angle units must be specified for all components, either by passing in a
  `~astropy.units.Quantity` object (e.g., ``10.5*u.degree``), by including them
  in the value (e.g., ``'+41d12m00s'``), or via the ``unit`` keyword.

|skycoord| and all other `~astropy.coordinates` objects also support
array coordinates. These work in the same way as single-value coordinates, but
they store multiple coordinates in a single object. When you are going
to apply the same operation to many different coordinates (say, from a
catalog), this is a better choice than a list of |skycoord| objects,
because it will be *much* faster than applying the operation to each
|skycoord| in a ``for`` loop. Like the underlying `~numpy.ndarray` instances
that contain the data, |skycoord| objects can be sliced, reshaped, etc.::

    >>> c = SkyCoord(ra=[10, 11, 12, 13]*u.degree, dec=[41, -5, 42, 0]*u.degree)
    >>> c  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec) in deg
        [(10., 41.), (11., -5.), (12., 42.), (13.,  0.)]>
    >>> c[1]  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec) in deg
        (11., -5.)>
    >>> c.reshape(2, 2)  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec) in deg
        [[(10., 41.), (11., -5.)],
         [(12., 42.), (13.,  0.)]]>

Coordinate Access
-----------------

Once you have a coordinate object you can access the components of that
coordinate (e.g., RA, Dec) to get string representations of the full
coordinate.

The component values are accessed using (typically lowercase) named attributes
that depend on the coordinate frame (e.g., ICRS, Galactic, etc.). For the
default, ICRS, the coordinate component names are ``ra`` and ``dec``::

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

Coordinates can be converted to strings using the
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

One convenient way to transform to a new coordinate frame is by accessing
the appropriately named attribute. For instance, to get the coordinate in
the `~astropy.coordinates.Galactic` frame use::

    >>> c_icrs = SkyCoord(ra=10.68458*u.degree, dec=41.26917*u.degree, frame='icrs')
    >>> c_icrs.galactic  # doctest: +FLOAT_CMP
    <SkyCoord (Galactic): (l, b) in deg
        (121.17424181, -21.57288557)>

For more control, you can use the `~astropy.coordinates.SkyCoord.transform_to`
method, which accepts a frame name, frame class, or frame instance::

    >>> c_fk5 = c_icrs.transform_to('fk5')  # c_icrs.fk5 does the same thing
    >>> c_fk5  # doctest: +FLOAT_CMP
    <SkyCoord (FK5: equinox=J2000.000): (ra, dec) in deg
        (10.68459154, 41.26917146)>

    >>> from astropy.coordinates import FK5
    >>> c_fk5.transform_to(FK5(equinox='J1975'))  # precess to a different equinox  # doctest: +FLOAT_CMP
    <SkyCoord (FK5: equinox=J1975.000): (ra, dec) in deg
        (10.34209135, 41.13232112)>

This form of `~astropy.coordinates.SkyCoord.transform_to` also makes it
possible to convert from celestial coordinates to
`~astropy.coordinates.AltAz` coordinates, allowing the use of |skycoord|
as a tool for planning observations. For a more complete example of
this, see :ref:`sphx_glr_generated_examples_coordinates_plot_obs-planning.py`.

Some coordinate frames such as `~astropy.coordinates.AltAz` require Earth
rotation information (UT1-UTC offset and/or polar motion) when transforming
to/from other frames. These Earth rotation values are automatically downloaded
from the International Earth Rotation and Reference Systems (IERS) service when
required. See :ref:`utils-iers` for details of this process.

Representation
--------------

So far we have been using a spherical coordinate representation in all of our
examples, and this is the default for the built-in frames. Frequently it is
convenient to initialize or work with a coordinate using a different
representation such as Cartesian or Cylindrical. This can be done by setting
the ``representation_type`` for either |skycoord| objects or low-level frame
coordinate objects::

    >>> c = SkyCoord(x=1, y=2, z=3, unit='kpc', representation_type='cartesian')
    >>> c  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (x, y, z) in kpc
        (1., 2., 3.)>
    >>> c.x, c.y, c.z  # doctest: +FLOAT_CMP
    (<Quantity 1. kpc>, <Quantity 2. kpc>, <Quantity 3. kpc>)

    >>> c.representation_type = 'cylindrical'
    >>> c  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (rho, phi, z) in (kpc, deg, kpc)
        (2.23606798, 63.43494882, 3.)>

For all of the details see :ref:`astropy-skycoord-representations`.

Distance
--------

|skycoord| and the individual frame classes also support specifying a distance
from the frame origin. The origin depends on the particular coordinate frame;
this can be, for example, centered on the earth, centered on the solar system
barycenter, etc. Two angles and a distance specify a unique point in 3D space,
which also allows converting the coordinates to a Cartesian representation::

    >>> c = SkyCoord(ra=10.68458*u.degree, dec=41.26917*u.degree, distance=770*u.kpc)
    >>> c.cartesian.x  # doctest: +FLOAT_CMP
    <Quantity 568.71286542 kpc>
    >>> c.cartesian.y  # doctest: +FLOAT_CMP
    <Quantity 107.3008974 kpc>
    >>> c.cartesian.z  # doctest: +FLOAT_CMP
    <Quantity 507.88994292 kpc>

With distances assigned, |skycoord| convenience methods are more powerful, as
they can make use of the 3D information. For example, to compute the physical,
3D separation between two points in space::

    >>> c1 = SkyCoord(ra=10*u.degree, dec=9*u.degree, distance=10*u.pc, frame='icrs')
    >>> c2 = SkyCoord(ra=11*u.degree, dec=10*u.degree, distance=11.5*u.pc, frame='icrs')
    >>> c1.separation_3d(c2)  # doctest: +FLOAT_CMP
    <Distance 1.52286024 pc>

Convenience Methods
-------------------

|skycoord| defines a number of convenience methods that support, for example,
computing on-sky (i.e., angular) and 3D separations between two coordinates::

    >>> c1 = SkyCoord(ra=10*u.degree, dec=9*u.degree, frame='icrs')
    >>> c2 = SkyCoord(ra=11*u.degree, dec=10*u.degree, frame='fk5')
    >>> c1.separation(c2)  # Differing frames handled correctly  # doctest: +FLOAT_CMP
    <Angle 1.40453359 deg>

Or cross-matching catalog coordinates (detailed in
:ref:`astropy-coordinates-matching`)::

    >>> target_c = SkyCoord(ra=10*u.degree, dec=9*u.degree, frame='icrs')
    >>> # read in coordinates from a catalog...
    >>> catalog_c = ... # doctest: +SKIP
    >>> idx, sep, _ = target_c.match_to_catalog_sky(catalog_c) # doctest: +SKIP

The `astropy.coordinates` sub-package also provides a quick way to get
coordinates for named objects, assuming you have an active internet
connection. The `~astropy.coordinates.SkyCoord.from_name` method of |skycoord|
uses `Sesame <http://cds.u-strasbg.fr/cgi-bin/Sesame>`_ to retrieve coordinates
for a particular named object::

    >>> SkyCoord.from_name("PSR J1012+5307")  # doctest: +REMOTE_DATA +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec) in deg
        (153.1393271, 53.117343)>

In some cases, the coordinates are embedded in the catalog name of the object.
For such object names, `~astropy.coordinates.SkyCoord.from_name` is able
to parse the coordinates from the name if given the ``parse=True`` option.
For slow connections, this may be much faster than a sesame query for the same
object name. It's worth noting, however, that the coordinates extracted in this
way may differ from the database coordinates by a few deci-arcseconds, so only
use this option if you do not need sub-arcsecond accuracy for your coordinates::

    >>> SkyCoord.from_name("CRTS SSS100805 J194428-420209", parse=True)  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec) in deg
        (296.11666667, -42.03583333)>


For sites (primarily observatories) on the Earth, `astropy.coordinates` provides
a quick way to get an `~astropy.coordinates.EarthLocation` - the
`~astropy.coordinates.EarthLocation.of_site` method::

    >>> from astropy.coordinates import EarthLocation
    >>> EarthLocation.of_site('Apache Point Observatory')  # doctest: +REMOTE_DATA +FLOAT_CMP
    <EarthLocation (-1463969.30185172, -5166673.34223433,  3434985.71204565) m>

To see the list of site names available, use
:func:`astropy.coordinates.EarthLocation.get_site_names`.

For arbitrary Earth addresses (e.g., not observatory sites), use the
`~astropy.coordinates.EarthLocation.of_address` classmethod. Any address passed
to this function uses Google maps to retrieve the latitude and longitude and can
also (optionally) query Google maps to get the height of the location. As with
Google maps, this works with fully specified addresses, location names, city
names, etc.:

.. doctest-skip::

    >>> EarthLocation.of_address('1002 Holy Grail Court, St. Louis, MO')
    <EarthLocation (-26726.98216371, -4997009.8604809, 3950271.16507911) m>
    >>> EarthLocation.of_address('1002 Holy Grail Court, St. Louis, MO',
    ...                          get_height=True)
    <EarthLocation (-26727.6272786, -4997130.47437768, 3950367.15622108) m>
    >>> EarthLocation.of_address('Danbury, CT')
    <EarthLocation ( 1364606.64511651, -4593292.9428273,  4195415.93695139) m>

.. note::
    `~astropy.coordinates.SkyCoord.from_name`,
    `~astropy.coordinates.EarthLocation.of_site`, and
    `~astropy.coordinates.EarthLocation.of_address` are for convenience, and
    hence are by design relatively low precision. If you need more precise coordinates for an
    object you should find the appropriate reference and input the coordinates
    manually, or use more specialized functionality like that in the `astroquery
    <http://www.astropy.org/astroquery/>`_ or `astroplan
    <https://astroplan.readthedocs.io/>`_ affiliated packages.

    Also note that these methods retrieve data from the internet to
    determine the celestial or Earth coordinates. The online data may be
    updated, so if you need to guarantee that your scripts are reproducible
    in the long term, see the :doc:`remote_methods` section.

This functionality can be combined to do more complicated tasks like computing
barycentric corrections to radial velocity observations (also a supported
high-level |skycoord| method - see :ref:`astropy-coordinates-rv-corrs`)::

    >>> from astropy.time import Time
    >>> obstime = Time('2017-2-14')
    >>> target = SkyCoord.from_name('M31')  # doctest: +REMOTE_DATA
    >>> keck = EarthLocation.of_site('Keck')  # doctest: +REMOTE_DATA
    >>> target.radial_velocity_correction(obstime=obstime, location=keck).to('km/s')  # doctest: +REMOTE_DATA +FLOAT_CMP
    <Quantity -22.359784554780255 km / s>

Velocities (Proper Motions and Radial Velocities)
-------------------------------------------------

In addition to positional coordinates, `~astropy.coordinates` supports storing
and transforming velocities. These are available both via the lower-level
:doc:`coordinate frame classes <frames>`, and (new in v3.0) via  |skycoord|
objects::

    >>> sc = SkyCoord(1*u.deg, 2*u.deg, radial_velocity=20*u.km/u.s)
    >>> sc  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec) in deg
        (1., 2.)
     (radial_velocity) in km / s
        (20.,)>

For more details on velocity support (and limitations), see the
:doc:`velocities` page.

.. _astropy-coordinates-overview:

Overview of `astropy.coordinates` Concepts
==========================================

.. note ::
    The `~astropy.coordinates` package from v0.4 onward builds from
    previous versions of  the package, and more detailed information and
    justification of the design is available in `APE (Astropy Proposal for Enhancement) 5 <https://github.com/astropy/astropy- APEs/blob/master/APE5.rst>`_.

Here we provide an overview of the package and associated framework.
This background information is not necessary for using `~astropy.coordinates`,
particularly if you use the |skycoord| high-level class, but it is helpful for
more advanced usage, particularly creating your own frame, transformations, or
representations. Another useful piece of background information are some
:ref:`astropy-coordinates-definitions` as they are used in
`~astropy.coordinates`.

`~astropy.coordinates` is built on a three-tiered system of objects:
representations, frames, and a high-level class. Representations
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
`~astropy.coordinates`. It allows representations, frames, and
transformations to be defined or extended separately, while still
preserving the high-level capabilities and ease-of-use of the |skycoord|
class.

.. topic:: Examples:

    See :ref:`sphx_glr_generated_examples_coordinates_plot_obs-planning.py` for
    an example of using the `~astropy.coordinates` functionality to prepare for
    an observing run.

Using `astropy.coordinates`
===========================

More detailed information on using the package is provided on separate pages,
listed below.

.. toctree::
   :maxdepth: 1

   angles
   skycoord
   transforming
   solarsystem
   formatting
   matchsep
   representations
   frames
   velocities
   apply_space_motion
   galactocentric
   remote_methods
   definitions
   inplace


In addition, another resource for the capabilities of this package is the
``astropy.coordinates.tests.test_api_ape5`` testing file. It showcases most of
the major capabilities of the package, and hence is a useful supplement to
this document. You can see it by either downloading a copy of the Astropy
source code, or typing the following in an IPython session::

    In [1]: from astropy.coordinates.tests import test_api_ape5
    In [2]: test_api_ape5??


.. note that if this section gets too long, it should be moved to a separate
   doc page - see the top of performance.inc.rst for the instructions on how to
   do that
.. include:: performance.inc.rst

.. _astropy-coordinates-seealso:

See Also
========

Some references that are particularly useful in understanding subtleties of the
coordinate systems implemented here include:

* `USNO Circular 179 <https://arxiv.org/abs/astro-ph/0602086>`_
    A useful guide to the IAU 2000/2003 work surrounding ICRS/IERS/CIRS and
    related problems in precision coordinate system work.
* `Standards Of Fundamental Astronomy <http://www.iausofa.org/>`_
    The definitive implementation of IAU-defined algorithms. The "SOFA Tools
    for Earth Attitude" document is particularly valuable for understanding
    the latest IAU standards in detail.
* `IERS Conventions (2010) <https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html>`_
    An exhaustive reference covering the ITRS, the IAU2000 celestial coordinates
    framework, and other related details of modern coordinate conventions.
* Meeus, J. "Astronomical Algorithms"
    A valuable text describing details of a wide range of coordinate-related
    problems and concepts.


Built-in Frame Classes
======================

.. automodapi:: astropy.coordinates.builtin_frames
    :skip: make_transform_graph_docs
    :no-inheritance-diagram:


.. _astropy-coordinates-api:

Reference/API
=============

.. automodapi:: astropy.coordinates
