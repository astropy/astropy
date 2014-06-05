.. include:: references.txt

.. _astropy-coordinates-high-level:

Using the SkyCoord High-level Class
-----------------------------------

The |SkyCoord| class provides a simple and flexible user interface for
celestial coordinate representation, manipulation, and transformation between
coordinate frames.  This is a high-level class that serves as a wrapper
around the low-level coordinate frame classes like `~astropy.coordinates.ICRS`
and `~astropy.coordinates.FK5` which do most of the heavy lifting.

The key distinctions between |SkyCoord| and the low-level classes
(:doc:`frames`) are as follows:

- The |SkyCoord| object can maintain the union of frame attributes for all
  built-in and user-defined coordinate frames in the
  ``~astropy.coordinates.frame_transform_graph``.  Individual frame classes hold
  only the required attributes (e.g. equinox, observation time or observer
  location) for that frame.  This means that a transformation from
  `~astropy.coordinates.FK4` (with equinox and observation time) to
  `~astropy.coordinates.ICRS` (with neither) and back to
  `~astropy.coordinates.FK4` via the low-level classes would not remember the
  original equinox and observation time.  Since the |SkyCoord| object stores
  all attributes, such a round-trip transformation will return to the same
  coordinate object.

- The |SkyCoord| class is more flexible with inputs to accomodate a wide
  variety of user preferences and available data formats.

- The |SkyCoord| class has a number of convenience methods that are useful
  in typical analysis.

- At present, |SkyCoord| objects can use only coordinate frames that have
  transformations defined in the ``astropy.coordinates.frame_transform_graph``
  transform graph object.

Creating SkyCoord objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The |SkyCoord| class accepts a wide variety of inputs for initialization.
At a minimum these must provide one or more celestial coordinate values
with unambiguous units.  Typically one also specifies the coordinate
frame, though this is not required.

Common patterns are shown below.  In this description the values in upper
case like ``COORD`` or ``FRAME`` represent inputs which are described in detail
in the `Initialization Syntax`_ section.  Elements in square brackets like
``[unit=UNIT]`` are optional.
::

  SkyCoord(COORD, [FRAME], keyword_args ...)
  SkyCoord(LON, LAT, [frame=FRAME], [unit=UNIT], keyword_args ...)
  SkyCoord([FRAME], <lon_attr>=LON, <lat_attr>=LAT, keyword_args ...)

The examples below illustrate common ways of initializing a |SkyCoord|
object.  For a complete description of the allowed syntax see the
full coordinates documentation.  First some imports::

  >>> from astropy.coordinates import SkyCoord  # High-level coordinates
  >>> from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
  >>> from astropy.coordinates import Angle, Latitude, Longitude  # Angles
  >>> import astropy.units as u
  >>> import numpy as np

The coordinate values and frame specification can now be provided using
positional and keyword arguments.  First we show positional arguments for
RA and Dec::

  >>> SkyCoord(10, 20, unit="deg")  # No frame (no transform to other frames)
  <SkyCoord (NoFrame): ra=10.0 deg, dec=20.0 deg>

  >>> SkyCoord([1, 2, 3], [-30, 45, 8], "icrs", unit="deg")
  <SkyCoord (ICRS): (ra, dec) in deg
      [(1.0, -30.0), (2.0, 45.0), (3.0, 8.0)]>

Notice that the first example above does not specify a frame.  This object
can be used for displaying the coordinates in different formats but cannot
be used in transformations, matching catalogs, nor in computing coordinate
separations.

String inputs in common formats are acceptable, and the frame can be supplied
as either a class type like `~astropy.coordinates.FK4` or the lower-case
version of the name as a string, e.g. ``"fk4"``::

  >>> coords = ["1:12:43.2 +1:12:43", "1 12 43.2 +1 12 43"]
  >>> sc = SkyCoord(coords, FK4, unit=(u.hourangle, u.deg), obstime="J1992.21")
  >>> sc = SkyCoord(coords, 'fk4', unit='hourangle,deg', obstime="J1992.21")

  >>> sc = SkyCoord("1h12m43.2s", "+1d12m43s", Galactic)  # Units from strings
  >>> sc = SkyCoord("1h12m43.2s +1d12m43s", Galactic)  # Units from string
  >>> sc = SkyCoord("galactic", l="1h12m43.2s", b="+1d12m43s")

Astropy `~astropy.units.Quantity`-type objects are acceptable and encouraged
as a form of input::

  >>> ra = Longitude([1, 2, 3], unit=u.deg)  # Could also use Angle
  >>> dec = np.array([4.5, 5.2, 6.3]) * u.deg  # Astropy Quantity
  >>> sc = SkyCoord(ra, dec, frame='icrs')
  >>> sc = SkyCoord(ICRS, ra=ra, dec=dec, obstime='2001-01-02T12:34:56')

Finally it is possible to initialize from a low-level coordinate frame object.

  >>> c = FK4(1 * u.deg, 2 * u.deg)
  >>> sc = SkyCoord(c, obstime='J2010.11', equinox='B1965')  # Override defaults

A key subtlety highlighted here is that when low-level objects are created they have
certain default attribute values.  For instance the `~astropy.coordinates.FK4`
frame uses ``equinox='B1950.0`` and ``obstime=equinox`` as defaults.  If
this object is used to initialize a |SkyCoord| it is possible to override
the low-level object attributes that were not explicitly set.  If the
coordinate above were created with
``c = FK4(1 * u.deg, 2 * u.deg, equinox='B1960')`` then creating a |SkyCoord|
with a different ``equinox`` would raise an exception.

Initialization Syntax
"""""""""""""""""""""""

The syntax for |SkyCoord| is given below::

  SkyCoord(COORD, [FRAME | frame=FRAME], [unit=UNIT], keyword_args ...)
  SkyCoord(LON, LAT, [FRAME | frame=FRAME], [unit=UNIT], keyword_args ...)
  SkyCoord([FRAME | frame=FRAME], <lon_name>=LON, <lat_name>=LAT, [unit=UNIT],
           keyword_args ...)

**LON**, **LAT**

For spherical coordinate frames, longitude and latitude values can be specified
as separate positional arguments.  The following options are available:

- Single angle value:

  - |Quantity| object
  - Plain numeric value with ``unit`` keyword specifying the unit
  - Angle string which is formatted for :ref:`angle-creation` of
    |Longitude| or |Latitude| objects

- List or |Quantity| array or numpy array of angle values
- |Angle|, |Longitude|, or |Latitude| object, which can be scalar or
  array-valued

**COORD**

This input form uses a single object to supply coordinate data.  For the case
of spherical coordinate frames, the coordinate can include one or more
longitude and latitude pairs in one of the following ways:

- Single coordinate string with a LON and LAT value separated by a space.  The
  respective values can be any string which is formatted for
  :ref:`angle-creation` of |Longitude| or |Latitude| objects, respectively.
- List or numpy array of coordinate strings
- List of (LON, LAT) tuples, where each LON and LAT are scalars (not arrays)
- ``N x 2`` numpy or |Quantity| array of values where the first column is
  longitude and the second column is latitude, e.g.
  ``[[270, -30], [355, +85]] * u.deg``

The input can also be more generalized objects that are not necessarily
represented in the standard spherical coordinates:

- Coordinate frame object, e.g. ``FK4(1*u.deg, 2*u.deg, obstime='J2012.2')``
- |SkyCoord| object (which just makes a copy of the object)
- `~astropy.coordinates.BaseRepresentation` subclass object like
  `~astropy.coordinates.SphericalRepresentation`,
  `~astropy.coordinates.CylindricalRepresentation`,  or
  `~astropy.coordinates.CartesianRepresentation`.

**FRAME**

This can be a `~astropy.coordinates.BaseCoordinateFrame` frame
class or the corresponding string alias.  The frame classes that are built in
to astropy are `~astropy.coordinates.ICRS`, `~astropy.coordinates.FK5`,
`~astropy.coordinates.FK4`, `~astropy.coordinates.FK4NoETerms`,
`~astropy.coordinates.Galactic`, and `~astropy.coordinates.AltAz`.
The string aliases are simply lower-case versions of the class name.

If the frame is not supplied then you will see a special ``NoFrame``
identifer.  This indicates that the frame is unspecified and operations
that require comparing coordinates (even within that object) are not allowed.

**unit=UNIT**

The unit specifier can be one of the following:

- `~astropy.units.Unit` object which is an angular unit that is equivalent to
  ``Unit('radian')``
- Single string with a valid angular unit name
- 2-tuple of `~astropy.units.Unit` objects or string unit names specifying the
  LON and LAT unit respectively, e.g. ``('hourangle', 'degree')``
- Single string with two unit names separated by a comma, e.g. ``'hourangle,degree'``

If only a single unit is provided then it applies to both LON and LAT.

**Other keyword arguments**

In lieu of positional arguments to specify the longitude and latitude, the
frame-specific names can be used as keyword arguments:

*ra*, *dec*: **LON**, **LAT** values, optional
    RA and Dec for frames where this is preferred representation, including
    `~astropy.coordinates.ICRS`, `~astropy.coordinates.FK5`,
    `~astropy.coordinates.FK4`, and `~astropy.coordinates.FK4NoETerms`.

*l*, *b*:  **LON**, **LAT** values, optional
    Galactic ``l`` and ``b`` for the `~astropy.coordinates.Galactic` frame.

The following keywords can be specified for any frame:

*distance*: valid `~astropy.coordinates.Distance` initializer, optional
    Distance from reference from center to source.

*obstime*: valid `~astropy.time.Time` initializer, optional
    Time of observation

*equinox*: valid `~astropy.time.Time` initializer, optional
    Coordinate frame equinox

If custom user-defined frames are included in the transform graph and they
have additional frame attributes, then those attributes can also be
set via corresponding keyword args in the |SkyCoord| initialization.

Array operations
^^^^^^^^^^^^^^^^^

It is possible to store arrays of coordinates in a |SkyCoord| object, and
manipulations done in this way will be orders of magnitude faster than
looping over a list of individual |SkyCoord| objects::

  >>> ra = np.random.uniform(0, 360, size=1000) * u.deg
  >>> dec = np.random.uniform(-90, 90, size=1000) * u.deg

  >>> sc_list = [SkyCoord(r, d, 'icrs') for r, d in zip(ra, dec)]
  >>> timeit sc_gal_list = [c.galactic for c in sc_list]  # doctest: +SKIP
  1 loops, best of 3: 7.66 s per loop

  >>> sc = SkyCoord(ra, dec, 'icrs')
  >>> timeit sc_gal = sc.galactic  # doctest: +SKIP
  100 loops, best of 3: 8.92 ms per loop

In addition to vectorized transformations, you can do the usual array
slicing, dicing, and selection::

  >>> north_mask = sc.dec > 0
  >>> sc_north = sc[north_mask]
  >>> len(sc_north)  # doctest: +SKIP
  504
  >>> sc[2:4]  # doctest: +SKIP
  <SkyCoord (ICRS): (ra, dec) in deg
      [(304.304015..., 6.900282...),
       (322.560148..., 34.872244...)]>
  >>> sc[2]  # doctest: +SKIP
  <SkyCoord (ICRS): ra=304.304015... deg, dec=6.900282... deg>


Attributes
^^^^^^^^^^^

The |SkyCoord| object has a number of useful attributes which come in handy.
By digging through these we'll learn a little bit about |SkyCoord| and how it
works.

To begin (if you don't know already) one of the most important tools for
learning about attributes and methods of objects is "TAB-discovery".  From
within IPython you can type an object name, the period, and then the <TAB> key
to see what's available.  This can often be faster than reading the
documentation::

  >>> sc = SkyCoord(1, 2, 'icrs', unit='deg', obstime='2013-01-02 14:25:36')
  >>> sc.<TAB>  # doctest: +SKIP
  sc.cartesian                 sc.has_data                  sc.representation
  sc.data                      sc.icrs                      sc.ra
  sc.dec                       sc.is_frame_attr_default     sc.realize_frame
  sc.distance                  sc.is_transformable_to       sc.represent_as
  sc.equinox                   sc.isscalar                  sc.separation
  sc.fk4                       sc.match_to_catalog_3d       sc.separation_3d
  sc.fk4noeterms               sc.match_to_catalog_sky      sc.shape
  sc.fk5                       sc.name                      sc.spherical
  sc.frame                     sc.obstime                   sc.time_attr_names
  sc.frame_attr_names          sc.position_angle            sc.to_string
  sc.from_name                 sc.preferred_attr_names      sc.transform_to
  sc.galactic                  sc.preferred_attr_units

Here we see a bunch of stuff there but much of it should be recognizable or
easily guessed.  The most obvious may be the longitude and latitude attributes
which are named ``ra`` and ``dec`` for the ``ICRS`` frame::

  >>> sc.ra
  <Longitude 1.0 deg>
  >>> sc.dec
  <Latitude 2.0 deg>

Next notice that all the built-in frame names ``icrs``, ``galactic``, ``fk5``
``fk4``, and ``fk4noeterms`` are there.  Through the magic of Python
properties, accessing these attributes calls the object
`~astropy.coordinates.SkyCoord.transform_to` method appropriately and returns a
new |SkyCoord| object in the requested frame::

  >>> sc_gal = sc.galactic
  >>> sc_gal
  <SkyCoord (Galactic): l=99.637943... deg, b=-58.709605... deg>

Other attributes you should recognize are ``distance``, ``equinox``,
``obstime``, ``shape``.

Digger deeper
"""""""""""""""
*[Casual users can skip this section]*

After transforming to Galactic the longitude and latitude values are now
labeled ``l`` and ``b``, following the normal convention for Galactic
coordinates.  How does the object know what to call its values?  The answer
lies in some less-obvious attributes::

  >>> sc_gal.preferred_attr_names
  OrderedDict([(u'l', u'lon'), (u'b', u'lat'), (u'distance', u'distance')])

  >>> sc_gal.preferred_attr_units
  OrderedDict([(u'l', Unit("deg")), (u'b', Unit("deg"))])

  >>> sc_gal.representation
  <class 'astropy.coordinates.representation.SphericalRepresentation'>

Together these tell the object that ``l`` and ``b`` are the longitude and
latitude, and that they should both be displayed in units of degrees as
a spherical-type coordinate (and not, e.g. a cartesian coordinate).
Furthermore the frame's ``preferred_attr_names`` attribute defines
the coordinate keyword arguments that |SkyCoord| will accept.

Another important attribute is ``frame_attr_names``, which defines the
additional attributes that are required to fully define the frame::

  >>> sc_fk4 = SkyCoord(1, 2, 'fk4', unit='deg')
  >>> sc_fk4.frame_attr_names
  {u'equinox': <Time object: scale='tai' format='byear_str' value=B1950.000>,
   u'obstime': None}

The key values correspond to the defaults if no explicit value is provide by
the user.  This example shows that the `~astropy.coordinates.FK4` frame has two
attributes ``equinox`` and ``obstime`` that are required to fully define the
frame.  These are actually a reference to the class attribute of the same
name::

  >>> FK4.frame_attr_names is sc_fk4.frame_attr_names
  True

Further trickery is happening here because many of these attributes are
actually owned by the underlying coordinate ``frame`` object which does much of
the real work.  This is the middle layer in the three-tiered system of objects:
representation (spherical, cartesian, etc.), frame (aka low-level frame class),
and |SkyCoord| (aka high-level class)::

  >>> sc.frame
  <ICRS Coordinate: ra=1.0 deg, dec=2.0 deg>

  >>> sc.preferred_attr_units is sc.frame.preferred_attr_units
  True

  >>> sc.frame.<TAB>  # doctest: +SKIP
  sc.frame.cartesian                 sc.frame.preferred_attr_units
  sc.frame.data                      sc.frame.ra
  sc.frame.dec                       sc.frame.realize_frame
  sc.frame.distance                  sc.frame.represent_as
  sc.frame.frame_attr_names          sc.frame.representation
  sc.frame.has_data                  sc.frame.separation
  sc.frame.is_frame_attr_default     sc.frame.separation_3d
  sc.frame.is_transformable_to       sc.frame.shape
  sc.frame.isscalar                  sc.frame.spherical
  sc.frame.name                      sc.frame.time_attr_names
  sc.frame.preferred_attr_names      sc.frame.transform_to

  >>> sc.frame.name
  'icrs'

The |SkyCoord| object exposes the ``frame`` object attributes as its own.  Though
it might seem a tad confusing at first, this a good thing because it makes
|SkyCoord| objects and `~astropy.coordinates.BaseCoordinateFrame` objects
behave very similarly and most routines can accept either one as input without
much bother (duck typing!).

The lowest layer in the stack is the abstract
`~astropy.coordinates.UnitSphericalRepresentation` object:

  >>> sc_gal.frame.data
  <UnitSphericalRepresentation lon=1.739010... rad, lat=-1.024675... rad>

Transformations
^^^^^^^^^^^^^^^^^

The topic of transformations is covered in detail in the section on
:ref:`astropy-coordinates-transforming`.

For completeness here we will give some simple examples.  Once you've defined
your coordinates and the reference frame, you can transform from that frame to
another frame.  You can do this a few different ways: if you just want the
default version of that frame, you can use attribute-style access (as mentioned
previously).  For more control, you can use the
`~astropy.coordinates.SkyCoord.transform_to` method, which accepts a frame
name, frame class, frame instance, or |SkyCoord|::

  >>> from astropy.coordinates import FK5
  >>> sc = SkyCoord(1, 2, 'icrs', unit='deg')
  >>> sc.galactic
  <SkyCoord (Galactic): l=99.637943... deg, b=-58.709605... deg>

  >>> sc.transform_to('fk5')  # Same as sc.fk5 and sc.transform_to(FK5)
  <SkyCoord (FK5): equinox=J2000.000, ra=1.00000655... deg, dec=2.00000243... deg>

  >>> sc.transform_to(FK5(equinox='J1975'))  # Transform to FK5 with a different equinox
  <SkyCoord (FK5): equinox=J1975.000, ra=0.679672... deg, dec=1.860830... deg>

Transforming to a |SkyCoord| instance is an easy way of ensuring that two
coordinates are in the exact same reference frame::

  >>> sc2 = SkyCoord(3, 4, 'fk4', unit='deg', obstime='J1978.123', equinox='B1960.0')
  >>> sc.transform_to(sc2)
  <SkyCoord (FK4): equinox=B1960.000, obstime=J1978.123, ra=0.487263... deg, dec=1.777316... deg>

Convenience methods
^^^^^^^^^^^^^^^^^^^^

A number of convenience methods are available, and you are encouraged to read
the available docstrings below:

- `~astropy.coordinates.SkyCoord.match_to_catalog_sky`,
- `~astropy.coordinates.SkyCoord.match_to_catalog_3d`,
- `~astropy.coordinates.SkyCoord.position_angle`,
- `~astropy.coordinates.SkyCoord.separation`,
- `~astropy.coordinates.SkyCoord.separation_3d`

Addition information and examples can be found in the section on
:ref:`astropy-coordinates-separations-matching`.

