.. include:: references.txt

.. _astropy-coordinates-high-level:

Using the SkyCoord High-level Class
-----------------------------------

The |SkyCoord| class provides a simple and flexible user interface for
celestial coordinate representation, manipulation, and transformation between
coordinate frames.  This is a high-level class that serves as a wrapper
around the low-level coordinate frame classes like `~astropy.coordinates.ICRS`
and `~astropy.coordinates.FK5` which actually do most of the heavy lifting.

The key distinctions between |SkyCoord| and the low-level classes (which are
derived from `~astropy.coordinates.BaseCoordinateFrame`) are as follows:

- The |SkyCoord| object can maintain the union of frame attributes for all
  built-in and user-defined coordinate frames in the
  `~astropy.coordinates.frame_transform_graph`.  Individual frame classes hold
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
  
Creating SkyCoord objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `SkyCoord` class accepts a wide variety of inputs for initialization.
At a minimum these must provide one or more celestial coordinate values
with unambiguous units.  Typically one also specifies the coordinate
frame, though this is not required.

Common patterns are shown below.  In this description the values in upper
case like ``COORD`` or ``FRAME`` represent inputs which are described in detail
in the `Initialization Syntax`_ section.  Elements in square brackets like
``[unit=UNIT]`` are optional.

  SkyCoord(COORD, [FRAME], keyword_args ...)
  SkyCoord(LON, LAT, [frame=FRAME], [unit=UNIT], keyword_args ...)
  SkyCoord([FRAME], <lon_attr>=LON, <lat_attr>=LAT, keyword_args ...)

The examples below illustrate common ways of initializing a `SkyCoord`
object.  For a complete description of the allowed syntax see the
full coordinates documentation.  First some imports::

  >>> from astropy.coordinates import SkyCoord  # High-level coordinates
  >>> from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
  >>> from astropy.coordinates import Angle, Latitude, Longitude  # Angles
  >>> import astropy.units as u

The coordinate values and frame specification can now be provided using
positional and keyword arguments.  First we show positional arguments for
RA and Dec::

  >>> SkyCoord(10, 20, unit="deg")  # No frame (no transform to other frames)
  <SkyCoord (NoFrame): ra=10.0 deg, dec=20.0 deg>

  >>> sc = SkyCoord([1, 2, 3], [-30, 45, 8], "icrs", unit="deg")
  <SkyCoord (ICRS): (ra, dec) in deg
      [(1.0, -30.0), (2.0, 45.0), (3.0, 8.0)]>

Notice that the first example above does not specify a frame.  This object
can be used for displaying the coordinates in different formats but cannot
be used in transformations, matching catalogs, nor in computing coordinate
separations.

String inputs in common formats are acceptable, and the frame can be supplied
as either a class type like `~astropy.coordinates.FK4` or the lower-case
version of the name as a string, e.g. `"fk4"`::

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

The syntax for |SkyCoord| is shown below::

  SkyCoord(COORD, [FRAME | frame=FRAME], [unit=UNIT], keyword_args ...)
  SkyCoord(LON, LAT, [FRAME | frame=FRAME], [unit=UNIT], keyword_args ...)
  SkyCoord([FRAME | frame=FRAME], <lon_name>=LON, <lat_name>=LAT, [unit=UNIT],
           keyword_args ...)

**LON**, **LAT**

Longitude and latitude values can be specified as separate positional
arguments.  The following options are available:

- Single angle value:

  - |Quantity| object
  - Plain numeric value with `unit` keyword specifying the unit
  - Angle string which is formatted for :ref:`angle-creation` of
    |Longitude| or |Latitude| objects

- List or |Quantity| array or numpy array of angle values
- |Angle|, |Longitude|, or |Latitude| object, which can be scalar or
  array-valued

**COORD**

A coordinate is an object that supplies one or more longitude and latitude
pairs in one of the following ways:

- Single coordinate string with a LON and LAT value separated by a space.  The
  respective values can be any string which is formatted for
  :ref:`angle-creation` of |Longitude| or |Latitude| objects, respectively.
- List or numpy array of coordinate strings
- List of (LON, LAT) tuples, where each LON and LAT are scalars (not arrays)
- ``N x 2`` numpy or |Quantity| array of values where the first column is
  longitude and the second column is latitude, e.g.
  ``[[270, -30], [355, +85]] * u.deg``
- Coordinate frame object, e.g. ``FK4(1*u.deg, 2*u.deg, obstime='J2012.2')``
- |SkyCoord| object (which just makes a copy of the object)

**FRAME**

This can be a `~astropy.coordinates.BaseCoordinateFrame` frame
class or the corresponding string alias.  The frame classes that are built in
to astropy are `~astropy.coordinates.ICRS`, `~astropy.coordinates.FK5`,
`~astropy.coordinates.FK4`, `~astropy.coordinates.FK4NoETerms`,
`~astropy.coordinates.Galactic`, and `~astropy.coordinates.AltAz`.
The string aliases are simply lower-case versions of the class name.

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
    `ICRS`, `FK5`, `FK4`, and `FK4NoETerms`.

*l*, *b*:  **LON**, **LAT** values, optional
    Galactic `l` and `b` for the `Galactic` frame.

The following keywords can be specified for any frame:

*distance*: valid `~astropy.coordinates.Distance` initializer, optional
    Distance from reference from center to source.

*obstime*: valid `~astropy.time.Time` initializer, optional
    Time of observation

*equinox*: valid `~astropy.time.Time` initializer, optional
    Coordinate frame equinox


Transformations
^^^^^^^^^^^^^^^^^

* mention transformation, but refer mostly to :ref:`astropy-coordinates-transforming`,

Convenience methods
^^^^^^^^^^^^^^^^^^^^

* matching/separation convinience methods (:ref:`astropy-coordinates-matching`)

