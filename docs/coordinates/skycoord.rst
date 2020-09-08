.. include:: references.txt

.. _astropy-coordinates-high-level:

Using the SkyCoord High-Level Class
***********************************

The |SkyCoord| class provides a simple and flexible user interface for
celestial coordinate representation, manipulation, and transformation between
coordinate frames. This is a high-level class that serves as a wrapper
around the low-level coordinate frame classes like `~astropy.coordinates.ICRS`
and `~astropy.coordinates.FK5` which do most of the heavy lifting.

The key distinctions between |SkyCoord| and the low-level classes
(:doc:`frames`) are as follows:

- The |SkyCoord| object can maintain the union of frame attributes for all
  built-in and user-defined coordinate frames in the
  ``astropy.coordinates.frame_transform_graph``. Individual frame classes hold
  only the required attributes (e.g., equinox, observation time, or observer
  location) for that frame. This means that a transformation from
  `~astropy.coordinates.FK4` (with equinox and observation time) to
  `~astropy.coordinates.ICRS` (with neither) and back to
  `~astropy.coordinates.FK4` via the low-level classes would not remember the
  original equinox and observation time. Since the |SkyCoord| object stores
  all attributes, such a round-trip transformation will return to the same
  coordinate object.

- The |SkyCoord| class is more flexible with inputs to accommodate a wide
  variety of user preferences and available data formats, whereas the frame
  classes expect to receive Quantity-like objects with angular units.

- The |SkyCoord| class has a number of convenience methods that are useful
  in typical analysis.

- At present, |SkyCoord| objects can use only coordinate frames that have
  transformations defined in the ``astropy.coordinates.frame_transform_graph``
  transform graph object.

Creating SkyCoord Objects
=========================

The |SkyCoord| class accepts a wide variety of inputs for initialization.
At a minimum, these must provide one or more celestial coordinate values
with unambiguous units. Typically you must also specify the coordinate
frame, though this is not required.

Common patterns are shown below. In this description the values in upper
case like ``COORD`` or ``FRAME`` represent inputs which are described in detail
in the `Initialization Syntax`_ section. Elements in square brackets like
``[unit=UNIT]`` are optional.
::

  SkyCoord(COORD, [FRAME], keyword_args ...)
  SkyCoord(LON, LAT, [frame=FRAME], [unit=UNIT], keyword_args ...)
  SkyCoord([FRAME], <lon_attr>=LON, <lat_attr>=LAT, keyword_args ...)

The examples below illustrate common ways of initializing a |SkyCoord| object.
These all reflect initializing using spherical coordinates, which is the
default for all built-in frames. In order to understand working with coordinates
using a different representation, such as Cartesian or cylindrical, see the
section on `Representations`_. First, some imports::

  >>> from astropy.coordinates import SkyCoord  # High-level coordinates
  >>> from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
  >>> from astropy.coordinates import Angle, Latitude, Longitude  # Angles
  >>> import astropy.units as u
  >>> import numpy as np

The coordinate values and frame specification can now be provided using
positional and keyword arguments. First we show positional arguments for
RA and Dec::

  >>> SkyCoord(10, 20, unit='deg')  # Defaults to ICRS  # doctest: +FLOAT_CMP
  <SkyCoord (ICRS): (ra, dec) in deg
      (10., 20.)>

  >>> SkyCoord([1, 2, 3], [-30, 45, 8], frame='icrs', unit='deg')  # doctest: +FLOAT_CMP
  <SkyCoord (ICRS): (ra, dec) in deg
      [(1., -30.), (2., 45.), (3.,   8.)]>

Notice that the first example above does not explicitly give a frame. In
this case, the default is taken to be the ICRS system (approximately
correct for "J2000" equatorial coordinates). It is always better to
explicitly specify the frame when it is known to be ICRS, however, as
anyone reading the code will be better able to understand the intent.

String inputs in common formats are acceptable, and the frame can be supplied
as either a class type like `~astropy.coordinates.FK4`, an instance of a
frame class, a `~astropy.coordinates.SkyCoord` instance (from which the frame
will be extracted), or the lowercase version of a frame name as a string, for
example, ``"fk4"``::

  >>> coords = ["1:12:43.2 +1:12:43", "1 12 43.2 +1 12 43"]
  >>> sc = SkyCoord(coords, frame=FK4, unit=(u.hourangle, u.deg), obstime="J1992.21")
  >>> sc = SkyCoord(coords, frame=FK4(obstime="J1992.21"), unit=(u.hourangle, u.deg))
  >>> sc = SkyCoord(coords, frame='fk4', unit='hourangle,deg', obstime="J1992.21")

  >>> sc = SkyCoord("1h12m43.2s", "+1d12m43s", frame=Galactic)  # Units from strings
  >>> sc = SkyCoord("1h12m43.2s +1d12m43s", frame=Galactic)  # Units from string
  >>> sc = SkyCoord(l="1h12m43.2s", b="+1d12m43s", frame='galactic')
  >>> sc = SkyCoord("1h12.72m +1d12.71m", frame='galactic')

Note that frame instances with data and `~astropy.coordinates.SkyCoord`
instances can only be passed as frames using the ``frame=`` keyword argument
and not as positional arguments.

For representations that have ``ra`` and ``dec`` attributes you can supply a
coordinate string in a number of other common formats. Examples include::

  >>> sc = SkyCoord("15h17+89d15")
  >>> sc = SkyCoord("275d11m15.6954s+17d59m59.876s")
  >>> sc = SkyCoord("8 00 -5 00.6", unit=(u.hour, u.deg))
  >>> sc = SkyCoord("J080000.00-050036.00", unit=(u.hour, u.deg))
  >>> sc = SkyCoord("J1874221.31+122328.03", unit=u.deg)

Astropy `~astropy.units.Quantity`-type objects are acceptable and encouraged
as a form of input::

  >>> ra = Longitude([1, 2, 3], unit=u.deg)  # Could also use Angle
  >>> dec = np.array([4.5, 5.2, 6.3]) * u.deg  # Astropy Quantity
  >>> sc = SkyCoord(ra, dec, frame='icrs')
  >>> sc = SkyCoord(ra=ra, dec=dec, frame=ICRS, obstime='2001-01-02T12:34:56')

Finally, it is possible to initialize from a low-level coordinate frame object.

  >>> c = FK4(1 * u.deg, 2 * u.deg)
  >>> sc = SkyCoord(c, obstime='J2010.11', equinox='B1965')  # Override defaults

A key subtlety highlighted here is that when low-level objects are created they
have certain default attribute values. For instance, the
`~astropy.coordinates.FK4` frame uses ``equinox='B1950.0`` and
``obstime=equinox`` as defaults. If this object is used to initialize a
|SkyCoord| it is possible to override the low-level object attributes that were
not explicitly set. If the coordinate above were created with
``c = FK4(1 * u.deg, 2 * u.deg, equinox='B1960')`` then creating a |SkyCoord|
with a different ``equinox`` would raise an exception.

Initialization Syntax
---------------------

For spherical representations, which are the most common and are the default
input format for all built-in frames, the syntax for |SkyCoord| is given
below::

  SkyCoord(COORD, [FRAME | frame=FRAME], [unit=UNIT], keyword_args ...)
  SkyCoord(LON, LAT, [DISTANCE], [FRAME | frame=FRAME], [unit=UNIT], keyword_args ...)
  SkyCoord([FRAME | frame=FRAME], <lon_name>=LON, <lat_name>=LAT, [unit=UNIT],
           keyword_args ...)

In the above description, elements in all capital letters (e.g., ``FRAME``)
describe a user input of that element type. Elements in square brackets are
optional. For nonspherical inputs, see the `Representations`_ section.


**LON**, **LAT**

Longitude and latitude value can be specified as separate positional arguments.
The following options are available for longitude and latitude:

- Single angle value:

  - |Quantity| object
  - Plain numeric value with ``unit`` keyword specifying the unit
  - Angle string which is formatted for :ref:`angle-creation` of
    |Longitude| or |Latitude| objects

- List or |Quantity| array, or NumPy array of angle values
- |Angle|, |Longitude|, or |Latitude| object, which can be scalar or
  array-valued

.. note::

    While |SkyCoord| is flexible with respect to specifying longitude and
    latitude component inputs, the frame classes expect to receive
    |Quantity|-like objects with angular units (i.e., |Angle| or |Quantity|).
    For example, when specifying components, the frame classes (e.g., ``ICRS``)
    must be created as

        >>> ICRS(0 * u.deg, 0 * u.deg) # doctest: +FLOAT_CMP
        <ICRS Coordinate: (ra, dec) in deg
            (0., 0.)>

    and other methods of flexible initialization (that work with |SkyCoord|)
    will not work

        >>> ICRS(0, 0, unit=u.deg) # doctest: +SKIP
        UnitTypeError: Longitude instances require units equivalent to 'rad', but no unit was given.

**DISTANCE**

The distance to the object from the frame center can be optionally specified:

- Single distance value:

  - |Quantity| or `~astropy.coordinates.Distance` object
  - Plain numeric value for a dimensionless distance
  - Plain numeric value with ``unit`` keyword specifying the unit

- List, or |Quantity|, or `~astropy.coordinates.Distance` array, or NumPy array
  of angle values

**COORD**

This input form uses a single object to supply coordinate data. For the case
of spherical coordinate frames, the coordinate can include one or more
longitude and latitude pairs in one of the following ways:

- Single coordinate string with a LON and LAT value separated by a space. The
  respective values can be any string which is formatted for
  :ref:`angle-creation` of |Longitude| or |Latitude| objects, respectively.
- List or NumPy array of such coordinate strings.
- List of (LON, LAT) tuples, where each LON and LAT are scalars (not arrays).
- ``N x 2`` NumPy or |Quantity| array of values where the first column is
  longitude and the second column is latitude, for example,
  ``[[270, -30], [355, +85]] * u.deg``.
- List of (LON, LAT, DISTANCE) tuples.
- ``N x 3`` NumPy or |Quantity| array of values where columns are
  longitude, latitude, and distance, respectively.

The input can also be more generalized objects that are not necessarily
represented in the standard spherical coordinates:

- Coordinate frame object (e.g., ``FK4(1*u.deg, 2*u.deg, obstime='J2012.2')``).
- |SkyCoord| object (which just makes a copy of the object).
- `~astropy.coordinates.BaseRepresentation` subclass object like
  `~astropy.coordinates.SphericalRepresentation`,
  `~astropy.coordinates.CylindricalRepresentation`, or
  `~astropy.coordinates.CartesianRepresentation`.

**FRAME**

This can be a `~astropy.coordinates.BaseCoordinateFrame` frame class, an
instance of such a class, or the corresponding string alias. The frame
classes that are built in to Astropy are `~astropy.coordinates.ICRS`,
`~astropy.coordinates.FK5`, `~astropy.coordinates.FK4`,
`~astropy.coordinates.FK4NoETerms`, `~astropy.coordinates.Galactic`, and
`~astropy.coordinates.AltAz`. The string aliases are lowercase versions of the
class name.

If the frame is not supplied then you will see a special ``ICRS``
identifier. This indicates that the frame is unspecified and operations
that require comparing coordinates (even within that object) are not allowed.

**unit=UNIT**

The unit specifier can be one of the following:

- `~astropy.units.Unit` object, which is an angular unit that is equivalent to
  ``Unit('radian')``.
- Single string with a valid angular unit name.
- 2-tuple of `~astropy.units.Unit` objects or string unit names specifying the
  LON and LAT unit, respectively (e.g., ``('hourangle', 'degree')``).
- Single string with two unit names separated by a comma (e.g.,
  ``'hourangle,degree'``).

If only a single unit is provided then it applies to both LON and LAT.

**Other keyword arguments**

In lieu of positional arguments to specify the longitude and latitude, the
frame-specific names can be used as keyword arguments:

*ra*, *dec*: **LON**, **LAT** values, optional
    RA and Dec for frames where these are representation, including [FIXME]
    `~astropy.coordinates.ICRS`, `~astropy.coordinates.FK5`,
    `~astropy.coordinates.FK4`, and `~astropy.coordinates.FK4NoETerms`.

*l*, *b*:  **LON**, **LAT** values, optional
    Galactic ``l`` and ``b`` for the `~astropy.coordinates.Galactic` frame.

The following keywords can be specified for any frame:

*distance*: valid `~astropy.coordinates.Distance` initializer, optional
    Distance from reference from center to source

*obstime*: valid `~astropy.time.Time` initializer, optional
    Time of observation

*equinox*: valid `~astropy.time.Time` initializer, optional
    Coordinate frame equinox

If custom user-defined frames are included in the transform graph and they
have additional frame attributes, then those attributes can also be
set via corresponding keyword arguments in the |SkyCoord| initialization.

.. _astropy-coordinates-array-operations:

Array Operations
================

It is possible to store arrays of coordinates in a |SkyCoord| object, and
manipulations done in this way will be orders of magnitude faster than
looping over a list of individual |SkyCoord| objects::

  >>> ra = np.linspace(0, 36000, 1001) * u.deg
  >>> dec = np.linspace(-90, 90, 1001) * u.deg

  >>> sc_list = [SkyCoord(r, d, frame='icrs') for r, d in zip(ra, dec)]  # doctest: +SKIP
  >>> timeit sc_gal_list = [c.galactic for c in sc_list]  # doctest: +SKIP
  1 loops, best of 3: 20.4 s per loop

  >>> sc = SkyCoord(ra, dec, frame='icrs')
  >>> timeit sc_gal = sc.galactic  # doctest: +SKIP
  100 loops, best of 3: 21.8 ms per loop

In addition to vectorized transformations, you can do the usual array slicing,
dicing, and selection using the same methods and attributes that you use for
`~numpy.ndarray` instances::

  >>> north_mask = sc.dec > 0
  >>> sc_north = sc[north_mask]
  >>> len(sc_north)
  500
  >>> sc[2:4]  # doctest: +FLOAT_CMP
  <SkyCoord (ICRS): (ra, dec) in deg
      [( 72., -89.64), (108., -89.46)]>
  >>> sc[500]  # doctest: +FLOAT_CMP
  <SkyCoord (ICRS): (ra, dec) in deg
      (0., 0.)>
  >>> sc[0:-1:100].reshape(2, 5)  # doctest: +FLOAT_CMP
  <SkyCoord (ICRS): (ra, dec) in deg
      [[(0., -90.), (0., -72.), (0., -54.), (0., -36.), (0., -18.)],
       [(0.,   0.), (0.,  18.), (0.,  36.), (0.,  54.), (0.,  72.)]]>

Note that similarly to the `~numpy.ndarray` methods, all but ``flatten`` try to
use new views of the data, with the data copied only if that is impossible
(as discussed, for example, in the documentation for NumPy
:func:`~numpy.reshape`).


Attributes
===========

The |SkyCoord| object has a number of useful attributes which come in handy.
By digging through these we will learn a little bit about |SkyCoord| and how it
works.

To begin, one of the most important tools for
learning about attributes and methods of objects is "TAB-discovery". From
within IPython you can type an object name, the period, and then the <TAB> key
to see what is available. This can often be faster than reading the
documentation::

  >>> sc = SkyCoord(1, 2, frame='icrs', unit='deg', obstime='2013-01-02 14:25:36')
  >>> sc.<TAB>  # doctest: +SKIP
  sc.T                                   sc.match_to_catalog_3d
  sc.altaz                               sc.match_to_catalog_sky
  sc.barycentrictrueecliptic             sc.name
  sc.cartesian                           sc.ndim
  sc.cirs                                sc.obsgeoloc
  sc.copy                                sc.obsgeovel
  sc.data                                sc.obstime
  sc.dec                                 sc.obswl
  sc.default_representation              sc.position_angle
  sc.diagonal                            sc.precessedgeocentric
  sc.distance                            sc.pressure
  sc.equinox                             sc.ra
  sc.fk4                                 sc.ravel
  sc.fk4noeterms                         sc.realize_frame
  sc.fk5                                 sc.relative_humidity
  sc.flatten                             sc.represent_as
  sc.frame                               sc.representation_component_names
  sc.frame_attributes                    sc.representation_component_units
  sc.frame_specific_representation_info  sc.representation_info
  sc.from_name                           sc.reshape
  sc.from_pixel                          sc.roll
  sc.galactic                            sc.search_around_3d
  sc.galactocentric                      sc.search_around_sky
  sc.galcen_distance                     sc.separation
  sc.gcrs                                sc.separation_3d
  sc.geocentrictrueecliptic              sc.shape
  sc.get_constellation                   sc.size
  sc.get_frame_attr_names                sc.skyoffset_frame
  sc.guess_from_table                    sc.spherical
  sc.has_data                            sc.spherical_offsets_to
  sc.hcrs                                sc.squeeze
  sc.heliocentrictrueecliptic            sc.supergalactic
  sc.icrs                                sc.swapaxes
  sc.info                                sc.take
  sc.is_equivalent_frame                 sc.temperature
  sc.is_frame_attr_default               sc.to_pixel
  sc.is_transformable_to                 sc.to_string
  sc.isscalar                            sc.transform_to
  sc.itrs                                sc.transpose
  sc.location                            sc.z_sun

Here we see many attributes and methods. The most recognizable may be the
longitude and latitude attributes which are named ``ra`` and ``dec`` for the
``ICRS`` frame::

  >>> sc.ra  # doctest: +FLOAT_CMP
  <Longitude 1. deg>
  >>> sc.dec  # doctest: +FLOAT_CMP
  <Latitude 2. deg>

Next, notice that all of the built-in frame names ``icrs``, ``galactic``,
``fk5``, ``fk4``, and ``fk4noeterms`` are there. Through the magic of Python
properties, accessing these attributes calls the object
`~astropy.coordinates.SkyCoord.transform_to` method appropriately and returns a
new |SkyCoord| object in the requested frame::

  >>> sc_gal = sc.galactic
  >>> sc_gal  # doctest: +FLOAT_CMP
  <SkyCoord (Galactic): (l, b) in deg
      (99.63785528, -58.70969293)>

Other attributes you may recognize are ``distance``, ``equinox``,
``obstime``, and ``shape``.

Digging Deeper
--------------
*[Casual users can skip this section]*

After transforming to Galactic, the longitude and latitude values are now
labeled ``l`` and ``b``, following the normal convention for Galactic
coordinates. How does the object know what to call its values? The answer
lies in some less obvious attributes::

  >>> sc_gal.representation_component_names
  OrderedDict([('l', 'lon'), ('b', 'lat'), ('distance', 'distance')])

  >>> sc_gal.representation_component_units
  OrderedDict([('l', Unit("deg")), ('b', Unit("deg"))])

  >>> sc_gal.representation_type
  <class 'astropy.coordinates.representation.SphericalRepresentation'>

Together these tell the object that ``l`` and ``b`` are the longitude and
latitude, and that they should both be displayed in units of degrees as
a spherical-type coordinate (and not, for example, a Cartesian coordinate).
Furthermore, the frame's ``representation_component_names`` attribute defines
the coordinate keyword arguments that |SkyCoord| will accept.

Another important attribute is ``frame_attr_names``, which defines the
additional attributes that are required to fully define the frame::

  >>> sc_fk4 = SkyCoord(1, 2, frame='fk4', unit='deg')
  >>> sc_fk4.get_frame_attr_names()
  OrderedDict([('equinox', <Time object: scale='tt' format='byear_str' value=B1950.000>), ('obstime', None)])

The key values correspond to the defaults if no explicit value is provided by
the user. This example shows that the `~astropy.coordinates.FK4` frame has two
attributes, ``equinox`` and ``obstime``, that are required to fully define the
frame.

Some trickery is happening here because many of these attributes are
actually owned by the underlying coordinate ``frame`` object which does much of
the real work. This is the middle layer in the three-tiered system of objects:
representation (spherical, Cartesian, etc.), frame (a.k.a. low-level frame
class), and |SkyCoord| (a.k.a. high-level class; see
:ref:`astropy-coordinates-overview` and
:ref:`astropy-coordinates-definitions`)::

  >>> sc.frame  # doctest: +FLOAT_CMP
  <ICRS Coordinate: (ra, dec) in deg
      (1., 2.)>

  >>> sc.has_data is sc.frame.has_data
  True

  >>> sc.frame.<TAB>  # doctest: +SKIP
  sc.frame.T                                   sc.frame.ra
  sc.frame.cartesian                           sc.frame.ravel
  sc.frame.copy                                sc.frame.realize_frame
  sc.frame.data                                sc.frame.represent_as
  sc.frame.dec                                 sc.frame.representation
  sc.frame.default_representation              sc.frame.representation_component_names
  sc.frame.diagonal                            sc.frame.representation_component_units
  sc.frame.distance                            sc.frame.representation_info
  sc.frame.flatten                             sc.frame.reshape
  sc.frame.frame_attributes                    sc.frame.separation
  sc.frame.frame_specific_representation_info  sc.frame.separation_3d
  sc.frame.get_frame_attr_names                sc.frame.shape
  sc.frame.has_data                            sc.frame.size
  sc.frame.is_equivalent_frame                 sc.frame.spherical
  sc.frame.is_frame_attr_default               sc.frame.squeeze
  sc.frame.is_transformable_to                 sc.frame.swapaxes
  sc.frame.isscalar                            sc.frame.take
  sc.frame.name                                sc.frame.transform_to
  sc.frame.ndim                                sc.frame.transpose

  >>> sc.frame.name
  'icrs'

The |SkyCoord| object exposes the ``frame`` object attributes as its own. Though
it might seem a tad confusing at first, this is a good thing because it makes
|SkyCoord| objects and `~astropy.coordinates.BaseCoordinateFrame` objects
behave very similarly and most routines can accept either one as input without
much bother (duck typing!).

The lowest layer in the stack is the abstract
`~astropy.coordinates.UnitSphericalRepresentation` object:

  >>> sc_gal.frame.data  # doctest: +FLOAT_CMP
  <UnitSphericalRepresentation (lon, lat) in rad
      (1.73900863, -1.02467744)>

Transformations
===============

The topic of transformations is covered in detail in the section on
:ref:`astropy-coordinates-transforming`.

For completeness, here we will give some examples. Once you have defined
your coordinates and the reference frame, you can transform from that frame to
another frame. You can do this in a few different ways: if you only want the
default version of that frame, you can use attribute-style access (as mentioned
previously). For more control, you can use the
`~astropy.coordinates.SkyCoord.transform_to` method, which accepts a frame
name, frame class, frame instance, or |SkyCoord|::

  >>> from astropy.coordinates import FK5
  >>> sc = SkyCoord(1, 2, frame='icrs', unit='deg')
  >>> sc.galactic  # doctest: +FLOAT_CMP
  <SkyCoord (Galactic): (l, b) in deg
      (99.63785528, -58.70969293)>

  >>> sc.transform_to('fk5')  # Same as sc.fk5 and sc.transform_to(FK5)  # doctest: +FLOAT_CMP
  <SkyCoord (FK5: equinox=J2000.000): (ra, dec) in deg
          (1.00000656, 2.00000243)>

  >>> sc.transform_to(FK5(equinox='J1975'))  # Transform to FK5 with a different equinox  # doctest: +FLOAT_CMP
  <SkyCoord (FK5: equinox=J1975.000): (ra, dec) in deg
          (0.67967282, 1.86083014)>

Transforming to a |SkyCoord| instance is a convenient way of ensuring that two
coordinates are in the exact same reference frame::

  >>> sc2 = SkyCoord(3, 4, frame='fk4', unit='deg', obstime='J1978.123', equinox='B1960.0')
  >>> sc.transform_to(sc2)  # doctest: +FLOAT_CMP
  <SkyCoord (FK4: equinox=B1960.000, obstime=J1978.123): (ra, dec) in deg
      (0.48726331, 1.77731617)>

.. _astropy-skycoord-representations:

Representations
===============

So far we have been using a spherical coordinate representation in all of the
examples, and this is the default for the built-in frames. Frequently it is
convenient to initialize or work with a coordinate using a different
representation such as Cartesian or cylindrical. In this section, we discuss
how to initialize an object using a different representation and how to
change the representation of an object. For more information about
representation objects themselves, see :ref:`astropy-coordinates-representations`.

Initialization
--------------

Most of what you need to know can be inferred from the examples below and
by extrapolating the previous documentation for spherical representations.
Initialization requires setting the ``representation_type`` keyword and
supplying the corresponding components for that representation::

    >>> c = SkyCoord(x=1, y=2, z=3, unit='kpc', representation_type='cartesian')
    >>> c  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (x, y, z) in kpc
        (1., 2., 3.)>
    >>> c.x, c.y, c.z  # doctest: +FLOAT_CMP
    (<Quantity 1. kpc>, <Quantity 2. kpc>, <Quantity 3. kpc>)

Other variations include::

    >>> SkyCoord(1, 2*u.deg, 3, representation_type='cylindrical')  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (rho, phi, z) in (, deg, )
        (1., 2., 3.)>

    >>> SkyCoord(rho=1*u.km, phi=2*u.deg, z=3*u.m, representation_type='cylindrical')  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (rho, phi, z) in (km, deg, m)
        (1., 2., 3.)>

    >>> SkyCoord(rho=1, phi=2, z=3, unit=(u.km, u.deg, u.m), representation_type='cylindrical')  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (rho, phi, z) in (km, deg, m)
        (1., 2., 3.)>

    >>> SkyCoord(1, 2, 3, unit=(None, u.deg, None), representation_type='cylindrical')  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (rho, phi, z) in (, deg, )
        (1., 2., 3.)>

In general terms, the allowed syntax is as follows::

  SkyCoord(COORD, [FRAME | frame=FRAME], [unit=UNIT], [representation_type=REPRESENTATION],
           keyword_args ...)
  SkyCoord(COMP1, COMP2, [COMP3], [FRAME | frame=FRAME], [unit=UNIT],
           [representation_type=REPRESENTATION], keyword_args ...)
  SkyCoord([FRAME | frame=FRAME], <comp1_name>=COMP1, <comp2_name>=COMP2,
           <comp3_name>=COMP3, [representation_type=REPRESENTATION], [unit=UNIT],
           keyword_args ...)

In this case, the ``keyword_args`` now includes the element
``representation_type=REPRESENTATION``. In the above description, elements in
all capital letters (e.g., ``FRAME``) describe a user input of that element
type. Elements in square brackets are optional.

**COMP1**, **COMP2**, **COMP3**

Component values can be specified as separate positional arguments or as
keyword arguments. In this formalism the exact type of allowed input depends
on the details of the representation. In general, the following input forms
are supported:

- Single value:

  - Component class object
  - Plain numeric value with ``unit`` keyword specifying the unit

- List or component class array, or NumPy array of values

Each representation component has a specified class (the "component class")
which is used to convert generic input data into a predefined object
class with a certain unit. These component classes are expected to be
subclasses of the `~astropy.units.Quantity` class.

**COORD**

This input form uses a single object to supply coordinate data. The coordinate
can specify one or more coordinate positions as follows:

- List of ``(COMP1, .., COMP<M>)`` tuples, where each component is a scalar (not
  array) and there are ``M`` components in the representation. Typically
  there are three components, but some
  (e.g., `~astropy.coordinates.UnitSphericalRepresentation`)
  can have fewer.
- ``N x M`` NumPy or |Quantity| array of values, where ``N`` is the number
  of coordinates and ``M`` is the number of components.

**REPRESENTATION**

The representation can be supplied either as a
`~astropy.coordinates.representation.BaseRepresentation` class (e.g.,
`~astropy.coordinates.CartesianRepresentation`) or as a string name
that is simply the class name in lowercase without the
``'representation'`` suffix (e.g., ``'cartesian'``).

The rest of the inputs for creating a |SkyCoord| object in the general case are
the same as for spherical.

Details
-------

The available set of representations is dynamic and may change depending on what
representation classes have been defined. The built-in representations are:

=====================  =======================================================
  Name                   Class
=====================  =======================================================
``spherical``          `~astropy.coordinates.SphericalRepresentation`
``unitspherical``      `~astropy.coordinates.UnitSphericalRepresentation`
``physicsspherical``   `~astropy.coordinates.PhysicsSphericalRepresentation`
``cartesian``          `~astropy.coordinates.CartesianRepresentation`
``cylindrical``        `~astropy.coordinates.CylindricalRepresentation`
=====================  =======================================================

Each frame knows about all of the available representations, but different
frames may use different names for the same components. A common example
is that the `~astropy.coordinates.Galactic` frame uses ``l`` and ``b``
instead of ``ra`` and ``dec`` for the ``lon`` and ``lat`` components of
the `~astropy.coordinates.SphericalRepresentation`.

For a particular frame, in order to see the full list of representations
and how it names all of the components, first make an instance of that frame
without any data, and then print the ``representation_info`` property::

    >>> ICRS().representation_info  # doctest: +SKIP
    {astropy.coordinates.representation.CartesianRepresentation:
      {'names': ('x', 'y', 'z'),
       'units': (None, None, None)},
     astropy.coordinates.representation.SphericalRepresentation:
      {'names': ('ra', 'dec', 'distance'),
       'units': (Unit("deg"), Unit("deg"), None)},
     astropy.coordinates.representation.UnitSphericalRepresentation:
      {'names': ('ra', 'dec'),
       'units': (Unit("deg"), Unit("deg"))},
     astropy.coordinates.representation.PhysicsSphericalRepresentation:
      {'names': ('phi', 'theta', 'r'),
       'units': (Unit("deg"), Unit("deg"), None)},
     astropy.coordinates.representation.CylindricalRepresentation:
      {'names': ('rho', 'phi', 'z'),
       'units': (None, Unit("deg"), None)}
    }

This is a bit messy but it shows that for each representation there is a
``dict`` with two keys:

- ``names``: defines how each component is named in that frame.
- ``units``: defines the units of each component when output, where ``None``
  means to not force a particular unit.

For a particular coordinate instance you can use the ``representation_type``
attribute in conjunction with the ``representation_component_names`` attribute
to figure out what keywords are accepted by a particular class object. The
former will be the representation class the system is expressed in (e.g.,
spherical for equatorial frames), and the latter will be a dictionary mapping
names for that frame to the component name on the representation class::

    >>> import astropy.units as u
    >>> icrs = ICRS(1*u.deg, 2*u.deg)
    >>> icrs.representation_type
    <class 'astropy.coordinates.representation.SphericalRepresentation'>
    >>> icrs.representation_component_names
    OrderedDict([('ra', 'lon'), ('dec', 'lat'), ('distance', 'distance')])

Changing Representation
--------------------------

The representation of the coordinate object can be changed, as shown
below. This actually does *nothing* to the object internal data which
stores the coordinate values, but it changes the external view of that
data in two ways:

- The object prints itself in accord with the new representation.
- The available attributes change to match those of the new representation
  (e.g., from ``ra, dec, distance`` to ``x, y, z``).

Setting the ``representation_type`` thus changes a *property* of the
object (how it appears) without changing the intrinsic object itself
which represents a point in 3D space::

    >>> c = SkyCoord(x=1, y=2, z=3, unit='kpc', representation_type='cartesian')
    >>> c  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (x, y, z) in kpc
        (1., 2., 3.)>

    >>> c.representation_type = 'cylindrical'
    >>> c  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (rho, phi, z) in (kpc, deg, kpc)
        (2.23606798, 63.43494882, 3.)>
    >>> c.phi.to(u.deg)  # doctest: +FLOAT_CMP
    <Angle 63.43494882 deg>
    >>> c.x
    Traceback (most recent call last):
    ...
    AttributeError: 'SkyCoord' object has no attribute 'x'

    >>> c.representation_type = 'spherical'
    >>> c  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec, distance) in (deg, deg, kpc)
        (63.43494882, 53.3007748, 3.74165739)>

    >>> c.representation_type = 'unitspherical'
    >>> c  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec) in deg
        (63.43494882, 53.3007748)>

You can also use any representation class to set the representation::

    >>> from astropy.coordinates import CartesianRepresentation
    >>> c.representation_type = CartesianRepresentation

Note that if all you want is a particular representation without changing the
state of the |SkyCoord| object, you should instead use the
``astropy.coordinates.SkyCoord.represent_as()`` method::

    >>> c.representation_type = 'spherical'
    >>> cart = c.represent_as(CartesianRepresentation)
    >>> cart  # doctest: +FLOAT_CMP
    <CartesianRepresentation (x, y, z) in kpc
        (1., 2., 3.)>
    >>> c.representation_type
    <class 'astropy.coordinates.representation.SphericalRepresentation'>


Example 1: Plotting random data in Aitoff projection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is an example of how to make a plot in the Aitoff projection using data
in a |SkyCoord| object. Here, a randomly generated data set will be used.

First we need to import the required packages. We use
`matplotlib <https://matplotlib.org/>`_ here for
plotting and `numpy <https://numpy.org/>`_  to get the value of pi and to
generate our random data.

    >>> from astropy import units as u
    >>> from astropy.coordinates import SkyCoord
    >>> import numpy as np

We now generate random data for visualization. For RA this is done in the range
of 0 and 360 degrees (``ra_random``), for DEC between -90 and +90 degrees
(``dec_random``). Finally, we multiply these values by degrees to get a
`~astropy.units.Quantity` with units of degrees.

    >>> ra_random = np.random.rand(100)*360.0 * u.degree
    >>> dec_random = (np.random.rand(100)*180.0-90.0) * u.degree

As the next step, those coordinates are transformed into an
`astropy.coordinates` |SkyCoord| object.

    >>> c = SkyCoord(ra=ra_random, dec=dec_random, frame='icrs')

Because matplotlib needs the coordinates in radians and between :math:`-\pi`
and :math:`\pi`, not 0 and :math:`2\pi`, we have to convert them.
For this purpose the `astropy.coordinates.Angle` object provides a special
method, which we use here to wrap at 180:

    >>> ra_rad = c.ra.wrap_at(180 * u.deg).radian
    >>> dec_rad = c.dec.radian

As a last step, we set up the plotting environment with matplotlib using the
Aitoff projection with a specific title, a grid, filled circles as markers with
a marker size of 2, and an alpha value of 0.3. We use a figure with an x-y ratio
that is well suited for such a projection and we move the title upwards from
its usual position to avoid overlap with the axis labels.

.. doctest-skip::

    >>> import matplotlib.pyplot as plt
    >>> plt.figure(figsize=(8,4.2))
    >>> plt.subplot(111, projection="aitoff")
    >>> plt.title("Aitoff projection of our random data")
    >>> plt.grid(True)
    >>> plt.plot(ra_rad, dec_rad, 'o', markersize=2, alpha=0.3)
    >>> plt.subplots_adjust(top=0.95,bottom=0.0)
    >>> plt.show()


.. plot::

    # This is an example how to make a plot in the Aitoff projection using data
    # in a SkyCoord object. Here a randomly generated data set will be used. The
    # final script can be found below.

    # First we need to import the required packages. We use
    # `matplotlib <https://matplotlib.org/>`_ here for
    # plotting and `numpy <https://numpy.org/>`_  to get the value of pi and to
    # generate our random data.
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    import matplotlib.pyplot as plt
    import numpy as np

    # We now generate random data for visualization. For RA this is done in the range
    # of 0 and 360 degrees (``ra_random``), for DEC between -90 and +90 degrees
    # (``dec_random``). Finally, we multiply these values by degrees to get a
    # `~astropy.units.Quantity` with units of degrees.
    ra_random = np.random.rand(100)*360.0 * u.degree
    dec_random = (np.random.rand(100)*180.0-90.0) * u.degree

    # As the next step, those coordinates are transformed into an astropy.coordinates
    # astropy.coordinates.SkyCoord object.
    c = SkyCoord(ra=ra_random, dec=dec_random, frame='icrs')

    # Because matplotlib needs the coordinates in radians and between :math:`-\pi`
    # and :math:`\pi`, not 0 and :math:`2\pi`, we have to convert them.
    # For this purpose the `astropy.coordinates.Angle` object provides a special method,
    # which we use here to wrap at 180:
    ra_rad = c.ra.wrap_at(180 * u.deg).radian
    dec_rad = c.dec.radian

    # As a last step we set up the plotting environment with matplotlib using the
    # Aitoff projection with a specific title, a grid, filled circles as markers with
    # a marker size of 2, and an alpha value of 0.3.
    plt.figure(figsize=(8,4.2))
    plt.subplot(111, projection="aitoff")
    plt.title("Aitoff projection of our random data", y=1.08)
    plt.grid(True)
    plt.plot(ra_rad, dec_rad, 'o', markersize=2, alpha=0.3)
    plt.subplots_adjust(top=0.95, bottom=0.0)
    plt.show()



Example 2: Plotting star positions in bulge and disk
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is a more realistic example of how to make a plot in the Aitoff projection
using data in a |SkyCoord| object. Here, a randomly generated data set
(multivariate normal distribution) for both stars in the bulge and in the disk
of a galaxy will be used. Both types will be plotted with different number
counts.

As in the last example, we first import the required packages.

    >>> from astropy import units as u
    >>> from astropy.coordinates import SkyCoord
    >>> import numpy as np

We now generate random data for visualization using
`numpy.random.Generator.multivariate_normal`.

    >>> disk = np.random.multivariate_normal(mean=[0,0,0], cov=np.diag([1,1,0.5]), size=5000)
    >>> bulge = np.random.multivariate_normal(mean=[0,0,0], cov=np.diag([1,1,1]), size=500)
    >>> galaxy = np.concatenate([disk, bulge])

As the next step, those coordinates are transformed into an
`astropy.coordinates` |SkyCoord| object.

    >>> c_gal = SkyCoord(galaxy, representation_type='cartesian', frame='galactic')
    >>> c_gal_icrs = c_gal.icrs

Again, as in the last example, we need to convert the coordinates in radians
and make sure they are between :math:`-\pi` and :math:`\pi`:

    >>> ra_rad = c_gal_icrs.ra.wrap_at(180 * u.deg).radian
    >>> dec_rad = c_gal_icrs.dec.radian

We use the same plotting setup as in the last example:

.. doctest-skip::

    >>> import matplotlib.pyplot as plt
    >>> plt.figure(figsize=(8,4.2))
    >>> plt.subplot(111, projection="aitoff")
    >>> plt.title("Aitoff projection of our random data")
    >>> plt.grid(True)
    >>> plt.plot(ra_rad, dec_rad, 'o', markersize=2, alpha=0.3)
    >>> plt.subplots_adjust(top=0.95,bottom=0.0)
    >>> plt.show()


.. plot::

    # This is more realistic example how to make a plot in the Aitoff projection
    # using data in a SkyCoord object.
    # Here a randomly generated data set (multivariate normal distribution)
    # for both stars in the bulge and in the disk of a galaxy
    # will be used. Both types will be plotted with different number counts. The
    # final script can be found below.

    # As in the last example, we first import the required packages.
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    import matplotlib.pyplot as plt
    import numpy as np

    # We now generate random data for visualization with
    # np.random.multivariate_normal.
    disk = np.random.multivariate_normal(mean=[0,0,0], cov=np.diag([1,1,0.5]), size=5000)
    bulge = np.random.multivariate_normal(mean=[0,0,0], cov=np.diag([1,1,1]), size=500)
    galaxy = np.concatenate([disk, bulge])

    # As the next step, those coordinates are transformed into an astropy.coordinates
    # astropy.coordinates.SkyCoord object.
    c_gal = SkyCoord(galaxy, representation_type='cartesian', frame='galactic')
    c_gal_icrs = c_gal.icrs

    # Again, as in the last example, we need to convert the coordinates in radians
    # and make sure they are between :math:`-\pi` and :math:`\pi`:
    ra_rad = c_gal_icrs.ra.wrap_at(180 * u.deg).radian
    dec_rad = c_gal_icrs.dec.radian

    # We use the same plotting setup as in the last example:
    plt.figure(figsize=(8,4.2))
    plt.subplot(111, projection="aitoff")
    plt.title("Aitoff projection of our random data", y=1.08)
    plt.grid(True)
    plt.plot(ra_rad, dec_rad, 'o', markersize=2, alpha=0.3)
    plt.subplots_adjust(top=0.95,bottom=0.0)
    plt.show()


Convenience Methods
====================

A number of convenience methods are available, and you are encouraged to read
the available docstrings below:

- `~astropy.coordinates.SkyCoord.match_to_catalog_sky`,
- `~astropy.coordinates.SkyCoord.match_to_catalog_3d`,
- `~astropy.coordinates.SkyCoord.position_angle`,
- `~astropy.coordinates.SkyCoord.separation`,
- `~astropy.coordinates.SkyCoord.separation_3d`
- `~astropy.coordinates.SkyCoord.apply_space_motion`

Additional information and examples can be found in the section on
:ref:`astropy-coordinates-separations-matching` and
:ref:`astropy-coordinates-apply-space-motion`.
