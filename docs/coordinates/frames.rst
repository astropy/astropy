.. include:: references.txt

.. We call EarthLocation.of_site here first to force the downloading
.. of sites.json so that future doctest output isn't clutted with
.. "Downloading ... [done]". This can be removed once we have a better
.. way of ignoring output lines based on pattern-matching, e.g.:
.. https://github.com/astropy/pytest-doctestplus/issues/11

.. testsetup::
    >>> from astropy.coordinates import EarthLocation
    >>> EarthLocation.of_site('greenwich') # doctest: +IGNORE_OUTPUT

Using and Designing Coordinate Frames
*************************************

In `astropy.coordinates`, as outlined in the
:ref:`astropy-coordinates-overview`, subclasses of |baseframe| ("frame
classes") define particular coordinate frames. They can (but do not
*have* to) contain representation objects storing the actual coordinate
data. The actual coordinate transformations are defined as functions
that transform representations between frame classes.  This approach
serves to separate high-level user functionality (see :doc:`skycoord`)
and details of how the coordinates are actually stored (see
:doc:`representations`) from the definition of frames and how they are
transformed.

Using Frame Objects
===================

Frames without Data
-------------------

Frame objects have two distinct (but related) uses.  The first is
storing the information needed to uniquely define a frame (e.g.,
equinox, observation time). This information is stored on the frame
objects as (read-only) Python attributes, which are set when the object
is first created::

    >>> from astropy.coordinates import ICRS, FK5
    >>> FK5(equinox='J1975')
    <FK5 Frame (equinox=J1975.000)>
    >>> ICRS()  # has no attributes
    <ICRS Frame>
    >>> FK5()  # uses default equinox
    <FK5 Frame (equinox=J2000.000)>

The specific names of attributes available for a particular frame (and
their default values)  are available as the class method
``get_frame_attr_names``::

    >>> FK5.get_frame_attr_names()
    OrderedDict([('equinox', <Time object: scale='utc' format='jyear_str' value=J2000.000>)])

You can access any of the attributes on a frame by using standard Python
attribute access.  Note that for cases like ``equinox``, which are time
inputs, if you pass in any unambiguous time string, it will be converted
into an `~astropy.time.Time` object with UTC scale (see
:ref:`astropy-time-inferring-input`)::

    >>> f = FK5(equinox='J1975')
    >>> f.equinox
    <Time object: scale='utc' format='jyear_str' value=J1975.000>
    >>> f = FK5(equinox='2011-05-15T12:13:14')
    >>> f.equinox
    <Time object: scale='utc' format='isot' value=2011-05-15T12:13:14.000>


Frames with Data
----------------

The second use for frame objects is to store actual realized coordinate
data for frames like those described above. In this use, it is similar
to the |skycoord| class, and in fact, the |skycoord| class internally
uses the frame classes as its implementation.  However, the frame
classes have fewer "convenience" features, thereby keeping the
implementation of frame classes simple.  As such, they are created
similarly to |skycoord| object.  The simplest way is to use
with keywords appropriate for the frame (e.g. ``ra`` and ``dec`` for
equatorial systems)::

    >>> from astropy import units as u
    >>> ICRS(ra=1.1*u.deg, dec=2.2*u.deg)  # doctest: +FLOAT_CMP
    <ICRS Coordinate: (ra, dec) in deg
        (1.1, 2.2)>
    >>> FK5(ra=1.1*u.deg, dec=2.2*u.deg, equinox='J1975')  # doctest: +FLOAT_CMP
    <FK5 Coordinate (equinox=J1975.000): (ra, dec) in deg
        (1.1, 2.2)>

These same attributes can be used to access the data in the frames, as
|Angle| objects (or |Angle| subclasses)::

    >>> coo = ICRS(ra=1.1*u.deg, dec=2.2*u.deg)
    >>> coo.ra  # doctest: +FLOAT_CMP
    <Longitude 1.1 deg>
    >>> coo.ra.value  # doctest: +FLOAT_CMP
    1.1
    >>> coo.ra.to(u.hourangle)  # doctest: +FLOAT_CMP
    <Longitude 0.07333333 hourangle>

You can use the ``representation_type`` attribute in conjunction
with the ``representation_component_names`` attribute to figure out what
keywords are accepted by a particular class object.  The former will be the
representation class the system is expressed in (e.g.,
spherical for equatorial frames), and the latter will be a dictionary
mapping names for that frame to the attribute name on the representation
class::

    >>> import astropy.units as u
    >>> icrs = ICRS(1*u.deg, 2*u.deg)
    >>> icrs.representation_type
    <class 'astropy.coordinates.representation.SphericalRepresentation'>
    >>> icrs.representation_component_names
    OrderedDict([('ra', 'lon'), ('dec', 'lat'), ('distance', 'distance')])

One can get the data in a different representation if needed::

    >>> icrs.represent_as('cartesian')  # doctest: +FLOAT_CMP
    <CartesianRepresentation (x, y, z) [dimensionless]
         (0.99923861, 0.01744177, 0.0348995)>

.. note::

    In previous versions of Astropy, both the frame attribute and the argument
    to frame classes that are now named ``representation_type`` used to be
    simply ``representation``. The name of this attribute/argument is confusing
    as it points to the representation *class*, not the object containing the
    underlying frame data (this is accessed via the frame attribute ``.data``).
    To clarify, we have renamed ``representation`` to ``representation_type``.
    In this version 3.0, we have only changed the references to this attribute
    in the documentation. In the next major version, we will issue a deprecation
    warning. In two major versions, we will remove the ``.representation``
    attribute and ``representation=`` argument.

The representation of the coordinate object can also be changed directly, as
shown below.  This actually does *nothing* to the object internal data which
stores the coordinate values, but it changes the external view of that data in
two ways: (1) the object prints itself in accord with the new representation,
and (2) the available attributes change to match those of the new
representation (e.g. from ``ra, dec, distance`` to ``x, y, z``).  Setting the
``representation_type`` thus changes a *property* of the object (how it appears)
without changing the intrinsic object itself which represents a point in 3d
space.::

    >>> from astropy.coordinates import CartesianRepresentation
    >>> icrs.representation_type = CartesianRepresentation
    >>> icrs  # doctest: +FLOAT_CMP
    <ICRS Coordinate: (x, y, z) [dimensionless]
        (0.99923861, 0.01744177, 0.0348995)>
    >>> icrs.x  # doctest: +FLOAT_CMP
    <Quantity 0.99923861>

The representation can also be set at the time of creating a coordinate
and affects the set of keywords used to supply the coordinate data.  For
example to create a coordinate with cartesian data do::

    >>> ICRS(x=1*u.kpc, y=2*u.kpc, z=3*u.kpc, representation_type='cartesian')  #  doctest: +FLOAT_CMP
    <ICRS Coordinate: (x, y, z) in kpc
        (1., 2., 3.)>

For more information about the use of representations in coordinates see the
:ref:`astropy-skycoord-representations` section, and for details about the
representations themselves see :ref:`astropy-coordinates-representations`.

There are two other ways to create frame classes with coordinates.  A
representation class can be passed in directly at creation, along with
any  frame attributes required::

    >>> from astropy.coordinates import SphericalRepresentation
    >>> rep = SphericalRepresentation(lon=1.1*u.deg, lat=2.2*u.deg, distance=3.3*u.kpc)
    >>> FK5(rep, equinox='J1975')  # doctest: +FLOAT_CMP
    <FK5 Coordinate (equinox=J1975.000): (ra, dec, distance) in (deg, deg, kpc)
        (1.1, 2.2, 3.3)>

A final way is to create a frame object from an already existing frame
(either one with or without data), using the ``realize_frame`` method. This
will yield a frame with the same attributes, but new data::

    >>> f1 = FK5(equinox='J1975')
    >>> f1
    <FK5 Frame (equinox=J1975.000)>
    >>> rep = SphericalRepresentation(lon=1.1*u.deg, lat=2.2*u.deg, distance=3.3*u.kpc)
    >>> f1.realize_frame(rep)  # doctest: +FLOAT_CMP
    <FK5 Coordinate (equinox=J1975.000): (ra, dec, distance) in (deg, deg, kpc)
        (1.1, 2.2, 3.3)>

You can check if a frame object has data using the ``has_data`` attribute, and
if it is preset, it can be accessed from the ``data`` attribute::

    >>> ICRS().has_data
    False
    >>> cooi = ICRS(ra=1.1*u.deg, dec=2.2*u.deg)
    >>> cooi.has_data
    True
    >>> cooi.data  # doctest: +FLOAT_CMP
    <UnitSphericalRepresentation (lon, lat) in deg
        (1.1, 2.2)>

All of the above methods can also accept array data (in the form of
class:`~astropy.units.Quantity`, or other Python sequences) to create arrays of
coordinates::

    >>> ICRS(ra=[1.5, 2.5]*u.deg, dec=[3.5, 4.5]*u.deg)  # doctest: +FLOAT_CMP
    <ICRS Coordinate: (ra, dec) in deg
        [(1.5, 3.5), (2.5, 4.5)]>

If you pass in mixed arrays and scalars, the arrays will be broadcast
over the scalars appropriately::

    >>> ICRS(ra=[1.5, 2.5]*u.deg, dec=[3.5, 4.5]*u.deg, distance=5*u.kpc)  # doctest: +FLOAT_CMP
    <ICRS Coordinate: (ra, dec, distance) in (deg, deg, kpc)
        [(1.5, 3.5, 5.), (2.5, 4.5, 5.)]>

Similar broadcasting happens if you transform to another frame.  E.g.::

    >>> import numpy as np
    >>> from astropy.coordinates import EarthLocation, AltAz
    >>> coo = ICRS(ra=180.*u.deg, dec=51.477811*u.deg)
    >>> lf = AltAz(location=EarthLocation.of_site('greenwich'),
    ...            obstime=['2012-03-21T00:00:00', '2012-06-21T00:00:00'])
    >>> lcoo = coo.transform_to(lf)  # this can load finals2000A.all # doctest: +IGNORE_OUTPUT
    >>> lcoo  # doctest: +FLOAT_CMP
    <AltAz Coordinate (obstime=['2012-03-21T00:00:00.000' '2012-06-21T00:00:00.000'], location=(3980608.9024681724, -102.47522910648239, 4966861.273100675) m, pressure=0.0 hPa, temperature=0.0 deg_C, relative_humidity=0.0, obswl=1.0 micron): (az, alt) in deg
        [( 94.71264944, 89.21424252), (307.69488825, 37.98077771)]>

Above, the shapes -- ``()`` for ``coo`` and ``(2,)`` for ``lf`` -- were
broadcast against each other.  If you wished to determine the positions for a
set of coordinates, you'd need to make sure that the shapes allowed this::

    >>> coo2 = ICRS(ra=[180., 225., 270.]*u.deg, dec=[51.5, 0., 51.5]*u.deg)
    >>> coo2.transform_to(lf)
    Traceback (most recent call last):
    ...
    ValueError: operands could not be broadcast together...
    >>> coo2.shape
    (3,)
    >>> lf.shape
    (2,)
    >>> lf2 = lf[:, np.newaxis]
    >>> lf2.shape
    (2, 1)
    >>> coo2.transform_to(lf2)  # doctest: +FLOAT_CMP
    <AltAz Coordinate (obstime=[['2012-03-21T00:00:00.000' '2012-03-21T00:00:00.000'
      '2012-03-21T00:00:00.000']
     ['2012-06-21T00:00:00.000' '2012-06-21T00:00:00.000'
      '2012-06-21T00:00:00.000']], location=(3980608.9024681724, -102.47522910648239, 4966861.273100675) m, pressure=0.0 hPa, temperature=0.0 deg_C, relative_humidity=0.0, obswl=1.0 micron): (az, alt) in deg
        [[( 93.09845202, 89.21613119), (126.85789652, 25.46600543),
          ( 51.37993229, 37.18532521)],
         [(307.71713699, 37.99437658), (231.37407871, 26.36768329),
          ( 85.42187335, 89.69297997)]]>

.. Note::
   One sees that frames without data have a ``shape`` that is determined by
   their frame attributes.  For frames with data the ``shape`` always is that
   of the data; any non-scalar attributes are broadcast to have matching shape
   (as can be seen for ``obstime`` in the last line above).

Transforming between Frames
===========================

To transform a frame object with data into another frame, use the
``transform_to`` method of an object, and provide it the frame you wish to
transform to.  This frame can either be a frame *class*, in which case
the default attributes will be used, or a frame object (with or without
data)::

    >>> cooi = ICRS(1.5*u.deg, 2.5*u.deg)
    >>> cooi.transform_to(FK5)  # doctest: +FLOAT_CMP
    <FK5 Coordinate (equinox=J2000.000): (ra, dec) in deg
        (1.50000661, 2.50000238)>
    >>> cooi.transform_to(FK5(equinox='J1975'))  # doctest: +FLOAT_CMP
    <FK5 Coordinate (equinox=J1975.000): (ra, dec) in deg
        (1.17960348, 2.36085321)>

The :ref:`astropy-coordinates-api` includes a list of all of the frames built
into `astropy.coordinates`, as well as the defined transformations between
them. Any transformation that has a valid path, even if it passes through
other frames, can be transformed to. To programmatically check for or
manipulate transformations, see the `~astropy.coordinates.TransformGraph`
documentation.


.. _astropy-coordinates-design:

Defining a New Frame
====================

Users can add new coordinate frames by creating new classes that are subclasses
of `~astropy.coordinates.BaseCoordinateFrame`.  Detailed instructions for
subclassing are in the docstrings for that class.  The key aspects are to
define the class attributes ``default_representation`` and
``frame_specific_representation_info`` along with frame attributes as
`~astropy.coordinates.Attribute` class instances (or subclasses like
`~astropy.coordinates.TimeAttribute`).  If these are
defined, there is often no need to define an ``__init__`` function, as the
initializer in `~astropy.coordinates.BaseCoordinateFrame` will probably behave
the way you want.  As an example::

  >>> from astropy.coordinates import BaseCoordinateFrame, Attribute, TimeAttribute, RepresentationMapping
  >>> import astropy.coordinates.representation as r
  >>> class MyFrame(BaseCoordinateFrame):
  ...     # Specify how coordinate values are represented when outputted
  ...      default_representation = r.SphericalRepresentation
  ...
  ...      # Specify overrides to the default names and units for all available
  ...      # representations (subclasses of BaseRepresentation).
  ...      frame_specific_representation_info = {
  ...          r.SphericalRepresentation: [RepresentationMapping(reprname='lon', framename='R', defaultunit=u.rad),
  ...                                      RepresentationMapping(reprname='lat', framename='D', defaultunit=u.rad),
  ...                                      RepresentationMapping(reprname='distance', framename='DIST', defaultunit=None)],
  ...          r.UnitSphericalRepresentation: [RepresentationMapping(reprname='lon', framename='R', defaultunit=u.rad),
  ...                                          RepresentationMapping(reprname='lat', framename='D', defaultunit=u.rad)],
  ...          r.CartesianRepresentation: [RepresentationMapping(reprname='x', framename='X'),
  ...                                      RepresentationMapping(reprname='y', framename='Y'),
  ...                                      RepresentationMapping(reprname='z', framename='Z')]
  ...      }
  ...
  ...      # Specify frame attributes required to fully specify the frame
  ...      location = Attribute(default=None)
  ...      equinox = TimeAttribute(default='B1950')
  ...      obstime = TimeAttribute(default=None, secondary_attribute='equinox')

  >>> c = MyFrame(R=10*u.deg, D=20*u.deg)
  >>> c  # doctest: +FLOAT_CMP
  <MyFrame Coordinate (location=None, equinox=B1950.000, obstime=B1950.000): (R, D) in rad
      (0.17453293, 0.34906585)>
  >>> c.equinox
  <Time object: scale='utc' format='byear_str' value=B1950.000>

If you also want to support velocity data in your coordinate frame, see the
velocities documentation at
:ref:`astropy-coordinate-custom-frame-with-velocities`.

You can also define arbitrary methods for any added functionality you
want your frame to have that's unique to that frame.  These methods will
be available in any |skycoord| that is created using your user-defined
frame.

For examples of defining frame classes, the first place to look is
probably the source code for the frames that are included in astropy
(available at ``astropy.coordinates.builtin_frames``).  These are not
"magic" in any way, and use all the same API and features available to
user-created frames.

.. topic:: Examples:

    See also :ref:`sphx_glr_generated_examples_coordinates_plot_sgr-coordinate-frame.py`
    for a more annotated example of defining a new coordinate frame.


Defining Transformations
========================

A frame may not be too useful without a way to transform coordinates
defined in it to or from other frames. Fortunately,
`astropy.coordinates` provides a framework to do just that.  The key
concept for these transformations is the frame transform graph,
available as ``astropy.coordinates.frame_transform_graph``, an instance of
the `~astropy.coordinates.TransformGraph` class.  This graph (in the
"graph theory" sense, not "plot"), stores all the transformations
between all of the builtin frames, as well as tools for finding shortest
paths through this graph to transform from any frame to any other.  All
of the power of this graph is available to user-created frames as well, meaning
that once you define even one transform from your frame to some frame in
the graph,  coordinates defined in your frame can be transformed to
*any* other frame in the graph.

The transforms themselves are represented as
`~astropy.coordinates.CoordinateTransform` objects or their subclasses. The useful
subclasses/types of transformations are:

* `~astropy.coordinates.FunctionTransform`

    A transform that is defined as a function that takes a frame object
    of one frame class and returns an object of another class.

* `~astropy.coordinates.AffineTransform`

    A transformation that includes a linear matrix operation and a translation
    (vector offset). These transformations are defined by a 3x3 matrix and a
    3-vector for the offset (supplied as a Cartesian representation). The
    transformation is applied to the Cartesian representation of one frame and
    transforms into the Cartesian representation of the target frame.

* `~astropy.coordinates.StaticMatrixTransform`
* `~astropy.coordinates.DynamicMatrixTransform`

    The matrix transforms are `~astropy.coordinates.AffineTransform`'s without
    a translation, i.e. a rotation. The static version is for
    the case where the matrix is independent of the frame attributes
    (e.g., the ICRS->FK5 transformation, because ICRS has no frame
    attributes).  The dynamic case is for transformations where the
    transformation matrix depends on the frame attributes of either the
    to or from frame.


Generally, it is not necessary to use these classes directly.  Instead,
use methods on ``frame_transform_graph`` that can be used as function
decorators.  Then just define functions that either do the actual
transformation (for FunctionTransform), or that compute the necessary
transformation matrices to transform. Then decorate the functions to
register these transformations with the frame transform graph::

    from astropy.coordinates import frame_transform_graph

    @frame_transform_graph.transform(DynamicMatrixTransform, ICRS, FK5)
    def icrs_to_fk5(icrscoord, fk5frame):
        ...

    @frame_transform_graph.transform(DynamicMatrixTransform, FK5, ICRS)
    def fk5_to_icrs(fk5coord, icrsframe):
        ...

If the transformation to your coordinate frame of interest is not
representable by a matrix operation, you can also specify a function to do
the actual transformation, and pass the
`~astropy.coordinates.FunctionTransform` class to the transform graph
decorator instead::

    @frame_transform_graph.transform(FunctionTransform, FK4NoETerms, FK4)
    def fk4_no_e_to_fk4(fk4noecoord, fk4frame):
        ...

Furthermore, the ``frame_transform_graph`` does some caching and
optimization to speed up transformations after the first attempt to go
from one frame to another, and shortcuts steps where relevant (for
example, combining multiple static matrix transforms into a single
matrix).  Hence, in general, it is better to define whatever are the
most natural transformations for a user-defined frame, rather than
worrying about optimizing or caching a transformation to speed up the
process.

For a demonstration of how to define transformation functions that also work for
transforming velocity components, see
:ref:`astropy-coordinate-transform-with-velocities`.
