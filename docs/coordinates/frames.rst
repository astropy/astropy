.. include:: references.txt

Using and Designing Coordinate Frames
-------------------------------------

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
^^^^^^^^^^^^^^^^^^^

Frame objects have two distinct (but related) uses.  The first is
storing the information needed to uniquely define a frame (e.g.,
equinox, observation time). This information is stored on the frame
objects as (read-only) Python attributes, which are set with the object
is first created::

    >>> from astropy.coordinates import ICRS, FK5
    >>> FK5(equinox='J1975')
    <FK5 Frame: equinox=J1975.000>
    >>> ICRS()  # has no attributes
    <ICRS Frame>
    >>> FK5()  # uses default equinox
    <FK5 Frame: equinox=J2000.000>

The specific names of attributes available for a particular frame (and
their default values)  are available as the class attribute
``frame_attr_names``::

    >>> FK5.frame_attr_names
    {'equinox': <Time object: scale='utc' format='jyear_str' value=J2000.000>}

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
^^^^^^^^^^^^^^^^

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
    >>> ICRS(ra=1.1*u.deg, dec=2.2*u.deg)
    <ICRS Coordinate: ra=1.1 deg, dec=2.2 deg>
    >>> FK5(ra=1.1*u.deg, dec=2.2*u.deg, equinox='J1975')
    <FK5 Coordinate: equinox=J1975.000, ra=1.1 deg, dec=2.2 deg>

These same attributes can be used to access the data in the frames, as
|Angle| objects (or |Angle| subclasses)::

    >>> coo = ICRS(ra=1.1*u.deg, dec=2.2*u.deg)
    >>> coo.ra
    <Longitude 1.1... deg>
    >>> coo.ra.value
    1.1...
    >>> coo.ra.to(u.hourangle)
    <Longitude 0.0733... hourangle>

You can use the ``representation`` attribute in conjunction
with the ``representation_names`` attribute to figure out what keywords
are accepted by a particular class object.  The former will be the
representation class the system is expressed in (e.g.,
spherical for equatorial frames), and the latter will be a dictionary
mapping names for that frame to the attribute name on the representation
class::

    >>> import astropy.units as u
    >>> icrs = ICRS(1*u.deg, 2*u.deg)
    >>> icrs.representation
    <class 'astropy.coordinates.representation.SphericalRepresentation'>
    >>> icrs.representation_names
    OrderedDict([(u'ra', u'lon'), (u'dec', u'lat'), (u'distance', u'distance')])

The representation of the coordinate object can be changed, as shown
below.  This actually does *nothing* to the object internal data which
stores the coordinate values, but it changes the external view of that
data in two ways: (1) the object prints itself in accord with the
new representation, and (2) the available attributes change to match
those of the new representation (e.g. from ``ra, dec, distance`` to
``x, y, z``).  Setting the ``representation`` thus changes a *property*
of the object (how it appears) without changing the intrinsic object
itself which represents a point in 3d space.
::

    >>> from astropy.coordinates import CartesianRepresentation
    >>> icrs.representation = CartesianRepresentation
    >>> icrs
    <ICRS Coordinate: x=0.999238... , y=0.017441... , z=0.034899... >
    >>> icrs.x
    <Quantity 0.999238...>

The representation can also be set at the time of creating a coordinate
and affects the set of keywords used to supply the coordinate data.  For
example to create a coordinate with cartesian data do::

    >>> ICRS(x=1*u.kpc, y=2*u.kpc, z=3*u.kpc, representation=CartesianRepresentation)
    <ICRS Coordinate: x=1.0 kpc, y=2.0 kpc, z=3.0 kpc>

For more information about the use of representations in coordinates see the
:ref:`astropy-skycoord-representations` section, and for details about the
representations themselves see :ref:`astropy-coordinates-representations`.

There are two other ways to create frame classes with coordinates.  A
representation class can be passed in directly at creation, along with
any  frame attributes required::

    >>> from astropy.coordinates import SphericalRepresentation
    >>> rep = SphericalRepresentation(lon=1.1*u.deg, lat=2.2*u.deg, distance=3.3*u.kpc)
    >>> FK5(rep, equinox='J1975')
    <FK5 Coordinate: equinox=J1975.000, ra=1.1 deg, dec=2.2 deg, distance=3.3 kpc>

A final way is to create a frame object from an already existing frame
(either one with or without data), using the ``realize_frame`` method. This
will yield a frame with the same attributes, but new data::

    >>> f1 = FK5(equinox='J1975')
    >>> f1
    <FK5 Frame: equinox=J1975.000>
    >>> rep = SphericalRepresentation(lon=1.1*u.deg, lat=2.2*u.deg, distance=3.3*u.kpc)
    >>> f1.realize_frame(rep)
    <FK5 Coordinate: equinox=J1975.000, ra=1.1 deg, dec=2.2 deg, distance=3.3 kpc>

You can check if a frame object has data using the ``has_data`` attribute, and
if it is preset, it can be accessed from the ``data`` attribute::

    >>> ICRS().has_data
    False
    >>> cooi = ICRS(ra=1.1*u.deg, dec=2.2*u.deg)
    >>> cooi.has_data
    True
    >>> cooi.data
    <UnitSphericalRepresentation lon=1.1 deg, lat=2.2 deg>

All of the above methods can also accept array data (or other Python
sequences) to create arrays of coordinates::

    >>> ICRS(ra=[1.5, 2.5]*u.deg, dec=[3.5, 4.5]*u.deg)
    <ICRS Coordinate: (ra, dec) in deg
        [(1.5, 3.5), (2.5, 4.5)]>

If you pass in mixed arrays and scalars, the arrays will be broadcast
over the scalars appropriately::

    >>> ICRS(ra=[1.5, 2.5]*u.deg, dec=[3.5, 4.5]*u.deg, distance=5*u.kpc)
    <ICRS Coordinate: (ra, dec, distance) in (deg, deg, kpc)
        [(1.5, 3.5, 5.0), (2.5, 4.5, 5.0)]>

An additional operation that may be useful is the ability to extract the
data in different representations.  E.g., to get the Cartesian form of
an ICRS coordinate::

    >>> from astropy.coordinates import CartesianRepresentation
    >>> cooi = ICRS(ra=0*u.deg, dec=45*u.deg, distance=10*u.pc)
    >>> cooi.represent_as(CartesianRepresentation)
    <CartesianRepresentation x=7.07106781... pc, y=0.0 pc, z=7.07106781... pc>

Transforming between Frames
===========================

To transform a frame object with data into another frame, use the
``transform_to`` method of an object, and provide it the frame you wish to
transform to.  This frame can either be a frame *class*, in which case
the default attributes will be used, or a frame object (with or without
data)::

    >>> cooi = ICRS(1.5*u.deg, 2.5*u.deg)
    >>> cooi.transform_to(FK5)
    <FK5 Coordinate: equinox=J2000.000, ra=1.500006605... deg, dec=2.500002382... deg>
    >>> cooi.transform_to(FK5(equinox='J1975'))
    <FK5 Coordinate: equinox=J1975.000, ra=1.179603481... deg, dec=2.360853208... deg>

The :ref:`astropy-coordinates-api` includes a list of all of the frames built
into `astropy.coordinates`, as well as the defined transformations between
them. Any transformation that has a valid path, even if it passes through
other frames, can be transformed to. To programmatically check for or
manipulate transformations, see the `~astropy.coordinates.TransformGraph`
documentation.


.. _astropy-coordinates-design:

Defining a New Frame
====================

Users can add new coordinate frames by creating new classes that
are subclasses of  `~astropy.coordinates.BaseCoordinateFrame`.  Detailed
instructions for subclassing are in the docstrings for that class.  The
key aspects are to define the class attributes ``frame_attr_names``and
``representation_attrs``. If these are defined, there is often no need to
define an ``__init__`` function, as the initializer in
`~astropy.coordinates.BaseCoordinateFrame` will probably behave the way
you want.

You can also define arbitrary methods for any added functionality you
want your frame to have that's unique to that frame.  These methods will
be available in any |skycoord| that is created using your user-defined
frame.

For examples of defining frame classes, the first place to look is
probably the source code for the frames that are included in astropy
(available at ``astropy.coordinates.builtin_frames``).  These are not
"magic" in any way, and use all the same API and features available to
user-created frames.  A more annotated example is also available in the
:ref:`sgr-example` documentation section.


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
of the power of this graph is available to user-created frames, meaning
that once you define even one transform from your frame to some frame in
the graph,  coordinates defined in your frame can be transformed to
*any* other frame in the graph.

The transforms themselves are represented as
`~astropy.coordinates.CoordinateTransform` objects or their subclasses. The useful
subclasses/types of transformations are:

* `~astropy.coordinates.FunctionTransform`

    A transform that is defined as a function that takes a frame object
    of one frame class and returns an object of another class.

* `~astropy.coordinates.StaticMatrixTransform`
* `~astropy.coordinates.DynamicMatrixTransform`

    These are both for transformations defined as a 3x3 matrix
    transforming the Cartesian representation of one frame into the
    target frame's Cartesian representation.  The static version is for
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
