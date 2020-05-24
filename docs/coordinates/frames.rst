.. include:: references.txt

.. We call EarthLocation.of_site here first to force the downloading
.. of sites.json so that future doctest output isn't cluttered with
.. "Downloading ... [done]". This can be removed once we have a better
.. way of ignoring output lines based on pattern-matching, e.g.:
.. https://github.com/astropy/pytest-doctestplus/issues/11

.. testsetup::
    >>> from astropy.coordinates import EarthLocation
    >>> EarthLocation.of_site('greenwich') # doctest: +IGNORE_OUTPUT +IGNORE_WARNINGS

Using and Designing Coordinate Frames
*************************************

In `astropy.coordinates`, as outlined in the
:ref:`astropy-coordinates-overview`, subclasses of |baseframe| ("frame
classes") define particular coordinate frames. They can (but do not
*have* to) contain representation objects storing the actual coordinate
data. The actual coordinate transformations are defined as functions
that transform representations between frame classes. This approach
serves to separate high-level user functionality (see :doc:`skycoord`)
and details of how the coordinates are actually stored (see
:doc:`representations`) from the definition of frames and how they are
transformed.

Using Frame Objects
===================

Frames without Data
-------------------

Frame objects have two distinct (but related) uses. The first is
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
their default values) are available as the class method
``get_frame_attr_names``::

    >>> FK5.get_frame_attr_names()
    OrderedDict([('equinox', <Time object: scale='tt' format='jyear_str' value=J2000.000>)])

You can access any of the attributes on a frame by using standard Python
attribute access. Note that for cases like ``equinox``, which are time
inputs, if you pass in any unambiguous time string, it will be converted
into an `~astropy.time.Time` object (see
:ref:`astropy-time-inferring-input`)::

    >>> f = FK5(equinox='J1975')
    >>> f.equinox
    <Time object: scale='tt' format='jyear_str' value=J1975.000>
    >>> f = FK5(equinox='2011-05-15T12:13:14')
    >>> f.equinox
    <Time object: scale='utc' format='isot' value=2011-05-15T12:13:14.000>


Frames with Data
----------------

The second use for frame objects is to store actual realized coordinate
data for frames like those described above. In this use, it is similar
to the |skycoord| class, and in fact, the |skycoord| class internally
uses the frame classes as its implementation. However, the frame
classes have fewer "convenience" features, thereby streamlining the
implementation of frame classes. As such, they are created
similarly to |skycoord| objects. One suggested way is to use
with keywords appropriate for the frame (e.g., ``ra`` and ``dec`` for
equatorial systems)::

    >>> from astropy import units as u
    >>> ICRS(ra=1.1*u.deg, dec=2.2*u.deg)  # doctest: +FLOAT_CMP
    <ICRS Coordinate: (ra, dec) in deg
        (1.1, 2.2)>
    >>> FK5(ra=1.1*u.deg, dec=2.2*u.deg, equinox='J1975')  # doctest: +FLOAT_CMP
    <FK5 Coordinate (equinox=J1975.000): (ra, dec) in deg
        (1.1, 2.2)>

These same attributes can be used to access the data in the frames as
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
keywords are accepted by a particular class object. The former will be the
representation class in which the system is expressed (e.g., spherical for
equatorial frames), and the latter will be a dictionary mapping names for that
frame to the attribute name on the representation class::

    >>> import astropy.units as u
    >>> icrs = ICRS(1*u.deg, 2*u.deg)
    >>> icrs.representation_type
    <class 'astropy.coordinates.representation.SphericalRepresentation'>
    >>> icrs.representation_component_names
    OrderedDict([('ra', 'lon'), ('dec', 'lat'), ('distance', 'distance')])

You can get the data in a different representation if needed::

    >>> icrs.represent_as('cartesian')  # doctest: +FLOAT_CMP
    <CartesianRepresentation (x, y, z) [dimensionless]
         (0.99923861, 0.01744177, 0.0348995)>

.. note::

    In previous versions of Astropy, both the frame attribute and the argument
    to frame classes that are now named ``representation_type`` used to be
    simply ``representation``. The name of this attribute/argument is confusing
    as it points to the representation *class*, not the object containing the
    underlying frame data (which is accessed via the frame attribute ``.data``).
    To clarify, we have renamed ``representation`` to ``representation_type``.
    In this version 3.0, we have only changed the references to this attribute
    in the documentation. In the next major version, we will issue a deprecation
    warning. In two major versions, we will remove the ``.representation``
    attribute and ``representation=`` argument.

The representation of the coordinate object can also be changed directly, as
shown below. This does *nothing* to the object internal data which stores the
coordinate values, but it changes the external view of that data in two ways:
(1) the object prints itself in accord with the new representation, and (2) the
available attributes change to match those of the new representation (e.g., from
``ra, dec, distance`` to ``x, y, z``). Setting the ``representation_type``
thus changes a *property* of the object (how it appears) without changing the
intrinsic object itself which represents a point in 3D space.::

    >>> from astropy.coordinates import CartesianRepresentation
    >>> icrs.representation_type = CartesianRepresentation
    >>> icrs  # doctest: +FLOAT_CMP
    <ICRS Coordinate: (x, y, z) [dimensionless]
        (0.99923861, 0.01744177, 0.0348995)>
    >>> icrs.x  # doctest: +FLOAT_CMP
    <Quantity 0.99923861>

The representation can also be set at the time of creating a coordinate
and affects the set of keywords used to supply the coordinate data. For
example, to create a coordinate with Cartesian data do::

    >>> ICRS(x=1*u.kpc, y=2*u.kpc, z=3*u.kpc, representation_type='cartesian')  #  doctest: +FLOAT_CMP
    <ICRS Coordinate: (x, y, z) in kpc
        (1., 2., 3.)>

For more information about the use of representations in coordinates see the
:ref:`astropy-skycoord-representations` section, and for details about the
representations themselves see :ref:`astropy-coordinates-representations`.

There are two other ways to create frame classes with coordinates. A
representation class can be passed in directly at creation, along with
any frame attributes required::

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
if it is present, it can be accessed from the ``data`` attribute::

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

Similar broadcasting happens if you transform to another frame. For example::

    >>> import numpy as np
    >>> from astropy.coordinates import EarthLocation, AltAz
    >>> coo = ICRS(ra=180.*u.deg, dec=51.477811*u.deg)
    >>> lf = AltAz(location=EarthLocation.of_site('greenwich'),
    ...            obstime=['2012-03-21T00:00:00', '2012-06-21T00:00:00'])
    >>> lcoo = coo.transform_to(lf)  # this can load finals2000A.all # doctest: +REMOTE_DATA +IGNORE_OUTPUT
    >>> lcoo  # doctest: +REMOTE_DATA +FLOAT_CMP
    <AltAz Coordinate (obstime=['2012-03-21T00:00:00.000' '2012-06-21T00:00:00.000'], location=(3980608.9024681724, -102.47522910648239, 4966861.273100675) m, pressure=0.0 hPa, temperature=0.0 deg_C, relative_humidity=0.0, obswl=1.0 micron): (az, alt) in deg
        [( 94.71264944, 89.21424252), (307.69488825, 37.98077771)]>

Above, the shapes — ``()`` for ``coo`` and ``(2,)`` for ``lf`` — were
broadcast against each other. If you wish to determine the positions for a
set of coordinates, you will need to make sure that the shapes allow this::

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
    >>> coo2.transform_to(lf2)  # doctest:  +REMOTE_DATA +FLOAT_CMP
    <AltAz Coordinate (obstime=[['2012-03-21T00:00:00.000' '2012-03-21T00:00:00.000'
      '2012-03-21T00:00:00.000']
     ['2012-06-21T00:00:00.000' '2012-06-21T00:00:00.000'
      '2012-06-21T00:00:00.000']], location=(3980608.90246817, -102.47522911, 4966861.27310068) m, pressure=0.0 hPa, temperature=0.0 deg_C, relative_humidity=0.0, obswl=1.0 micron): (az, alt) in deg
        [[( 93.09845185, 89.21613128), (126.85789663, 25.4660055 ),
          ( 51.37993234, 37.18532527)],
         [(307.71713698, 37.99437658), (231.3740787 , 26.36768329),
          ( 85.4218724 , 89.69297998)]]>

.. Note::
   Frames without data have a ``shape`` that is determined by their frame
   attributes. For frames with data, the ``shape`` always is that of the data;
   any non-scalar attributes are broadcast to have matching shapes
   (as can be seen for ``obstime`` in the last line above).

Transforming between Frames
===========================

To transform a frame object with data into another frame, use the
``transform_to`` method of an object, and provide it the frame you wish to
transform to. This frame can either be a frame *class*, in which case
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
other frames, can be transformed too. To programmatically check for or
manipulate transformations, see the `~astropy.coordinates.TransformGraph`
documentation.


.. _astropy-coordinates-design:

Defining a New Frame
====================

Implementing a new frame class that connects to the ``astropy.coordinates``
infrastructure can be done by subclassing
`~astropy.coordinates.BaseCoordinateFrame`. Some guidance and examples are given
below, but detailed instructions for creating new frames are given in the
docstring of `~astropy.coordinates.BaseCoordinateFrame`.

All frame classes must specify a default representation for the coordinate
positions by, at minimum, defining a ``default_representation`` class attribute
(see :ref:`astropy-coordinates-representations` for more information about the
supported ``Representation`` objects). For example, to create a new frame that,
by default, expects to receive its coordinate data in spherical coordinates, we
would create a subclass as follows::

    >>> from astropy.coordinates import BaseCoordinateFrame
    >>> import astropy.coordinates.representation as r
    >>> class MyFrame1(BaseCoordinateFrame):
    ...     # Specify how coordinate values are represented when outputted
    ...     default_representation = r.SphericalRepresentation

Already, this is a valid frame class::

    >>> fr = MyFrame1(1*u.deg, 2*u.deg)
    >>> fr # doctest: +FLOAT_CMP
    <MyFrame1 Coordinate: (lon, lat) in deg
        (1., 2.)>
    >>> fr.lon # doctest: +FLOAT_CMP
    <Longitude 1. deg>

However, as we have defined it above, (1) the coordinate component names will be
the same as used in the specified ``default_representation`` (in this case,
``lon``, ``lat``, and ``distance`` for longitude, latitude, and distance,
respectively), (2) this frame does not have any additional attributes or
metadata, (3) this frame does not support transformations to any other
coordinate frame, and (4) this frame does not support velocity data. We can
address each of these points by seeing some other ways of customizing frame
subclasses.

Customizing Frame Component Names
---------------------------------

First, as mentioned in the point (1) :ref:`above <astropy-coordinates-design>`,
some frame classes have special names for their components. For example, the
`~astropy.coordinates.ICRS` frame and other equatorial frame classes often use
"Right Ascension" or "RA" in place of longitude, and "Declination" or "Dec." in
place of latitude. These component name overrides, which change the frame
component name defaults taken from the ``Representation`` classes, are defined
by specifying a set of `~astropy.coordinates.RepresentationMapping` instances
(one per component) as a part of defining an additional class attribute on a
frame class: ``frame_specific_representation_info``. This attribute must be a
dictionary, and the keys should be either ``Representation`` or ``Differential``
classes (see below for a discussion about customizing behavior for velocity
components, which is done with the ``Differential`` classes). Using our example
frame implemented above, we can customize it to use the names "R" and "D" instead
of "lon" and "lat"::

    >>> from astropy.coordinates import RepresentationMapping
    >>> class MyFrame2(BaseCoordinateFrame):
    ...     # Specify how coordinate values are represented when outputted
    ...     default_representation = r.SphericalRepresentation
    ...
    ...     # Override component names (e.g., "ra" instead of "lon")
    ...     frame_specific_representation_info = {
    ...         r.SphericalRepresentation: [RepresentationMapping('lon', 'R'),
    ...                                     RepresentationMapping('lat', 'D')]
    ...     }

With this frame, we can now use the names ``R`` and ``D`` to access the frame
data::

    >>> fr = MyFrame2(3*u.deg, 4*u.deg)
    >>> fr # doctest: +FLOAT_CMP
    <MyFrame2 Coordinate: (R, D) in deg
        (3., 4.)>
    >>> fr.R # doctest: +FLOAT_CMP
    <Longitude 3. deg>

We can specify name mappings for any ``Representation`` class in
``astropy.coordinates`` to change the default component names. For example, the
`~astropy.coordinates.Galactic` frame uses the standard longitude and latitude
names "l" and "b" when used with a
`~astropy.coordinates.SphericalRepresentation`, but uses the component names
"x", "y", and "z" when the representation is changed to a
`~astropy.coordinates.CartesianRepresentation`. With our example above, we could
add an additional set of mappings to override the Cartesian component names to
be "a", "b", and "c" instead of the default "x", "y", and "z"::

    >>> class MyFrame3(BaseCoordinateFrame):
    ...     # Specify how coordinate values are represented when outputted
    ...     default_representation = r.SphericalRepresentation
    ...
    ...     # Override component names (e.g., "ra" instead of "lon")
    ...     frame_specific_representation_info = {
    ...         r.SphericalRepresentation: [RepresentationMapping('lon', 'R'),
    ...                                     RepresentationMapping('lat', 'D')],
    ...         r.CartesianRepresentation: [RepresentationMapping('x', 'a'),
    ...                                     RepresentationMapping('y', 'b'),
    ...                                     RepresentationMapping('z', 'c')]
    ...     }

For any `~astropy.coordinates.RepresentationMapping`, you can also specify a
default unit for the component by setting the ``defaultunit`` keyword argument.


Defining Frame Attributes
-------------------------

Second, as indicated by the point (2) in the :ref:`introduction above
<astropy-coordinates-design>`, it is often useful for coordinate frames to allow
specifying frame "attributes" that may specify additional data or parameters
needed in order to fully specify transformations between a given frame and some
other frame. For example, the `~astropy.coordinates.FK5` frame allows specifying
an ``equinox`` that helps define the transformation between
`~astropy.coordinates.FK5` and the `~astropy.coordinates.ICRS` frame. Frame
attributes are defined by creating class attributes that are instances of
`~astropy.coordinates.Attribute` or its subclasses (e.g.,
`~astropy.coordinates.TimeAttribute`, `~astropy.coordinates.QuantityAttribute`,
etc.). If attributes are defined using these classes, there is often no need to
define an ``__init__`` function, as the initializer in
`~astropy.coordinates.BaseCoordinateFrame` will probably behave in the way you
want. Let us now modify the above toy frame class implementation to add two frame
attributes::

    >>> from astropy.coordinates import TimeAttribute, QuantityAttribute
    >>> class MyFrame4(BaseCoordinateFrame):
    ...     # Specify how coordinate values are represented when outputted
    ...     default_representation = r.SphericalRepresentation
    ...
    ...     # Override component names (e.g., "ra" instead of "lon")
    ...     frame_specific_representation_info = {
    ...         r.SphericalRepresentation: [RepresentationMapping('lon', 'R'),
    ...                                     RepresentationMapping('lat', 'D')],
    ...         r.CartesianRepresentation: [RepresentationMapping('x', 'a'),
    ...                                     RepresentationMapping('y', 'b'),
    ...                                     RepresentationMapping('z', 'c')]
    ...     }
    ...
    ...     # Specify frame attributes required to fully specify the frame
    ...     time = TimeAttribute(default='B1950')
    ...     orientation = QuantityAttribute(default=42*u.deg)

Without specifying an initializer, defining these attributes tells the
`~astropy.coordinates.BaseCoordinateFrame` what to expect in terms of additional
arguments passed in to our subclass initializer. For example, when defining a
frame instance with our subclass, we can now optionally specify values for these
attributes::

    >>> fr = MyFrame4(R=1*u.deg, D=2*u.deg, orientation=21*u.deg)
    >>> fr # doctest: +FLOAT_CMP
    <MyFrame4 Coordinate (time=B1950.000, orientation=21.0 deg): (R, D) in deg
        (1., 2.)>

Note that we specified both frame attributes with default values, so they are
optional arguments to the frame initializer. Note also that the frame attributes
now appear in the ``repr`` of the frame instance above. As a bonus, for most of
the ``Attribute`` subclasses, even without defining an initializer, attributes
specified as arguments will be validated. For example, arguments passed in to
`~astropy.coordinates.QuantityAttribute` attributes will be checked that they
have valid and compatible units with the expected attribute units. Using our
frame example above, which expects an ``orientation`` with angular units,
passing in a time results in an error::

    >>> MyFrame4(R=1*u.deg, D=2*u.deg, orientation=55*u.microyear) # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ...
    UnitConversionError: 'uyr' (time) and 'deg' (angle) are not convertible

When defining frame attributes, you do not always have to specify a default
value as long as the ``Attribute`` subclass is able to validate the input. For
example, with the above frame, if the ``orientation`` does not require a default
value but we still want to enforce it to have angular units, we could instead
define it as::

    orientation = QuantityAttribute(unit=u.deg)

In the above case, if ``orientation`` is not specified when a new frame instance
is created, its value will be `None`: Note that it is up to the frame
classes and transformation function implementations to define how to handle a
`None` value. In most cases `None` should signify a special case like "use a
different frame attribute for this value" or similar.

Customizing Display of Attributes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While the default `repr` for coordinate frames is suitable for most cases, you
may want to customize how frame attributes are displayed in certain cases. To
do this you can define a method named ``_astropy_repr_in_frame``. This method
should be defined on the object that is set to the frame attribute itself,
**not** the `~astropy.coordinates.Attribute` descriptor.

For example, you could have an object ``Spam`` which you have as an attribute
of your frame::

  >>> class Spam:
  ...     def _astropy_repr_in_frame(self):
  ...         return "<A can of Spam>"

If your frame has this class as an attribute::

  >>> from astropy.coordinates import Attribute
  >>> class Egg(BaseCoordinateFrame):
  ...     can = Attribute(default=Spam())

When it is displayed by the frame it will use the result of
``_astropy_repr_in_frame``::

  >>> Egg()
  <Egg Frame (can=<A can of Spam>)>


Defining Transformations between Frames
---------------------------------------

As indicated by the point (3) in the :ref:`introduction above
<astropy-coordinates-design>`, a frame class on its own is likely not very
useful until transformations are defined between it and other coordinate frame
classes. The key concept for defining transformations in ``astropy.coordinates``
is the "frame transform graph" (in the "graph theory" sense, not "plot"), which
stores all of the transformations between the built-in frames, as well as tools
for finding the shortest paths through this graph to transform from any frame to
any other by composing the transformations. The power behind this concept is
available to user-created frames as well, meaning that once you define even one
transform from your frame to any frame in the graph, coordinates defined in your
frame can be transformed to *any* other frame in the graph. The "frame transform
graph" is available in code as ``astropy.coordinates.frame_transform_graph``,
which is an instance of the `~astropy.coordinates.TransformGraph` class.

The transformations themselves are represented as
`~astropy.coordinates.CoordinateTransform` objects or their subclasses. The
useful subclasses/types of transformations are:

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

    The matrix transforms are `~astropy.coordinates.AffineTransform`
    transformations without a translation (i.e., only a rotation). The static
    version is for the case where the matrix is independent of the frame
    attributes (e.g., the ICRS->FK5 transformation, because ICRS has no frame
    attributes). The dynamic case is for transformations where the
    transformation matrix depends on the frame attributes of either the
    to or from frame.

Generally, it is not necessary to use these classes directly. Instead,
use methods on the ``frame_transform_graph`` that can be used as function
decorators. Define functions that either do the actual
transformation (for `~astropy.coordinates.FunctionTransform`), or that compute
the necessary transformation matrices to transform. Then decorate the functions
to register these transformations with the frame transform graph::

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
matrix). Hence, in general, it is better to define whatever are the
most natural transformations for a user-defined frame, rather than
worrying about optimizing or caching a transformation to speed up the
process.

For a demonstration of how to define transformation functions that also work for
transforming velocity components, see
:ref:`astropy-coordinate-transform-with-velocities`.


Supporting Velocity Data in Frames
----------------------------------

As alluded to by point (4) in the :ref:`introduction above
<astropy-coordinates-design>`, the examples we have seen above mostly deal with
customizing frame behavior for positional information. (For some context about
how velocities are handled in ``astropy.coordinates``, it may be useful to read
the overview: :ref:`astropy-coordinate-custom-frame-with-velocities`.)

When defining a frame class, it is also possible to set a
``default_differential`` (analogous to ``default_representation``), and to
customize how velocity data components are named. Expanding on our custom frame
example above, we can use `~astropy.coordinates.RepresentationMapping` to
override ``Differential`` component names. The default ``Differential``
components are typically named after the corresponding ``Representation``
component, preceded by ``d_``. So, for example, the longitude ``Differential``
component is, by default, ``d_lon``. However, there are some defaults to be
aware of. Here, if we set the default ``Differential`` class to also be
Spherical, it will implement a set of default "nicer" names for the velocity
components, mapping ``pm_R`` to ``d_lon``, ``pm_D`` to ``d_lat``, and
``radial_velocity`` to ``d_distance`` (taking the previously overridden
longitude and latitude component names)::

    >>> class MyFrame4WithVelocity(BaseCoordinateFrame):
    ...     # Specify how coordinate values are represented when outputted
    ...     default_representation = r.SphericalRepresentation
    ...     default_differential = r.SphericalDifferential
    ...
    ...     # Override component names (e.g., "ra" instead of "lon")
    ...     frame_specific_representation_info = {
    ...         r.SphericalRepresentation: [RepresentationMapping('lon', 'R'),
    ...                                     RepresentationMapping('lat', 'D')],
    ...         r.CartesianRepresentation: [RepresentationMapping('x', 'a'),
    ...                                     RepresentationMapping('y', 'b'),
    ...                                     RepresentationMapping('z', 'c')]
    ...     }
    >>> fr = MyFrame4WithVelocity(R=1*u.deg, D=2*u.deg,
    ...                           pm_R=3*u.mas/u.yr, pm_D=4*u.mas/u.yr)
    >>> fr # doctest: +FLOAT_CMP
    <MyFrame4WithVelocity Coordinate: (R, D) in deg
        (1., 2.)
    (pm_R, pm_D) in mas / yr
        (3., 4.)>

If you want to override the default "nicer" names, you can specify a new key in
the ``frame_specific_representation_info`` for any of the ``Differential``
classes, for example::

    >>> class MyFrame4WithVelocity2(BaseCoordinateFrame):
    ...     # Specify how coordinate values are represented when outputted
    ...     default_representation = r.SphericalRepresentation
    ...     default_differential = r.SphericalDifferential
    ...
    ...     # Override component names (e.g., "ra" instead of "lon")
    ...     frame_specific_representation_info = {
    ...         r.SphericalRepresentation: [RepresentationMapping('lon', 'R'),
    ...                                     RepresentationMapping('lat', 'D')],
    ...         r.CartesianRepresentation: [RepresentationMapping('x', 'a'),
    ...                                     RepresentationMapping('y', 'b'),
    ...                                     RepresentationMapping('z', 'c')],
    ...         r.SphericalDifferential: [RepresentationMapping('d_lon', 'pm1'),
    ...                                   RepresentationMapping('d_lat', 'pm2'),
    ...                                   RepresentationMapping('d_distance', 'rv')]
    ...     }
    >>> fr = MyFrame4WithVelocity2(R=1*u.deg, D=2*u.deg,
    ...                           pm1=3*u.mas/u.yr, pm2=4*u.mas/u.yr)
    >>> fr # doctest: +FLOAT_CMP
    <MyFrame4WithVelocity2 Coordinate: (R, D) in deg
        (1., 2.)
    (pm1, pm2) in mas / yr
        (3., 4.)>


Final Notes
-----------

You can also define arbitrary methods for any added functionality you
want your frame to have that is unique to that frame. These methods will
be available in any |skycoord| that is created using your user-defined
frame.

For examples of defining frame classes, the first place to look is
at the source code for the frames that are included in ``astropy``
(available at ``astropy.coordinates.builtin_frames``). These are not
special-cased, but rather use all of the same API and features available to
user-created frames.

.. topic:: Examples:

    See also :ref:`sphx_glr_generated_examples_coordinates_plot_sgr-coordinate-frame.py`
    for a more annotated example of defining a new coordinate frame.
