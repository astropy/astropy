.. include:: references.txt

Using and Designing Coordinate Frames
-------------------------------------

.. todo:: rewrite for APE5: @eteq

The key concept is that a
registry of all the transformations is used to determine which
coordinates can convert to others.  When you ask for a transformation,
the registry (or "transformation graph") is searched for the shortest
path from your starting coordinate to your target, and it applies all
of the transformations in that path in series.  This allows only the
simplest transformations to be defined, and the package will
automatically determine how to combine those transformations to get
from one system to another.

TODO: need a better explanation here...

Using Frame Objects
-------------------

TODO: We might want to first introduce SkyCoordinate with many of the
same example as below.

Creating new coordinate objects is of course crucial to using
`~astropy.coordinates`.  The typical way to create a new coordinate object
is to directly initialize your preferred coordinate system using standard
python class creation, using the name of the class representing that
system and a number for the two angles.  For example::

    >>> from astropy.coordinates import ICRS, FK4, Galactic
    >>> import astropy.units as u
    >>> ICRS(187.70592*u.degree, 12.39112*u.degree)
    <ICRS Coordinate: ra=187.70592 deg, dec=12.39112 deg>
    >>> FK4(187.07317*u.degree, 12.66715*u.degree)
    <FK4 Coordinate: equinox=B1950.000, obstime=B1950.000, ra=187.07317 deg, dec=12.66715 deg>
    >>> Galactic(283.77763*u.degree, 74.49108*u.degree)
    <Galactic Coordinate: l=283.77763 deg, b=74.49108 deg>

Note that if you do not provide units explicitly, this will fail::

    >>> ICRS(23, 1)
    Traceback (most recent call last):
        ...
    UnitsError: No unit was given - must be some kind of angle

While the above example uses python numerical types, you can also provide
strings to create coordinates.  Strings will be interpreted using the
`~astropy.coordinates.Angle` class' parsing scheme, and has a guiding
principal of being able to interpret any *unambiguous* string specifying an
angle. For the exact rules for how each string is parsed, see the
`~astropy.coordinates.Angle` documentation.  Some examples::

    >>> ICRS("3h36m29.7888s", "-41d08m15.162342s")
    <SkyCoord (ICRS): ra=54.12412 deg, dec=-41.137545095 deg>
    >>> ICRS("14.12412 hours", "-41:08:15.162342 degrees")
    <SkyCoord (ICRS): ra=211.8618 deg, dec=-41.137545095 deg>
    >>> ICRS("14.12412", "-41:08:15.162342")
    Traceback (most recent call last):
        ...
    UnitsError: No unit specified

It's also possible to create coordinates using lists or `numpy` arrays.  The
same unit rules apply as for scalar angles.::

    >>> ICRS([187.70592, 123.45678]*u.degree, [12.39112, 9.87654]*u.degree)  # doctest: +SKIP
    <ICRS Coordinate: (ra, dec) in deg
        [(187.70592, 12.39112), (123.45677999999998, 9.87654)]>
    >>> ICRS([187.70592*u.degree, 8.23*u.hourangle], [12.39112*u.degree, 1.2*u.radian])  # doctest: +SKIP
    <SkyCoord (ICRS): (ra, dec) in deg
        [(187.70592, 12.39112), (123.44999999999999, 68.75493541569878)]>
    >>> ICRS([187.70592, 123.45678], [12.39112, 9.87654])
    Traceback (most recent call last):
        ...
    UnitsError: No unit was given - must be some kind of angle

.. warning::
    If you try to create an angle using a tuple for each angle instead of a
    list or `numpy` array, it will be interpreted as ``(hours, minutes,
    seconds)`` or ``(degrees, arcmin, arcsec)``.  So if you actually want
    multiple coordinates from a tuple, convert it to a list or array.

TODO: update for SkyCoordinate?
One final way to create coordinates is to copy them from an already
existing coordinate object::


.. _astropy-coordinates-design

Defining a New Frame
====================

New coordinate systems can easily be added by users by subclassing
the `~astropy.coordinates.BaseCoordinateFrame` object.  Detailed
instructions for subclassing are in the docstrings for that class.
See the built-in frame classes (e.g., `~astropy.coordinates.ICRS`) for
examples of how to subclass the base class.

To define transformations to and from this coordinate, the easiest method is
to define a function that accepts an object in one coordinate system and
returns the other. For example, to transform from ICRS to FK5 coordinates,
the transformation operator is a precession matrix. We just need to define
functions that compute the necessary matrices to transform from ICRS to FK5
and vice versa, then decorate the functions to register these transformations
with the global transform graph::

    @frame_transform_graph.transform(DynamicMatrixTransform, ICRS, FK5)
    def icrs_to_fk5(icrscoord, fk5frame):
        ...

    @frame_transform_graph.transform(DynamicMatrixTransform, FK5, ICRS)
    def fk5_to_icrs(fk5coord, icrsframe):
        ...

If the transformation to your coordinate frame of interest is not
representable by a matrix operation, you can also specify a function to
do the actual transformation, and pass the `FunctionTransform` class to
the transform graph decorator instead::

    @frame_transform_graph.transform(FunctionTransform, FK4NoETerms, FK4)
    def fk4_no_e_to_fk4(fk4noecoord, fk4frame):
        ...


Defining Transformations
========================
