Designing Coordinate Systems
----------------------------

TODO: need a better explanation here...

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