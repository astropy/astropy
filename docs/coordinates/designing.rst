Designing Coordinate Systems
----------------------------

New coordinate systems can easily be added by users by simply subclassing the
`~astropy.coordinates.coordsystems.SphericalCoordinatesBase` object.
Detailed instructions for subclassing are in the docstrings for that class.  If
defining a latitude/longitude style of coordinate system, the
`_initialize_latlong` method and `_init_docstring_param_templ` attribute are
helpful for automated processing of the inputs.

To define transformations to and from this coordinate, the easiest method is to
define a function that accepts an object in one coordinate system and returns
the other.  Decorate this function with
`~astropy.coordinates.transformations.transform_function` function decorator,
supplying the information to determine which coordinates the function transforms
to or from.  This will register the transformation, allowing any other
coordinate object to use this converter.  You can also use the
`~astropy.coordinates.transformations.static_transform_matrix` and
`~astropy.coordinates.transformations.dynamic_transform_matrix` decorators to
specify the transformation in terms of 3 x 3 cartesian coordinate transformation
matrices (typically rotations).