Matching Catalogs/Finding Nearest Coordinates
---------------------------------------------

`astropy.coordinates` supports leverages the coordinate framework to make it
straightforward to find the closest coordinates in a catalog to a desired set
of other coordinates.  For example, assuming `ra1`/`dec1` and `ra2`/`dec2` are
arrays loaded from some file ::

    >>> #assume ra1/dec1 and ra/dec2 are arrays loaded from some file
    >>> c = coord.ICRSCoordinates(ra1, dec1, unit=(u.degree, u.degree))
    >>> catalog = coord.ICRSCoordinates(ra2, dec2, unit=(u.degree, u.degree))
    >>> idx, d2d, d3d = c1.match_to_catalog(catalog)
    >>> from astropy.coordinates import match_coordinates
    >>> idx, d2d, d3d = match_coordinates(c1, catalog)  # same thing

Now `idx` are indicies into `catalog` that are the closest objects to each of
the coordinates in `c`, `d2d` are the on-sky distances from between them, and
`d3d` are the 3-dimensional distances.

.. warning::
    The `~astropy.coordinates.coordsystems.SphericalCoordinatesBase.match_to_catalog` 
    method and the  `~astropy.coordinates.matching.match_coordinates` function
    both use the three-dimensional distance between the catalog and the 
    coordinates to determine the closest match.  If no `distance` is present 
    for a spherical coordinate, this is the same as the on-sky closest 
    neighbor.  However, if `distance` is present, then it finds the nearest 
    neighbor in 3D space,  which is not necessarily the same thing. If you want
    the on-sky closest object and your coordinate objects have `distance`, you
    should make a copy of the coordinates with the `distance` removed before
    matching.
