.. doctest_skip

Matching Catalogs/Finding Nearest Coordinates
---------------------------------------------

`astropy.coordinates` supports leverages the coordinate framework to make it
straightforward to find the closest coordinates in a catalog to a desired set
of other coordinates.  For example, assuming `ra1`/`dec1` and `ra2`/`dec2` are
arrays loaded from some file ::

    >>> from astropy.coordinates import ICRS
    >>> from astropy import units as u
    >>> #assume ra1/dec1 and ra/dec2 are arrays loaded from some file
    >>> c = ICRS(ra1, dec1, unit=(u.degree, u.degree))
    >>> catalog = ICRS(ra2, dec2, unit=(u.degree, u.degree))
    >>> idx, d2d, d3d = c.match_to_catalog_sky(catalog)
    >>> from astropy.coordinates import match_coordinates_sky
    >>> idx, d2d, d3d = match_coordinates_sky(c1, catalog)  # same thing

You can also find the nearest 3d matches, different from the above when
the coordinates have distances ::

    >>> #assume ra1/dec1 and ra/dec2 are arrays loaded from some file
    >>> c = ICRS(ra1, dec1, unit=(u.degree, u.degree))
    >>> catalog = ICRS(ra2, dec2, unit=(u.degree, u.degree))
    >>> idx, d2d, d3d = c.match_to_catalog_3d(catalog)

Now `idx` are indices into `catalog` that are the closest objects to each of
the coordinates in `c`, `d2d` are the on-sky distances between them, and
`d3d` are the 3-dimensional distances.  Because coordinate objects support
indexing, `idx` enables easy access to the matched set of coordinates in 
the catalog::

    >>> matches = catalog[idx]
    >>> (matches.separation_3d(c) == d3d).all()
    True
    >>> dra = (matches.ra - c.ra).arcmin
    >>> ddec = (matches.dec - c.dec).arcmin
