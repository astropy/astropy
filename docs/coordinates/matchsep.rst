.. _astropy-coordinates-matching:

Matching Catalogs/Finding Nearest Coordinates
---------------------------------------------

.. todo:: update for SkyCoordinate

`astropy.coordinates` supports leverages the coordinate framework to make it
straightforward to find the closest coordinates in a catalog to a desired set
of other coordinates. For example, assuming ``ra1``/``dec1`` and
``ra2``/``dec2`` are arrays loaded from some file ::

    >>> from astropy.coordinates import SkyCoordinate
    >>> from astropy import units as u
    >>> #assume ra1/dec1 and ra/dec2 are arrays loaded from some file
    >>> c = SkyCoordinate(ra=ra1*u.degree, dec=dec1*u.degree)
    >>> catalog = SkyCoordinate(ra=ra2*u.degree, dec=dec2*u.degree)
    >>> idx, d2d, d3d = c.match_to_catalog_sky(catalog)
    >>> from astropy.coordinates import match_coordinates_sky
    >>> idx, d2d, d3d = match_coordinates_sky(c, catalog)  # same thing

You can also find the nearest 3d matches, different from the above when
the coordinates have distances ::

    >>> #assume ra1/dec1 and ra/dec2 are arrays loaded from some file
    >>> c = SkyCoordinate(ra=ra1*u.degree, dec=dec1*u.degree)
    >>> catalog = SkyCoordinate(ra=ra2*u.degree, dec=dec2*u.degree)
    >>> idx, d2d, d3d = c.match_to_catalog_3d(catalog)

Now ``idx`` are indices into ``catalog`` that are the closest objects to each
of the coordinates in ``c``, ``d2d`` are the on-sky distances between them, and
``d3d`` are the 3-dimensional distances.  Because coordinate objects support
indexing, ``idx`` enables easy access to the matched set of coordinates in
the catalog::

    >>> matches = catalog[idx]
    >>> (matches.separation_3d(c) == d3d).all()
    True
    >>> dra = (matches.ra - c.ra).arcmin
    >>> ddec = (matches.dec - c.dec).arcmin

Separations
===========

The on-sky separation is easily computed with the
:meth:`astropy.coordinates.BaseCoordinateFrame.separation` or
:meth:`astropy.coordinates.SkyCoord.separation` methods,
which computes the great-circle distance (*not* the small-angle
approximation)::

    >>> from astropy import units as u
    >>> from astropy.coordinates import ICRS
    >>> c1 = ICRS('5h23m34.5s', '-69d45m22s')
    >>> c2 = ICRS('0h52m44.8s', '-72d49m43s')
    >>> sep = c1.separation(c2)
    >>> sep
    <Angle 20.746118437707... deg>


The returned object is an `~astropy.coordinates.Angle` instance, so it
is straightforward to access the angle in any of several equivalent angular
units::

    >>> sep.radian
    0.36208807374669...
    >>> sep.hour
    1.38307456251383...
    >>> sep.arcminute
    1244.76710626244...
    >>> sep.arcsecond
    74686.0263757469...

The :meth:`astropy.coordinates.BaseCoordinateFrame.separation_3d` or
:meth:`astropy.coordinates.SkyCoord.separation_3d` methods can similarly
be used to get 3D distance between two coordinates, rather than on-sky::

    >>> from astropy.coordinates import SkyCoord
    >>> c1 = SkyCoord('5h23m34.5s', '-69d45m22s', distance=70*u.kpc, frame='icrs')
    >>> c2 = SkyCoord('0h52m44.8s', '-72d49m43s', distance=80*u.kpc, frame='icrs')
    >>> sep = c1.separation_3d(c2)
    >>> sep
    <Distance 28.743988... kpc>
