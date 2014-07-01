.. include:: references.txt

.. _astropy-coordinates-separations-matching:

Separations, Catalog Matching, and Related Functionality
--------------------------------------------------------

`astropy.coordinates` contains commonly-used tools for comparing or
matching coordinate objects.  Of particular importance are those for
determining separations between coordinates and those for matching a
coordinate (or coordinates) to a catalog.  These are mainly implemented
as methods on the coordinate objects.

Separations
===========

The on-sky separation is easily computed with the
:meth:`astropy.coordinates.BaseCoordinateFrame.separation` or
:meth:`astropy.coordinates.SkyCoord.separation` methods,
which computes the great-circle distance (*not* the small-angle
approximation)::

    >>> from astropy import units as u
    >>> from astropy.coordinates import SkyCoord
    >>> c1 = SkyCoord('5h23m34.5s', '-69d45m22s', frame='icrs')
    >>> c2 = SkyCoord('0h52m44.8s', '-72d49m43s', frame='fk5')
    >>> sep = c1.separation(c2)
    >>> sep  # doctest: +FLOAT_CMP
    <Angle 20.74611447604398 deg>

The returned object is an `~astropy.coordinates.Angle` instance, so it
is straightforward to access the angle in any of several equivalent angular
units::

    >>> sep.radian  # doctest: +FLOAT_CMP
    0.36208800460262575
    >>> sep.hour  # doctest: +FLOAT_CMP
    1.3830742984029323
    >>> sep.arcminute  # doctest: +FLOAT_CMP
    1244.7668685626388
    >>> sep.arcsecond  # doctest: +FLOAT_CMP
    74686.01211375833

Also note that the two input coordinates were not in the same frame -
one is  automatically converted to match the other, ensuring that even
though they are  in different frames, the separation is determined
consistently.  This does mean, however, that a |skycoord| without a
frame cannot be compared in this manner::

    >>> c1 = SkyCoord('5h23m34.5s', '-69d45m22s')
    >>> c2 = SkyCoord('0h52m44.8s', '-72d49m43s')
    >>> sep = c1.separation(c2)  # doctest: +SKIP
    ValueError: Cannot transform to/from this SkyCoord because the frame was not specified at creation.

In addition to the on-sky separation described above,
:meth:`astropy.coordinates.BaseCoordinateFrame.separation_3d` or
:meth:`astropy.coordinates.SkyCoord.separation_3d` methods will
determine the 3D distance between two coordinates that have ``distance``
defined::

    >>> from astropy.coordinates import SkyCoord
    >>> c1 = SkyCoord('5h23m34.5s', '-69d45m22s', distance=70*u.kpc, frame='icrs')
    >>> c2 = SkyCoord('0h52m44.8s', '-72d49m43s', distance=80*u.kpc, frame='icrs')
    >>> sep = c1.separation_3d(c2)
    >>> sep  # doctest: +FLOAT_CMP
    <Distance 28.743988157814094 kpc>

.. _astropy-coordinates-matching:

Matching Catalogs
=================

`~astropy.coordinates` supports leverages the coordinate framework to make it
straightforward to find the closest coordinates in a catalog to a desired set
of other coordinates. For example, assuming ``ra1``/``dec1`` and
``ra2``/``dec2`` are numpy arrays loaded from some file::

    >>> from astropy.coordinates import SkyCoord
    >>> from astropy import units as u
    >>> c = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)  # doctest: +SKIP
    >>> catalog = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)  # doctest: +SKIP
    >>> idx, d2d, d3d = c.match_to_catalog_sky(catalog)  # doctest: +SKIP

You can also find the nearest 3d matches, different from the on-sky
separation shown above only when the coordinates were initialized with
a ``distance``::

    >>> c = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)  # doctest: +SKIP
    >>> catalog = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)  # doctest: +SKIP
    >>> idx, d2d, d3d = c.match_to_catalog_3d(catalog)  # doctest: +SKIP

Now ``idx`` are indices into ``catalog`` that are the closest objects to each
of the coordinates in ``c``, ``d2d`` are the on-sky distances between them, and
``d3d`` are the 3-dimensional distances.  Because coordinate objects support
indexing, ``idx`` enables easy access to the matched set of coordinates in
the catalog::

    >>> matches = catalog[idx]  # doctest: +SKIP
    >>> (matches.separation_3d(c) == d3d).all()  # doctest: +SKIP
    True
    >>> dra = (matches.ra - c.ra).arcmin  # doctest: +SKIP
    >>> ddec = (matches.dec - c.dec).arcmin  # doctest: +SKIP

This functionality can also be accessed from the
:func:`~astropy.coordinates.match_coordinates_sky` and
:func:`~astropy.coordinates.match_coordinates_3d` functions. These
will work on either |skycoord| objects *or* the lower-level frame classes::

    >>> from astropy.coordinates import match_coordinates_sky
    >>> idx, d2d, d3d = match_coordinates_sky(c, catalog)  # doctest: +SKIP
    >>> idx, d2d, d3d = match_coordinates_sky(c.frame, catalog.frame)  # doctest: +SKIP

