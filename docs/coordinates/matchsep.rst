.. include:: references.txt


.. _astropy-coordinates-separations-matching:

Separations, Offsets, Catalog Matching, and Related Functionality
*****************************************************************

`astropy.coordinates` contains commonly-used tools for comparing or
matching coordinate objects. Of particular importance are those for
determining separations between coordinates and those for matching a
coordinate (or coordinates) to a catalog. These are mainly implemented
as methods on the coordinate objects.

Separations
===========

The on-sky separation is easily computed with the
:meth:`astropy.coordinates.BaseCoordinateFrame.separation` or
:meth:`astropy.coordinates.SkyCoord.separation` methods,
which computes the great-circle distance (*not* the small-angle approximation)::

    >>> import numpy as np
    >>> from astropy import units as u
    >>> from astropy.coordinates import SkyCoord
    >>> c1 = SkyCoord('5h23m34.5s', '-69d45m22s', frame='icrs')
    >>> c2 = SkyCoord('0h52m44.8s', '-72d49m43s', frame='fk5')
    >>> sep = c1.separation(c2)
    >>> sep  # doctest: +FLOAT_CMP
    <Angle 20.74611448 deg>

The returned object is an `~astropy.coordinates.Angle` instance, so it
is possible to access the angle in any of several equivalent angular
units::

    >>> sep.radian  # doctest: +FLOAT_CMP
    0.36208800460262563
    >>> sep.hour  # doctest: +FLOAT_CMP
    1.3830742984029318
    >>> sep.arcminute  # doctest: +FLOAT_CMP
    1244.7668685626384
    >>> sep.arcsecond  # doctest: +FLOAT_CMP
    74686.0121137583

Also note that the two input coordinates were not in the same frame —
one is automatically converted to match the other, ensuring that even
though they are in different frames, the separation is determined
consistently.

In addition to the on-sky separation described above,
:meth:`astropy.coordinates.BaseCoordinateFrame.separation_3d` or
:meth:`astropy.coordinates.SkyCoord.separation_3d` methods will
determine the 3D distance between two coordinates that have ``distance``
defined::

    >>> from astropy import units as u
    >>> from astropy.coordinates import SkyCoord
    >>> c1 = SkyCoord('5h23m34.5s', '-69d45m22s', distance=70*u.kpc, frame='icrs')
    >>> c2 = SkyCoord('0h52m44.8s', '-72d49m43s', distance=80*u.kpc, frame='icrs')
    >>> sep = c1.separation_3d(c2)
    >>> sep  # doctest: +FLOAT_CMP
    <Distance 28.74398816 kpc>


Offsets
=======

Closely related to angular separations are offsets between coordinates. The key
distinction for offsets is generally the concept of a "from" and "to" coordinate
rather than the single scalar angular offset of a separation.
`~astropy.coordinates` contains conveniences to compute some of the common
offsets encountered in astronomy.

The first piece of such functionality is the
:meth:`~astropy.coordinates.SkyCoord.position_angle` method. This method
computes the position angle between one
|skycoord| instance and another (passed as the argument) following the
astronomy convention (positive angles East of North)::

    >>> from astropy import units as u
    >>> from astropy.coordinates import SkyCoord
    >>> c1 = SkyCoord(1*u.deg, 1*u.deg, frame='icrs')
    >>> c2 = SkyCoord(2*u.deg, 2*u.deg, frame='icrs')
    >>> c1.position_angle(c2).to(u.deg)  # doctest: +FLOAT_CMP
    <Angle 44.97818294 deg>

The combination of :meth:`~astropy.coordinates.SkyCoord.separation` and
:meth:`~astropy.coordinates.SkyCoord.position_angle` thus give a set of
directional offsets. To do the inverse operation — determining the new
"destination" coordinate given a separation and position angle — the
:meth:`~astropy.coordinates.SkyCoord.directional_offset_by` method is provided::

    >>> from astropy import units as u
    >>> from astropy.coordinates import SkyCoord
    >>> c1 = SkyCoord(1*u.deg, 1*u.deg, frame='icrs')
    >>> position_angle = 45 * u.deg
    >>> separation = 1.414 * u.deg
    >>> c1.directional_offset_by(position_angle, separation)  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec) in deg
        (2.0004075, 1.99964588)>

This technique is also useful for computing the midpoint (or indeed any point)
between two coordinates in a way that accounts for spherical geometry
(i.e., instead of averaging the RAs/Decs separately)::

    >>> from astropy import units as u
    >>> from astropy.coordinates import SkyCoord
    >>> coord1 = SkyCoord(0*u.deg, 0*u.deg, frame='icrs')
    >>> coord2 = SkyCoord(1*u.deg, 1*u.deg, frame='icrs')
    >>> pa = coord1.position_angle(coord2)
    >>> sep = coord1.separation(coord2)
    >>> coord1.directional_offset_by(pa, sep/2)  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec) in deg
        (0.49996192, 0.50001904)>



There is also a :meth:`~astropy.coordinates.SkyCoord.spherical_offsets_to`
method for computing angular offsets (e.g., small shifts like you might give a
telescope operator to move from a bright star to a fainter target)::

    >>> from astropy import units as u
    >>> from astropy.coordinates import SkyCoord
    >>> bright_star = SkyCoord('8h50m59.75s', '+11d39m22.15s', frame='icrs')
    >>> faint_galaxy = SkyCoord('8h50m47.92s', '+11d39m32.74s', frame='icrs')
    >>> dra, ddec = bright_star.spherical_offsets_to(faint_galaxy)
    >>> dra.to(u.arcsec)  # doctest: +FLOAT_CMP
    <Angle -173.78873354 arcsec>
    >>> ddec.to(u.arcsec)  # doctest: +FLOAT_CMP
    <Angle 10.60510342 arcsec>


.. _astropy-skyoffset-frames:

"Sky Offset" Frames
-------------------

To extend the concept of spherical offsets, `~astropy.coordinates` has
a frame class :class:`~astropy.coordinates.builtin_frames.skyoffset.SkyOffsetFrame`
which creates distinct frames that are centered on a specific point.
These are known as "sky offset frames," as they are a convenient way to create
a frame centered on an arbitrary position on the sky suitable for computing
positional offsets (e.g., for astrometry)::

    >>> from astropy import units as u
    >>> from astropy.coordinates import SkyOffsetFrame, ICRS, SkyCoord
    >>> center = ICRS(10*u.deg, 45*u.deg)
    >>> center.transform_to(SkyOffsetFrame(origin=center)) # doctest: +FLOAT_CMP
    <SkyOffsetICRS Coordinate (rotation=0.0 deg, origin=<ICRS Coordinate: (ra, dec) in deg
        (10., 45.)>): (lon, lat) in deg
        (0., 0.)>
    >>> target = ICRS(11*u.deg, 46*u.deg)
    >>> target.transform_to(SkyOffsetFrame(origin=center))  # doctest: +FLOAT_CMP
    <SkyOffsetICRS Coordinate (rotation=0.0 deg, origin=<ICRS Coordinate: (ra, dec) in deg
        (10., 45.)>): (lon, lat) in deg
        (0.69474685, 1.00428706)>


Alternatively, the convenience method
:meth:`~astropy.coordinates.SkyCoord.skyoffset_frame` lets you create a sky
offset frame from an existing |SkyCoord|::

    >>> center = SkyCoord(10*u.deg, 45*u.deg)
    >>> aframe = center.skyoffset_frame()
    >>> target.transform_to(aframe)  # doctest: +FLOAT_CMP
    <SkyOffsetICRS Coordinate (rotation=0.0 deg, origin=<ICRS Coordinate: (ra, dec) in deg
        (10., 45.)>): (lon, lat) in deg
        (0.69474685, 1.00428706)>
    >>> other = SkyCoord(9*u.deg, 44*u.deg, frame='fk5')
    >>> other.transform_to(aframe)  # doctest: +FLOAT_CMP
    <SkyCoord (SkyOffsetICRS: rotation=0.0 deg, origin=<ICRS Coordinate: (ra, dec) in deg
        (10., 45.)>): (lon, lat) in deg
        (-0.71943945, -0.99556216)>

.. note ::

    While sky offset frames *appear* to be all the same class, this not the
    case: the sky offset frame for each different type of frame for ``origin`` is
    actually a distinct class. E.g., ``SkyOffsetFrame(origin=ICRS(...))``
    yields an object of class ``SkyOffsetICRS``, *not* ``SkyOffsetFrame``.
    While this is not important for most uses of this class, it is important for
    things like type-checking, because something like
    ``SkyOffsetFrame(origin=ICRS(...)).__class__ is SkyOffsetFrame`` will
    *not* be ``True``, as it would be for most classes.

This same frame is also useful as a tool for defining frames that are relative
to a specific, known object useful for hierarchical physical systems like galaxy
groups. For example, objects around M31 are sometimes shown in a coordinate
frame aligned with standard ICRA RA/Dec, but on M31::

    >>> m31 = SkyCoord(10.6847083*u.deg, 41.26875*u.deg, frame='icrs')
    >>> ngc147 = SkyCoord(8.3005*u.deg, 48.5087389*u.deg, frame='icrs')
    >>> ngc147_inm31 = ngc147.transform_to(m31.skyoffset_frame())
    >>> xi, eta = ngc147_inm31.lon, ngc147_inm31.lat
    >>> xi  # doctest: +FLOAT_CMP
    <Longitude -1.59206948 deg>
    >>> eta  # doctest: +FLOAT_CMP
    <Latitude 7.26183757 deg>

.. note::

    Currently, distance information in the ``origin`` of a
    :class:`~astropy.coordinates.builtin_frames.skyoffset.SkyOffsetFrame` is not
    used to compute any part of the transform. The ``origin`` is only used for
    on-sky rotation. This may change in the future, however.


.. _astropy-coordinates-matching:

Matching Catalogs
=================

`~astropy.coordinates` leverages the coordinate framework to make it
possible to find the closest coordinates in a catalog to a desired set
of other coordinates. For example, assuming ``ra1``/``dec1`` and
``ra2``/``dec2`` are NumPy arrays loaded from some file:

.. testsetup::
    >>> ra1 = [5.3517]
    >>> dec1 = [-5.2328]
    >>> distance1 = 1344
    >>> ra2 = [6.459]
    >>> dec2 = [-16.4258]
    >>> distance2 = 8.611

.. doctest-requires:: scipy

    >>> from astropy.coordinates import SkyCoord
    >>> from astropy import units as u
    >>> c = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)
    >>> catalog = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)
    >>> idx, d2d, d3d = c.match_to_catalog_sky(catalog)

The distances returned ``d3d`` are 3-dimensional distances.
Unless both source (``c``) and catalog (``catalog``) coordinates have
associated distances, this quantity assumes that all sources are at a distance
of 1 (dimensionless).

You can also find the nearest 3D matches, different from the on-sky
separation shown above only when the coordinates were initialized with
a ``distance``:

.. doctest-requires:: scipy

    >>> c = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree, distance=distance1*u.kpc)
    >>> catalog = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree, distance=distance2*u.kpc)
    >>> idx, d2d, d3d = c.match_to_catalog_3d(catalog)

Now ``idx`` are indices into ``catalog`` that are the closest objects to each
of the coordinates in ``c``, ``d2d`` are the on-sky distances between them, and
``d3d`` are the 3-dimensional distances. Because coordinate objects support
indexing, ``idx`` enables easy access to the matched set of coordinates in
the catalog:

.. doctest-requires:: scipy

    >>> matches = catalog[idx]
    >>> (matches.separation_3d(c) == d3d).all()
    True
    >>> dra, ddec = c.spherical_offsets_to(matches)

This functionality can also be accessed from the
:func:`~astropy.coordinates.match_coordinates_sky` and
:func:`~astropy.coordinates.match_coordinates_3d` functions. These
will work on either |skycoord| objects *or* the lower-level frame classes:

.. doctest-requires:: scipy

    >>> from astropy.coordinates import match_coordinates_sky
    >>> idx, d2d, d3d = match_coordinates_sky(c, catalog)
    >>> idx, d2d, d3d = match_coordinates_sky(c.frame, catalog.frame)

It is possible to impose a separation constraint (e.g., the maximum separation to be
considered a match) by creating a boolean mask with ``d2d`` or ``d3d``. For example:

.. doctest-requires:: scipy

    >>> max_sep = 1.0 * u.arcsec
    >>> idx, d2d, d3d = c.match_to_catalog_3d(catalog)
    >>> sep_constraint = d2d < max_sep
    >>> c_matches = c[sep_constraint]
    >>> catalog_matches = catalog[idx[sep_constraint]]

Now, ``c_matches`` and ``catalog_matches`` are the matched sources in ``c``
and ``catalog``, respectively, which are separated by less than 1 arcsecond.

.. _astropy-searching-coordinates:

Searching around Coordinates
============================

Closely related functionality can be used to search for *all* coordinates within
a certain distance (either 3D distance or on-sky) of another set of coordinates.
The ``search_around_*`` methods (and functions) provide this functionality,
with an interface very similar to ``match_coordinates_*``:

..  doctest-requires:: scipy

    >>> import numpy as np
    >>> idxc, idxcatalog, d2d, d3d = catalog.search_around_sky(c, 1*u.deg)
    >>> np.all(d2d < 1*u.deg)
    True

.. doctest-requires:: scipy

    >>> idxc, idxcatalog, d2d, d3d = catalog.search_around_3d(c, 1*u.kpc)
    >>> np.all(d3d < 1*u.kpc)
    True

The key difference for these methods is that there can be multiple (or no)
matches in ``catalog`` around any locations in ``c``. Hence, indices into both
``c`` and ``catalog`` are returned instead of just indices into ``catalog``.
These can then be indexed back into the two |skycoord| objects, or, for that
matter, any array with the same order:

..  doctest-requires:: scipy

    >>> np.all(c[idxc].separation(catalog[idxcatalog]) == d2d)
    True
    >>> np.all(c[idxc].separation_3d(catalog[idxcatalog]) == d3d)
    True
    >>> print(catalog_objectnames[idxcatalog]) #doctest: +SKIP
    ['NGC 1234' 'NGC 4567' ...]

Note, though, that this dual-indexing means that ``search_around_*`` does not
work well if one of the coordinates is a scalar, because the returned index
would not make sense for a scalar::

    >>> scalarc = SkyCoord(ra=1*u.deg, dec=2*u.deg, distance=distance1*u.kpc)
    >>> idxscalarc, idxcatalog, d2d, d3d = catalog.search_around_sky(scalarc, 1*u.deg) # doctest: +SKIP
    ValueError: One of the inputs to search_around_sky is a scalar.

As a result (and because the ``search_around_*`` algorithm is inefficient in
the scalar case), the best approach for this scenario is to instead
use the ``separation*`` methods:

..  doctest-requires:: scipy

    >>> d2d = scalarc.separation(catalog)
    >>> catalogmsk = d2d < 1*u.deg
    >>> d3d = scalarc.separation_3d(catalog)
    >>> catalog3dmsk = d3d < 1*u.kpc

The resulting ``catalogmsk`` or ``catalog3dmsk`` variables are boolean arrays
rather than arrays of indices, but in practice they usually can be used in
the same way as ``idxcatalog`` from the above examples. If you definitely do
need indices instead of boolean masks, you can do:

..  doctest-requires:: scipy

    >>> idxcatalog = np.where(catalogmsk)[0]
    >>> idxcatalog3d = np.where(catalog3dmsk)[0]
