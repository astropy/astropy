.. include:: references.txt

.. _astropy-coordinates-separations-matching:

Separations, Offsets, Catalog Matching, and Related Functionality
*****************************************************************

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

    >>> import numpy as np
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
consistently.

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


Offsets
=======

Closely related to angular separations are offsets between coordinates. The key
distinction for offsets is generally the concept of a "from" and "to" coordinate
rather than the single scalar angular offset of a separation.
`~astropy.coordinates` contains conveniences to compute some of the common
offsets encountered in astronomy.

The first piece of such functionality is the
:meth:`~astropy.coordinates.SkyCoord.position_angle` method, which gives the
conventional astronomy position angle (positive angles East of North) from one
the |skycoord| it is called on to another given as the argument::

    >>> from astropy.coordinates import SkyCoord
    >>> c1 = SkyCoord(1*u.deg, 1*u.deg, frame='icrs')
    >>> c2 = SkyCoord(2*u.deg, 2*u.deg, frame='icrs')
    >>> c1.position_angle(c2).to(u.deg)  # doctest: +FLOAT_CMP
    <Angle 44.97818294 deg>

The combination of :meth:`~astropy.coordinates.SkyCoord.separation` and
:meth:`~astropy.coordinates.SkyCoord.position_angle` thus give a set of
directional offsets. To do the inverse operation - determining the new
"destination" coordinate given a separation and position angle - the
:meth:`~astropy.coordinates.SkyCoord.directional_offset_by` method is provided::

    >>> from astropy.coordinates import SkyCoord
    >>> c1 = SkyCoord(1*u.deg, 1*u.deg, frame='icrs')
    >>> position_angle = 45 * u.deg
    >>> separation = 1.414 * u.deg
    >>> c1.directional_offset_by(position_angle, separation)  # doctest: +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec) in deg
    (2.0004075, 1.99964588)>

There is also a :meth:`~astropy.coordinates.SkyCoord.spherical_offsets_to` method
for computing angular offsets (e.g., small shifts like you might give a
telescope operator to move from a bright star to a fainter target.)::

    >>> from astropy.coordinates import SkyCoord
    >>> bright_star = SkyCoord('8h50m59.75s', '+11d39m22.15s', frame='icrs')
    >>> faint_galaxy = SkyCoord('8h50m47.92s', '+11d39m32.74s', frame='icrs')
    >>> dra, ddec = bright_star.spherical_offsets_to(faint_galaxy)
    >>> dra.to(u.arcsec)  # doctest: +FLOAT_CMP
    <Angle -173.78873354064126 arcsec>
    >>> ddec.to(u.arcsec)  # doctest: +FLOAT_CMP
    <Angle 10.605103417374696 arcsec>


.. _astropy-skyoffset-frames:

"Sky Offset" Frames
-------------------

To extend the concept of spherical offsets, `~astropy.coordinates` has
a frame class :class:`~astropy.coordinates.builtin_frames.skyoffset.SkyOffsetFrame`
which creates distinct frames that are centered on a specific point.
These are known as "sky offset frames", as they are a convenient way to create
a frame centered on an arbitrary position on the sky, suitable for computing
positional offsets (e.g., for astrometry)::

    >>> from astropy.coordinates import SkyOffsetFrame, ICRS
    >>> center = ICRS(10*u.deg, 45*u.deg)
    >>> center.transform_to(SkyOffsetFrame(origin=center))  # doctest: +SKIP
    <SkyOffsetICRS Coordinate (rotation=0.0 deg, origin=<ICRS Coordinate: (ra, dec) in deg
        (10.0, 45.0)>): (lon, lat) in deg
        (0.0, 0.0)>
    >>> target = ICRS(11*u.deg, 46*u.deg)
    >>> target.transform_to(SkyOffsetFrame(origin=center))  # doctest: +FLOAT_CMP
    <SkyOffsetICRS Coordinate (rotation=0.0 deg, origin=<ICRS Coordinate: (ra, dec) in deg
        ( 10.,  45.)>): (lon, lat) in deg
        ( 0.69474685,  1.00428706)>


Alternatively, the convenience method
:meth:`~astropy.coordinates.SkyCoord.skyoffset_frame` lets you create an skyoffset
frame from an already-existing |SkyCoord|::

    >>> center = SkyCoord(10*u.deg, 45*u.deg)
    >>> aframe = center.skyoffset_frame()
    >>> target.transform_to(aframe)  # doctest: +FLOAT_CMP
    <SkyOffsetICRS Coordinate (rotation=0.0 deg, origin=<ICRS Coordinate: (ra, dec) in deg
        ( 10.,  45.)>): (lon, lat) in deg
        ( 0.69474685,  1.00428706)>
    >>> other = SkyCoord(9*u.deg, 44*u.deg, frame='fk5')
    >>> other.transform_to(aframe)  # doctest: +FLOAT_CMP
    <SkyCoord (SkyOffsetICRS: rotation=0.0 deg, origin=<ICRS Coordinate: (ra, dec) in deg
        ( 10.,  45.)>): (lon, lat) in deg
        (-0.71943945, -0.99556216)>

.. note ::
    While sky offset frames *appear* to be all the same class, this not the
    case: the sky offset frame for each different type of frame for ``origin`` is
    actually a distinct class.  E.g., ``SkyOffsetFrame(origin=ICRS(...))``
    yields an object of class ``SkyOffsetICRS``, *not* ``SkyOffsetFrame``.
    While this is not important for most uses of this class, it is important for
    things like type-checking, because something like
    ``SkyOffsetFrame(origin=ICRS(...)).__class__ is SkyOffsetFrame`` will
    *not* be ``True``, as it would be for most classes.

This same frame is also useful as a tool for defining frames that are relative
to a specific known object, useful for hierarchical physical systems like galaxy
groups.  For example, objects around M31 are sometimes shown in a coordinate
frame aligned with standard ICRA RA/Dec, but on M31::

    >>> m31 = SkyCoord(10.6847083*u.deg, 41.26875*u.deg, frame='icrs')
    >>> ngc147 = SkyCoord(8.3005*u.deg, 48.5087389*u.deg, frame='icrs')
    >>> ngc147_inm31 = ngc147.transform_to(m31.skyoffset_frame())
    >>> xi, eta = ngc147_inm31.lon, ngc147_inm31.lat
    >>> xi  # doctest: +FLOAT_CMP
    <Longitude -1.5920694752086249 deg>
    >>> eta  # doctest: +FLOAT_CMP
    <Latitude 7.261837574183891 deg>

.. note::

    Currently, distance information in the ``origin`` of a
    :class:`~astropy.coordinates.builtin_frames.skyoffset.SkyOffsetFrame` is not
    used to compute any part of the transform. The ``origin`` is only used for
    on-sky rotation.  This may change in the future, however.


.. _astropy-coordinates-matching:

Matching Catalogs
=================

`~astropy.coordinates` leverages the coordinate framework to make it
straightforward to find the closest coordinates in a catalog to a desired set
of other coordinates. For example, assuming ``ra1``/``dec1`` and
``ra2``/``dec2`` are numpy arrays loaded from some file::

    >>> from astropy.coordinates import SkyCoord
    >>> from astropy import units as u
    >>> c = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)  # doctest: +SKIP
    >>> catalog = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)  # doctest: +SKIP
    >>> idx, d2d, d3d = c.match_to_catalog_sky(catalog)  # doctest: +SKIP

The 3-dimensional distances returned ``d3d`` are 3-dimensional distances.
Unless both source (``c``) and catalog (``catalog``) coordinates have
associated distances, this quantity assumes that all sources are at a distance
of 1 (dimensionless).

You can also find the nearest 3d matches, different from the on-sky
separation shown above only when the coordinates were initialized with
a ``distance``::

    >>> c = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree, distance=distance1*u.kpc)  # doctest: +SKIP
    >>> catalog = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree, distance=distance2*u.kpc)  # doctest: +SKIP
    >>> idx, d2d, d3d = c.match_to_catalog_3d(catalog)  # doctest: +SKIP

Now ``idx`` are indices into ``catalog`` that are the closest objects to each
of the coordinates in ``c``, ``d2d`` are the on-sky distances between them, and
``d3d`` are the 3-dimensional distances.  Because coordinate objects support
indexing, ``idx`` enables easy access to the matched set of coordinates in
the catalog::

    >>> matches = catalog[idx]  # doctest: +SKIP
    >>> (matches.separation_3d(c) == d3d).all()  # doctest: +SKIP
    True
    >>> dra, ddec = c.spherical_offsets_to(matches)  # doctest: +SKIP

This functionality can also be accessed from the
:func:`~astropy.coordinates.match_coordinates_sky` and
:func:`~astropy.coordinates.match_coordinates_3d` functions. These
will work on either |skycoord| objects *or* the lower-level frame classes::

    >>> from astropy.coordinates import match_coordinates_sky
    >>> idx, d2d, d3d = match_coordinates_sky(c, catalog)  # doctest: +SKIP
    >>> idx, d2d, d3d = match_coordinates_sky(c.frame, catalog.frame)  # doctest: +SKIP

It is possible to impose a separation constraint (e.g., the maximum separation to be
considered a match) by creating a boolean mask with ``d2d`` or ``d3d``. For example,::

    >>> max_sep = 1.0 * u.arcsec # doctest: +SKIP
    >>> idx, d2d, d3d = c.match_to_catalog_3d(catalog) # doctest: +SKIP
    >>> sep_constraint = d2d < max_sep # doctest: +SKIP
    >>> c_matches = c[sep_constraint] # doctest: +SKIP
    >>> catalog_matches = catalog[idx[sep_constraint]] # doctest: +SKIP

Now, ``c_matches`` and ``catalog_matches`` are the matched sources in ``c``
and ``catalog``, respectively, which are separated by less than 1 arcsecond.

.. _astropy-searching-coordinates:

Searching Around Coordinates
============================

Closely-related functionality can be used to search for *all* coordinates within
a certain distance (either 3D distance or on-sky) of another set of coordinates.
The ``search_around_*`` methods (and functions) provide this functionality,
with an interface very similar to ``match_coordinates_*``::

    >>> idxc, idxcatalog, d2d, d3d = catalog.search_around_sky(c, 1*u.deg)  # doctest: +SKIP
    >>> np.all(d2d < 1*u.deg)  # doctest: +SKIP
    True
    >>> idxc, idxcatalog, d2d, d3d = catalog.search_around_3d(c, 1*u.kpc)  # doctest: +SKIP
    >>> np.all(d3d < 1*u.kpc)  # doctest: +SKIP
    True

The key difference for these methods is that there can be multiple (or no)
matches in ``catalog`` around any locations in ``c``.  Hence, indices into both
``c`` and ``catalog`` are returned instead of just indices into ``catalog``.
These can then be indexed back into the two |skycoord| objects, or, for that
matter, any array with the same order::

    >>> np.all(c[idxc].separation(catalog[idxcatalog]) == d2d)  # doctest: +SKIP
    True
    >>> np.all(c[idxc].separation_3d(catalog[idxcatalog]) == d3d)  # doctest: +SKIP
    True
    >>> print catalog_objectnames[idxcatalog]  # doctest: +SKIP
    ['NGC 1234' 'NGC 4567' ...]

Note, though, that this dual-indexing means that ``search_around_*`` does not
work well if one of the coordinates is a scalar, because the returned index
would not make sense for a scalar::

    >>> scalarc = SkyCoord(1*u.deg, 2*u.deg)  # doctest: +SKIP
    >>> idxscalarc, idxcatalog, d2d, d3d = catalog.search_around_sky(scalarc, 1*u.deg)  # THIS DOESN'T ACTUALLY WORK  # doctest: +SKIP
    >>> scalarc[idxscalarc]  # doctest: +SKIP
    IndexError: 0-d arrays can't be indexed

As a result (and because the ``search_around_*`` algorithm is inefficient in
the scalar case, anyway), the best approach for this scenario is to instead
use the ``separation*`` methods::

    >>> d2d = scalarc.separation(catalog)  # doctest: +SKIP
    >>> catalogmsk = d2d < 1*u.deg  # doctest: +SKIP
    >>> d3d = scalarc.separation_3d(catalog)  # doctest: +SKIP
    >>> catalog3dmsk = d3d < 1*u.kpc  # doctest: +SKIP

The resulting ``catalogmsk`` or ``catalog3dmsk`` variables are boolean arrays
rather than arrays of indices, but in practice they usually can be used in
the same way as ``idxcatalog`` from the above examples.  If you definitely do
need indices instead of boolean masks, you can do:

    >>> idxcatalog = np.where(catalogmsk)[0]  # doctest: +SKIP
    >>> idxcatalog3d = np.where(catalog3dmsk)[0]  # doctest: +SKIP
