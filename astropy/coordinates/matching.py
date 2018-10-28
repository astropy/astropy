# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains functions for matching coordinate catalogs.
"""

import numpy as np

from .representation import UnitSphericalRepresentation
from .. import units as u
from . import Angle

__all__ = ['match_coordinates_3d', 'match_coordinates_sky', 'search_around_3d',
           'search_around_sky']


def match_coordinates_3d(matchcoord, catalogcoord, nthneighbor=1, storekdtree='kdtree_3d'):
    """
    Finds the nearest 3-dimensional matches of a coordinate or coordinates in
    a set of catalog coordinates.

    This finds the 3-dimensional closest neighbor, which is only different
    from the on-sky distance if ``distance`` is set in either ``matchcoord``
    or ``catalogcoord``.

    Parameters
    ----------
    matchcoord : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
        The coordinate(s) to match to the catalog.
    catalogcoord : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
        The base catalog in which to search for matches. Typically this will
        be a coordinate object that is an array (i.e.,
        ``catalogcoord.isscalar == False``)
    nthneighbor : int, optional
        Which closest neighbor to search for.  Typically ``1`` is desired here,
        as that is correct for matching one set of coordinates to another.
        The next likely use case is ``2``, for matching a coordinate catalog
        against *itself* (``1`` is inappropriate because each point will find
        itself as the closest match).
    storekdtree : bool or str, optional
        If a string, will store the KD-Tree used for the computation
        in the ``catalogcoord``, as in ``catalogcoord.cache`` with the
        provided name.  This dramatically speeds up subsequent calls with the
        same catalog. If False, the KD-Tree is discarded after use.

    Returns
    -------
    idx : integer array
        Indices into ``catalogcoord`` to get the matched points for each
        ``matchcoord``. Shape matches ``matchcoord``.
    sep2d : `~astropy.coordinates.Angle`
        The on-sky separation between the closest match for each ``matchcoord``
        and the ``matchcoord``. Shape matches ``matchcoord``.
    dist3d : `~astropy.units.Quantity`
        The 3D distance between the closest match for each ``matchcoord`` and
        the ``matchcoord``. Shape matches ``matchcoord``.

    Notes
    -----
    This function requires `SciPy <https://www.scipy.org/>`_ to be installed
    or it will fail.
    """
    if catalogcoord.isscalar or len(catalogcoord) < 1:
        raise ValueError('The catalog for coordinate matching cannot be a '
                         'scalar or length-0.')

    kdt = _get_cartesian_kdtree(catalogcoord, storekdtree)

    # make sure coordinate systems match
    matchcoord = matchcoord.transform_to(catalogcoord)

    # make sure units match
    catunit = catalogcoord.cartesian.x.unit
    matchxyz = matchcoord.cartesian.xyz.to(catunit)

    matchflatxyz = matchxyz.reshape((3, np.prod(matchxyz.shape) // 3))
    dist, idx = kdt.query(matchflatxyz.T, nthneighbor)

    if nthneighbor > 1:  # query gives 1D arrays if k=1, 2D arrays otherwise
        dist = dist[:, -1]
        idx = idx[:, -1]

    sep2d = catalogcoord[idx].separation(matchcoord)
    return idx.reshape(matchxyz.shape[1:]), sep2d, dist.reshape(matchxyz.shape[1:]) * catunit


def match_coordinates_sky(matchcoord, catalogcoord, nthneighbor=1, storekdtree='kdtree_sky'):
    """
    Finds the nearest on-sky matches of a coordinate or coordinates in
    a set of catalog coordinates.

    This finds the on-sky closest neighbor, which is only different from the
    3-dimensional match if ``distance`` is set in either ``matchcoord``
    or ``catalogcoord``.

    Parameters
    ----------
    matchcoord : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
        The coordinate(s) to match to the catalog.
    catalogcoord : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
        The base catalog in which to search for matches. Typically this will
        be a coordinate object that is an array (i.e.,
        ``catalogcoord.isscalar == False``)
    nthneighbor : int, optional
        Which closest neighbor to search for.  Typically ``1`` is desired here,
        as that is correct for matching one set of coordinates to another.
        The next likely use case is ``2``, for matching a coordinate catalog
        against *itself* (``1`` is inappropriate because each point will find
        itself as the closest match).
    storekdtree : bool or str, optional
        If a string, will store the KD-Tree used for the computation
        in the ``catalogcoord`` in ``catalogcoord.cache`` with the
        provided name.  This dramatically speeds up subsequent calls with the
        same catalog. If False, the KD-Tree is discarded after use.

    Returns
    -------
    idx : integer array
        Indices into ``catalogcoord`` to get the matched points for each
        ``matchcoord``. Shape matches ``matchcoord``.
    sep2d : `~astropy.coordinates.Angle`
        The on-sky separation between the closest match for each
        ``matchcoord`` and the ``matchcoord``. Shape matches ``matchcoord``.
    dist3d : `~astropy.units.Quantity`
        The 3D distance between the closest match for each ``matchcoord`` and
        the ``matchcoord``. Shape matches ``matchcoord``.  If either
        ``matchcoord`` or ``catalogcoord`` don't have a distance, this is the 3D
        distance on the unit sphere, rather than a true distance.

    Notes
    -----
    This function requires `SciPy <https://www.scipy.org/>`_ to be installed
    or it will fail.
    """
    if catalogcoord.isscalar or len(catalogcoord) < 1:
        raise ValueError('The catalog for coordinate matching cannot be a '
                         'scalar or length-0.')

    # send to catalog frame
    newmatch = matchcoord.transform_to(catalogcoord)

    # strip out distance info
    match_urepr = newmatch.data.represent_as(UnitSphericalRepresentation)
    newmatch_u = newmatch.realize_frame(match_urepr)

    cat_urepr = catalogcoord.data.represent_as(UnitSphericalRepresentation)
    newcat_u = catalogcoord.realize_frame(cat_urepr)

    # Check for a stored KD-tree on the passed-in coordinate. Normally it will
    # have a distinct name from the "3D" one, so it's safe to use even though
    # it's based on UnitSphericalRepresentation.
    storekdtree = catalogcoord.cache.get(storekdtree, storekdtree)

    idx, sep2d, sep3d = match_coordinates_3d(newmatch_u, newcat_u, nthneighbor, storekdtree)
    # sep3d is *wrong* above, because the distance information was removed,
    # unless one of the catalogs doesn't have a real distance
    if not (isinstance(catalogcoord.data, UnitSphericalRepresentation) or
            isinstance(newmatch.data, UnitSphericalRepresentation)):
        sep3d = catalogcoord[idx].separation_3d(newmatch)

    # update the kdtree on the actual passed-in coordinate
    if isinstance(storekdtree, str):
        catalogcoord.cache[storekdtree] = newcat_u.cache[storekdtree]
    elif storekdtree is True:
        # the old backwards-compatible name
        catalogcoord.cache['kdtree'] = newcat_u.cache['kdtree']

    return idx, sep2d, sep3d


def search_around_3d(coords1, coords2, distlimit, storekdtree='kdtree_3d'):
    """
    Searches for pairs of points that are at least as close as a specified
    distance in 3D space.

    This is intended for use on coordinate objects with arrays of coordinates,
    not scalars.  For scalar coordinates, it is better to use the
    ``separation_3d`` methods.

    Parameters
    ----------
    coords1 : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
        The first set of coordinates, which will be searched for matches from
        ``coords2`` within ``seplimit``.  Cannot be a scalar coordinate.
    coords2 : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
        The second set of coordinates, which will be searched for matches from
        ``coords1`` within ``seplimit``.  Cannot be a scalar coordinate.
    distlimit : `~astropy.units.Quantity` with distance units
        The physical radius to search within.
    storekdtree : bool or str, optional
        If a string, will store the KD-Tree used in the search with the name
        ``storekdtree`` in ``coords2.cache``. This speeds up subsequent calls
        to this function. If False, the KD-Trees are not saved.

    Returns
    -------
    idx1 : integer array
        Indices into ``coords1`` that matches to the corresponding element of
        ``idx2``. Shape matches ``idx2``.
    idx2 : integer array
        Indices into ``coords2`` that matches to the corresponding element of
        ``idx1``. Shape matches ``idx1``.
    sep2d : `~astropy.coordinates.Angle`
        The on-sky separation between the coordinates. Shape matches ``idx1``
        and ``idx2``.
    dist3d : `~astropy.units.Quantity`
        The 3D distance between the coordinates. Shape matches ``idx1`` and
        ``idx2``. The unit is that of ``coords1``.

    Notes
    -----
    This function requires `SciPy <https://www.scipy.org/>`_ (>=0.12.0)
    to be installed or it will fail.

    If you are using this function to search in a catalog for matches around
    specific points, the convention is for ``coords2`` to be the catalog, and
    ``coords1`` are the points to search around.  While these operations are
    mathematically the same if ``coords1`` and ``coords2`` are flipped, some of
    the optimizations may work better if this convention is obeyed.

    In the current implementation, the return values are always sorted in the
    same order as the ``coords1`` (so ``idx1`` is in ascending order).  This is
    considered an implementation detail, though, so it could change in a future
    release.
    """
    if not distlimit.isscalar:
        raise ValueError('distlimit must be a scalar in search_around_3d')

    if coords1.isscalar or coords2.isscalar:
        raise ValueError('One of the inputs to search_around_3d is a scalar. '
                         'search_around_3d is intended for use with array '
                         'coordinates, not scalars.  Instead, use '
                         '``coord1.separation_3d(coord2) < distlimit`` to find '
                         'the coordinates near a scalar coordinate.')

    if len(coords1) == 0 or len(coords2) == 0:
        # Empty array input: return empty match
        return (np.array([], dtype=int), np.array([], dtype=int),
                Angle([], u.deg),
                u.Quantity([], coords1.distance.unit))

    kdt2 = _get_cartesian_kdtree(coords2, storekdtree)
    cunit = coords2.cartesian.x.unit

    # we convert coord1 to match coord2's frame.  We do it this way
    # so that if the conversion does happen, the KD tree of coord2 at least gets
    # saved. (by convention, coord2 is the "catalog" if that makes sense)
    coords1 = coords1.transform_to(coords2)

    kdt1 = _get_cartesian_kdtree(coords1, storekdtree, forceunit=cunit)

    # this is the *cartesian* 3D distance that corresponds to the given angle
    d = distlimit.to_value(cunit)

    idxs1 = []
    idxs2 = []
    for i, matches in enumerate(kdt1.query_ball_tree(kdt2, d)):
        for match in matches:
            idxs1.append(i)
            idxs2.append(match)
    idxs1 = np.array(idxs1, dtype=int)
    idxs2 = np.array(idxs2, dtype=int)

    if idxs1.size == 0:
        d2ds = Angle([], u.deg)
        d3ds = u.Quantity([], coords1.distance.unit)
    else:
        d2ds = coords1[idxs1].separation(coords2[idxs2])
        d3ds = coords1[idxs1].separation_3d(coords2[idxs2])

    return idxs1, idxs2, d2ds, d3ds


def search_around_sky(coords1, coords2, seplimit, storekdtree='kdtree_sky'):
    """
    Searches for pairs of points that have an angular separation at least as
    close as a specified angle.

    This is intended for use on coordinate objects with arrays of coordinates,
    not scalars.  For scalar coordinates, it is better to use the ``separation``
    methods.

    Parameters
    ----------
    coords1 : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
        The first set of coordinates, which will be searched for matches from
        ``coords2`` within ``seplimit``. Cannot be a scalar coordinate.
    coords2 : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
        The second set of coordinates, which will be searched for matches from
        ``coords1`` within ``seplimit``. Cannot be a scalar coordinate.
    seplimit : `~astropy.units.Quantity` with angle units
        The on-sky separation to search within.
    storekdtree : bool or str, optional
        If a string, will store the KD-Tree used in the search with the name
        ``storekdtree`` in ``coords2.cache``. This speeds up subsequent calls
        to this function. If False, the KD-Trees are not saved.

    Returns
    -------
    idx1 : integer array
        Indices into ``coords1`` that matches to the corresponding element of
        ``idx2``. Shape matches ``idx2``.
    idx2 : integer array
        Indices into ``coords2`` that matches to the corresponding element of
        ``idx1``. Shape matches ``idx1``.
    sep2d : `~astropy.coordinates.Angle`
        The on-sky separation between the coordinates. Shape matches ``idx1``
        and ``idx2``.
    dist3d : `~astropy.units.Quantity`
        The 3D distance between the coordinates. Shape matches ``idx1``
        and ``idx2``; the unit is that of ``coords1``.
        If either ``coords1`` or ``coords2`` don't have a distance,
        this is the 3D distance on the unit sphere, rather than a
        physical distance.

    Notes
    -----
    This function requires `SciPy <https://www.scipy.org/>`_ (>=0.12.0)
    to be installed or it will fail.

    In the current implementation, the return values are always sorted in the
    same order as the ``coords1`` (so ``idx1`` is in ascending order).  This is
    considered an implementation detail, though, so it could change in a future
    release.
    """
    if not seplimit.isscalar:
        raise ValueError('seplimit must be a scalar in search_around_sky')

    if coords1.isscalar or coords2.isscalar:
        raise ValueError('One of the inputs to search_around_sky is a scalar. '
                         'search_around_sky is intended for use with array '
                         'coordinates, not scalars.  Instead, use '
                         '``coord1.separation(coord2) < seplimit`` to find the '
                         'coordinates near a scalar coordinate.')

    if len(coords1) == 0 or len(coords2) == 0:
        # Empty array input: return empty match
        if coords2.distance.unit == u.dimensionless_unscaled:
            distunit = u.dimensionless_unscaled
        else:
            distunit = coords1.distance.unit
        return (np.array([], dtype=int), np.array([], dtype=int),
                Angle([], u.deg),
                u.Quantity([], distunit))

    # we convert coord1 to match coord2's frame.  We do it this way
    # so that if the conversion does happen, the KD tree of coord2 at least gets
    # saved. (by convention, coord2 is the "catalog" if that makes sense)
    coords1 = coords1.transform_to(coords2)

    # strip out distance info
    urepr1 = coords1.data.represent_as(UnitSphericalRepresentation)
    ucoords1 = coords1.realize_frame(urepr1)

    kdt1 = _get_cartesian_kdtree(ucoords1, storekdtree)

    if storekdtree and coords2.cache.get(storekdtree):
        # just use the stored KD-Tree
        kdt2 = coords2.cache[storekdtree]
    else:
        # strip out distance info
        urepr2 = coords2.data.represent_as(UnitSphericalRepresentation)
        ucoords2 = coords2.realize_frame(urepr2)

        kdt2 = _get_cartesian_kdtree(ucoords2, storekdtree)
        if storekdtree:
            # save the KD-Tree in coords2, *not* ucoords2
            coords2.cache['kdtree' if storekdtree is True else storekdtree] = kdt2

    # this is the *cartesian* 3D distance that corresponds to the given angle
    r = (2 * np.sin(Angle(seplimit) / 2.0)).value

    idxs1 = []
    idxs2 = []
    for i, matches in enumerate(kdt1.query_ball_tree(kdt2, r)):
        for match in matches:
            idxs1.append(i)
            idxs2.append(match)
    idxs1 = np.array(idxs1, dtype=int)
    idxs2 = np.array(idxs2, dtype=int)

    if idxs1.size == 0:
        if coords2.distance.unit == u.dimensionless_unscaled:
            distunit = u.dimensionless_unscaled
        else:
            distunit = coords1.distance.unit
        d2ds = Angle([], u.deg)
        d3ds = u.Quantity([], distunit)
    else:
        d2ds = coords1[idxs1].separation(coords2[idxs2])
        try:
            d3ds = coords1[idxs1].separation_3d(coords2[idxs2])
        except ValueError:
            # they don't have distances, so we just fall back on the cartesian
            # distance, computed from d2ds
            d3ds = 2 * np.sin(d2ds / 2.0)

    return idxs1, idxs2, d2ds, d3ds


def _get_cartesian_kdtree(coord, attrname_or_kdt='kdtree', forceunit=None):
    """
    This is a utility function to retrieve (and build/cache, if necessary)
    a 3D cartesian KD-Tree from various sorts of astropy coordinate objects.

    Parameters
    ----------
    coord : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
        The coordinates to build the KD-Tree for.
    attrname_or_kdt : bool or str or KDTree
        If a string, will store the KD-Tree used for the computation in the
        ``coord``, in ``coord.cache`` with the provided name. If given as a
        KD-Tree, it will just be used directly.
    forceunit : unit or None
        If a unit, the cartesian coordinates will convert to that unit before
        being put in the KD-Tree.  If None, whatever unit it's already in
        will be used

    Returns
    -------
    kdt : `~scipy.spatial.cKDTree` or `~scipy.spatial.KDTree`
        The KD-Tree representing the 3D cartesian representation of the input
        coordinates.
    """
    from warnings import warn

    # without scipy this will immediately fail
    from scipy import spatial
    try:
        KDTree = spatial.cKDTree
    except Exception:
        warn('C-based KD tree not found, falling back on (much slower) '
             'python implementation')
        KDTree = spatial.KDTree

    if attrname_or_kdt is True:  # backwards compatibility for pre v0.4
        attrname_or_kdt = 'kdtree'

    # figure out where any cached KDTree might be
    if isinstance(attrname_or_kdt, str):
        kdt = coord.cache.get(attrname_or_kdt, None)
        if kdt is not None and not isinstance(kdt, KDTree):
            raise TypeError('The `attrname_or_kdt` "{0}" is not a scipy KD tree!'.format(attrname_or_kdt))
    elif isinstance(attrname_or_kdt, KDTree):
        kdt = attrname_or_kdt
        attrname_or_kdt = None
    elif not attrname_or_kdt:
        kdt = None
    else:
        raise TypeError('Invalid `attrname_or_kdt` argument for KD-Tree:' +
                         str(attrname_or_kdt))

    if kdt is None:
        # need to build the cartesian KD-tree for the catalog
        if forceunit is None:
            cartxyz = coord.cartesian.xyz
        else:
            cartxyz = coord.cartesian.xyz.to(forceunit)
        flatxyz = cartxyz.reshape((3, np.prod(cartxyz.shape) // 3))
        try:
            # Set compact_nodes=False, balanced_tree=False to use
            # "sliding midpoint" rule, which is much faster than standard for
            # many common use cases
            kdt = KDTree(flatxyz.value.T, compact_nodes=False, balanced_tree=False)
        except TypeError:
            # Python implementation does not take compact_nodes and balanced_tree
            # as arguments.  However, it uses sliding midpoint rule by default
            kdt = KDTree(flatxyz.value.T)

    if attrname_or_kdt:
        # cache the kdtree in `coord`
        coord.cache[attrname_or_kdt] = kdt

    return kdt
