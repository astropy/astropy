# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains functions for matching coordinate catalogs.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from ..extern import six
from .representation import UnitSphericalRepresentation

__all__ = ['match_coordinates_3d', 'match_coordinates_sky']


def match_coordinates_3d(matchcoord, catalogcoord, nthneighbor=1, storekdtree='_kdtree_3d'):
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
        in the ``catalogcoord``, as an attribute in ``catalogcoord`` with the
        provided name.  This dramatically speeds up subsequent calls with the
        same catalog. If False, the KD-Tree is discarded after use.

    Returns
    -------
    idx : integer array
        Indecies into ``catalogcoord`` to get the matched points for each
        ``matchcoord``. Shape matches ``matchcoord``.
    sep2d : `~astropy.coordinates.Angle`
        The on-sky separation between the closest match for each ``matchcoord``
        and the ``matchcoord``. Shape matches ``matchcoord``.
    dist3d : `~astropy.units.Quantity`
        The 3D distance between the closest match for each ``matchcoord`` and
        the ``matchcoord``. Shape matches ``matchcoord``.

    Notes
    -----
    This function requires `SciPy <http://www.scipy.org>`_ to be installed
    or it will fail.
    """
    from warnings import warn

    #without scipy this will immediately fail
    from scipy import spatial
    try:
        KDTree = spatial.cKDTree
    except:
        warn('C-base KD tree not found, falling back on (much slower) '
             'python implementation')
        KDTree = spatial.KDTree

    if storekdtree is True:  # backwards compatibility for pre v0.4
        storekdtree = '_kdtree'

    # figure out where any cached KDTree might be
    if isinstance(storekdtree, six.string_types):
        kdt = getattr(catalogcoord, storekdtree, None)
        if kdt is not None and not isinstance(kdt, KDTree):
            raise ValueError('Invalid `storekdtree` string:' + storekdtree)
    elif isinstance(storekdtree, KDTree):
        kdt = storekdtree
        storekdtree = None
    elif not storekdtree:
        kdt = None
    else:
        raise ValueError('Invalid `storekdtree` argument:' +
                          str(storekdtree))

    if kdt is None:
        #need to build the cartesian KD-tree for the catalog
        cartxyz = catalogcoord.cartesian.xyz
        flatxyz = cartxyz.reshape((3, np.prod(cartxyz.shape) // 3))
        kdt = KDTree(flatxyz.value.T)

    #make sure coordinate systems match
    matchcoord = matchcoord.transform_to(catalogcoord)

    #make sure units match
    catunit = catalogcoord.cartesian.x.unit
    matchxyz = matchcoord.cartesian.xyz.to(catunit)

    matchflatxyz = matchxyz.reshape((3, np.prod(matchxyz.shape) // 3))
    dist, idx = kdt.query(matchflatxyz.T, nthneighbor)

    if nthneighbor > 1:  # query gives 1D arrays if k=1, 2D arrays otherwise
        dist = dist[:, -1]
        idx = idx[:, -1]

    if storekdtree:
        #cache the kdtree in `catalogcoord`
        setattr(catalogcoord, storekdtree, kdt)

    sep2d = catalogcoord[idx].separation(matchcoord)
    return idx.reshape(matchxyz.shape[1:]), sep2d, dist.reshape(matchxyz.shape[1:]) * catunit


def match_coordinates_sky(matchcoord, catalogcoord, nthneighbor=1, storekdtree='_kdtree_sky'):
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
        in the ``catalogcoord``, as an attrbute in ``catalogcoord`` with the
        provided name.  This dramatically speeds up subsequent calls with the
        same catalog. If False, the KD-Tree is discarded after use.

    Returns
    -------
    idx : integer array
        Indecies into ``catalogcoord`` to get the matched points for each
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
    This function requires `SciPy <http://www.scipy.org>`_ to be installed
    or it will fail.
    """

    # send to catalog frame
    newmatch = matchcoord.transform_to(catalogcoord)

    #strip out distance info
    match_urepr = newmatch.data.represent_as(UnitSphericalRepresentation)
    newmatch_u = newmatch.realize_frame(match_urepr)

    cat_urepr = catalogcoord.data.represent_as(UnitSphericalRepresentation)
    newcat_u = catalogcoord.realize_frame(cat_urepr)

    idx, sep2d, sep3d = match_coordinates_3d(newmatch_u, newcat_u, nthneighbor, storekdtree)
    # sep3d is *wrong* above, because the distance information was removed,
    # unless one of the catalogs doesn't have a real distance
    if not (isinstance(catalogcoord.data, UnitSphericalRepresentation) or
            isinstance(newmatch.data, UnitSphericalRepresentation)):
        sep3d = catalogcoord[idx].separation_3d(newmatch)

    #update the kdtree on the actual passed-in coordinate
    if storekdtree:
        setattr(catalogcoord, storekdtree, getattr(newcat_u, storekdtree))

    return idx, sep2d, sep3d
