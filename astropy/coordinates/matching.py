# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains functions for matching coordinate catalogs.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from .. import units as u

__all__ = ['match_coordinates']

def match_coordinates(matchcoord, catalogcoord, nthneighbor=1, storekdtree=True):
    """
    Matches a coordinate or coordinates to a set of catalog coordinates 
    
    Parameters
    ----------
    matchcoord : `~astropy.coordinates.SphericalCoordinatesBase`
        The coordinates to match to the catalog.
    catalogcoord : `~astropy.coordinates.SphericalCoordinatesBase`
        The base catalog in which to search for matches.
    nthneighbor : int, optional
        Which closest neighbor to search for.  Typically ``1`` is desired here,
        as that is correct for matching one set of coordinates to another.
        The next likely use case is ``2``, for matching a coordinate catalog
        against *itself* (``1`` is inappropriate because each point will find 
        itself as the closest match).
    storekdtree : bool or str, optional
        If True or a string, will store the KD-Tree used for the computation
        in the `catalogcoord`.  This dramatically speeds up subsequent calls
        with the same catalog. If a str, it specifies the attribute name for
        `catalogcoord` that should store the KD-tree.

    Returns
    -------
    idx : integer array
        Indecies into `catalogcoord` to get the matched points for each 
        `matchcoord`. Shape matches `matchcoord`.
    sep2d : `~astropy.units.quantity.Angle` 
        The on-sky separation between the closest match for each `matchcoord` and 
        the `matchcoord`. Shape matches `matchcoord`.
    dist3d : `~astropy.units.quantity.Quantity` 
        The 3D distance between the closest match for each `matchcoord` and 
        the `matchcoord`. Shape matches `matchcoord`.

    Notes
    -----
    This function requires `scipy` to be installed or it will fail.
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

    if storekdtree is True:
        storekdtree = '_kdtree'
    
    # figure out where any cached KDTree might be
    if isinstance(storekdtree, basestring):
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
        cart = catalogcoord.cartesian
        flatxyz = cart.reshape((3, np.prod(cart.shape) // 3))
        kdt = KDTree(flatxyz.value.T)

    #make sure coordinate systems match
    matchcoord = matchcoord.transform_to(catalogcoord.__class__)

    #make sure units match
    catunit = catalogcoord.cartesian.unit
    cart = matchcoord.cartesian.to(catunit) 

    flatxyz = cart.reshape((3, np.prod(cart.shape) // 3))
    dist, idx = kdt.query(flatxyz.T, nthneighbor)

    if nthneighbor > 1:  # query gives 1D arrays if k=1, 2D arrays otherwise
        dist = dist[-1]
        idx = idx[-1]

    if storekdtree:
        #cache the kdtree
        setattr(catalogcoord, storekdtree, kdt)

    #TODO: switch to this once __getitem__ on coordinates is working
    #sep2d = catalogcoord[idx].separation(matchcoord)
    from .angle_utilities import angular_separation
    sep2d = angular_separation(catalogcoord.lonangle[idx], 
                               catalogcoord.latangle[idx], 
                               matchcoord.lonangle, 
                               matchcoord.latangle)

    return idx.reshape(cart.shape[1:]), sep2d, dist.reshape(cart.shape[1:]) * catunit



