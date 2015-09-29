# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains the transformation functions for getting to/from ecliptic systems.
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np

from ... import units as u
from ..baseframe import frame_transform_graph
from ..transformations import FunctionTransform, DynamicMatrixTransform
from ..angles import rotation_matrix
from ..representation import CartesianRepresentation
from ... import _erfa as erfa

from .icrs import ICRS
from .gcrs import GCRS
from .ecliptic import GeocentricTrueEcliptic, BarycentricTrueEcliptic, HeliocentricTrueEcliptic
from .utils import cartrepr_from_matmul


def _ecliptic_rotation_matrix(equinox):
    rnpb = erfa.pnm06a(equinox.jd1, equinox.jd2)
    obl = erfa.obl06(equinox.jd1, equinox.jd2)*u.radian
    return np.asarray(np.dot(rotation_matrix(obl, 'x'), rnpb))


@frame_transform_graph.transform(FunctionTransform, GCRS, GeocentricTrueEcliptic)
def gcrs_to_geoecliptic(from_coo, to_frame):
    if np.any(from_coo.obstime != to_frame.equinox):
        # if they GCRS obstime and ecliptic equinox are not the same, we first
        # have to move to a GCRS where they are.
        frameattrs = from_coo.get_frame_attr_names()
        frameattrs['obstime'] = to_frame.equinox
        from_coo = from_coo.transform_to(GCRS(**frameattrs))

    rmat = _ecliptic_rotation_matrix(to_frame.equinox)
    newrepr = cartrepr_from_matmul(rmat, from_coo)
    return to_frame.realize_frame(newrepr)


@frame_transform_graph.transform(FunctionTransform, GeocentricTrueEcliptic, GCRS)
def geoecliptic_to_gcrs(from_coo, to_frame):
    rmat = _ecliptic_rotation_matrix(to_frame.equinox)
    newrepr = cartrepr_from_matmul(rmat, from_coo, transpose=True)

    if np.all(from_coo.equinox == to_frame.obstime):
        return to_frame.realize_frame(newrepr)
    else:
        # if the GCRS obstime and ecliptic equinox don't match, need to move
        # to one where they do
        frameattrs = to_frame.get_frame_attr_names()
        frameattrs['obstime'] = from_coo.equinox
        return GCRS(newrepr, **frameattrs).transform_to(to_frame)


@frame_transform_graph.transform(DynamicMatrixTransform, ICRS, BarycentricTrueEcliptic)
def icrs_to_baryecliptic(from_coo, to_frame):
    return _ecliptic_rotation_matrix(to_frame.equinox)


@frame_transform_graph.transform(DynamicMatrixTransform, BarycentricTrueEcliptic, ICRS)
def baryecliptic_to_icrs(from_coo, to_frame):
    return icrs_to_baryecliptic(to_frame, from_coo).T


@frame_transform_graph.transform(FunctionTransform, ICRS, HeliocentricTrueEcliptic)
def icrs_to_helioecliptic(from_coo, to_frame):
    pvh, pvb = erfa.epv00(to_frame.equinox.jd1, to_frame.equinox.jd2)
    delta_bary_to_helio = pvh[..., 0, :] - pvb[..., 0, :]

    #first offset to heliocentric
    heliocart = from_coo.cartesian.xyz + delta_bary_to_helio * u.au

    # now compute the matrix to precess to the right orientation
    rmat = _ecliptic_rotation_matrix(to_frame.equinox)

    # it's not really ICRS because of the offset, but this is digestible by cartrepr_from_matmul
    newrepr = cartrepr_from_matmul(rmat, ICRS(CartesianRepresentation(heliocart)))
    return to_frame.realize_frame(newrepr)


@frame_transform_graph.transform(FunctionTransform, HeliocentricTrueEcliptic, ICRS)
def helioecliptic_to_icrs(from_coo, to_frame):

    # first un-precess from ecliptic to ICRS orientation
    rmat = _ecliptic_rotation_matrix(from_coo.equinox)
    intermed_repr = cartrepr_from_matmul(rmat, from_coo, transpose=True)

    # now offset back to barycentric, which is the correct center for ICRS
    pvh, pvb = erfa.epv00(from_coo.equinox.jd1, from_coo.equinox.jd2)
    delta_bary_to_helio = pvh[..., 0, :] - pvb[..., 0, :]
    newrepr = CartesianRepresentation(intermed_repr.cartesian.xyz - delta_bary_to_helio*u.au)

    return to_frame.realize_frame(newrepr)
