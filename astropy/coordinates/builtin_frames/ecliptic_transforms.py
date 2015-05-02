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
from ... import _erfa as erfa

from .icrs import ICRS
from .gcrs import GCRS
from .ecliptic import GeocentricEcliptic, HeliocentricEcliptic
from .utils import cartrepr_from_matmul


def _gcrs_to_geoecliptic_matrix(equinox):
    rnpb = erfa.pnm06a(equinox.jd1, equinox.jd2)
    obl = erfa.obl06(equinox.jd1, equinox.jd2)*u.radian
    return np.dot(rotation_matrix(obl, 'x'), rnpb)


@frame_transform_graph.transform(FunctionTransform, GCRS, GeocentricEcliptic)
def gcrs_to_geoecliptic(from_coo, to_frame):
    if np.any(from_coo.obstime != to_frame.equinox):
        # if they GCRS obstime and ecliptic equinox are not the same, we first
        # have to move to a GCRS where they are.
        frameattrs = from_coo.get_frame_attr_names()
        frameattrs['obstime'] = to_frame.equinox
        from_coo = from_coo.transform_to(GCRS(**frameattrs))

    rmat = _gcrs_to_geoecliptic_matrix(to_frame.equinox)
    newrepr = cartrepr_from_matmul(rmat, from_coo)
    return to_frame.realize_frame(newrepr)


@frame_transform_graph.transform(FunctionTransform, GeocentricEcliptic, GCRS)
def geoecliptic_to_gcrs(from_coo, to_frame):
    rmat = _gcrs_to_geoecliptic_matrix(to_frame.equinox)
    newrepr = cartrepr_from_matmul(rmat, from_coo, transpose=True)

    if np.all(from_coo.equinox == to_frame.obstime):
        return to_frame.realize_frame(newrepr)
    else:
        # if the GCRS obstime and ecliptic equinox don't match, need to move
        # to one where they do
        frameattrs = to_frame.get_frame_attr_names()
        frameattrs['obstime'] = from_coo.equinox
        return GCRS(newrepr, **frameattrs).transform_to(to_frame)



@frame_transform_graph.transform(DynamicMatrixTransform, ICRS, HeliocentricEcliptic)
def icrs_to_helioecliptic(from_coo, to_frame):
    rnpb = erfa.pnm06a(to_frame.equinox.jd1, to_frame.equinox.jd2)

    obl = erfa.obl06(to_frame.equinox.jd1, to_frame.equinox.jd2)*u.radian
    return np.dot(rotation_matrix(obl, 'x'), rnpb)


@frame_transform_graph.transform(DynamicMatrixTransform, HeliocentricEcliptic, ICRS)
def helioecliptic_to_icrs(from_coo, to_frame):
    return icrs_to_helioecliptic(to_frame, from_coo).T
