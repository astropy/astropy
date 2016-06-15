# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains the transformation functions for getting to/from ecliptic systems.
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np

from ... import units as u
from ...utils.compat import NUMPY_LT_1_10
from ..baseframe import frame_transform_graph
from ..transformations import FunctionTransform, DynamicMatrixTransform
from ..angles import rotation_matrix
from ..representation import CartesianRepresentation
from ... import _erfa as erfa

from .icrs import ICRS
from .gcrs import GCRS
from .ecliptic import GeocentricTrueEcliptic, BarycentricTrueEcliptic, HeliocentricTrueEcliptic
from .utils import cartrepr_from_matmul, get_jd12
from ..errors import UnitsError


def _ecliptic_rotation_matrix(equinox):
    jd1, jd2 = get_jd12(equinox, 'tt')
    rnpb = erfa.pnm06a(jd1, jd2)
    obl = erfa.obl06(jd1, jd2)*u.radian
    """
    The following code is the equivalent of looping over obl and rnpb,
    creating a rotation matrix from obl and then taking
    the dot product of the resulting matrices, finally combining
    into a new array.
    """
    try:
        rmat = np.array([rotation_matrix(this_obl, 'x') for this_obl in obl])
        if NUMPY_LT_1_10:
            result = np.einsum('...ij,...jk->...ik', rmat, rnpb)
        else:
            result = np.matmul(rmat, rnpb)
    except:
        # must be a scalar obliquity
        result = np.asarray(np.dot(rotation_matrix(obl, 'x'), rnpb))
    return result


@frame_transform_graph.transform(FunctionTransform, GCRS, GeocentricTrueEcliptic)
def gcrs_to_geoecliptic(gcrs_coo, to_frame):
    # first get us to a 0 pos/vel GCRS at the target equinox
    gcrs_coo2 = gcrs_coo.transform_to(GCRS(obstime=to_frame.equinox))

    rmat = _ecliptic_rotation_matrix(to_frame.equinox)
    newrepr = cartrepr_from_matmul(rmat, gcrs_coo2)
    return to_frame.realize_frame(newrepr)


@frame_transform_graph.transform(FunctionTransform, GeocentricTrueEcliptic, GCRS)
def geoecliptic_to_gcrs(from_coo, gcrs_frame):
    rmat = _ecliptic_rotation_matrix(from_coo.equinox)
    newrepr = cartrepr_from_matmul(rmat, from_coo, transpose=True)
    gcrs = GCRS(newrepr, obstime=from_coo.equinox)

    # now do any needed offsets (no-op if same obstime and 0 pos/vel)
    return gcrs.transform_to(gcrs_frame)


@frame_transform_graph.transform(DynamicMatrixTransform, ICRS, BarycentricTrueEcliptic)
def icrs_to_baryecliptic(from_coo, to_frame):
    return _ecliptic_rotation_matrix(to_frame.equinox)


@frame_transform_graph.transform(DynamicMatrixTransform, BarycentricTrueEcliptic, ICRS)
def baryecliptic_to_icrs(from_coo, to_frame):
    return icrs_to_baryecliptic(to_frame, from_coo).T


_NEED_ORIGIN_HINT = ("The input {0} coordinates do not have length units. This "
                     "probably means you created coordinates with lat/lon but "
                     "no distance.  Heliocentric<->ICRS transforms cannot "
                     "function in this case because there is an origin shift.")


@frame_transform_graph.transform(FunctionTransform, ICRS, HeliocentricTrueEcliptic)
def icrs_to_helioecliptic(from_coo, to_frame):
    if not u.m.is_equivalent(from_coo.cartesian.x.unit):
        raise UnitsError(_NEED_ORIGIN_HINT.format(from_coo.__class__.__name__))

    pvh, pvb = erfa.epv00(*get_jd12(to_frame.equinox, 'tdb'))
    delta_bary_to_helio = pvh[..., 0, :] - pvb[..., 0, :]

    # first offset to heliocentric
    heliocart = (from_coo.cartesian.xyz).T + delta_bary_to_helio * u.au

    # now compute the matrix to precess to the right orientation
    rmat = _ecliptic_rotation_matrix(to_frame.equinox)

    # it's not really ICRS because of the offset, but this is digestible by cartrepr_from_matmul
    newrepr = cartrepr_from_matmul(rmat, ICRS(CartesianRepresentation(heliocart.T)))
    return to_frame.realize_frame(newrepr)


@frame_transform_graph.transform(FunctionTransform, HeliocentricTrueEcliptic, ICRS)
def helioecliptic_to_icrs(from_coo, to_frame):
    if not u.m.is_equivalent(from_coo.cartesian.x.unit):
        raise UnitsError(_NEED_ORIGIN_HINT.format(from_coo.__class__.__name__))

    # first un-precess from ecliptic to ICRS orientation
    rmat = _ecliptic_rotation_matrix(from_coo.equinox)
    intermed_repr = cartrepr_from_matmul(rmat, from_coo, transpose=True)

    # now offset back to barycentric, which is the correct center for ICRS
    pvh, pvb = erfa.epv00(*get_jd12(from_coo.equinox, 'tdb'))
    delta_bary_to_helio = pvh[..., 0, :] - pvb[..., 0, :]
    newrepr = CartesianRepresentation((intermed_repr.to_cartesian().xyz.T - delta_bary_to_helio*u.au).T)

    return to_frame.realize_frame(newrepr)
