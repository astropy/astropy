# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains the transformation functions for getting to/from ecliptic systems.
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ... import units as u
from ..baseframe import frame_transform_graph
from ..transformations import FunctionTransformWithFiniteDifference, DynamicMatrixTransform
from ..matrix_utilities import (rotation_matrix,
                                matrix_product, matrix_transpose)
from ..representation import CartesianRepresentation
from ... import _erfa as erfa

from .icrs import ICRS
from .gcrs import GCRS
from .ecliptic import GeocentricTrueEcliptic, BarycentricTrueEcliptic, HeliocentricTrueEcliptic
from .utils import get_jd12
from ..errors import UnitsError


def _ecliptic_rotation_matrix(equinox):
    # This code calls pmat06 from ERFA, which retrieves the precession
    # matrix (including frame bias) according to the IAU 2006 model, but
    # leaves out the nutation. This matches what ERFA does in the ecm06
    # function and also brings the results closer to what other libraries
    # give (see https://github.com/astropy/astropy/pull/6508). However,
    # notice that this makes the name "TrueEcliptic" misleading, and might
    # be changed in the future (discussion in the same pull request)
    jd1, jd2 = get_jd12(equinox, 'tt')
    rbp = erfa.pmat06(jd1, jd2)
    obl = erfa.obl06(jd1, jd2)*u.radian
    return matrix_product(rotation_matrix(obl, 'x'), rbp)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 GCRS, GeocentricTrueEcliptic,
                                 finite_difference_frameattr_name='equinox')
def gcrs_to_geoecliptic(gcrs_coo, to_frame):
    # first get us to a 0 pos/vel GCRS at the target equinox
    gcrs_coo2 = gcrs_coo.transform_to(GCRS(obstime=to_frame.equinox))

    rmat = _ecliptic_rotation_matrix(to_frame.equinox)
    newrepr = gcrs_coo2.cartesian.transform(rmat)
    return to_frame.realize_frame(newrepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, GeocentricTrueEcliptic, GCRS)
def geoecliptic_to_gcrs(from_coo, gcrs_frame):
    rmat = _ecliptic_rotation_matrix(from_coo.equinox)
    newrepr = from_coo.cartesian.transform(matrix_transpose(rmat))
    gcrs = GCRS(newrepr, obstime=from_coo.equinox)

    # now do any needed offsets (no-op if same obstime and 0 pos/vel)
    return gcrs.transform_to(gcrs_frame)


@frame_transform_graph.transform(DynamicMatrixTransform, ICRS, BarycentricTrueEcliptic)
def icrs_to_baryecliptic(from_coo, to_frame):
    return _ecliptic_rotation_matrix(to_frame.equinox)


@frame_transform_graph.transform(DynamicMatrixTransform, BarycentricTrueEcliptic, ICRS)
def baryecliptic_to_icrs(from_coo, to_frame):
    return matrix_transpose(icrs_to_baryecliptic(to_frame, from_coo))


_NEED_ORIGIN_HINT = ("The input {0} coordinates do not have length units. This "
                     "probably means you created coordinates with lat/lon but "
                     "no distance.  Heliocentric<->ICRS transforms cannot "
                     "function in this case because there is an origin shift.")


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 ICRS, HeliocentricTrueEcliptic,
                                 finite_difference_frameattr_name='equinox')
def icrs_to_helioecliptic(from_coo, to_frame):
    if not u.m.is_equivalent(from_coo.cartesian.x.unit):
        raise UnitsError(_NEED_ORIGIN_HINT.format(from_coo.__class__.__name__))

    # get barycentric sun coordinate
    # this goes here to avoid circular import errors
    from ..solar_system import get_body_barycentric
    bary_sun_pos = get_body_barycentric('sun', to_frame.obstime)

    # offset to heliocentric
    heliocart = from_coo.cartesian - bary_sun_pos

    # now compute the matrix to precess to the right orientation
    rmat = _ecliptic_rotation_matrix(to_frame.equinox)

    newrepr = heliocart.transform(rmat)
    return to_frame.realize_frame(newrepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliocentricTrueEcliptic, ICRS,
                                 finite_difference_frameattr_name='equinox')
def helioecliptic_to_icrs(from_coo, to_frame):
    if not u.m.is_equivalent(from_coo.cartesian.x.unit):
        raise UnitsError(_NEED_ORIGIN_HINT.format(from_coo.__class__.__name__))

    # first un-precess from ecliptic to ICRS orientation
    rmat = _ecliptic_rotation_matrix(from_coo.equinox)
    intermed_repr = from_coo.cartesian.transform(matrix_transpose(rmat))

    # now offset back to barycentric, which is the correct center for ICRS

    # this goes here to avoid circular import errors
    from ..solar_system import get_body_barycentric

    # get barycentric sun coordinate
    bary_sun_pos = get_body_barycentric('sun', from_coo.obstime)

    newrepr = intermed_repr + bary_sun_pos
    return to_frame.realize_frame(newrepr)
