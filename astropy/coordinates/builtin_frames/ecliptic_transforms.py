# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains the transformation functions for getting to/from ecliptic systems.
"""
import erfa

from astropy import units as u
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.transformations import (
    FunctionTransformWithFiniteDifference, DynamicMatrixTransform,
    AffineTransform,
)
from astropy.coordinates.matrix_utilities import (rotation_matrix,
                                                  matrix_product,
                                                  matrix_transpose)

from .icrs import ICRS
from .gcrs import GCRS
from .ecliptic import (GeocentricMeanEcliptic, BarycentricMeanEcliptic, HeliocentricMeanEcliptic,
                       GeocentricTrueEcliptic, BarycentricTrueEcliptic, HeliocentricTrueEcliptic,
                       HeliocentricEclipticIAU76, CustomBarycentricEcliptic)
from .utils import get_jd12, get_offset_sun_from_barycenter, EQUINOX_J2000
from astropy.coordinates.errors import UnitsError


def _mean_ecliptic_rotation_matrix(equinox):
    # This code just calls ecm06, which uses the precession matrix according to the
    # IAU 2006 model, but leaves out nutation. This brings the results closer to what
    # other libraries give (see https://github.com/astropy/astropy/pull/6508).
    return erfa.ecm06(*get_jd12(equinox, 'tt'))


def _true_ecliptic_rotation_matrix(equinox):
    # This code calls the same routines as done in pnm06a from ERFA, which
    # retrieves the precession matrix (including frame bias) according to
    # the IAU 2006 model, and including the nutation.
    # This family of systems is less popular
    # (see https://github.com/astropy/astropy/pull/6508).
    jd1, jd2 = get_jd12(equinox, 'tt')
    # Here, we call the three routines from erfa.pnm06a separately,
    # so that we can keep the nutation for calculating the true obliquity
    # (which is a fairly expensive operation); see gh-11000.
    # pnm06a: Fukushima-Williams angles for frame bias and precession.
    # (ERFA names short for F-W's gamma_bar, phi_bar, psi_bar and epsilon_A).
    gamb, phib, psib, epsa = erfa.pfw06(jd1, jd2)
    # pnm06a: Nutation components (in longitude and obliquity).
    dpsi, deps = erfa.nut06a(jd1, jd2)
    # pnm06a: Equinox based nutation x precession x bias matrix.
    rnpb = erfa.fw2m(gamb, phib, psib+dpsi, epsa+deps)
    # calculate the true obliquity of the ecliptic
    obl = erfa.obl06(jd1, jd2)+deps
    return matrix_product(rotation_matrix(obl << u.radian, 'x'), rnpb)


def _obliquity_only_rotation_matrix(obl=erfa.obl80(EQUINOX_J2000.jd1, EQUINOX_J2000.jd2) * u.radian):
    # This code only accounts for the obliquity,
    # which can be passed explicitly.
    # The default value is the IAU 1980 value for J2000,
    # which is computed using obl80 from ERFA:
    #
    # obl = erfa.obl80(EQUINOX_J2000.jd1, EQUINOX_J2000.jd2) * u.radian
    return rotation_matrix(obl, "x")


# MeanEcliptic frames


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 GCRS, GeocentricMeanEcliptic,
                                 finite_difference_frameattr_name='equinox')
def gcrs_to_geoecliptic(gcrs_coo, to_frame):
    # first get us to a 0 pos/vel GCRS at the target equinox
    gcrs_coo2 = gcrs_coo.transform_to(GCRS(obstime=to_frame.obstime))

    rmat = _mean_ecliptic_rotation_matrix(to_frame.equinox)
    newrepr = gcrs_coo2.cartesian.transform(rmat)
    return to_frame.realize_frame(newrepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, GeocentricMeanEcliptic, GCRS)
def geoecliptic_to_gcrs(from_coo, gcrs_frame):
    rmat = _mean_ecliptic_rotation_matrix(from_coo.equinox)
    newrepr = from_coo.cartesian.transform(matrix_transpose(rmat))
    gcrs = GCRS(newrepr, obstime=from_coo.obstime)

    # now do any needed offsets (no-op if same obstime and 0 pos/vel)
    return gcrs.transform_to(gcrs_frame)


@frame_transform_graph.transform(DynamicMatrixTransform, ICRS, BarycentricMeanEcliptic)
def icrs_to_baryecliptic(from_coo, to_frame):
    return _mean_ecliptic_rotation_matrix(to_frame.equinox)


@frame_transform_graph.transform(DynamicMatrixTransform, BarycentricMeanEcliptic, ICRS)
def baryecliptic_to_icrs(from_coo, to_frame):
    return matrix_transpose(icrs_to_baryecliptic(to_frame, from_coo))


_NEED_ORIGIN_HINT = ("The input {0} coordinates do not have length units. This "
                     "probably means you created coordinates with lat/lon but "
                     "no distance.  Heliocentric<->ICRS transforms cannot "
                     "function in this case because there is an origin shift.")


@frame_transform_graph.transform(AffineTransform,
                                 ICRS, HeliocentricMeanEcliptic)
def icrs_to_helioecliptic(from_coo, to_frame):
    if not u.m.is_equivalent(from_coo.cartesian.x.unit):
        raise UnitsError(_NEED_ORIGIN_HINT.format(from_coo.__class__.__name__))

    # get the offset of the barycenter from the Sun
    ssb_from_sun = get_offset_sun_from_barycenter(to_frame.obstime, reverse=True,
                                                  include_velocity=bool(from_coo.data.differentials))

    # now compute the matrix to precess to the right orientation
    rmat = _mean_ecliptic_rotation_matrix(to_frame.equinox)

    return rmat, ssb_from_sun.transform(rmat)


@frame_transform_graph.transform(AffineTransform,
                                 HeliocentricMeanEcliptic, ICRS)
def helioecliptic_to_icrs(from_coo, to_frame):
    if not u.m.is_equivalent(from_coo.cartesian.x.unit):
        raise UnitsError(_NEED_ORIGIN_HINT.format(from_coo.__class__.__name__))

    # first un-precess from ecliptic to ICRS orientation
    rmat = _mean_ecliptic_rotation_matrix(from_coo.equinox)

    # now offset back to barycentric, which is the correct center for ICRS
    sun_from_ssb = get_offset_sun_from_barycenter(from_coo.obstime,
                                                  include_velocity=bool(from_coo.data.differentials))

    return matrix_transpose(rmat), sun_from_ssb


# TrueEcliptic frames


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 GCRS, GeocentricTrueEcliptic,
                                 finite_difference_frameattr_name='equinox')
def gcrs_to_true_geoecliptic(gcrs_coo, to_frame):
    # first get us to a 0 pos/vel GCRS at the target equinox
    gcrs_coo2 = gcrs_coo.transform_to(GCRS(obstime=to_frame.obstime))

    rmat = _true_ecliptic_rotation_matrix(to_frame.equinox)
    newrepr = gcrs_coo2.cartesian.transform(rmat)
    return to_frame.realize_frame(newrepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, GeocentricTrueEcliptic, GCRS)
def true_geoecliptic_to_gcrs(from_coo, gcrs_frame):
    rmat = _true_ecliptic_rotation_matrix(from_coo.equinox)
    newrepr = from_coo.cartesian.transform(matrix_transpose(rmat))
    gcrs = GCRS(newrepr, obstime=from_coo.obstime)

    # now do any needed offsets (no-op if same obstime and 0 pos/vel)
    return gcrs.transform_to(gcrs_frame)


@frame_transform_graph.transform(DynamicMatrixTransform, ICRS, BarycentricTrueEcliptic)
def icrs_to_true_baryecliptic(from_coo, to_frame):
    return _true_ecliptic_rotation_matrix(to_frame.equinox)


@frame_transform_graph.transform(DynamicMatrixTransform, BarycentricTrueEcliptic, ICRS)
def true_baryecliptic_to_icrs(from_coo, to_frame):
    return matrix_transpose(icrs_to_true_baryecliptic(to_frame, from_coo))


@frame_transform_graph.transform(AffineTransform,
                                 ICRS, HeliocentricTrueEcliptic)
def icrs_to_true_helioecliptic(from_coo, to_frame):
    if not u.m.is_equivalent(from_coo.cartesian.x.unit):
        raise UnitsError(_NEED_ORIGIN_HINT.format(from_coo.__class__.__name__))

    # get the offset of the barycenter from the Sun
    ssb_from_sun = get_offset_sun_from_barycenter(to_frame.obstime, reverse=True,
                                                  include_velocity=bool(from_coo.data.differentials))

    # now compute the matrix to precess to the right orientation
    rmat = _true_ecliptic_rotation_matrix(to_frame.equinox)

    return rmat, ssb_from_sun.transform(rmat)


@frame_transform_graph.transform(AffineTransform,
                                 HeliocentricTrueEcliptic, ICRS)
def true_helioecliptic_to_icrs(from_coo, to_frame):
    if not u.m.is_equivalent(from_coo.cartesian.x.unit):
        raise UnitsError(_NEED_ORIGIN_HINT.format(from_coo.__class__.__name__))

    # first un-precess from ecliptic to ICRS orientation
    rmat = _true_ecliptic_rotation_matrix(from_coo.equinox)

    # now offset back to barycentric, which is the correct center for ICRS
    sun_from_ssb = get_offset_sun_from_barycenter(from_coo.obstime,
                                                  include_velocity=bool(from_coo.data.differentials))

    return matrix_transpose(rmat), sun_from_ssb


# Other ecliptic frames


@frame_transform_graph.transform(AffineTransform,
                                 HeliocentricEclipticIAU76, ICRS)
def ecliptic_to_iau76_icrs(from_coo, to_frame):
    # first un-precess from ecliptic to ICRS orientation
    rmat = _obliquity_only_rotation_matrix()

    # now offset back to barycentric, which is the correct center for ICRS
    sun_from_ssb = get_offset_sun_from_barycenter(from_coo.obstime,
                                                  include_velocity=bool(from_coo.data.differentials))

    return matrix_transpose(rmat), sun_from_ssb


@frame_transform_graph.transform(AffineTransform,
                                 ICRS, HeliocentricEclipticIAU76)
def icrs_to_iau76_ecliptic(from_coo, to_frame):
    # get the offset of the barycenter from the Sun
    ssb_from_sun = get_offset_sun_from_barycenter(to_frame.obstime, reverse=True,
                                                  include_velocity=bool(from_coo.data.differentials))

    # now compute the matrix to precess to the right orientation
    rmat = _obliquity_only_rotation_matrix()

    return rmat, ssb_from_sun.transform(rmat)


@frame_transform_graph.transform(DynamicMatrixTransform,
                                 ICRS, CustomBarycentricEcliptic)
def icrs_to_custombaryecliptic(from_coo, to_frame):
    return _obliquity_only_rotation_matrix(to_frame.obliquity)


@frame_transform_graph.transform(DynamicMatrixTransform,
                                 CustomBarycentricEcliptic, ICRS)
def custombaryecliptic_to_icrs(from_coo, to_frame):
    return icrs_to_custombaryecliptic(to_frame, from_coo).T


# Create loopback transformations
frame_transform_graph._add_merged_transform(GeocentricMeanEcliptic, ICRS, GeocentricMeanEcliptic)
frame_transform_graph._add_merged_transform(GeocentricTrueEcliptic, ICRS, GeocentricTrueEcliptic)
frame_transform_graph._add_merged_transform(HeliocentricMeanEcliptic, ICRS, HeliocentricMeanEcliptic)
frame_transform_graph._add_merged_transform(HeliocentricTrueEcliptic, ICRS, HeliocentricTrueEcliptic)
frame_transform_graph._add_merged_transform(HeliocentricEclipticIAU76, ICRS, HeliocentricEclipticIAU76)
frame_transform_graph._add_merged_transform(BarycentricMeanEcliptic, ICRS, BarycentricMeanEcliptic)
frame_transform_graph._add_merged_transform(BarycentricTrueEcliptic, ICRS, BarycentricTrueEcliptic)
frame_transform_graph._add_merged_transform(CustomBarycentricEcliptic, ICRS, CustomBarycentricEcliptic)
