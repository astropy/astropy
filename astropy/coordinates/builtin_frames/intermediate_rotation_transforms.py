# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains the transformation functions for getting to/from ITRS, TEME, GCRS, and CIRS.
These are distinct from the ICRS and AltAz functions because they are just
rotations without aberration corrections or offsets.
"""

import numpy as np
import erfa

import astropy.units as u
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.transformations import FunctionTransformWithFiniteDifference
from astropy.coordinates.matrix_utilities import matrix_transpose


from .gcrs import GCRS, PrecessedGeocentric
from .cirs import CIRS
from .itrs import ITRS
from .equatorial import TEME, TETE
from .utils import get_polar_motion, get_jd12, EARTH_CENTER

# # first define helper functions


def teme_to_itrs_mat(time):
    # Sidereal time, rotates from ITRS to mean equinox
    # Use 1982 model for consistency with Vallado et al (2006)
    # http://www.celestrak.com/publications/aiaa/2006-6753/AIAA-2006-6753.pdf
    gst = erfa.gmst82(*get_jd12(time, 'ut1'))

    # Polar Motion
    # Do not include TIO locator s' because it is not used in Vallado 2006
    xp, yp = get_polar_motion(time)
    pmmat = erfa.pom00(xp, yp, 0)

    # rotation matrix
    # c2tcio expects a GCRS->CIRS matrix as it's first argument.
    # Here, we just set that to an I-matrix, because we're already
    # in TEME and the difference between TEME and CIRS is just the
    # rotation by the sidereal time rather than the Earth Rotation Angle
    return erfa.c2tcio(np.eye(3), gst, pmmat)


def gcrs_to_cirs_mat(time):
    # celestial-to-intermediate matrix
    return erfa.c2i06a(*get_jd12(time, 'tt'))


def cirs_to_itrs_mat(time):
    # compute the polar motion p-matrix
    xp, yp = get_polar_motion(time)
    sp = erfa.sp00(*get_jd12(time, 'tt'))
    pmmat = erfa.pom00(xp, yp, sp)

    # now determine the Earth Rotation Angle for the input obstime
    # era00 accepts UT1, so we convert if need be
    era = erfa.era00(*get_jd12(time, 'ut1'))

    # c2tcio expects a GCRS->CIRS matrix, but we just set that to an I-matrix
    # because we're already in CIRS
    return erfa.c2tcio(np.eye(3), era, pmmat)


def tete_to_itrs_mat(time):
    # compute the polar motion p-matrix
    xp, yp = get_polar_motion(time)
    sp = erfa.sp00(*get_jd12(time, 'tt'))
    pmmat = erfa.pom00(xp, yp, sp)

    # now determine the greenwich apparent siderial time for the input obstime
    # we use the 2006A model for consistency with RBPN matrix use in GCRS <-> TETE
    ujd1, ujd2 = get_jd12(time, 'ut1')
    jd1, jd2 = get_jd12(time, 'tt')
    gast = erfa.gst06a(ujd1, ujd2, jd1, jd2)

    # c2tcio expects a GCRS->CIRS matrix, but we just set that to an I-matrix
    # because we're already in CIRS equivalent frame
    return erfa.c2tcio(np.eye(3), gast, pmmat)


def gcrs_precession_mat(equinox):
    gamb, phib, psib, epsa = erfa.pfw06(*get_jd12(equinox, 'tt'))
    return erfa.fw2m(gamb, phib, psib, epsa)


# now the actual transforms

@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, GCRS, TETE)
def gcrs_to_tete(gcrs_coo, tete_frame):
    # first apply GCRS self transform to correct time and observatory position/velocity
    gcrs_coo2 = gcrs_coo.transform_to(GCRS(obstime=tete_frame.obstime,
                                           obsgeoloc=tete_frame.obsgeoloc,
                                           obsgeovel=tete_frame.obsgeovel))

    jd1, jd2 = get_jd12(tete_frame.obstime, 'tt')
    # Classical NPB matrix, IAU 2006/2000A
    # (same as in builtin_frames.utils.get_cip).
    rbpn = erfa.pnm06a(jd1, jd2)

    # These rotations are defined at the geocenter, but can be applied to
    # topocentric positions as well, assuming rigid Earth. See p57 of
    # https://www.usno.navy.mil/USNO/astronomical-applications/publications/Circular_179.pdf
    crepr = gcrs_coo2.cartesian.transform(rbpn)
    return tete_frame.realize_frame(crepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, TETE, GCRS)
def tete_to_gcrs(tete_coo, gcrs_frame):
    # compute the pn matrix, and then multiply by its transpose

    jd1, jd2 = get_jd12(tete_coo.obstime, 'tt')
    # Classical NPB matrix, IAU 2006/2000A
    # (same as in builtin_frames.utils.get_cip).
    rbpn = erfa.pnm06a(jd1, jd2)
    newrepr = tete_coo.cartesian.transform(matrix_transpose(rbpn))

    gcrs = GCRS(newrepr, obstime=tete_coo.obstime, obsgeoloc=tete_coo.obsgeoloc,
                obsgeovel=tete_coo.obsgeovel)

    # now do any needed offsets (no-op if same obstime and 0 pos/vel)
    return gcrs.transform_to(gcrs_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, TETE, ITRS)
def tete_to_itrs(tete_coo, itrs_frame):
    # first get us to TETE at the target obstime, and geocentric position
    tete_coo2 = tete_coo.transform_to(TETE(obstime=itrs_frame.obstime,
                                           obsgeoloc=[0, 0, 0]*u.km,
                                           obsgeovel=[0, 0, 0]*u.km/u.s))

    # now get the pmatrix
    pmat = tete_to_itrs_mat(itrs_frame.obstime)
    crepr = tete_coo2.cartesian.transform(pmat)
    return itrs_frame.realize_frame(crepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ITRS, TETE)
def itrs_to_tete(itrs_coo, tete_frame):
    # compute the pmatrix, and then multiply by its transpose
    pmat = tete_to_itrs_mat(itrs_coo.obstime)
    newrepr = itrs_coo.cartesian.transform(matrix_transpose(pmat))
    tete = TETE(newrepr, obstime=itrs_coo.obstime)

    # now do any needed offsets (no-op if same obstime)
    return tete.transform_to(tete_frame)


def get_location_gcrs(location, obstime, gcrs_to_cirs_matrix):
    """Create a GCRS frame at given location and obstime.

    Helper function that avoids location.get_gcrs (which would
    trigger infinite recursion), and uses the already calculated
    GCRS to CIRS matrix to calculate obsgeoloc and obsgeovel
    required for GCRS.
    """
    # TODO: ideally, GCRS would just use a location too;
    # See gh-10996.
    gcrs_cart = (location.get_itrs(obstime, include_velocity=True)
                 .transform_to(CIRS(obstime=obstime))
                 .cartesian
                 .transform(matrix_transpose(gcrs_to_cirs_matrix)))
    return GCRS(obstime=obstime,
                obsgeoloc=gcrs_cart.without_differentials(),
                obsgeovel=gcrs_cart.differentials['s'].to_cartesian())


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, GCRS, CIRS)
def gcrs_to_cirs(gcrs_coo, cirs_frame):
    # first get the pmatrix
    pmat = gcrs_to_cirs_mat(cirs_frame.obstime)
    # Get GCRS coordinates for the target observer location and time.
    loc_gcrs = get_location_gcrs(cirs_frame.location, cirs_frame.obstime, pmat)
    gcrs_coo2 = gcrs_coo.transform_to(loc_gcrs)
    # Now we are relative to the correct observer, do the transform to CIRS.
    crepr = gcrs_coo2.cartesian.transform(pmat)
    return cirs_frame.realize_frame(crepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, CIRS, GCRS)
def cirs_to_gcrs(cirs_coo, gcrs_frame):
    # compute the pmatrix, and then multiply by its transpose
    pmat = gcrs_to_cirs_mat(cirs_coo.obstime)
    newrepr = cirs_coo.cartesian.transform(matrix_transpose(pmat))
    # Transform to GCRS but still at the CIRS location and obstime.
    loc_gcrs = get_location_gcrs(cirs_coo.location, cirs_coo.obstime, pmat)
    gcrs = loc_gcrs.realize_frame(newrepr)
    # now do any needed offsets (no-op if same obstime and pos/vel)
    return gcrs.transform_to(gcrs_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, CIRS, ITRS)
def cirs_to_itrs(cirs_coo, itrs_frame):
    # first get us to geocentric CIRS at the target obstime
    cirs_coo2 = cirs_coo.transform_to(CIRS(obstime=itrs_frame.obstime,
                                           location=EARTH_CENTER))

    # now get the pmatrix
    pmat = cirs_to_itrs_mat(itrs_frame.obstime)
    crepr = cirs_coo2.cartesian.transform(pmat)
    return itrs_frame.realize_frame(crepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ITRS, CIRS)
def itrs_to_cirs(itrs_coo, cirs_frame):
    # compute the pmatrix, and then multiply by its transpose
    pmat = cirs_to_itrs_mat(itrs_coo.obstime)
    newrepr = itrs_coo.cartesian.transform(matrix_transpose(pmat))
    cirs = CIRS(newrepr, obstime=itrs_coo.obstime)

    # now do any needed offsets (no-op if same obstime)
    return cirs.transform_to(cirs_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ITRS, ITRS)
def itrs_to_itrs(from_coo, to_frame):
    # this self-transform goes through CIRS right now, which implicitly also
    # goes back to ICRS
    return from_coo.transform_to(CIRS()).transform_to(to_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, GCRS, PrecessedGeocentric)
def gcrs_to_precessedgeo(from_coo, to_frame):
    # first get us to GCRS with the right attributes (might be a no-op)
    gcrs_coo = from_coo.transform_to(GCRS(obstime=to_frame.obstime,
                                          obsgeoloc=to_frame.obsgeoloc,
                                          obsgeovel=to_frame.obsgeovel))

    # now precess to the requested equinox
    pmat = gcrs_precession_mat(to_frame.equinox)
    crepr = gcrs_coo.cartesian.transform(pmat)
    return to_frame.realize_frame(crepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, PrecessedGeocentric, GCRS)
def precessedgeo_to_gcrs(from_coo, to_frame):
    # first un-precess
    pmat = gcrs_precession_mat(from_coo.equinox)
    crepr = from_coo.cartesian.transform(matrix_transpose(pmat))
    gcrs_coo = GCRS(crepr,
                    obstime=to_frame.obstime,
                    obsgeoloc=to_frame.obsgeoloc,
                    obsgeovel=to_frame.obsgeovel)

    # then move to the GCRS that's actually desired
    return gcrs_coo.transform_to(to_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, TEME, ITRS)
def teme_to_itrs(teme_coo, itrs_frame):
    # first get us to TEME at the target obstime
    # TODO: self transform?
    teme_coo2 = teme_coo.transform_to(TEME(obstime=itrs_frame.obstime))

    # now get the pmatrix
    pmat = teme_to_itrs_mat(itrs_frame.obstime)
    crepr = teme_coo2.cartesian.transform(pmat)
    return itrs_frame.realize_frame(crepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ITRS, TEME)
def itrs_to_teme(itrs_coo, teme_frame):
    # compute the pmatrix, and then multiply by its transpose
    pmat = teme_to_itrs_mat(itrs_coo.obstime)
    newrepr = itrs_coo.cartesian.transform(matrix_transpose(pmat))
    teme = TEME(newrepr, obstime=itrs_coo.obstime)

    # now do any needed offsets (no-op if same obstime)
    return teme.transform_to(teme_frame)
