# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains the transformation functions for getting to/from ITRS, TEME, GCRS, and CIRS.
These are distinct from the ICRS and AltAz functions because they are just
rotations without aberration corrections or offsets.
"""

import erfa
import numpy as np

from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.matrix_utilities import matrix_transpose
from astropy.coordinates.transformations import FunctionTransformWithFiniteDifference

from .cirs import CIRS
from .equatorial import TEME, TETE
from .gcrs import GCRS, PrecessedGeocentric
from .icrs import ICRS
from .itrs import ITRS
from .utils import get_jd12, get_polar_motion

# # first define helper functions


def teme_to_itrs_mat(time):
    # Sidereal time, rotates from ITRS to mean equinox
    # Use 1982 model for consistency with Vallado et al (2006)
    # https://celestrak.org/publications/aiaa/2006-6753/AIAA-2006-6753.pdf
    gst = erfa.gmst82(*get_jd12(time, "ut1"))

    # Polar Motion
    # Do not include TIO locator s' because it is not used in Vallado 2006
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
    return erfa.c2i06a(*get_jd12(time, "tt"))


def cirs_to_itrs_mat(time):
    # compute the polar motion p-matrix
    xp, yp = get_polar_motion(time)
    sp = erfa.sp00(*get_jd12(time, "tt"))
    pmmat = erfa.pom00(xp, yp, sp)

    # now determine the Earth Rotation Angle for the input obstime
    # era00 accepts UT1, so we convert if need be
    era = erfa.era00(*get_jd12(time, "ut1"))

    # c2tcio expects a GCRS->CIRS matrix, but we just set that to an I-matrix
    # because we're already in CIRS
    return erfa.c2tcio(np.eye(3), era, pmmat)


def tete_to_itrs_mat(time, rbpn=None):
    """Compute the polar motion p-matrix at the given time.

    If the nutation-precession matrix is already known, it should be passed in,
    as this is by far the most expensive calculation.
    """
    xp, yp = get_polar_motion(time)
    sp = erfa.sp00(*get_jd12(time, "tt"))
    pmmat = erfa.pom00(xp, yp, sp)

    # now determine the greenwich apparent sidereal time for the input obstime
    # we use the 2006A model for consistency with RBPN matrix use in GCRS <-> TETE
    ujd1, ujd2 = get_jd12(time, "ut1")
    jd1, jd2 = get_jd12(time, "tt")
    if rbpn is None:
        # erfa.gst06a calls pnm06a to calculate rbpn and then gst06. Use it in
        # favour of getting rbpn with erfa.pnm06a to avoid a possibly large array.
        gast = erfa.gst06a(ujd1, ujd2, jd1, jd2)
    else:
        gast = erfa.gst06(ujd1, ujd2, jd1, jd2, rbpn)

    # c2tcio expects a GCRS->CIRS matrix, but we just set that to an I-matrix
    # because we're already in CIRS equivalent frame
    return erfa.c2tcio(np.eye(3), gast, pmmat)


def gcrs_precession_mat(equinox):
    gamb, phib, psib, epsa = erfa.pfw06(*get_jd12(equinox, "tt"))
    return erfa.fw2m(gamb, phib, psib, epsa)


def get_location_gcrs(location, obstime, ref_to_itrs, gcrs_to_ref):
    """Create a GCRS frame at the location and obstime.

    The reference frame z axis must point to the Celestial Intermediate Pole
    (as is the case for CIRS and TETE).

    This function is here to avoid location.get_gcrs(obstime), which would
    recalculate matrices that are already available below (and return a GCRS
    coordinate, rather than a frame with obsgeoloc and obsgeovel).  Instead,
    it uses the private method that allows passing in the matrices.

    """
    obsgeoloc, obsgeovel = location._get_gcrs_posvel(obstime, ref_to_itrs, gcrs_to_ref)
    return GCRS(obstime=obstime, obsgeoloc=obsgeoloc, obsgeovel=obsgeovel)


# now the actual transforms


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, GCRS, TETE)
def gcrs_to_tete(gcrs_coo, tete_frame):
    # Classical NPB matrix, IAU 2006/2000A
    # (same as in builtin_frames.utils.get_cip).
    rbpn = erfa.pnm06a(*get_jd12(tete_frame.obstime, "tt"))
    # Get GCRS coordinates for the target observer location and time.
    loc_gcrs = get_location_gcrs(
        tete_frame.location,
        tete_frame.obstime,
        tete_to_itrs_mat(tete_frame.obstime, rbpn=rbpn),
        rbpn,
    )
    gcrs_coo2 = gcrs_coo.transform_to(loc_gcrs)
    # Now we are relative to the correct observer, do the transform to TETE.
    # These rotations are defined at the geocenter, but can be applied to
    # topocentric positions as well, assuming rigid Earth. See p57 of
    # https://www.usno.navy.mil/USNO/astronomical-applications/publications/Circular_179.pdf
    crepr = gcrs_coo2.cartesian.transform(rbpn)
    return tete_frame.realize_frame(crepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, TETE, GCRS)
def tete_to_gcrs(tete_coo, gcrs_frame):
    # Compute the pn matrix, and then multiply by its transpose.
    rbpn = erfa.pnm06a(*get_jd12(tete_coo.obstime, "tt"))
    newrepr = tete_coo.cartesian.transform(matrix_transpose(rbpn))
    # We now have a GCRS vector for the input location and obstime.
    # Turn it into a GCRS frame instance.
    loc_gcrs = get_location_gcrs(
        tete_coo.location,
        tete_coo.obstime,
        tete_to_itrs_mat(tete_coo.obstime, rbpn=rbpn),
        rbpn,
    )
    gcrs = loc_gcrs.realize_frame(newrepr)
    # Finally, do any needed offsets (no-op if same obstime and location)
    return gcrs.transform_to(gcrs_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, TETE, ITRS)
def tete_to_itrs(tete_coo, itrs_frame):
    # first get us to TETE at the target obstime, and location (no-op if same)
    tete_coo2 = tete_coo.transform_to(
        TETE(obstime=itrs_frame.obstime, location=itrs_frame.location)
    )

    # now get the pmatrix
    pmat = tete_to_itrs_mat(itrs_frame.obstime)
    crepr = tete_coo2.cartesian.transform(pmat)
    return itrs_frame.realize_frame(crepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ITRS, TETE)
def itrs_to_tete(itrs_coo, tete_frame):
    # compute the pmatrix, and then multiply by its transpose
    pmat = tete_to_itrs_mat(itrs_coo.obstime)
    newrepr = itrs_coo.cartesian.transform(matrix_transpose(pmat))
    tete = TETE(newrepr, obstime=itrs_coo.obstime, location=itrs_coo.location)

    # now do any needed offsets (no-op if same obstime and location)
    return tete.transform_to(tete_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, GCRS, CIRS)
def gcrs_to_cirs(gcrs_coo, cirs_frame):
    # first get the pmatrix
    pmat = gcrs_to_cirs_mat(cirs_frame.obstime)
    # Get GCRS coordinates for the target observer location and time.
    loc_gcrs = get_location_gcrs(
        cirs_frame.location,
        cirs_frame.obstime,
        cirs_to_itrs_mat(cirs_frame.obstime),
        pmat,
    )
    gcrs_coo2 = gcrs_coo.transform_to(loc_gcrs)
    # Now we are relative to the correct observer, do the transform to CIRS.
    crepr = gcrs_coo2.cartesian.transform(pmat)
    return cirs_frame.realize_frame(crepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, CIRS, GCRS)
def cirs_to_gcrs(cirs_coo, gcrs_frame):
    # Compute the pmatrix, and then multiply by its transpose,
    pmat = gcrs_to_cirs_mat(cirs_coo.obstime)
    newrepr = cirs_coo.cartesian.transform(matrix_transpose(pmat))
    # We now have a GCRS vector for the input location and obstime.
    # Turn it into a GCRS frame instance.
    loc_gcrs = get_location_gcrs(
        cirs_coo.location, cirs_coo.obstime, cirs_to_itrs_mat(cirs_coo.obstime), pmat
    )
    gcrs = loc_gcrs.realize_frame(newrepr)
    # Finally, do any needed offsets (no-op if same obstime and location)
    return gcrs.transform_to(gcrs_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, CIRS, ITRS)
def cirs_to_itrs(cirs_coo, itrs_frame):
    # first get us to CIRS at the target obstime, and location (no-op if same)
    cirs_coo2 = cirs_coo.transform_to(
        CIRS(obstime=itrs_frame.obstime, location=itrs_frame.location)
    )

    # now get the pmatrix
    pmat = cirs_to_itrs_mat(itrs_frame.obstime)
    crepr = cirs_coo2.cartesian.transform(pmat)
    return itrs_frame.realize_frame(crepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ITRS, CIRS)
def itrs_to_cirs(itrs_coo, cirs_frame):
    # compute the pmatrix, and then multiply by its transpose
    pmat = cirs_to_itrs_mat(itrs_coo.obstime)
    newrepr = itrs_coo.cartesian.transform(matrix_transpose(pmat))
    cirs = CIRS(newrepr, obstime=itrs_coo.obstime, location=itrs_coo.location)

    # now do any needed offsets (no-op if same obstime and location)
    return cirs.transform_to(cirs_frame)


# TODO: implement GCRS<->CIRS if there's call for it.  The thing that's awkward
# is that they both have obstimes, so an extra set of transformations are necessary.
# so unless there's a specific need for that, better to just have it go through the above
# two steps anyway


@frame_transform_graph.transform(
    FunctionTransformWithFiniteDifference, GCRS, PrecessedGeocentric
)
def gcrs_to_precessedgeo(from_coo, to_frame):
    # first get us to GCRS with the right attributes (might be a no-op)
    gcrs_coo = from_coo.transform_to(
        GCRS(
            obstime=to_frame.obstime,
            obsgeoloc=to_frame.obsgeoloc,
            obsgeovel=to_frame.obsgeovel,
        )
    )

    # now precess to the requested equinox
    pmat = gcrs_precession_mat(to_frame.equinox)
    crepr = gcrs_coo.cartesian.transform(pmat)
    return to_frame.realize_frame(crepr)


@frame_transform_graph.transform(
    FunctionTransformWithFiniteDifference, PrecessedGeocentric, GCRS
)
def precessedgeo_to_gcrs(from_coo, to_frame):
    # first un-precess
    pmat = gcrs_precession_mat(from_coo.equinox)
    crepr = from_coo.cartesian.transform(matrix_transpose(pmat))
    gcrs_coo = GCRS(
        crepr,
        obstime=from_coo.obstime,
        obsgeoloc=from_coo.obsgeoloc,
        obsgeovel=from_coo.obsgeovel,
    )

    # then move to the GCRS that's actually desired
    return gcrs_coo.transform_to(to_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, TEME, ITRS)
def teme_to_itrs(teme_coo, itrs_frame):
    # use the pmatrix to transform to ITRS in the source obstime
    pmat = teme_to_itrs_mat(teme_coo.obstime)
    crepr = teme_coo.cartesian.transform(pmat)
    itrs = ITRS(crepr, obstime=teme_coo.obstime)

    # transform the ITRS coordinate to the target obstime
    return itrs.transform_to(itrs_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ITRS, TEME)
def itrs_to_teme(itrs_coo, teme_frame):
    # transform the ITRS coordinate to the target obstime
    itrs_coo2 = itrs_coo.transform_to(ITRS(obstime=teme_frame.obstime))

    # compute the pmatrix, and then multiply by its transpose
    pmat = teme_to_itrs_mat(teme_frame.obstime)
    newrepr = itrs_coo2.cartesian.transform(matrix_transpose(pmat))
    return teme_frame.realize_frame(newrepr)


# Create loopback transformations
frame_transform_graph._add_merged_transform(ITRS, CIRS, ITRS)
frame_transform_graph._add_merged_transform(
    PrecessedGeocentric, GCRS, PrecessedGeocentric
)
frame_transform_graph._add_merged_transform(TEME, ITRS, TEME)
frame_transform_graph._add_merged_transform(TETE, ICRS, TETE)
