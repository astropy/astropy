# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains the transformation functions for getting to/from ITRS, GCRS, and CIRS.
These are distinct from the ICRS and AltAz functions because they are just
rotations without aberration corrections or offsets.
"""

import numpy as np

from ..baseframe import frame_transform_graph
from ..transformations import FunctionTransformWithFiniteDifference
from ..matrix_utilities import matrix_transpose
from ... import _erfa as erfa

from .gcrs import GCRS, PrecessedGeocentric
from .cirs import CIRS
from .itrs import ITRS
from .utils import get_polar_motion, get_jd12

# # first define helper functions


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


def gcrs_precession_mat(equinox):
    gamb, phib, psib, epsa = erfa.pfw06(*get_jd12(equinox, 'tt'))
    return erfa.fw2m(gamb, phib, psib, epsa)


# now the actual transforms

@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, GCRS, CIRS)
def gcrs_to_cirs(gcrs_coo, cirs_frame):
    # first get us to a 0 pos/vel GCRS at the target obstime
    gcrs_coo2 = gcrs_coo.transform_to(GCRS(obstime=cirs_frame.obstime))

    # now get the pmatrix
    pmat = gcrs_to_cirs_mat(cirs_frame.obstime)
    crepr = gcrs_coo2.cartesian.transform(pmat)
    return cirs_frame.realize_frame(crepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, CIRS, GCRS)
def cirs_to_gcrs(cirs_coo, gcrs_frame):
    # compute the pmatrix, and then multiply by its transpose
    pmat = gcrs_to_cirs_mat(cirs_coo.obstime)
    newrepr = cirs_coo.cartesian.transform(matrix_transpose(pmat))
    gcrs = GCRS(newrepr, obstime=cirs_coo.obstime)

    # now do any needed offsets (no-op if same obstime and 0 pos/vel)
    return gcrs.transform_to(gcrs_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, CIRS, ITRS)
def cirs_to_itrs(cirs_coo, itrs_frame):
    # first get us to CIRS at the target obstime
    cirs_coo2 = cirs_coo.transform_to(CIRS(obstime=itrs_frame.obstime))

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
    return from_coo.transform_to(CIRS).transform_to(to_frame)

# TODO: implement GCRS<->CIRS if there's call for it.  The thing that's awkward
# is that they both have obstimes, so an extra set of transformations are necessary.
# so unless there's a specific need for that, better to just have it go through the above
# two steps anyway


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
    gcrs_coo = GCRS(crepr, obstime=to_frame.obstime,
                           obsgeoloc=to_frame.obsgeoloc,
                           obsgeovel=to_frame.obsgeovel)

    # then move to the GCRS that's actually desired
    return gcrs_coo.transform_to(to_frame)
