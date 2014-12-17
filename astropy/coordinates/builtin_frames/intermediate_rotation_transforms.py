# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains the transformation functions for getting to/from ITRS, GCRS, and CIRS.
These are distinct from the ICRS and AltAz functions because they are just
rotations without aberration corrections or offsets.
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np

from ..baseframe import frame_transform_graph
from ..transformations import FunctionTransform
from ..representation import CartesianRepresentation
from ... import erfa

from .gcrs import GCRS
from .cirs import CIRS
from .itrs import ITRS
from .utils import get_polar_motion

# first define helper functions
def gcrs_to_itrs_mat(time):
    #first compute the celestial-to-intermediate matrix
    c2imat = erfa.c2i06a(time.jd1, time.jd2)

    #now compute the polar motion p-matrix
    xp, yp = get_polar_motion(time)
    sp = erfa.sp00(time.jd1, time.jd2)
    pmmat = erfa.pom00(xp, yp, sp)

    #now determine the Earth Rotation Angle for the input obstime
    era = erfa.era00(time.jd1, time.jd2)

    return erfa.c2tcio(c2imat, era, pmmat)


def cirs_to_itrs_mat(time):
    #compute the polar motion p-matrix
    xp, yp = get_polar_motion(time)
    sp = erfa.sp00(time.jd1, time.jd2)
    pmmat = erfa.pom00(xp, yp, sp)

    #now determine the Earth Rotation Angle for the input obstime
    era = erfa.era00(time.jd1, time.jd2)

    #c2tcio expects a GCRS->CIRS matrix, but we just set that to an I-matrix
    #because we're already in CIRS
    return erfa.c2tcio(np.eye(3), era, pmmat)


def cartrepr_from_matmul(pmat, coo):
    xyz = coo.cartesian.xyz
    newxyz = np.dot(pmat, xyz.reshape(3,-1)).reshape(xyz.shape)
    return CartesianRepresentation(newxyz)


# now the actual transforms

@frame_transform_graph.transform(FunctionTransform, GCRS, ITRS)
def gcrs_to_itrs(gcrs_coo, itrs_frame):
    # first get us to a 0 pos/vel GCRS at the target obstime
    gcrs_coo2 = gcrs_coo.transform_to(GCRS(obstime=itrs_frame.obstime))

    #now get the pmatrix
    pmat = gcrs_to_itrs_mat(itrs_frame.obstime)
    crepr = cartrepr_from_matmul(pmat, gcrs_coo2)
    return itrs_frame.realize_frame(crepr)

@frame_transform_graph.transform(FunctionTransform, ITRS, GCRS)
def itrs_to_gcrs(itrs_coo, gcrs_frame):
    #compute the pmatrix, and then multiply by its transpose
    pmat = gcrs_to_itrs_mat(itrs_coo.obstime)
    newrepr = cartrepr_from_matmul(pmat.T, itrs_coo)
    gcrs = GCRS(newrepr, obstime=itrs_coo.obstime)

    #now do any needed offsets (no-op if same obstime and 0 pos/vel)
    return gcrs.transform_to(gcrs_frame)



@frame_transform_graph.transform(FunctionTransform, CIRS, ITRS)
def cirs_to_itrs(cirs_coo, itrs_frame):
    # first get us to CIRS at the target obstime
    cirs_coo2 = cirs_coo.transform_to(CIRS(obstime=itrs_frame.obstime))

    #now get the pmatrix
    pmat = cirs_to_itrs_mat(itrs_frame.obstime)
    crepr = cartrepr_from_matmul(pmat, cirs_coo2)
    return itrs_frame.realize_frame(crepr)


@frame_transform_graph.transform(FunctionTransform, ITRS, CIRS)
def itrs_to_cirs(itrs_coo, cirs_frame):
    #compute the pmatrix, and then multiply by its transpose
    pmat = cirs_to_itrs_mat(itrs_coo.obstime)
    newrepr = cartrepr_from_matmul(pmat.T, itrs_coo)
    cirs = CIRS(newrepr, obstime=itrs_coo.obstime)

    #now do any needed offsets (no-op if same obstime)
    return cirs.transform_to(cirs_frame)


#TODO: implement GCRS<->CIRS if there's call for it.  The thing that's awkward
#is that they both have obstimes, so an extra set of transformations are necessary.
#so unless there's a specific need for that, better to just have it go through te above
#two steps anyway
