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
from ... import _erfa as erfa

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


def cartrepr_from_matmul(pmat, coo, transpose=False):
    if pmat.shape[-2:] != (3, 3):
        raise ValueError("tried to do matrix multiplication with an array that "
                         "doesn't end in 3x3")
    xyz = coo.cartesian.xyz.T
    # these expression are the same as iterating over the first dimension of
    # pmat and xyz and doing matrix multiplication on each in turn.  resulting
    # dimension is <coo shape> x 3
    pmat = pmat.reshape(pmat.size//9, 3, 3)
    if transpose:
        pmat = pmat.transpose(0, 2, 1)
    newxyz = np.sum(pmat * xyz.reshape(xyz.size//3, 1, 3), axis=-1)

    return CartesianRepresentation(newxyz.T)


# now the actual transforms

# the priority for the GCRS<->ITRS trasnsforms are higher (=less traveled) to
#make GCRS<->ICRS<->CIRS the preferred route over GCRS<->ITRS<->CIRS
@frame_transform_graph.transform(FunctionTransform, GCRS, ITRS, priority=1.01)
def gcrs_to_itrs(gcrs_coo, itrs_frame):
    # first get us to a 0 pos/vel GCRS at the target obstime
    gcrs_coo2 = gcrs_coo.transform_to(GCRS(obstime=itrs_frame.obstime))

    #now get the pmatrix
    pmat = gcrs_to_itrs_mat(itrs_frame.obstime)
    crepr = cartrepr_from_matmul(pmat, gcrs_coo2)
    return itrs_frame.realize_frame(crepr)


@frame_transform_graph.transform(FunctionTransform, ITRS, GCRS, priority=1.01)
def itrs_to_gcrs(itrs_coo, gcrs_frame):
    #compute the pmatrix, and then multiply by its transpose
    pmat = gcrs_to_itrs_mat(itrs_coo.obstime)
    newrepr = cartrepr_from_matmul(pmat, itrs_coo, transpose=True)
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
    newrepr = cartrepr_from_matmul(pmat, itrs_coo, transpose=True)
    cirs = CIRS(newrepr, obstime=itrs_coo.obstime)

    #now do any needed offsets (no-op if same obstime)
    return cirs.transform_to(cirs_frame)


@frame_transform_graph.transform(FunctionTransform, ITRS, ITRS)
def itrs_to_itrs(from_coo, to_frame):
    # this self-transform goes through CIRS right now, which implicitly also
    # goes back to ICRS
    return from_coo.transform_to(CIRS).transform_to(to_frame)



#TODO: implement GCRS<->CIRS if there's call for it.  The thing that's awkward
#is that they both have obstimes, so an extra set of transformations are necessary.
#so unless there's a specific need for that, better to just have it go through the above
#two steps anyway
