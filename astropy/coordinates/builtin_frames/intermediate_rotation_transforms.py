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

from ... import units as u
from ..baseframe import frame_transform_graph
from ..transformations import FunctionTransform, DynamicMatrixTransform
from ..representation import UnitSphericalRepresentation, CartesianRepresentation
from ... import erfa

from .gcrs import GCRS
from .cirs import CIRS
from .itrs import ITRS
from .utils import get_polar_motion


@frame_transform_graph.transform(DynamicMatrixTransform, GCRS, ITRS)
def gcrs_to_itrs(gcrs_coo, itrs_frame):
    #first compute the celestial-to-intermediate matrix
    c2imat = erfa.c2i06a(gcrs_coo.obstime.jd1, gcrs_coo.obstime.jd2)

    #now compute the polar motion p-matrix
    xp, yp = get_polar_motion(gcrs_coo.obstime)
    sp = erfa.sp00(gcrs_coo.obstime.jd1, gcrs_coo.obstime.jd2)
    pmmat = erfa.pom00(xp, yp, sp)

    #now determine the Earth Rotation Angle for the input obstime
    era = erfa.era00(gcrs_coo.obstime.jd1, gcrs_coo.obstime.jd2)

    return erfa.c2tcio(c2imat, era, pmmat)


@frame_transform_graph.transform(DynamicMatrixTransform, ITRS, GCRS)
def itrs_to_gcrs(itrs_coo, gcrs_frame):
    return gcrs_to_itrs(gcrs_frame, itrs_coo).T


@frame_transform_graph.transform(DynamicMatrixTransform, CIRS, ITRS)
def cirs_to_itrs(cirs_coo, itrs_frame):
    #compute the polar motion p-matrix
    xp, yp = get_polar_motion(cirs_coo.obstime)
    sp = erfa.sp00(cirs_coo.obstime.jd1, cirs_coo.obstime.jd2)
    pmmat = erfa.pom00(xp, yp, sp)

    #now determine the Earth Rotation Angle for the input obstime
    era = erfa.era00(cirs_coo.obstime.jd1, cirs_coo.obstime.jd2)

    #c2tcio expects a GCRS->CIRS matrix, but we just set that to an I-matrix
    return erfa.c2tcio(np.eye(3), era, pmmat)


@frame_transform_graph.transform(DynamicMatrixTransform, ITRS, CIRS)
def itrs_to_cirs(itrs_coo, cirs_frame):
    return cirs_to_itrs(cirs_frame, itrs_coo).T

#TODO: implement GCRS<->CIRS if there's call for it.  The thing that's awkward
#is that they both have obstimes, so an extra set of transformations are necessary.
#so unless there's a specific need for that, better to just have it go through te above
#two steps anyway
