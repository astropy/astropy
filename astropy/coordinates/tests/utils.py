# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst


import numpy as np

from ... import units as u
from ...utils import NumpyRNGContext


def randomly_sample_sphere(ntosample, randomseed=12345):
    """
    Generates a set of spherical coordinates uniformly distributed over the
    sphere in a way that gives the same answer for the same seed.  Also
    generates a random distance vector on [0, 1] (no units)

    This simply returns (lon, lat, r) instead of a representation to avoid
    failures due to the representation module.
    """
    with NumpyRNGContext(randomseed):
        lat = np.arcsin(np.random.rand(ntosample)*2-1)
        lon = np.random.rand(ntosample)*np.pi*2
        r = np.random.rand(ntosample)

    return lon*u.rad, lat*u.rad, r
