# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst


import numpy as np

from astropy import units as u
from astropy.utils import NumpyRNGContext
from astropy.utils.decorators import deprecated


# TODO: remove this function in v5.0. I think we can have a fairly fast
# deprecation cycle here because it is not meant to be public API.
@deprecated(since='v4.3',
            message='This function has been deprecated in favor of the '
                    'public-facing utilities in '
                    'astropy.coordinates.angle_utilities',
            alternative='Use uniform_spherical_random_surface() from '
                        'astropy.coordinates.angle_utilities instead.')
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
