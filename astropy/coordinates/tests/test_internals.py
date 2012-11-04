# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import print_function

import numpy as np
from ...tests.helper import pytest

from .. import angle_utilities

distance_funcs = {'small_angle': angle_utilities.small_angle_dist,
                  'sphere': angle_utilities.sphere_dist,
                  'haversine': angle_utilities.haversine_dist,
                  'haversine_atan': angle_utilities.haversine_dist_atan,
                  'vicenty': angle_utilities.vicenty_dist,
                 }
# lat1, long1, lat2, long2 in degrees
coords = [(0, 0, 1, 0),
          (0, 0, 0, 10),
          (0, 0, 0, 90),
          (0, 0, 0, 180),
          (45, 0, -45, 0),
          (60, 0, -30, 0),
          (0, 0, 0, 0),
          (0, 0, 1. / 60., 1. / 60.)
         ]
correct_seps = [1, 10, 90, 180, 90, 90, 0, 0.023570225877234643]
correctness_margin = 1e-98


@pytest.mark.parametrize(('coords', 'corrsep'), zip(coords, correct_seps))
def test_2dseparations(coords, corrsep):
    """
    A variety of tests to examine how close the various distance estimators are
    """
    from math import fabs, radians, degrees

    lat1, long1, lat2, long2 = coords

    print('({0},{1}) - ({2},{3})'.format(lat1, long1, lat2, long2))
    print('Correct separation', corrsep)

    seps = {}
    for n, dfunc in distance_funcs.iteritems():
        sep = degrees(dfunc(radians(lat1), radians(long1), radians(lat2), radians(long2)))
        seps[n] = sep - corrsep
    print(seps)
    assert all([fabs(sepval) < 1e-9 for sepval in seps.values()])

