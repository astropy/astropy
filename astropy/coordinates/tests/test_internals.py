# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import print_function

from pytest import raises

from ... import units as u


def test_2dseparations():
    """
    A variety of tests to examine how close the various distance estimators are
    """
    from math import fabs

    from .. import angles

    distances = {'small_angle': angles.AngularSeparation._small_angle_dist,
                 'sphere': angles.AngularSeparation._sphere_dist,
                 'haversine': angles.AngularSeparation._haversine_dist,
                 'haversine_atan': angles.AngularSeparation._haversine_dist_atan,
                 'vicenty': angles.AngularSeparation._vicenty_dist,
                }

             # lat1, long1, lat2, long2
    coords = [(0, 0, 1, 0),
              (0, 0, 0, 10),
              (0, 0, 0, 90),
              (0, 0, 0, 180),
              (45, 0, -45, 0),
              (0, 0, 0, 0),
              (0, 0, 1 / 60., .1 / 60.)
             ]
    coorect_seps = [1, 10, 90, 180, 90, 0, 2 ** 0.5 / 60.]

    for (lat1, long1, lat2, long2), corrsep in zip(coords, coorect_seps):
        print('({0},{1}) - ({2},{3})'.format(lat1, long1, lat2, long2))
        print('Correct separation', corrsep)
        for n, dfunc in distances.iteritems():
            sep = dfunc(lat1, long1, lat2, long2)
            print(n, sep - corrsep)
            assert fabs(sep - corrsep) < 1e-10
