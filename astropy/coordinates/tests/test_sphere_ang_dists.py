# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import print_function

from ...tests.helper import pytest

from .. import angle_utilities

distance_funcs = [angle_utilities.small_angle_sphere_dist,
                  angle_utilities.simple_sphere_dist,
                  angle_utilities.haversine_sphere_dist,
                  angle_utilities.haversine_atan_sphere_dist,
                  angle_utilities.vincenty_sphere_dist,
                 ]

# lat1, lon1, lat2, lon2 in degrees
coords = [(0, 0, 1, 0),
          (0, 0, 0, 10),
          (0, 0, 0, 90),
          (0, 0, 0, 180),
          (45, 0, -45, 0),
          (60, 0, -30, 0),
          (-15, -135, 15, 45),
          (-89, 100, 89, -80),
          (0, 0, 0, 0),
          (0, 0, 1. / 60., 1. / 60.)
         ]
correct_seps = [1, 10, 90, 180, 90, 90, 180, 180, 0, 0.023570225877234643]
correctness_margin = 2e-10

# set this to a numer to run the timing tests and trigger a failure so stdout is
# read back in pytest - the number gives the number of iterations
dotiming = False

paramsets = []
for coord, coorsep in zip(coords, correct_seps):
    for df in distance_funcs:
        paramsets.append((coord, coorsep, df))


@pytest.mark.parametrize(('coord', 'correctsep', 'dfunc'), paramsets)
def test_2dseparations(coord, correctsep, dfunc):
    """
    A variety of tests to examine how close the various sphereical
    distance/great circle measurements are from the expectations
    """
    from time import time
    from math import fabs, radians, degrees

    lat1, lon1, lat2, lon2 = coord

    print('distance function', dfunc)
    print('({0},{1}) - ({2},{3})'.format(lon1, lat1, lon2, lat2))
    print('Correct separation', correctsep)

    inputs = (radians(lon1), radians(lat1), radians(lon2), radians(lat2))
    sep = degrees(dfunc(*inputs))
    print('Reported:', sep)
    print('Deviation from correct:', sep - correctsep)

    if dotiming:
        starttime = time()
        for i in range(int(dotiming)):
            dfunc(*inputs)
        endtime = time()
        dt = endtime - starttime
        print('{0} in {1} sec'.format(int(dotiming), dt))
        print('{0} sec per execution'.format(dt / float(int(dotiming))))
        assert False  # Triggers py.test failures so we can see the timings

    #a few cases are known to fail because of bad approximations - let them fail
    if dfunc is angle_utilities.small_angle_sphere_dist:
        if fabs(lat2 - lat1) > 1 and fabs(lon2 - lon1) > 1:  # radians
            pytest.xfail('Small angle approximation fails for large angles')

    assert fabs(sep - correctsep) < correctness_margin


def test_fk5_seps():
    """
    This tests if `separation` works for FK5Coordinate objects.

    This is a regression test for github issue #891
    """
    from astropy.coordinates import FK5Coordinates

    a = FK5Coordinates(1., 1., unit=('deg', 'deg'))
    b = FK5Coordinates(2., 2., unit=('deg', 'deg'))
    a.separation(b)

def test_angsep():
    """
    Tests that the angular separation object also behaves correctly.
    """
    from math import fabs

    from ... import units as u
    from ..angles import AngularSeparation

    for (lat1, lon1, lat2, lon2), corrsep in zip(coords, correct_seps):
        angsep = AngularSeparation(lon1, lat1, lon2, lat2, u.deg)
        assert fabs(angsep.degrees - corrsep) < correctness_margin
