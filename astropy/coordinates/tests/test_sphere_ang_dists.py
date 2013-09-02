# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import print_function

# lon1, lat1, lon2, lat2 in degrees
coords = [(1, 0, 0, 0),
          (0, 1, 0, 0),
          (0, 0, 1, 0),
          (0, 0, 0, 1),
          (0, 0, 10, 0),
          (0, 0, 90, 0),
          (0, 0, 180, 0),
          (0, 45, 0, -45),
          (0, 60, 0, -30),
          (-135, -15, 45, 15),
          (100, -89, -80, 89),
          (0, 0, 0, 0),
          (0, 0, 1. / 60., 1. / 60.)]
correct_seps = [1, 1, 1, 1, 10, 90, 180, 90, 90, 180, 180, 0,
                0.023570225877234643]
correctness_margin = 2e-10


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
    from ..angle_utilities import angular_separation

    for (lon1, lat1, lon2, lat2), corrsep in zip(coords, correct_seps):
        angsep = angular_separation(lon1, lat1, lon2, lat2)
        assert fabs(angsep.degree - corrsep) < correctness_margin
