# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

"""
Tests for the projected separation stuff
"""

import numpy as np

from ...tests.helper import pytest, assert_quantity_allclose as assert_allclose
from ...extern.six.moves import zip
from ... import units as u
from ..builtin_frames import ICRS, FK5, Galactic
from .. import Angle, Distance

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


def test_angsep():
    """
    Tests that the angular separation object also behaves correctly.
    """
    from ..angle_utilities import angular_separation

    # check it both works with floats in radians, Quantities, or Angles
    for conv in (np.deg2rad,
                 lambda x: u.Quantity(x, "deg"),
                 lambda x: Angle(x, "deg")):
        for (lon1, lat1, lon2, lat2), corrsep in zip(coords, correct_seps):
            angsep = angular_separation(conv(lon1), conv(lat1),
                                        conv(lon2), conv(lat2))
            assert np.fabs(angsep - conv(corrsep)) < conv(correctness_margin)


def test_fk5_seps():
    """
    This tests if `separation` works for FK5 objects.

    This is a regression test for github issue #891
    """
    a = FK5(1.*u.deg, 1.*u.deg)
    b = FK5(2.*u.deg, 2.*u.deg)
    a.separation(b)


def test_proj_separations():
    """
    Test angular separation functionality
    """
    c1 = ICRS(ra=0*u.deg, dec=0*u.deg)
    c2 = ICRS(ra=0*u.deg, dec=1*u.deg)

    sep = c2.separation(c1)
    #returns an Angle object
    assert isinstance(sep, Angle)

    assert sep.degree == 1
    assert_allclose(sep.arcminute, 60.)

    # these operations have ambiguous interpretations for points on a sphere
    with pytest.raises(TypeError):
        c1 + c2
    with pytest.raises(TypeError):
        c1 - c2

    ngp = Galactic(l=0*u.degree, b=90*u.degree)
    ncp = ICRS(ra=0*u.degree, dec=90*u.degree)

    # if there is a defined conversion between the relevant coordinate systems,
    # it will be automatically performed to get the right angular separation
    assert_allclose(ncp.separation(ngp.transform_to(ICRS)).degree,
                    ncp.separation(ngp).degree)

    # distance from the north galactic pole to celestial pole
    assert_allclose(ncp.separation(ngp.transform_to(ICRS)).degree,
                    62.87174758503201)


def test_3d_separations():
    """
    Test 3D separation functionality
    """
    c1 = ICRS(ra=1*u.deg, dec=1*u.deg, distance=9*u.kpc)
    c2 = ICRS(ra=1*u.deg, dec=1*u.deg, distance=10*u.kpc)

    sep3d = c2.separation_3d(c1)

    assert isinstance(sep3d, Distance)
    assert_allclose(sep3d - 1*u.kpc, 0*u.kpc, atol=1e-12*u.kpc)
