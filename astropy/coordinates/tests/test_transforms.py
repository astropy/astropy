# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import print_function

import numpy as np
from numpy import testing as npytest


def test_sphere_cart():
    """
    Tests the spherical <-> cartesian transform functions
    """
    from ..transformations import spherical_to_cartesian, cartesian_to_spherical

    x, y ,z = spherical_to_cartesian(1, 0, 0)
    npytest.assert_almost_equal(x, 1)
    npytest.assert_almost_equal(y, 0)
    npytest.assert_almost_equal(z, 0)


    x, y ,z = spherical_to_cartesian(0, 1, 1)
    npytest.assert_almost_equal(x, 0)
    npytest.assert_almost_equal(y, 0)
    npytest.assert_almost_equal(z, 0)

    x, y ,z = spherical_to_cartesian(5, 0, np.arcsin(4. / 5.))
    npytest.assert_almost_equal(x, 3)
    npytest.assert_almost_equal(y, 4)
    npytest.assert_almost_equal(z, 0)

    r, lat, lng = cartesian_to_spherical(0, 1, 0)
    npytest.assert_almost_equal(r, 1)
    npytest.assert_almost_equal(lat, 0)
    npytest.assert_almost_equal(lng, np.pi / 2)


    #test round-tripping
    #TODO: add NumpRNG context after merging w/master
    x, y, z = np.random.randn(3, 5)

    r, lat, lng = cartesian_to_spherical(x, y, z)
    x2, y2, z2 = spherical_to_cartesian(r, lat, lng)

    npytest.assert_allclose(x, x2)
    npytest.assert_allclose(y, y2)
    npytest.assert_allclose(z, z2)