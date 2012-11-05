# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import print_function

import numpy as np
from numpy import testing as npytest

from ... import units as u
from .. import transformations as t
from ..builtin_systems import ICRSCoordinates
from ...tests.helper import pytest


#Coordinates just for these tests. TODO: expunge
class TestCoo1(ICRSCoordinates):
    pass


class TestCoo2(ICRSCoordinates):
    pass



def test_transform_classes():
    """
    Tests the class-based/OO syntax for creating transforms
    """
    t.FunctionTransform(TestCoo1, TestCoo2,
        lambda c: TestCoo2(c.ra.r, c.dec.r, unit=u.radian))

    c1 = TestCoo1(1, 2, unit=u.radian)
    c2 = c1.transform_to(TestCoo2)
    assert c2.ra.r == 1
    assert c2.dec.r == 2

    def matfunc(coo):
        return [[1, 0, 0],
                [0, coo.ra.d, 0],
                [0, 0, 1]]
    t.DynamicMatrixTransform(TestCoo1, TestCoo2, matfunc)

    c3 = TestCoo1(1, 2, unit=u.degree)
    c4 = c3.transform_to(TestCoo2)

    assert c4.ra.radian == 1
    assert c4.dec.radian == 2


def test_transform_decos():
    """
    Tests the decorator syntax for creating transforms
    """
    c1 = TestCoo1(1, 2, unit=u.degree)

    @t.transform_function(TestCoo1, TestCoo2)
    def trans(coo1):
        return TestCoo2(coo1.ra.radians, coo1.dec.radians * 2, unit=u.radian)

    c2 = c1.transform_to(TestCoo2)
    assert c2.ra.d == 1
    assert c2.dec.d == 4

    c3 = TestCoo1(x=1, y=1, z=2, unit=u.pc)

    @t.static_transform_matrix(TestCoo1, TestCoo2)
    def matrix():
        return [[2, 0, 0],
                [0, 1, 0],
                [0, 0, 1]]

    c4 = c3.transform_to(TestCoo2)

    assert c4.x == 2
    assert c4.y == 1
    assert c4.z == 2

def test_coo_alias():
    """
    Tests the shortname/attribute-style accessing of transforms
    """
    t.coordinate_alias('coo2', TestCoo2)

    t.FunctionTransform(TestCoo1, TestCoo2, lambda c: TestCoo2(c.ra, c.dec))

    c1 = TestCoo1(1, 2, unit=u.degree)
    assert c1.coo2.ra.d == c1.ra.d
    assert c1.coo2.dec.d == c1.dec.d


def test_sphere_cart():
    """
    Tests the spherical <-> cartesian transform functions
    """
    from ..coordsystems import spherical_to_cartesian, cartesian_to_spherical

    x, y, z = spherical_to_cartesian(1, 0, 0)
    npytest.assert_almost_equal(x, 1)
    npytest.assert_almost_equal(y, 0)
    npytest.assert_almost_equal(z, 0)

    x, y, z = spherical_to_cartesian(0, 1, 1)
    npytest.assert_almost_equal(x, 0)
    npytest.assert_almost_equal(y, 0)
    npytest.assert_almost_equal(z, 0)

    x, y, z = spherical_to_cartesian(5, 0, np.arcsin(4. / 5.))
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
