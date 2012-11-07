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
        lambda c: TestCoo2(c.ra.radians, c.dec.radians, unit=u.radian))

    c1 = TestCoo1(1, 2, unit=u.radian)
    c1._make_cart()
    c2 = c1.transform_to(TestCoo2)
    npytest.assert_almost_equal(c2.ra.radians, 1)
    npytest.assert_almost_equal(c2.dec.radians, 2)

    def matfunc(coo):
        return [[1, 0, 0],
                [0, coo.ra.degrees, 0],
                [0, 0, 1]]
    t.DynamicMatrixTransform(TestCoo1, TestCoo2, matfunc)

    c3 = TestCoo1(1, 2, unit=u.degree)
    c3._make_cart()
    c4 = c3.transform_to(TestCoo2)

    npytest.assert_almost_equal(c4.ra.degrees, 1)
    npytest.assert_almost_equal(c4.ra.degrees, 1)


def test_transform_decos():
    """
    Tests the decorator syntax for creating transforms
    """
    c1 = TestCoo1(1, 2, unit=u.degree)

    @t.transform_function(TestCoo1, TestCoo2)
    def trans(coo1):
        return TestCoo2(coo1.ra.radians, coo1.dec.radians * 2, unit=u.radian)

    c1._make_cart()
    c2 = c1.transform_to(TestCoo2)
    npytest.assert_almost_equal(c2.ra.degrees, 1)
    npytest.assert_almost_equal(c2.dec.degrees, 4)

    c3 = TestCoo1(x=1, y=1, z=2, unit=u.pc)

    @t.static_transform_matrix(TestCoo1, TestCoo2)
    def matrix():
        return [[2, 0, 0],
                [0, 1, 0],
                [0, 0, 1]]

    c3._make_cart()
    c4 = c3.transform_to(TestCoo2)

    npytest.assert_almost_equal(c4.x, 2)
    npytest.assert_almost_equal(c4.y, 1)
    npytest.assert_almost_equal(c4.z, 2)

def test_coo_alias():
    """
    Tests the shortname/attribute-style accessing of transforms
    """
    t.coordinate_alias('coo2', TestCoo2)

    t.FunctionTransform(TestCoo1, TestCoo2, lambda c: TestCoo2(c.ra, c.dec))

    c1 = TestCoo1(1, 2, unit=u.degree)
    assert c1.coo2.ra.degrees == c1.ra.degrees
    assert c1.coo2.dec.degrees == c1.dec.degrees

def test_shortest_path():
    class FakeTransform(object):
        def __init__(self, pri):
            self.priority = pri

    g = t.TransformGraph()

    #cheating by adding graph elements directly that are not classes - the
    #graphing algorithm still works fine with integers - it just isn't a valid
    #TransformGraph

    #the graph looks is a down-going diamond graph with the lower-right slightly
    #heavier and a cycle from the bottom to the top
    #also, a pair of nodes isolated from 1

    g._graph[1][2] = FakeTransform(1)
    g._graph[1][3] = FakeTransform(1)
    g._graph[2][4] = FakeTransform(1)
    g._graph[3][4] = FakeTransform(2)
    g._graph[4][1] = FakeTransform(5)

    g._graph[5][6] = FakeTransform(1)

    path, d = g.find_shortest_path(1, 2)
    assert path == [2]
    assert d == 1
    path, d = g.find_shortest_path(1, 3)
    assert path == [3]
    assert d == 1
    path, d = g.find_shortest_path(1, 4)
    print('Cached paths:', g._shortestpaths)
    assert path == [2, 4]
    assert d == 2

    #unreachable
    path, d = g.find_shortest_path(1, 5)
    assert path is None
    assert d == float('inf')

    path, d = g.find_shortest_path(5, 6)
    assert path == [6]
    assert d == 1


def test_sphere_cart():
    """
    Tests the spherical <-> cartesian transform functions
    """
    from ..distances import spherical_to_cartesian, cartesian_to_spherical

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


def test_icrs_gal():
    """
    This tests the ICRS to galactic (and vice versa) conversion
    based on some known-good coordinates.

    Implicitly, this tests the ICRS<->FK5 conversion as well
    """
    raise NotImplementedError


def test_ICRS_FK5():
    """
    This tests the FK5 <-> ICRS conversion
    based on some known-good coordinates
    """
    raise NotImplementedError


def test_ICRS_FK4():
    """
    This tests the FK4 <-> ICRS conversion
    based on some known-good coordinates
    """
    raise NotImplementedError
