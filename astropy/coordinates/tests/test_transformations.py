# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import print_function, division

import numpy as np
from numpy import testing as npt

from ... import units as u
from ..distances import Distance
from .. import transformations as t
from ..builtin_systems import ICRSCoordinates, FK5Coordinates, FK4Coordinates, FK4NoETermCoordinates
from ..builtin_systems import GalacticCoordinates
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
        lambda c: TestCoo2(c.ra.radians, c.dec.radians, unit=(u.radian, u.radian)))

    c1 = TestCoo1(1, 2, unit=(u.radian, u.radian))
    c1._make_cart()
    c2 = c1.transform_to(TestCoo2)
    npt.assert_almost_equal(c2.ra.radians, 1)
    npt.assert_almost_equal(c2.dec.radians, 2)

    def matfunc(coo):
        return [[1, 0, 0],
                [0, coo.ra.degrees, 0],
                [0, 0, 1]]
    t.DynamicMatrixTransform(TestCoo1, TestCoo2, matfunc)

    c3 = TestCoo1(1, 2, unit=(u.degree, u.degree))
    c3._make_cart()
    c4 = c3.transform_to(TestCoo2)

    npt.assert_almost_equal(c4.ra.degrees, 1)
    npt.assert_almost_equal(c4.ra.degrees, 1)


def test_transform_decos():
    """
    Tests the decorator syntax for creating transforms
    """
    c1 = TestCoo1(1, 2, unit=(u.degree, u.degree))

    @t.transform_function(TestCoo1, TestCoo2)
    def trans(coo1):
        return TestCoo2(coo1.ra.radians, coo1.dec.radians * 2, unit=(u.radian, u.radian))

    c1._make_cart()
    c2 = c1.transform_to(TestCoo2)
    npt.assert_almost_equal(c2.ra.degrees, 1)
    npt.assert_almost_equal(c2.dec.degrees, 4)

    c3 = TestCoo1(x=1, y=1, z=2, unit=u.pc)

    @t.static_transform_matrix(TestCoo1, TestCoo2)
    def matrix():
        return [[2, 0, 0],
                [0, 1, 0],
                [0, 0, 1]]

    c3._make_cart()
    c4 = c3.transform_to(TestCoo2)

    npt.assert_almost_equal(c4.x, 2)
    npt.assert_almost_equal(c4.y, 1)
    npt.assert_almost_equal(c4.z, 2)

def test_coo_alias():
    """
    Tests the shortname/attribute-style accessing of transforms
    """
    t.coordinate_alias('coo2', TestCoo2)

    t.FunctionTransform(TestCoo1, TestCoo2, lambda c: TestCoo2(c.ra, c.dec))

    c1 = TestCoo1(1, 2, unit=(u.degree, u.degree))
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
    assert path == [1, 2]
    assert d == 1
    path, d = g.find_shortest_path(1, 3)
    assert path == [1, 3]
    assert d == 1
    path, d = g.find_shortest_path(1, 4)
    print('Cached paths:', g._shortestpaths)
    assert path == [1, 2, 4]
    assert d == 2

    #unreachable
    path, d = g.find_shortest_path(1, 5)
    assert path is None
    assert d == float('inf')

    path, d = g.find_shortest_path(5, 6)
    assert path == [5, 6]
    assert d == 1


def test_sphere_cart():
    """
    Tests the spherical <-> cartesian transform functions
    """
    from ...tests.compat import assert_allclose
    from ...utils import NumpyRNGContext
    from ..distances import spherical_to_cartesian, cartesian_to_spherical


    x, y, z = spherical_to_cartesian(1, 0, 0)
    npt.assert_almost_equal(x, 1)
    npt.assert_almost_equal(y, 0)
    npt.assert_almost_equal(z, 0)

    x, y, z = spherical_to_cartesian(0, 1, 1)
    npt.assert_almost_equal(x, 0)
    npt.assert_almost_equal(y, 0)
    npt.assert_almost_equal(z, 0)

    x, y, z = spherical_to_cartesian(5, 0, np.arcsin(4. / 5.))
    npt.assert_almost_equal(x, 3)
    npt.assert_almost_equal(y, 4)
    npt.assert_almost_equal(z, 0)

    r, lat, lon = cartesian_to_spherical(0, 1, 0)
    npt.assert_almost_equal(r, 1)
    npt.assert_almost_equal(lat, 0)
    npt.assert_almost_equal(lon, np.pi / 2)

    #test round-tripping
    with NumpyRNGContext(13579):
        x, y, z = np.random.randn(3, 5)

    r, lat, lon = cartesian_to_spherical(x, y, z)
    x2, y2, z2 = spherical_to_cartesian(r, lat, lon)

    assert_allclose(x, x2)
    assert_allclose(y, y2)
    assert_allclose(z, z2)


m31_sys = [(ICRSCoordinates, 'icrs'), (FK5Coordinates, 'fk5'), (FK4Coordinates, 'fk4'), (GalacticCoordinates, 'galactic')]
m31_coo = [(10.6847929, 41.2690650), (10.6847929, 41.2690650), (10.0004738, 40.9952444), (121.1744050, -21.5729360)]
m31_dist = Distance(770, u.kpc)
convert_precision = 1 / 3600.  # 1 arcsec
roundtrip_precision = 1e-4
dist_precision = 1e-9

m31_params =[]
for i in range(len(m31_sys)):
    for j in range(len(m31_sys)):
        if i < j:
            m31_params.append((m31_sys[i], m31_sys[j], m31_coo[i], m31_coo[j]))

@pytest.mark.parametrize(('fromsys', 'tosys', 'fromcoo', 'tocoo'), m31_params)
def test_m31_coord_transforms(fromsys, tosys, fromcoo, tocoo):
    """
    This tests a variety of coordinate conversions for the Chandra point-source
    catalog location of M31 from NED.
    """
    from math import fabs

    from ...time import Time

    coo1 = fromsys[0](fromcoo[0], fromcoo[1], unit=(u.degree, u.degree), distance=m31_dist)
    coo2 = coo1.transform_to(tosys[0])
    if tosys[0] is FK4Coordinates:
        coo2_prec = coo2.precess_to(Time('B1950', scale='utc'))
        assert fabs(coo2_prec.lonangle.degrees - tocoo[0]) < convert_precision  # <1 arcsec
        assert fabs(coo2_prec.latangle.degrees - tocoo[1]) < convert_precision
    else:
        assert fabs(coo2.lonangle.degrees - tocoo[0]) < convert_precision  # <1 arcsec
        assert fabs(coo2.latangle.degrees - tocoo[1]) < convert_precision
    assert fabs(coo2.distance.kpc - m31_dist.kpc) < dist_precision

    if fromsys[1] is not None:
        coo1_2 = getattr(coo2, fromsys[1])  # implicit `transform_to` call.

        #check round-tripping
        assert fabs(coo1_2.lonangle.degrees - fromcoo[0]) < roundtrip_precision
        assert fabs(coo1_2.latangle.degrees - fromcoo[1]) < roundtrip_precision
        assert fabs(coo1_2.distance.kpc - m31_dist.kpc) < dist_precision

        if tosys[1] is not None:
            coo2_2 = getattr(coo1_2, tosys[1])
            assert fabs(coo2_2.lonangle.degrees - coo2.lonangle.degrees) < roundtrip_precision
            assert fabs(coo2_2.latangle.degrees - coo2.latangle.degrees) < roundtrip_precision
            assert fabs(coo2_2.distance.kpc - m31_dist.kpc) < dist_precision


def test_precession():
    """
    Ensures that FK4 and FK5 coordinates precess their equinoxes
    """
    from ...time import Time

    j2000 = Time('J2000', scale='utc')
    b1950 = Time('B1950', scale='utc')
    j1975 = Time('J1975', scale='utc')
    b1975 = Time('B1975', scale='utc')

    fk4 = FK4Coordinates(1, 2, unit=(u.radian, u.radian))
    assert fk4.equinox.byear == b1950.byear
    fk4_2 = fk4.precess_to(b1975)
    assert fk4_2.equinox.byear == b1975.byear

    fk5 = FK5Coordinates(1, 2, unit=(u.radian, u.radian))
    assert fk5.equinox.jyear == j2000.jyear
    fk5_2 = fk5.precess_to(j1975)
    assert fk5_2.equinox.jyear == j1975.jyear


def test_alias_transform():
    """
    Tests the use of aliases to do trasnforms and also transforming from
    a system to itself.  Also checks that `dir` correctly includes
    valid transforms
    """
    c = ICRSCoordinates(12.34, 56.78, unit=(u.hour, u.degree))
    assert isinstance(c.galactic, GalacticCoordinates)
    assert isinstance(c.icrs, ICRSCoordinates)

    d = dir(c)
    assert 'galactic' in d
    assert 'fk4' in d
    assert 'fk5' in d


def test_transform_path_pri():
    """
    This checks that the transformation path prioritization works by
    making sure the ICRS -> Gal transformation always goes through FK5
    and not FK4.
    """
    t.master_transform_graph.invalidate_cache()
    tpath, td = t.master_transform_graph.find_shortest_path(ICRSCoordinates, GalacticCoordinates)
    assert tpath == [ICRSCoordinates, FK5Coordinates, GalacticCoordinates]
    assert td == 2

    #but direct from FK4 to Galactic should still be possible
    tpath, td = t.master_transform_graph.find_shortest_path(FK4Coordinates, GalacticCoordinates)
    assert tpath == [FK4Coordinates, FK4NoETermCoordinates, GalacticCoordinates]
    assert td == 1.01 + 1.02


def test_obstime():
    """
    Checks to make sure observation time survives transforms, and that it's
    accounted for at least in FK4 <-> ICRS transformations
    """
    from ...time import Time

    b1950 = Time('B1950', scale='utc')
    j1975 = Time('J1975', scale='utc')

    fk4_50 = FK4Coordinates(1, 2, unit=(u.degree, u.degree), obstime=b1950)
    fk4_75 = FK4Coordinates(1, 2, unit=(u.degree, u.degree), obstime=j1975)

    icrs_50 = fk4_50.transform_to(ICRSCoordinates)
    icrs_75 = fk4_75.transform_to(ICRSCoordinates)

    assert icrs_50.obstime == fk4_50.obstime
    assert icrs_75.obstime == fk4_75.obstime

    # now check that the resulting coordinates are *different* - they should be,
    # because the obstime is different
    assert (icrs_50.ra.degrees != icrs_75.ra.degrees and
            icrs_50.dec.degrees != icrs_75.dec.degrees)
