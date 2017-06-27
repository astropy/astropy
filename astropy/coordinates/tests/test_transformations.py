# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
import pytest

from ... import units as u
from .. import transformations as t
from ..builtin_frames import ICRS, FK5, FK4, FK4NoETerms, Galactic, AltAz
from .. import representation as r
from ..baseframe import frame_transform_graph
from ...tests.helper import (assert_quantity_allclose as assert_allclose,
                             quantity_allclose, catch_warnings)
from ...time import Time


# Coordinates just for these tests.
class TCoo1(ICRS):
    pass


class TCoo2(ICRS):
    pass


class TCoo3(ICRS):
    pass


def test_transform_classes():
    """
    Tests the class-based/OO syntax for creating transforms
    """

    tfun = lambda c, f: f.__class__(ra=c.ra, dec=c.dec)
    trans1 = t.FunctionTransform(tfun, TCoo1, TCoo2,
                        register_graph=frame_transform_graph)

    c1 = TCoo1(ra=1*u.radian, dec=0.5*u.radian)
    c2 = c1.transform_to(TCoo2)
    assert_allclose(c2.ra.radian, 1)
    assert_allclose(c2.dec.radian, 0.5)

    def matfunc(coo, fr):
        return [[1, 0, 0],
                [0, coo.ra.degree, 0],
                [0, 0, 1]]
    trans2 = t.DynamicMatrixTransform(matfunc, TCoo1, TCoo2)
    trans2.register(frame_transform_graph)

    c3 = TCoo1(ra=1*u.deg, dec=2*u.deg)
    c4 = c3.transform_to(TCoo2)

    assert_allclose(c4.ra.degree, 1)
    assert_allclose(c4.ra.degree, 1)

    # be sure to unregister the second one - no need for trans1 because it
    # already got unregistered when trans2 was created.
    trans2.unregister(frame_transform_graph)


def test_transform_decos():
    """
    Tests the decorator syntax for creating transforms
    """
    c1 = TCoo1(ra=1*u.deg, dec=2*u.deg)

    @frame_transform_graph.transform(t.FunctionTransform, TCoo1, TCoo2)
    def trans(coo1, f):
        return TCoo2(ra=coo1.ra, dec=coo1.dec * 2)

    c2 = c1.transform_to(TCoo2)
    assert_allclose(c2.ra.degree, 1)
    assert_allclose(c2.dec.degree, 4)

    c3 = TCoo1(r.CartesianRepresentation(x=1*u.pc, y=1*u.pc, z=2*u.pc))

    @frame_transform_graph.transform(t.StaticMatrixTransform, TCoo1, TCoo2)
    def matrix():
        return [[2, 0, 0],
                [0, 1, 0],
                [0, 0, 1]]

    c4 = c3.transform_to(TCoo2)

    assert_allclose(c4.cartesian.x, 2*u.pc)
    assert_allclose(c4.cartesian.y, 1*u.pc)
    assert_allclose(c4.cartesian.z, 2*u.pc)


def test_shortest_path():
    class FakeTransform(object):
        def __init__(self, pri):
            self.priority = pri

    g = t.TransformGraph()

    # cheating by adding graph elements directly that are not classes - the
    # graphing algorithm still works fine with integers - it just isn't a valid
    # TransformGraph

    # the graph looks is a down-going diamond graph with the lower-right slightly
    # heavier and a cycle from the bottom to the top
    # also, a pair of nodes isolated from 1

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

    # unreachable
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
    from ...utils import NumpyRNGContext
    from .. import spherical_to_cartesian, cartesian_to_spherical

    x, y, z = spherical_to_cartesian(1, 0, 0)
    assert_allclose(x, 1)
    assert_allclose(y, 0)
    assert_allclose(z, 0)

    x, y, z = spherical_to_cartesian(0, 1, 1)
    assert_allclose(x, 0)
    assert_allclose(y, 0)
    assert_allclose(z, 0)

    x, y, z = spherical_to_cartesian(5, 0, np.arcsin(4. / 5.))
    assert_allclose(x, 3)
    assert_allclose(y, 4)
    assert_allclose(z, 0)

    r, lat, lon = cartesian_to_spherical(0, 1, 0)
    assert_allclose(r, 1)
    assert_allclose(lat, 0 * u.deg)
    assert_allclose(lon, np.pi / 2 * u.rad)

    # test round-tripping
    with NumpyRNGContext(13579):
        x, y, z = np.random.randn(3, 5)

    r, lat, lon = cartesian_to_spherical(x, y, z)
    x2, y2, z2 = spherical_to_cartesian(r, lat, lon)

    assert_allclose(x, x2)
    assert_allclose(y, y2)
    assert_allclose(z, z2)


def test_transform_path_pri():
    """
    This checks that the transformation path prioritization works by
    making sure the ICRS -> Gal transformation always goes through FK5
    and not FK4.
    """
    frame_transform_graph.invalidate_cache()
    tpath, td = frame_transform_graph.find_shortest_path(ICRS, Galactic)
    assert tpath == [ICRS, FK5, Galactic]
    assert td == 2

    # but direct from FK4 to Galactic should still be possible
    tpath, td = frame_transform_graph.find_shortest_path(FK4, Galactic)
    assert tpath == [FK4, FK4NoETerms, Galactic]
    assert td == 2


def test_obstime():
    """
    Checks to make sure observation time is
    accounted for at least in FK4 <-> ICRS transformations
    """
    b1950 = Time('B1950', scale='utc')
    j1975 = Time('J1975', scale='utc')

    fk4_50 = FK4(ra=1*u.deg, dec=2*u.deg, obstime=b1950)
    fk4_75 = FK4(ra=1*u.deg, dec=2*u.deg, obstime=j1975)

    icrs_50 = fk4_50.transform_to(ICRS)
    icrs_75 = fk4_75.transform_to(ICRS)

    # now check that the resulting coordinates are *different* - they should be,
    # because the obstime is different
    assert icrs_50.ra.degree != icrs_75.ra.degree
    assert icrs_50.dec.degree != icrs_75.dec.degree

# ------------------------------------------------------------------------------
# Affine transform tests and helpers:

# just acting as a namespace


class transfunc(object):
    rep = r.CartesianRepresentation(np.arange(3)*u.pc)
    dif = r.CartesianDifferential(*np.arange(3, 6)*u.pc/u.Myr)
    rep0 = r.CartesianRepresentation(np.zeros(3)*u.pc)

    @classmethod
    def both(cls, coo, fr):
        # exchange x <-> z and offset
        M = np.array([[0., 0., 1.],
                      [0., 1., 0.],
                      [1., 0., 0.]])
        return M, cls.rep.with_differentials(cls.dif)

    @classmethod
    def just_matrix(cls, coo, fr):
        # exchange x <-> z and offset
        M = np.array([[0., 0., 1.],
                      [0., 1., 0.],
                      [1., 0., 0.]])
        return M, None

    @classmethod
    def no_matrix(cls, coo, fr):
        return None, cls.rep.with_differentials(cls.dif)

    @classmethod
    def no_pos(cls, coo, fr):
        return None, cls.rep0.with_differentials(cls.dif)

    @classmethod
    def no_vel(cls, coo, fr):
        return None, cls.rep


@pytest.mark.parametrize('transfunc', [transfunc.both, transfunc.no_matrix,
                                       transfunc.no_pos, transfunc.no_vel,
                                       transfunc.just_matrix])
@pytest.mark.parametrize('rep', [
    r.CartesianRepresentation(5, 6, 7, unit=u.pc),
    r.CartesianRepresentation(5, 6, 7, unit=u.pc,
                              differentials=r.CartesianDifferential(8, 9, 10,
                                                                    unit=u.pc/u.Myr)),
    r.CartesianRepresentation(5, 6, 7, unit=u.pc,
                              differentials=r.CartesianDifferential(8, 9, 10,
                                                                    unit=u.pc/u.Myr))
     .represent_as(r.CylindricalRepresentation, r.CylindricalDifferential)
])
def test_affine_transform_succeed(transfunc, rep):
    c = TCoo1(rep)

    # compute expected output
    M, offset = transfunc(c, TCoo2)

    _rep = rep.to_cartesian()
    diffs = dict([(k, diff.represent_as(r.CartesianDifferential, rep))
                  for k, diff in rep.differentials.items()])
    expected_rep = _rep.with_differentials(diffs)

    if M is not None:
        expected_rep = expected_rep.transform(M)

    expected_pos = expected_rep.without_differentials()
    if offset is not None:
        expected_pos = expected_pos + offset.without_differentials()

    expected_vel = None
    if c.data.differentials:
        expected_vel = expected_rep.differentials['s']

        if offset and offset.differentials:
            expected_vel = (expected_vel + offset.differentials['s'])

    # register and do the transformation and check against expected
    trans = t.AffineTransform(transfunc, TCoo1, TCoo2)
    trans.register(frame_transform_graph)

    c2 = c.transform_to(TCoo2)

    assert quantity_allclose(c2.data.to_cartesian().xyz,
                             expected_pos.to_cartesian().xyz)

    if expected_vel is not None:
        diff = c2.data.differentials['s'].to_cartesian(base=c2.data)
        assert quantity_allclose(diff.xyz, expected_vel.d_xyz)

    trans.unregister(frame_transform_graph)


# these should fail
def transfunc_invalid_matrix(coo, fr):
    return np.eye(4), None

# Leaving this open in case we want to add more functions to check for failures


@pytest.mark.parametrize('transfunc', [transfunc_invalid_matrix])
def test_affine_transform_fail(transfunc):
    diff = r.CartesianDifferential(8, 9, 10, unit=u.pc/u.Myr)
    rep = r.CartesianRepresentation(5, 6, 7, unit=u.pc, differentials=diff)
    c = TCoo1(rep)

    # register and do the transformation and check against expected
    trans = t.AffineTransform(transfunc, TCoo1, TCoo2)
    trans.register(frame_transform_graph)

    with pytest.raises(ValueError):
        c2 = c.transform_to(TCoo2)

    trans.unregister(frame_transform_graph)


def test_too_many_differentials():
    dif1 = r.CartesianDifferential(*np.arange(3, 6)*u.pc/u.Myr)
    dif2 = r.CartesianDifferential(*np.arange(3, 6)*u.pc/u.Myr**2)
    rep = r.CartesianRepresentation(np.arange(3)*u.pc,
                                    differentials={'s': dif1, 's2': dif2})

    with pytest.raises(ValueError):
        c = TCoo1(rep)

    # register and do the transformation and check against expected
    trans = t.AffineTransform(transfunc.both, TCoo1, TCoo2)
    trans.register(frame_transform_graph)

    # Check that if frame somehow gets through to transformation, multiple
    # differentials are caught
    c = TCoo1(rep.without_differentials())
    c._data = c._data.with_differentials({'s': dif1, 's2': dif2})
    with pytest.raises(ValueError):
        c2 = c.transform_to(TCoo2)

    trans.unregister(frame_transform_graph)

# A matrix transform of a unit spherical with differentials should work


@pytest.mark.parametrize('rep', [
    r.UnitSphericalRepresentation(lon=15*u.degree, lat=-11*u.degree,
        differentials=r.SphericalDifferential(d_lon=15*u.mas/u.yr,
                                              d_lat=11*u.mas/u.yr,
                                              d_distance=-110*u.km/u.s)),
    r.UnitSphericalRepresentation(lon=15*u.degree, lat=-11*u.degree,
        differentials={'s': r.RadialDifferential(d_distance=-110*u.km/u.s)}),
    r.SphericalRepresentation(lon=15*u.degree, lat=-11*u.degree,
                              distance=150*u.pc,
        differentials={'s': r.RadialDifferential(d_distance=-110*u.km/u.s)})
])
def test_unit_spherical_with_differentials(rep):

    c = TCoo1(rep)

    # register and do the transformation and check against expected
    trans = t.AffineTransform(transfunc.just_matrix, TCoo1, TCoo2)
    trans.register(frame_transform_graph)
    c2 = c.transform_to(TCoo2)

    assert 's' in rep.differentials
    assert isinstance(c2.data.differentials['s'],
                      rep.differentials['s'].__class__)

    if isinstance(rep.differentials['s'], r.RadialDifferential):
        assert c2.data.differentials['s'] is rep.differentials['s']

    trans.unregister(frame_transform_graph)

    # should fail if we have to do offsets
    trans = t.AffineTransform(transfunc.both, TCoo1, TCoo2)
    trans.register(frame_transform_graph)

    with pytest.raises(TypeError):
        c.transform_to(TCoo2)

    trans.unregister(frame_transform_graph)


def test_vel_transformation_obstime_err():
    # TODO: replace after a final decision on PR #6280
    from ..sites import get_builtin_sites

    diff = r.CartesianDifferential([.1, .2, .3]*u.km/u.s)
    rep = r.CartesianRepresentation([1, 2, 3]*u.au, differentials=diff)

    loc = get_builtin_sites()['example_site']

    aaf = AltAz(obstime='J2010', location=loc)
    aaf2 = AltAz(obstime=aaf.obstime + 3*u.day, location=loc)
    aaf3 = AltAz(obstime=aaf.obstime + np.arange(3)*u.day, location=loc)
    aaf4 = AltAz(obstime=aaf.obstime, location=loc)

    aa = aaf.realize_frame(rep)

    with pytest.raises(NotImplementedError) as exc:
        aa.transform_to(aaf2)
    assert 'cannot transform' in exc.value.args[0]

    with pytest.raises(NotImplementedError) as exc:
        aa.transform_to(aaf3)
    assert 'cannot transform' in exc.value.args[0]

    aa.transform_to(aaf4)

    aa.transform_to(ICRS())


def test_function_transform_with_differentials():
    tfun = lambda c, f: f.__class__(ra=c.ra, dec=c.dec)
    ftrans = t.FunctionTransform(tfun, TCoo3, TCoo2,
                                 register_graph=frame_transform_graph)

    t3 = TCoo3(ra=1*u.deg, dec=2*u.deg, pm_ra_cosdec=1*u.marcsec/u.yr,
               pm_dec=1*u.marcsec/u.yr,)

    with catch_warnings() as w:
        t2 = t3.transform_to(TCoo2)
        assert len(w) == 1
        assert 'they have been dropped' in str(w[0].message)
