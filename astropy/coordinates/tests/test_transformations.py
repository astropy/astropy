# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from numpy import testing as npt

from ... import units as u
from ..distances import Distance
from .. import transformations as t
from ..builtin_frames import ICRS, FK5, FK4, FK4NoETerms, Galactic, \
                             Galactocentric, CIRS, GCRS, AltAz, ITRS
from .. import representation as r
from ..baseframe import frame_transform_graph
from ...tests.helper import pytest



#Coordinates just for these tests.
class TestCoo1(ICRS):
    pass


class TestCoo2(ICRS):
    pass


def test_transform_classes():
    """
    Tests the class-based/OO syntax for creating transforms
    """

    tfun = lambda c, f: f.__class__(ra=c.ra, dec=c.dec)
    trans1 = t.FunctionTransform(tfun, TestCoo1, TestCoo2,
                        register_graph=frame_transform_graph)

    c1 = TestCoo1(ra=1*u.radian, dec=0.5*u.radian)
    c2 = c1.transform_to(TestCoo2)
    npt.assert_allclose(c2.ra.radian, 1)
    npt.assert_allclose(c2.dec.radian, 0.5)


    def matfunc(coo, fr):
        return [[1, 0, 0],
                [0, coo.ra.degree, 0],
                [0, 0, 1]]
    trans2 = t.DynamicMatrixTransform(matfunc, TestCoo1, TestCoo2)
    trans2.register(frame_transform_graph)

    c3 = TestCoo1(ra=1*u.deg, dec=2*u.deg)
    c4 = c3.transform_to(TestCoo2)

    npt.assert_allclose(c4.ra.degree, 1)
    npt.assert_allclose(c4.ra.degree, 1)

    # be sure to unregister the second one - no need for trans1 because it
    # already got unregistered when trans2 was created.
    trans2.unregister(frame_transform_graph)


def test_transform_decos():
    """
    Tests the decorator syntax for creating transforms
    """
    c1 = TestCoo1(ra=1*u.deg, dec=2*u.deg)

    @frame_transform_graph.transform(t.FunctionTransform, TestCoo1, TestCoo2)
    def trans(coo1, f):
        return TestCoo2(ra=coo1.ra, dec=coo1.dec * 2)

    c2 = c1.transform_to(TestCoo2)
    npt.assert_allclose(c2.ra.degree, 1)
    npt.assert_allclose(c2.dec.degree, 4)

    c3 = TestCoo1(r.CartesianRepresentation(x=1*u.pc, y=1*u.pc, z=2*u.pc))

    @frame_transform_graph.transform(t.StaticMatrixTransform, TestCoo1, TestCoo2)
    def matrix():
        return [[2, 0, 0],
                [0, 1, 0],
                [0, 0, 1]]

    c4 = c3.transform_to(TestCoo2)

    npt.assert_allclose(c4.cartesian.x.value, 2)
    npt.assert_allclose(c4.cartesian.y.value, 1)
    npt.assert_allclose(c4.cartesian.z.value, 2)


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
    from numpy.testing.utils import assert_allclose
    from ...utils import NumpyRNGContext
    from .. import spherical_to_cartesian, cartesian_to_spherical

    x, y, z = spherical_to_cartesian(1, 0, 0)
    npt.assert_allclose(x, 1)
    npt.assert_allclose(y, 0)
    npt.assert_allclose(z, 0)

    x, y, z = spherical_to_cartesian(0, 1, 1)
    npt.assert_allclose(x, 0)
    npt.assert_allclose(y, 0)
    npt.assert_allclose(z, 0)

    x, y, z = spherical_to_cartesian(5, 0, np.arcsin(4. / 5.))
    npt.assert_allclose(x, 3)
    npt.assert_allclose(y, 4)
    npt.assert_allclose(z, 0)

    r, lat, lon = cartesian_to_spherical(0, 1, 0)
    npt.assert_allclose(r, 1)
    npt.assert_allclose(lat, 0)
    npt.assert_allclose(lon, np.pi / 2)

    #test round-tripping
    with NumpyRNGContext(13579):
        x, y, z = np.random.randn(3, 5)

    r, lat, lon = cartesian_to_spherical(x, y, z)
    x2, y2, z2 = spherical_to_cartesian(r, lat, lon)

    assert_allclose(x, x2)
    assert_allclose(y, y2)
    assert_allclose(z, z2)


m31_sys = [ICRS, FK5, FK4, Galactic]
m31_coo = [(10.6847929, 41.2690650), (10.6847929, 41.2690650), (10.0004738, 40.9952444), (121.1744050, -21.5729360)]
m31_dist = Distance(770, u.kpc)
convert_precision = 1 * u.arcsec
roundtrip_precision = 1e-4 * u.degree
dist_precision = 1e-9 * u.kpc

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

    from ...time import Time

    coo1 = fromsys(ra=fromcoo[0]*u.deg, dec=fromcoo[1]*u.deg, distance=m31_dist)
    coo2 = coo1.transform_to(tosys)
    if tosys is FK4:
        coo2_prec = coo2.transform_to(FK4(equinox=Time('B1950', scale='utc')))
        assert (coo2_prec.spherical.lon - tocoo[0]*u.deg) < convert_precision  # <1 arcsec
        assert (coo2_prec.spherical.lat - tocoo[1]*u.deg) < convert_precision
    else:
        assert (coo2.spherical.lon - tocoo[0]*u.deg) < convert_precision  # <1 arcsec
        assert (coo2.spherical.lat - tocoo[1]*u.deg) < convert_precision
    assert coo1.distance.unit == u.kpc
    assert coo2.distance.unit == u.kpc
    assert m31_dist.unit == u.kpc
    assert (coo2.distance - m31_dist) < dist_precision

    #check round-tripping
    coo1_2 = coo2.transform_to(fromsys)
    assert (coo1_2.spherical.lon - fromcoo[0]*u.deg) < roundtrip_precision
    assert (coo1_2.spherical.lat - fromcoo[1]*u.deg) < roundtrip_precision
    assert (coo1_2.distance - m31_dist) < dist_precision


def test_precession():
    """
    Ensures that FK4 and FK5 coordinates precess their equinoxes
    """
    from ...time import Time

    j2000 = Time('J2000', scale='utc')
    b1950 = Time('B1950', scale='utc')
    j1975 = Time('J1975', scale='utc')
    b1975 = Time('B1975', scale='utc')

    fk4 = FK4(ra=1*u.radian, dec=0.5*u.radian)
    assert fk4.equinox.byear == b1950.byear
    fk4_2 = fk4.transform_to(FK4(equinox=b1975))
    assert fk4_2.equinox.byear == b1975.byear

    fk5 = FK5(ra=1*u.radian, dec=0.5*u.radian)
    assert fk5.equinox.jyear == j2000.jyear
    fk5_2 = fk5.transform_to(FK4(equinox=j1975))
    assert fk5_2.equinox.jyear == j1975.jyear


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

    #but direct from FK4 to Galactic should still be possible
    tpath, td = frame_transform_graph.find_shortest_path(FK4, Galactic)
    assert tpath == [FK4, FK4NoETerms, Galactic]
    assert td == 2


def test_obstime():
    """
    Checks to make sure observation time is
    accounted for at least in FK4 <-> ICRS transformations
    """
    from ...time import Time

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


def test_fk5_galactic():
    """
    Check that FK5 -> Galactic gives the same as FK5 -> FK4 -> Galactic.
    """

    fk5 = FK5(ra=1*u.deg, dec=2*u.deg)

    direct = fk5.transform_to(Galactic)
    indirect = fk5.transform_to(FK4).transform_to(Galactic)

    assert direct.separation(indirect).degree < 1.e-10

    direct = fk5.transform_to(Galactic)
    indirect = fk5.transform_to(FK4NoETerms).transform_to(Galactic)

    assert direct.separation(indirect).degree < 1.e-10


def test_galactocentric():
    # when z_sun=0, transformation should be very similar to Galactic
    icrs_coord = ICRS(ra=np.linspace(0,360,10)*u.deg,
                      dec=np.linspace(-90,90,10)*u.deg,
                      distance=1.*u.kpc)

    g_xyz = icrs_coord.transform_to(Galactic).cartesian.xyz
    gc_xyz = icrs_coord.transform_to(Galactocentric(z_sun=0*u.kpc)).cartesian.xyz
    diff = np.array(np.abs(g_xyz - gc_xyz))

    assert np.allclose(diff[0], 8.3, atol=1E-5)
    assert np.allclose(diff[1:], 0, atol=1E-5)

    # generate some test coordinates
    g = Galactic(l=[0,0,45,315]*u.deg, b=[-45,45,0,0]*u.deg,
                 distance=[np.sqrt(2)]*4*u.kpc)
    xyz = g.transform_to(Galactocentric(galcen_distance=1.*u.kpc, z_sun=0.*u.pc)).cartesian.xyz
    true_xyz = np.array([[0,0,-1.],[0,0,1],[0,1,0],[0,-1,0]]).T*u.kpc
    assert np.allclose(xyz.to(u.kpc).value, true_xyz.to(u.kpc).value, atol=1E-5)

    # check that ND arrays work

    # from Galactocentric to Galactic
    x = np.linspace(-10., 10., 100) * u.kpc
    y = np.linspace(-10., 10., 100) * u.kpc
    z = np.zeros_like(x)

    g1 = Galactocentric(x=x, y=y, z=z)
    g2 = Galactocentric(x=x.reshape(100,1,1), y=y.reshape(100,1,1), z=z.reshape(100,1,1))

    g1t = g1.transform_to(Galactic)
    g2t = g2.transform_to(Galactic)

    np.testing.assert_almost_equal(g1t.cartesian.xyz.value, g2t.cartesian.xyz.value[:,:,0,0])

    # from Galactic to Galactocentric
    l = np.linspace(15, 30., 100) * u.deg
    b = np.linspace(-10., 10., 100) * u.deg
    d = np.ones_like(l.value) * u.kpc

    g1 = Galactic(l=l, b=b, distance=d)
    g2 = Galactic(l=l.reshape(100,1,1), b=b.reshape(100,1,1), distance=d.reshape(100,1,1))

    g1t = g1.transform_to(Galactocentric)
    g2t = g2.transform_to(Galactocentric)

    np.testing.assert_almost_equal(g1t.cartesian.xyz.value, g2t.cartesian.xyz.value[:,:,0,0])

def test_icrs_cirs():
    """
    Check a few cases of ICRS<->CIRS for consistency.

    Also includes the CIRS<->CIRS transforms at differnt times, as those go
    through ICRS
    """
    from ...time import Time
    from ...utils import NumpyRNGContext

    with NumpyRNGContext(12345):
        inod = ICRS(ra=np.random.rand(100)*360*u.deg,
                     dec=(np.random.rand(100)*180-90)*u.deg)
        iwd = ICRS(ra=inod.ra, dec=inod.dec, distance=np.random.rand(100)*u.pc)

    cframe1 = CIRS()
    cirsnod = inod.transform_to(cframe1)  #uses the default time
    #first do a round-tripping test
    inod2 = cirsnod.transform_to(ICRS)
    npt.assert_allclose(inod.ra, inod2.ra)
    npt.assert_allclose(inod.dec, inod2.dec)

    #now check that a different time yields different answers
    cframe2 = CIRS(obstime=Time('J2005', scale='utc'))
    cirsnod2 = inod.transform_to(cframe2)
    assert not np.allclose(cirsnod.ra, cirsnod2.ra, rtol=1e-8)
    assert not np.allclose(cirsnod.dec, cirsnod2.dec, rtol=1e-8)

    # parallax effects should be included, so with and w/o distance should be different
    cirswd = iwd.transform_to(cframe1)
    assert not np.allclose(cirswd.ra, cirsnod.ra, rtol=1e-8)
    assert not np.allclose(cirswd.dec, cirsnod.dec, rtol=1e-8)
    # and the distance should transform at least somehow
    assert not np.allclose(cirswd.distance, iwd.distance, rtol=1e-8)

    #now check that the cirs self-transform works as expected
    cirsnod3 = cirsnod.transform_to(cframe1)  # should be a no-op
    npt.assert_allclose(cirsnod.ra, cirsnod3.ra)
    npt.assert_allclose(cirsnod.dec, cirsnod3.dec)

    cirsnod4 = cirsnod.transform_to(cframe2)  # should be different
    assert not np.allclose(cirsnod4.ra, cirsnod.ra, rtol=1e-8)
    assert not np.allclose(cirsnod4.dec, cirsnod.dec, rtol=1e-8)

    cirsnod5 = cirsnod4.transform_to(cframe1)  # should be back to the same
    npt.assert_allclose(cirsnod.ra, cirsnod5.ra)
    npt.assert_allclose(cirsnod.dec, cirsnod5.dec)


def test_icrs_gcrs():
    """
    Check ICRS<->GCRS for consistency
    """
    from ...time import Time
    from ...utils import NumpyRNGContext

    with NumpyRNGContext(12345):
        inod = ICRS(ra=np.random.rand(100)*360*u.deg,
                     dec=(np.random.rand(100)*180-90)*u.deg)
        iwd = ICRS(ra=inod.ra, dec=inod.dec, distance=np.random.rand(100)*u.pc)

    gframe1 = GCRS()
    gcrsnod = inod.transform_to(gframe1)  #uses the default time
    #first do a round-tripping test
    inod2 = gcrsnod.transform_to(ICRS)
    npt.assert_allclose(inod.ra, inod2.ra)
    npt.assert_allclose(inod.dec, inod2.dec)

    #now check that a different time yields different answers
    gframe2 = GCRS(obstime=Time('J2005', scale='utc'))
    gcrsnod2 = inod.transform_to(gframe2)
    assert not np.allclose(gcrsnod.ra, gcrsnod2.ra, rtol=1e-8)
    assert not np.allclose(gcrsnod.dec, gcrsnod2.dec, rtol=1e-8)

    # parallax effects should be included, so with and w/o distance should be different
    gcrswd = iwd.transform_to(gframe1)
    assert not np.allclose(gcrswd.ra, gcrsnod.ra, rtol=1e-8)
    assert not np.allclose(gcrswd.dec, gcrsnod.dec, rtol=1e-8)
    # and the distance should transform at least somehow
    assert not np.allclose(gcrswd.distance, iwd.distance, rtol=1e-8)

    #now check that the cirs self-transform works as expected
    gcrsnod3 = gcrsnod.transform_to(gframe1)  # should be a no-op
    npt.assert_allclose(gcrsnod.ra, gcrsnod3.ra)
    npt.assert_allclose(gcrsnod.dec, gcrsnod3.dec)

    gcrsnod4 = gcrsnod.transform_to(gframe2)  # should be different
    assert not np.allclose(gcrsnod4.ra, gcrsnod.ra, rtol=1e-8)
    assert not np.allclose(gcrsnod4.dec, gcrsnod.dec, rtol=1e-8)

    gcrsnod5 = gcrsnod4.transform_to(gframe1)  # should be back to the same
    npt.assert_allclose(gcrsnod.ra, gcrsnod5.ra)
    npt.assert_allclose(gcrsnod.dec, gcrsnod5.dec)

    #also make sure that a GCRS with a different geoloc/geovel gets a different answer
    # roughly a moon-like frame
    gframe3 = GCRS(obsgeoloc=[385000., 0, 0]*u.km, obsgeovel=[1, 0, 0]*u.km/u.s)
    gcrsnod6 = inod.transform_to(gframe3)  # should be different
    assert not np.allclose(gcrsnod.ra, gcrsnod6.ra, rtol=1e-8)
    assert not np.allclose(gcrsnod.dec, gcrsnod6.dec, rtol=1e-8)
    inodviag3 = gcrsnod6.transform_to(ICRS)  # and now back to the original
    npt.assert_allclose(inod.ra, inodviag3.ra)
    npt.assert_allclose(inod.dec, inodviag3.dec)


def test_cirs_to_altaz():
    """
    Check the basic CIRS<->AltAz transforms.  More thorough checks implicitly
    happen in `test_iau_fullstack`
    """
    from ...time import Time
    from ...utils import NumpyRNGContext
    from .. import EarthLocation

    with NumpyRNGContext(12345):
        cirs = CIRS(ra=np.random.rand(100)*360*u.deg,
                    dec=(np.random.rand(100)*180-90)*u.deg,
                    obstime='J2000')
        cartrepr = r.SphericalRepresentation(lon=cirs.ra, lat=cirs.dec,
                                             distance=np.random.rand(len(cirs))
                                            ).to_cartesian()
        cirscart = CIRS(cartrepr, obstime=cirs.obstime,
                                  representation=r.CartesianRepresentation)
    loc = EarthLocation(lat=0*u.deg, lon=0*u.deg, height=0*u.m)
    altazframe = AltAz(location=loc, obstime=Time('J2005'))

    cirs2 = cirs.transform_to(altazframe).transform_to(cirs)
    cirs3 = cirscart.transform_to(altazframe).transform_to(cirs)

    #check round-tripping
    npt.assert_allclose(cirs.ra.deg, cirs2.ra.deg)
    npt.assert_allclose(cirs.dec.deg, cirs2.dec.deg)
    npt.assert_allclose(cirs.ra.deg, cirs3.ra.deg)
    npt.assert_allclose(cirs.dec.deg, cirs3.dec.deg)

def test_gcrs_itrs():
    """
    Check basic GCRS<->ITRS transforms for round-tripping.
    """
    from ...utils import NumpyRNGContext

    with NumpyRNGContext(12345):
        gcrs = GCRS(ra=np.random.rand(100)*360*u.deg,
                    dec=(np.random.rand(100)*180-90)*u.deg,
                    obstime='J2000')
        gcrs6 = GCRS(ra=np.random.rand(100)*360*u.deg,
                     dec=(np.random.rand(100)*180-90)*u.deg,
                     obstime='J2006')

    gcrs2 = gcrs.transform_to(ITRS).transform_to(gcrs)
    gcrs6_2 = gcrs6.transform_to(ITRS).transform_to(gcrs)

    npt.assert_allclose(gcrs.ra.deg, gcrs2.ra.deg)
    npt.assert_allclose(gcrs.dec.deg, gcrs2.dec.deg)
    assert not np.allclose(gcrs.ra.deg, gcrs6_2.ra.deg)
    assert not np.allclose(gcrs.dec.deg, gcrs6_2.dec.deg)

    #also try with the cartesian representation
    gcrsc = gcrs.realize_frame(gcrs.data)
    gcrsc.representation = r.CartesianRepresentation
    gcrsc2 = gcrsc.transform_to(ITRS).transform_to(gcrsc)
    npt.assert_allclose(gcrsc.spherical.lon.deg, gcrsc2.ra.deg)
    npt.assert_allclose(gcrsc.spherical.lat.deg, gcrsc2.dec.deg)




def test_cirs_itrs():
    """
    Check basic CIRS<->ITRS transforms for round-tripping.
    """
    from ...utils import NumpyRNGContext

    with NumpyRNGContext(12345):
        cirs = CIRS(ra=np.random.rand(100)*360*u.deg,
                    dec=(np.random.rand(100)*180-90)*u.deg,
                    obstime='J2000')
        cirs6 = CIRS(ra=np.random.rand(100)*360*u.deg,
                     dec=(np.random.rand(100)*180-90)*u.deg,
                     obstime='J2006')

    cirs2 = cirs.transform_to(ITRS).transform_to(cirs)
    cirs6_2 = cirs6.transform_to(ITRS).transform_to(cirs) # different obstime

    #just check round-tripping
    npt.assert_allclose(cirs.ra.deg, cirs2.ra.deg)
    npt.assert_allclose(cirs.dec.deg, cirs2.dec.deg)
    assert not np.allclose(cirs.ra.deg, cirs6_2.ra.deg)
    assert not np.allclose(cirs.dec.deg, cirs6_2.dec.deg)



def test_gcrs_cirs():
    """
    Check GCRS<->CIRS transforms for round-tripping.  More complicated than the
    above two because it's multi-hop
    """
    from ...utils import NumpyRNGContext

    with NumpyRNGContext(12345):
        gcrs = GCRS(ra=np.random.rand(100)*360*u.deg,
                    dec=(np.random.rand(100)*180-90)*u.deg,
                    obstime='J2000')
        gcrs6 = GCRS(ra=np.random.rand(100)*360*u.deg,
                     dec=(np.random.rand(100)*180-90)*u.deg,
                     obstime='J2006')

    gcrs2 = gcrs.transform_to(CIRS).transform_to(gcrs)
    gcrs6_2 = gcrs6.transform_to(CIRS).transform_to(gcrs)

    npt.assert_allclose(gcrs.ra.deg, gcrs2.ra.deg)
    npt.assert_allclose(gcrs.dec.deg, gcrs2.dec.deg)
    assert not np.allclose(gcrs.ra.deg, gcrs6_2.ra.deg)
    assert not np.allclose(gcrs.dec.deg, gcrs6_2.dec.deg)

    #now try explicit intermediate pathways and ensure they're all consistent
    gcrs3 = gcrs.transform_to(ITRS).transform_to(CIRS).transform_to(ITRS).transform_to(gcrs)
    npt.assert_allclose(gcrs.ra.deg, gcrs3.ra.deg)
    npt.assert_allclose(gcrs.dec.deg, gcrs3.dec.deg)

    gcrs4 = gcrs.transform_to(ICRS).transform_to(CIRS).transform_to(ICRS).transform_to(gcrs)
    npt.assert_allclose(gcrs.ra.deg, gcrs4.ra.deg)
    npt.assert_allclose(gcrs.dec.deg, gcrs4.dec.deg)
