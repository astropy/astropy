# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from ... import units as u
from ..distances import Distance
from .. import transformations as t
from ..builtin_frames import ICRS, FK5, FK4, FK4NoETerms, Galactic, \
                             Supergalactic, Galactocentric, CIRS, GCRS, AltAz, \
                             ITRS, PrecessedGeocentric
from .. import representation as r
from ..baseframe import frame_transform_graph
from ...tests.helper import (pytest, quantity_allclose as allclose,
                             assert_quantity_allclose as assert_allclose)
from .utils import randomly_sample_sphere
from ...time import Time


#Coordinates just for these tests.
class TCoo1(ICRS):
    pass


class TCoo2(ICRS):
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
    diff = np.abs(g_xyz - gc_xyz)

    assert allclose(diff[0], 8.3*u.kpc, atol=1E-5*u.kpc)
    assert allclose(diff[1:], 0*u.kpc, atol=1E-5*u.kpc)

    # generate some test coordinates
    g = Galactic(l=[0,0,45,315]*u.deg, b=[-45,45,0,0]*u.deg,
                 distance=[np.sqrt(2)]*4*u.kpc)
    xyz = g.transform_to(Galactocentric(galcen_distance=1.*u.kpc, z_sun=0.*u.pc)).cartesian.xyz
    true_xyz = np.array([[0,0,-1.],[0,0,1],[0,1,0],[0,-1,0]]).T*u.kpc
    assert allclose(xyz.to(u.kpc), true_xyz.to(u.kpc), atol=1E-5*u.kpc)

    # check that ND arrays work

    # from Galactocentric to Galactic
    x = np.linspace(-10., 10., 100) * u.kpc
    y = np.linspace(-10., 10., 100) * u.kpc
    z = np.zeros_like(x)

    g1 = Galactocentric(x=x, y=y, z=z)
    g2 = Galactocentric(x=x.reshape(100,1,1), y=y.reshape(100,1,1), z=z.reshape(100,1,1))

    g1t = g1.transform_to(Galactic)
    g2t = g2.transform_to(Galactic)

    assert_allclose(g1t.cartesian.xyz, g2t.cartesian.xyz[:,:,0,0])

    # from Galactic to Galactocentric
    l = np.linspace(15, 30., 100) * u.deg
    b = np.linspace(-10., 10., 100) * u.deg
    d = np.ones_like(l.value) * u.kpc

    g1 = Galactic(l=l, b=b, distance=d)
    g2 = Galactic(l=l.reshape(100,1,1), b=b.reshape(100,1,1), distance=d.reshape(100,1,1))

    g1t = g1.transform_to(Galactocentric)
    g2t = g2.transform_to(Galactocentric)

    np.testing.assert_almost_equal(g1t.cartesian.xyz.value, g2t.cartesian.xyz.value[:,:,0,0])


def test_supergalactic():
    """
    Check Galactic<->Supergalactic and Galactic<->ICRS conversion.
    """
    # Check supergalactic North pole.
    npole = Galactic(l=47.37*u.degree, b=+6.32*u.degree)
    assert allclose(npole.transform_to(Supergalactic).sgb.deg, +90, atol=1e-9)

    # Check the origin of supergalactic longitude.
    lon0 = Supergalactic(sgl=0*u.degree, sgb=0*u.degree)
    lon0_gal = lon0.transform_to(Galactic)
    assert allclose(lon0_gal.l.deg, 137.37, atol=1e-9)
    assert allclose(lon0_gal.b.deg, 0, atol=1e-9)

    # Test Galactic<->ICRS with some positions that appear in Foley et al. 2008
    # (http://adsabs.harvard.edu/abs/2008A%26A...484..143F)

    # GRB 021219
    supergalactic = Supergalactic(sgl=29.91*u.degree, sgb=+73.72*u.degree)
    icrs = ICRS('18h50m27s +31d57m17s')
    assert supergalactic.separation(icrs) < 0.005 * u.degree

    # GRB 030320
    supergalactic = Supergalactic(sgl=-174.44*u.degree, sgb=+46.17*u.degree)
    icrs = ICRS('17h51m36s -25d18m52s')
    assert supergalactic.separation(icrs) < 0.005 * u.degree


def test_icrs_cirs():
    """
    Check a few cases of ICRS<->CIRS for consistency.

    Also includes the CIRS<->CIRS transforms at different times, as those go
    through ICRS
    """
    ra, dec, dist = randomly_sample_sphere(200)
    inod = ICRS(ra=ra, dec=dec)
    iwd = ICRS(ra=ra, dec=dec, distance=dist*u.pc)

    cframe1 = CIRS()
    cirsnod = inod.transform_to(cframe1)  #uses the default time
    #first do a round-tripping test
    inod2 = cirsnod.transform_to(ICRS)
    assert_allclose(inod.ra, inod2.ra)
    assert_allclose(inod.dec, inod2.dec)

    #now check that a different time yields different answers
    cframe2 = CIRS(obstime=Time('J2005', scale='utc'))
    cirsnod2 = inod.transform_to(cframe2)
    assert not allclose(cirsnod.ra, cirsnod2.ra, rtol=1e-8)
    assert not allclose(cirsnod.dec, cirsnod2.dec, rtol=1e-8)

    # parallax effects should be included, so with and w/o distance should be different
    cirswd = iwd.transform_to(cframe1)
    assert not allclose(cirswd.ra, cirsnod.ra, rtol=1e-8)
    assert not allclose(cirswd.dec, cirsnod.dec, rtol=1e-8)
    # and the distance should transform at least somehow
    assert not allclose(cirswd.distance, iwd.distance, rtol=1e-8)

    #now check that the cirs self-transform works as expected
    cirsnod3 = cirsnod.transform_to(cframe1)  # should be a no-op
    assert_allclose(cirsnod.ra, cirsnod3.ra)
    assert_allclose(cirsnod.dec, cirsnod3.dec)

    cirsnod4 = cirsnod.transform_to(cframe2)  # should be different
    assert not allclose(cirsnod4.ra, cirsnod.ra, rtol=1e-8)
    assert not allclose(cirsnod4.dec, cirsnod.dec, rtol=1e-8)

    cirsnod5 = cirsnod4.transform_to(cframe1)  # should be back to the same
    assert_allclose(cirsnod.ra, cirsnod5.ra)
    assert_allclose(cirsnod.dec, cirsnod5.dec)


def test_icrs_gcrs():
    """
    Check ICRS<->GCRS for consistency
    """
    ra, dec, dist = randomly_sample_sphere(200)
    inod = ICRS(ra=ra, dec=dec)
    iwd = ICRS(ra=ra, dec=dec, distance=dist*u.pc)

    gframe1 = GCRS()
    gcrsnod = inod.transform_to(gframe1)  #uses the default time
    #first do a round-tripping test
    inod2 = gcrsnod.transform_to(ICRS)
    assert_allclose(inod.ra, inod2.ra)
    assert_allclose(inod.dec, inod2.dec)

    #now check that a different time yields different answers
    gframe2 = GCRS(obstime=Time('J2005', scale='utc'))
    gcrsnod2 = inod.transform_to(gframe2)
    assert not allclose(gcrsnod.ra, gcrsnod2.ra, rtol=1e-8, atol=1e-10*u.deg)
    assert not allclose(gcrsnod.dec, gcrsnod2.dec, rtol=1e-8, atol=1e-10*u.deg)

    # parallax effects should be included, so with and w/o distance should be different
    gcrswd = iwd.transform_to(gframe1)
    assert not allclose(gcrswd.ra, gcrsnod.ra, rtol=1e-8, atol=1e-10*u.deg)
    assert not allclose(gcrswd.dec, gcrsnod.dec, rtol=1e-8, atol=1e-10*u.deg)
    # and the distance should transform at least somehow
    assert not allclose(gcrswd.distance, iwd.distance, rtol=1e-8,
                        atol=1e-10*u.pc)

    #now check that the cirs self-transform works as expected
    gcrsnod3 = gcrsnod.transform_to(gframe1)  # should be a no-op
    assert_allclose(gcrsnod.ra, gcrsnod3.ra)
    assert_allclose(gcrsnod.dec, gcrsnod3.dec)

    gcrsnod4 = gcrsnod.transform_to(gframe2)  # should be different
    assert not allclose(gcrsnod4.ra, gcrsnod.ra, rtol=1e-8, atol=1e-10*u.deg)
    assert not allclose(gcrsnod4.dec, gcrsnod.dec, rtol=1e-8, atol=1e-10*u.deg)

    gcrsnod5 = gcrsnod4.transform_to(gframe1)  # should be back to the same
    assert_allclose(gcrsnod.ra, gcrsnod5.ra, rtol=1e-8, atol=1e-10*u.deg)
    assert_allclose(gcrsnod.dec, gcrsnod5.dec, rtol=1e-8, atol=1e-10*u.deg)

    #also make sure that a GCRS with a different geoloc/geovel gets a different answer
    # roughly a moon-like frame
    gframe3 = GCRS(obsgeoloc=[385000., 0, 0]*u.km, obsgeovel=[1, 0, 0]*u.km/u.s)
    gcrsnod6 = inod.transform_to(gframe3)  # should be different
    assert not allclose(gcrsnod.ra, gcrsnod6.ra, rtol=1e-8, atol=1e-10*u.deg)
    assert not allclose(gcrsnod.dec, gcrsnod6.dec, rtol=1e-8, atol=1e-10*u.deg)
    inodviag3 = gcrsnod6.transform_to(ICRS)  # and now back to the original
    assert_allclose(inod.ra, inodviag3.ra)
    assert_allclose(inod.dec, inodviag3.dec)


def test_cirs_to_altaz():
    """
    Check the basic CIRS<->AltAz transforms.  More thorough checks implicitly
    happen in `test_iau_fullstack`
    """
    from .. import EarthLocation

    ra, dec, dist = randomly_sample_sphere(200)
    cirs = CIRS(ra=ra, dec=dec, obstime='J2000')
    crepr = r.SphericalRepresentation(lon=ra, lat=dec, distance=dist)
    cirscart = CIRS(crepr, obstime=cirs.obstime, representation=r.CartesianRepresentation)

    loc = EarthLocation(lat=0*u.deg, lon=0*u.deg, height=0*u.m)
    altazframe = AltAz(location=loc, obstime=Time('J2005'))

    cirs2 = cirs.transform_to(altazframe).transform_to(cirs)
    cirs3 = cirscart.transform_to(altazframe).transform_to(cirs)

    #check round-tripping
    assert_allclose(cirs.ra, cirs2.ra)
    assert_allclose(cirs.dec, cirs2.dec)
    assert_allclose(cirs.ra, cirs3.ra)
    assert_allclose(cirs.dec, cirs3.dec)


def test_gcrs_itrs():
    """
    Check basic GCRS<->ITRS transforms for round-tripping.
    """
    ra, dec, _ = randomly_sample_sphere(200)
    gcrs = GCRS(ra=ra, dec=dec, obstime='J2000')
    gcrs6 = GCRS(ra=ra, dec=dec, obstime='J2006')

    gcrs2 = gcrs.transform_to(ITRS).transform_to(gcrs)
    gcrs6_2 = gcrs6.transform_to(ITRS).transform_to(gcrs)

    assert_allclose(gcrs.ra, gcrs2.ra)
    assert_allclose(gcrs.dec, gcrs2.dec)
    assert not allclose(gcrs.ra, gcrs6_2.ra)
    assert not allclose(gcrs.dec, gcrs6_2.dec)

    #also try with the cartesian representation
    gcrsc = gcrs.realize_frame(gcrs.data)
    gcrsc.representation = r.CartesianRepresentation
    gcrsc2 = gcrsc.transform_to(ITRS).transform_to(gcrsc)
    assert_allclose(gcrsc.spherical.lon.deg, gcrsc2.ra.deg)
    assert_allclose(gcrsc.spherical.lat, gcrsc2.dec)


def test_cirs_itrs():
    """
    Check basic CIRS<->ITRS transforms for round-tripping.
    """
    ra, dec, _ = randomly_sample_sphere(200)
    cirs = CIRS(ra=ra, dec=dec, obstime='J2000')
    cirs6 = CIRS(ra=ra, dec=dec, obstime='J2006')

    cirs2 = cirs.transform_to(ITRS).transform_to(cirs)
    cirs6_2 = cirs6.transform_to(ITRS).transform_to(cirs) # different obstime

    #just check round-tripping
    assert_allclose(cirs.ra, cirs2.ra)
    assert_allclose(cirs.dec, cirs2.dec)
    assert not allclose(cirs.ra, cirs6_2.ra)
    assert not allclose(cirs.dec, cirs6_2.dec)


def test_gcrs_cirs():
    """
    Check GCRS<->CIRS transforms for round-tripping.  More complicated than the
    above two because it's multi-hop
    """
    ra, dec, _ = randomly_sample_sphere(200)
    gcrs = GCRS(ra=ra, dec=dec, obstime='J2000')
    gcrs6 = GCRS(ra=ra, dec=dec, obstime='J2006')

    gcrs2 = gcrs.transform_to(CIRS).transform_to(gcrs)
    gcrs6_2 = gcrs6.transform_to(CIRS).transform_to(gcrs)

    assert_allclose(gcrs.ra, gcrs2.ra)
    assert_allclose(gcrs.dec, gcrs2.dec)
    assert not allclose(gcrs.ra, gcrs6_2.ra)
    assert not allclose(gcrs.dec, gcrs6_2.dec)

    #now try explicit intermediate pathways and ensure they're all consistent
    gcrs3 = gcrs.transform_to(ITRS).transform_to(CIRS).transform_to(ITRS).transform_to(gcrs)
    assert_allclose(gcrs.ra, gcrs3.ra)
    assert_allclose(gcrs.dec, gcrs3.dec)

    gcrs4 = gcrs.transform_to(ICRS).transform_to(CIRS).transform_to(ICRS).transform_to(gcrs)
    assert_allclose(gcrs.ra, gcrs4.ra)
    assert_allclose(gcrs.dec, gcrs4.dec)


def test_gcrs_altaz():
    """
    Check GCRS<->AltAz transforms for round-tripping.  Has multiple paths
    """
    from .. import EarthLocation

    ra, dec, _ = randomly_sample_sphere(1)
    gcrs = GCRS(ra=ra[0], dec=dec[0], obstime='J2000')

    # check array times sure N-d arrays work
    times = Time(np.linspace(2456293.25, 2456657.25, 51) * u.day,
                 format='jd', scale='utc')

    loc = EarthLocation(lon=10 * u.deg, lat=80. * u.deg)
    aaframe = AltAz(obstime=times, location=loc)

    aa1 = gcrs.transform_to(aaframe)
    aa2 = gcrs.transform_to(ICRS).transform_to(CIRS).transform_to(aaframe)
    aa3 = gcrs.transform_to(ITRS).transform_to(CIRS).transform_to(aaframe)

    # make sure they're all consistent
    assert_allclose(aa1.alt, aa2.alt)
    assert_allclose(aa1.az, aa2.az)
    assert_allclose(aa1.alt, aa3.alt)
    assert_allclose(aa1.az, aa3.az)


def test_precessed_geocentric():
    assert PrecessedGeocentric().equinox.jd == Time('J2000', scale='utc').jd

    gcrs_coo = GCRS(180*u.deg, 2*u.deg, distance=10000*u.km)
    pgeo_coo = gcrs_coo.transform_to(PrecessedGeocentric)
    assert np.abs(gcrs_coo.ra - pgeo_coo.ra) > 10*u.marcsec
    assert np.abs(gcrs_coo.dec - pgeo_coo.dec) > 10*u.marcsec
    assert_allclose(gcrs_coo.distance, pgeo_coo.distance)

    gcrs_roundtrip = pgeo_coo.transform_to(GCRS)
    assert_allclose(gcrs_coo.ra, gcrs_roundtrip.ra)
    assert_allclose(gcrs_coo.dec, gcrs_roundtrip.dec)
    assert_allclose(gcrs_coo.distance, gcrs_roundtrip.distance)

    pgeo_coo2 = gcrs_coo.transform_to(PrecessedGeocentric(equinox='B1850'))
    assert np.abs(gcrs_coo.ra - pgeo_coo2.ra) > 1.5*u.deg
    assert np.abs(gcrs_coo.dec - pgeo_coo2.dec) > 0.5*u.deg
    assert_allclose(gcrs_coo.distance, pgeo_coo2.distance)

    gcrs2_roundtrip = pgeo_coo2.transform_to(GCRS)
    assert_allclose(gcrs_coo.ra, gcrs2_roundtrip.ra)
    assert_allclose(gcrs_coo.dec, gcrs2_roundtrip.dec)
    assert_allclose(gcrs_coo.distance, gcrs2_roundtrip.distance)
