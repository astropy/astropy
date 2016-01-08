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
                             ITRS, PrecessedGeocentric, Astrometric
from .. import SkyCoord
from .. import representation as r
from ..baseframe import frame_transform_graph
from ...tests.helper import (pytest, quantity_allclose as allclose,
                             assert_quantity_allclose as assert_allclose)
from .utils import randomly_sample_sphere
from ...time import Time

# used below in the next parametrized test
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

def test_astrometric():
    #Setup
    input_ra = np.linspace(0,360,10)
    input_dec = np.linspace(-90,90,10)
    input_ra_rad = np.deg2rad(input_ra)
    input_dec_rad = np.deg2rad(input_dec)
    icrs_coord = ICRS(ra = input_ra*u.deg,
                      dec = input_dec*u.deg, 
                      distance=1.*u.kpc)
    #RA rotations

    for ra in np.linspace(0,360,24):
        # expected rotation
        expected = ICRS(ra=np.linspace(0-ra,360-ra,10)*u.deg,
                        dec=np.linspace(-90,90,10)*u.deg,
                        distance=1.*u.kpc)  
        expected_xyz = expected.cartesian.xyz

        # actual transformation to the frame
        astrometric_frame = Astrometric(origin_ra=ra*u.deg,
                                        origin_dec=0.0*u.deg,
                                        origin_distance=0.*u.kpc)
        actual = icrs_coord.transform_to(astrometric_frame)
        actual_xyz = actual.cartesian.xyz
        assert allclose(actual_xyz.to(u.kpc), expected_xyz.to(u.kpc), atol=1E-5*u.kpc)

    #Dec rotations
    #Done in xyz space because dec must be [-90,90]

    for dec in np.linspace(-90,90,13):
        # expected rotation
        dec_rad = -np.deg2rad(dec)
        expected_x = (-np.sin(input_dec_rad) * np.sin(dec_rad) + 
                       np.cos(input_ra_rad) * np.cos(input_dec_rad) * np.cos(dec_rad))
        expected_y = (np.sin(input_ra_rad) * np.cos(input_dec_rad))  
        expected_z = (np.sin(input_dec_rad) * np.cos(dec_rad) + 
                      np.sin(dec_rad) * np.cos(input_ra_rad) * np.cos(input_dec_rad))
        expected = SkyCoord(x=expected_x,
                            y=expected_y,
                            z=expected_z, unit='kpc', representation='cartesian')    
        expected_xyz = expected.cartesian.xyz
    
        # actual transformation to the frame
        astrometric_frame = Astrometric(origin_ra=0.0*u.deg,
                                        origin_dec=dec*u.deg,
                                        origin_distance=0.*u.kpc)
        actual = icrs_coord.transform_to(astrometric_frame)
        actual_xyz = actual.cartesian.xyz
            
        assert allclose(actual_xyz.to(u.kpc), expected_xyz.to(u.kpc), atol=1E-5*u.kpc)

    #Both rotations
    for ra in np.linspace(0,360,24):
        for dec in np.linspace(-90,90,13):
            # expected rotation
            dec_rad = -np.deg2rad(dec)
            ra_rad = np.deg2rad(ra)
            expected_x = (-np.sin(input_dec_rad) * np.sin(dec_rad) + 
                           np.cos(input_ra_rad) * np.cos(input_dec_rad) * np.cos(dec_rad) * np.cos(ra_rad) + 
                           np.sin(input_ra_rad) * np.cos(input_dec_rad) * np.cos(dec_rad) * np.sin(ra_rad))
            expected_y = (np.sin(input_ra_rad) * np.cos(input_dec_rad) * np.cos(ra_rad) - 
                          np.cos(input_ra_rad) * np.cos(input_dec_rad) * np.sin(ra_rad))
            expected_z = (np.sin(input_dec_rad) * np.cos(dec_rad) + 
                          np.sin(dec_rad) * np.cos(ra_rad) * np.cos(input_ra_rad) * np.cos(input_dec_rad) +
                          np.sin(dec_rad) * np.sin(ra_rad) * np.sin(input_ra_rad) * np.cos(input_dec_rad))
            expected = SkyCoord(x=expected_x,
                                y=expected_y,
                                z=expected_z, unit='kpc', representation='cartesian')    
            expected_xyz = expected.cartesian.xyz
        
            # actual transformation to the frame
            astrometric_frame = Astrometric(origin_ra=ra*u.deg,
                                            origin_dec=dec*u.deg,
                                            origin_distance=0.*u.kpc)
            actual = icrs_coord.transform_to(astrometric_frame)
            actual_xyz = actual.cartesian.xyz
                
            assert allclose(actual_xyz.to(u.kpc), expected_xyz.to(u.kpc), atol=1E-5*u.kpc)
