# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np

from astropy import units as u
from astropy.coordinates import galactocentric_frame_defaults
from astropy.coordinates.distances import Distance
from astropy.coordinates.builtin_frames import (
    ICRS, FK5, FK4, FK4NoETerms, Galactic,
    Supergalactic, Galactocentric, HCRS, GCRS, LSR)
from astropy.coordinates import SkyCoord
from astropy.tests.helper import assert_quantity_allclose as assert_allclose
from astropy.coordinates import EarthLocation, CartesianRepresentation
from astropy.time import Time
from astropy.units import allclose

# used below in the next parametrized test
m31_sys = [ICRS, FK5, FK4, Galactic]
m31_coo = [(10.6847929, 41.2690650), (10.6847929, 41.2690650),
           (10.0004738, 40.9952444), (121.1744050, -21.5729360)]
m31_dist = Distance(770, u.kpc)
convert_precision = 1 * u.arcsec
roundtrip_precision = 1e-4 * u.degree
dist_precision = 1e-9 * u.kpc

m31_params = []
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
    coo2 = coo1.transform_to(tosys())
    if tosys is FK4:
        coo2_prec = coo2.transform_to(FK4(equinox=Time('B1950')))
        assert (coo2_prec.spherical.lon - tocoo[0]*u.deg) < convert_precision  # <1 arcsec
        assert (coo2_prec.spherical.lat - tocoo[1]*u.deg) < convert_precision
    else:
        assert (coo2.spherical.lon - tocoo[0]*u.deg) < convert_precision  # <1 arcsec
        assert (coo2.spherical.lat - tocoo[1]*u.deg) < convert_precision
    assert coo1.distance.unit == u.kpc
    assert coo2.distance.unit == u.kpc
    assert m31_dist.unit == u.kpc
    assert (coo2.distance - m31_dist) < dist_precision

    # check round-tripping
    coo1_2 = coo2.transform_to(fromsys())
    assert (coo1_2.spherical.lon - fromcoo[0]*u.deg) < roundtrip_precision
    assert (coo1_2.spherical.lat - fromcoo[1]*u.deg) < roundtrip_precision
    assert (coo1_2.distance - m31_dist) < dist_precision


def test_precession():
    """
    Ensures that FK4 and FK5 coordinates precess their equinoxes
    """
    j2000 = Time('J2000')
    b1950 = Time('B1950')
    j1975 = Time('J1975')
    b1975 = Time('B1975')

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

    direct = fk5.transform_to(Galactic())
    indirect = fk5.transform_to(FK4()).transform_to(Galactic())

    assert direct.separation(indirect).degree < 1.e-10

    direct = fk5.transform_to(Galactic())
    indirect = fk5.transform_to(FK4NoETerms()).transform_to(Galactic())

    assert direct.separation(indirect).degree < 1.e-10


def test_galactocentric():
    # when z_sun=0, transformation should be very similar to Galactic
    icrs_coord = ICRS(ra=np.linspace(0, 360, 10)*u.deg,
                      dec=np.linspace(-90, 90, 10)*u.deg,
                      distance=1.*u.kpc)

    g_xyz = icrs_coord.transform_to(Galactic()).cartesian.xyz
    with galactocentric_frame_defaults.set('pre-v4.0'):
        gc_xyz = icrs_coord.transform_to(Galactocentric(z_sun=0*u.kpc)).cartesian.xyz
    diff = np.abs(g_xyz - gc_xyz)

    assert allclose(diff[0], 8.3*u.kpc, atol=1E-5*u.kpc)
    assert allclose(diff[1:], 0*u.kpc, atol=1E-5*u.kpc)

    # generate some test coordinates
    g = Galactic(l=[0, 0, 45, 315]*u.deg, b=[-45, 45, 0, 0]*u.deg,
                 distance=[np.sqrt(2)]*4*u.kpc)
    with galactocentric_frame_defaults.set('pre-v4.0'):
        xyz = g.transform_to(Galactocentric(galcen_distance=1.*u.kpc, z_sun=0.*u.pc)).cartesian.xyz
    true_xyz = np.array([[0, 0, -1.], [0, 0, 1], [0, 1, 0], [0, -1, 0]]).T*u.kpc
    assert allclose(xyz.to(u.kpc), true_xyz.to(u.kpc), atol=1E-5*u.kpc)

    # check that ND arrays work

    # from Galactocentric to Galactic
    x = np.linspace(-10., 10., 100) * u.kpc
    y = np.linspace(-10., 10., 100) * u.kpc
    z = np.zeros_like(x)

    # from Galactic to Galactocentric
    l = np.linspace(15, 30., 100) * u.deg
    b = np.linspace(-10., 10., 100) * u.deg
    d = np.ones_like(l.value) * u.kpc

    with galactocentric_frame_defaults.set('latest'):
        g1 = Galactocentric(x=x, y=y, z=z)
        g2 = Galactocentric(x=x.reshape(100, 1, 1), y=y.reshape(100, 1, 1),
                            z=z.reshape(100, 1, 1))

        g1t = g1.transform_to(Galactic())
        g2t = g2.transform_to(Galactic())

        assert_allclose(g1t.cartesian.xyz, g2t.cartesian.xyz[:, :, 0, 0])

        g1 = Galactic(l=l, b=b, distance=d)
        g2 = Galactic(l=l.reshape(100, 1, 1), b=b.reshape(100, 1, 1),
                      distance=d.reshape(100, 1, 1))

        g1t = g1.transform_to(Galactocentric())
        g2t = g2.transform_to(Galactocentric())

        np.testing.assert_almost_equal(g1t.cartesian.xyz.value,
                                       g2t.cartesian.xyz.value[:, :, 0, 0])


def test_supergalactic():
    """
    Check Galactic<->Supergalactic and Galactic<->ICRS conversion.
    """
    # Check supergalactic North pole.
    npole = Galactic(l=47.37*u.degree, b=+6.32*u.degree)
    assert allclose(npole.transform_to(Supergalactic()).sgb.deg, +90, atol=1e-9)

    # Check the origin of supergalactic longitude.
    lon0 = Supergalactic(sgl=0*u.degree, sgb=0*u.degree)
    lon0_gal = lon0.transform_to(Galactic())
    assert allclose(lon0_gal.l.deg, 137.37, atol=1e-9)
    assert allclose(lon0_gal.b.deg, 0, atol=1e-9)

    # Test Galactic<->ICRS with some positions that appear in Foley et al. 2008
    # (https://ui.adsabs.harvard.edu/abs/2008A%26A...484..143F)

    # GRB 021219
    supergalactic = Supergalactic(sgl=29.91*u.degree, sgb=+73.72*u.degree)
    icrs = SkyCoord('18h50m27s +31d57m17s')
    assert supergalactic.separation(icrs) < 0.005 * u.degree

    # GRB 030320
    supergalactic = Supergalactic(sgl=-174.44*u.degree, sgb=+46.17*u.degree)
    icrs = SkyCoord('17h51m36s -25d18m52s')
    assert supergalactic.separation(icrs) < 0.005 * u.degree


class TestHCRS():
    """
    Check HCRS<->ICRS coordinate conversions.

    Uses ICRS Solar positions predicted by get_body_barycentric; with `t1` and
    `tarr` as defined below, the ICRS Solar positions were predicted using, e.g.
    coord.ICRS(coord.get_body_barycentric(tarr, 'sun')).
    """

    def setup(self):
        self.t1 = Time("2013-02-02T23:00")
        self.t2 = Time("2013-08-02T23:00")
        self.tarr = Time(["2013-02-02T23:00", "2013-08-02T23:00"])

        self.sun_icrs_scalar = ICRS(ra=244.52984668*u.deg,
                                    dec=-22.36943723*u.deg,
                                    distance=406615.66347377*u.km)
        # array of positions corresponds to times in `tarr`
        self.sun_icrs_arr = ICRS(ra=[244.52989062, 271.40976248]*u.deg,
                                 dec=[-22.36943605, -25.07431079]*u.deg,
                                 distance=[406615.66347377, 375484.13558956]*u.km)

        # corresponding HCRS positions
        self.sun_hcrs_t1 = HCRS(CartesianRepresentation([0.0, 0.0, 0.0] * u.km),
                                obstime=self.t1)
        twod_rep = CartesianRepresentation([[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]] * u.km)
        self.sun_hcrs_tarr = HCRS(twod_rep, obstime=self.tarr)
        self.tolerance = 5*u.km

    def test_from_hcrs(self):
        # test scalar transform
        transformed = self.sun_hcrs_t1.transform_to(ICRS())
        separation = transformed.separation_3d(self.sun_icrs_scalar)
        assert_allclose(separation, 0*u.km, atol=self.tolerance)

        # test non-scalar positions and times
        transformed = self.sun_hcrs_tarr.transform_to(ICRS())
        separation = transformed.separation_3d(self.sun_icrs_arr)
        assert_allclose(separation, 0*u.km, atol=self.tolerance)

    def test_from_icrs(self):
        # scalar positions
        transformed = self.sun_icrs_scalar.transform_to(HCRS(obstime=self.t1))
        separation = transformed.separation_3d(self.sun_hcrs_t1)
        assert_allclose(separation, 0*u.km, atol=self.tolerance)
        # nonscalar positions
        transformed = self.sun_icrs_arr.transform_to(HCRS(obstime=self.tarr))
        separation = transformed.separation_3d(self.sun_hcrs_tarr)
        assert_allclose(separation, 0*u.km, atol=self.tolerance)


class TestHelioBaryCentric():
    """
    Check GCRS<->Heliocentric and Barycentric coordinate conversions.

    Uses the WHT observing site (information grabbed from data/sites.json).
    """

    def setup(self):
        wht = EarthLocation(342.12*u.deg, 28.758333333333333*u.deg, 2327*u.m)
        self.obstime = Time("2013-02-02T23:00")
        self.wht_itrs = wht.get_itrs(obstime=self.obstime)

    def test_heliocentric(self):
        gcrs = self.wht_itrs.transform_to(GCRS(obstime=self.obstime))
        helio = gcrs.transform_to(HCRS(obstime=self.obstime))
        # Check it doesn't change from previous times.
        previous = [-1.02597256e+11, 9.71725820e+10, 4.21268419e+10] * u.m
        assert_allclose(helio.cartesian.xyz, previous)

        # And that it agrees with SLALIB to within 14km
        helio_slalib = [-0.685820296, 0.6495585893, 0.2816005464] * u.au
        assert np.sqrt(((helio.cartesian.xyz -
                         helio_slalib)**2).sum()) < 14. * u.km

    def test_barycentric(self):
        gcrs = self.wht_itrs.transform_to(GCRS(obstime=self.obstime))
        bary = gcrs.transform_to(ICRS())
        previous = [-1.02758958e+11, 9.68331109e+10, 4.19720938e+10] * u.m
        assert_allclose(bary.cartesian.xyz, previous)

        # And that it agrees with SLALIB answer to within 14km
        bary_slalib = [-0.6869012079, 0.6472893646, 0.2805661191] * u.au
        assert np.sqrt(((bary.cartesian.xyz -
                         bary_slalib)**2).sum()) < 14. * u.km


def test_lsr_sanity():

    # random numbers, but zero velocity in ICRS frame
    icrs = ICRS(ra=15.1241*u.deg, dec=17.5143*u.deg, distance=150.12*u.pc,
                pm_ra_cosdec=0*u.mas/u.yr, pm_dec=0*u.mas/u.yr,
                radial_velocity=0*u.km/u.s)
    lsr = icrs.transform_to(LSR())

    lsr_diff = lsr.data.differentials['s']
    cart_lsr_vel = lsr_diff.represent_as(CartesianRepresentation, base=lsr.data)
    lsr_vel = ICRS(cart_lsr_vel)
    gal_lsr = lsr_vel.transform_to(Galactic()).cartesian.xyz
    assert allclose(gal_lsr.to(u.km/u.s, u.dimensionless_angles()),
                    lsr.v_bary.d_xyz)

    # moving with LSR velocity
    lsr = LSR(ra=15.1241*u.deg, dec=17.5143*u.deg, distance=150.12*u.pc,
              pm_ra_cosdec=0*u.mas/u.yr, pm_dec=0*u.mas/u.yr,
              radial_velocity=0*u.km/u.s)
    icrs = lsr.transform_to(ICRS())

    icrs_diff = icrs.data.differentials['s']
    cart_vel = icrs_diff.represent_as(CartesianRepresentation, base=icrs.data)
    vel = ICRS(cart_vel)
    gal_icrs = vel.transform_to(Galactic()).cartesian.xyz
    assert allclose(gal_icrs.to(u.km/u.s, u.dimensionless_angles()),
                    -lsr.v_bary.d_xyz)


def test_hcrs_icrs_differentials():
    # Regression to ensure that we can transform velocities from HCRS to LSR.
    # Numbers taken from the original issue, gh-6835.
    hcrs = HCRS(ra=8.67*u.deg, dec=53.09*u.deg, distance=117*u.pc,
                pm_ra_cosdec=4.8*u.mas/u.yr, pm_dec=-15.16*u.mas/u.yr,
                radial_velocity=23.42*u.km/u.s)
    icrs = hcrs.transform_to(ICRS())

    # The position and velocity should not change much
    assert allclose(hcrs.cartesian.xyz, icrs.cartesian.xyz, rtol=1e-8)
    assert allclose(hcrs.velocity.d_xyz, icrs.velocity.d_xyz, rtol=1e-2)

    hcrs2 = icrs.transform_to(HCRS())

    # The values should round trip
    assert allclose(hcrs.cartesian.xyz, hcrs2.cartesian.xyz, rtol=1e-12)
    assert allclose(hcrs.velocity.d_xyz, hcrs2.velocity.d_xyz, rtol=1e-12)
