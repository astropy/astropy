# Licensed under a 3-clause BSD style license - see LICENSE.rst
from datetime import datetime

import pytest
import numpy as np
from numpy.testing import assert_array_equal

from astropy import _erfa as erfa
from astropy.time import Time
from astropy.tests.helper import catch_warnings


def test_erfa_wrapper():
    """
    Runs a set of tests that mostly make sure vectorization is
    working as expected
    """

    jd = np.linspace(2456855.5, 2456855.5+1.0/24.0/60.0, 60*2+1)
    ra = np.linspace(0.0, np.pi*2.0, 5)
    dec = np.linspace(-np.pi/2.0, np.pi/2.0, 4)

    aob, zob, hob, dob, rob, eo = erfa.atco13(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, jd, 0.0, 0.0, 0.0, np.pi/4.0, 0.0, 0.0, 0.0, 1014.0, 0.0, 0.0, 0.5)
    assert aob.shape == (121,)

    aob, zob, hob, dob, rob, eo = erfa.atco13(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, jd[0], 0.0, 0.0, 0.0, np.pi/4.0, 0.0, 0.0, 0.0, 1014.0, 0.0, 0.0, 0.5)
    assert aob.shape == ()

    aob, zob, hob, dob, rob, eo = erfa.atco13(ra[:, None, None], dec[None, :, None], 0.0, 0.0, 0.0, 0.0, jd[None, None, :], 0.0, 0.0, 0.0, np.pi/4.0, 0.0, 0.0, 0.0, 1014.0, 0.0, 0.0, 0.5)
    (aob.shape) == (5, 4, 121)

    iy, im, id, ihmsf = erfa.d2dtf("UTC", 3, jd, 0.0)
    assert iy.shape == (121,)
    assert ihmsf.shape == (121,)
    assert ihmsf.dtype == erfa.dt_hmsf

    iy, im, id, ihmsf = erfa.d2dtf("UTC", 3, jd[0], 0.0)
    assert iy.shape == ()
    assert ihmsf.shape == ()
    assert ihmsf.dtype == erfa.dt_hmsf


def test_angle_ops():

    sign, idmsf = erfa.a2af(6, -np.pi)
    assert sign == b'-'
    assert idmsf.item() == (180, 0, 0, 0)

    sign, ihmsf = erfa.a2tf(6, np.pi)
    assert sign == b'+'
    assert ihmsf.item() == (12, 0, 0, 0)

    rad = erfa.af2a('-', 180, 0, 0.0)
    np.testing.assert_allclose(rad, -np.pi)

    rad = erfa.tf2a('+', 12, 0, 0.0)
    np.testing.assert_allclose(rad, np.pi)

    rad = erfa.anp(3.*np.pi)
    np.testing.assert_allclose(rad, np.pi)

    rad = erfa.anpm(3.*np.pi)
    np.testing.assert_allclose(rad, -np.pi)

    sign, ihmsf = erfa.d2tf(1, -1.5)
    assert sign == b'-'
    assert ihmsf.item() == (36, 0, 0, 0)

    days = erfa.tf2d('+', 3, 0, 0.0)
    np.testing.assert_allclose(days, 0.125)


def test_spherical_cartesian():

    theta, phi = erfa.c2s([0.0, np.sqrt(2.0), np.sqrt(2.0)])
    np.testing.assert_allclose(theta, np.pi/2.0)
    np.testing.assert_allclose(phi, np.pi/4.0)

    theta, phi, r = erfa.p2s([0.0, np.sqrt(2.0), np.sqrt(2.0)])
    np.testing.assert_allclose(theta, np.pi/2.0)
    np.testing.assert_allclose(phi, np.pi/4.0)
    np.testing.assert_allclose(r, 2.0)

    pv = np.array(([0.0, np.sqrt(2.0), np.sqrt(2.0)], [1.0, 0.0, 0.0]),
                  dtype=erfa.dt_pv)
    theta, phi, r, td, pd, rd = erfa.pv2s(pv)
    np.testing.assert_allclose(theta, np.pi/2.0)
    np.testing.assert_allclose(phi, np.pi/4.0)
    np.testing.assert_allclose(r, 2.0)
    np.testing.assert_allclose(td, -np.sqrt(2.0)/2.0)
    np.testing.assert_allclose(pd, 0.0)
    np.testing.assert_allclose(rd, 0.0)

    c = erfa.s2c(np.pi/2.0, np.pi/4.0)
    np.testing.assert_allclose(c, [0.0, np.sqrt(2.0)/2.0, np.sqrt(2.0)/2.0], atol=1e-14)

    c = erfa.s2p(np.pi/2.0, np.pi/4.0, 1.0)
    np.testing.assert_allclose(c, [0.0, np.sqrt(2.0)/2.0, np.sqrt(2.0)/2.0], atol=1e-14)

    pv = erfa.s2pv(np.pi/2.0, np.pi/4.0, 2.0, np.sqrt(2.0)/2.0, 0.0, 0.0)
    np.testing.assert_allclose(pv['p'], [0.0, np.sqrt(2.0), np.sqrt(2.0)], atol=1e-14)
    np.testing.assert_allclose(pv['v'],  [-1.0, 0.0, 0.0], atol=1e-14)


def test_errwarn_reporting():
    """
    Test that the ERFA error reporting mechanism works as it should
    """

    # no warning
    erfa.dat(1990, 1, 1, 0.5)

    # check warning is raised for a scalar
    with catch_warnings() as w:
        erfa.dat(100, 1, 1, 0.5)
        assert len(w) == 1
        assert w[0].category == erfa.ErfaWarning
        assert '1 of "dubious year (Note 1)"' in str(w[0].message)

    # and that the count is right for a vector.
    with catch_warnings() as w:
        erfa.dat([100, 200, 1990], 1, 1, 0.5)
        assert len(w) == 1
        assert w[0].category == erfa.ErfaWarning
        assert '2 of "dubious year (Note 1)"' in str(w[0].message)

    try:
        erfa.dat(1990, [1, 34, 2], [1, 1, 43], 0.5)
    except erfa.ErfaError as e:
        if '1 of "bad day (Note 3)", 1 of "bad month"' not in e.args[0]:
            assert False, 'Raised the correct type of error, but wrong message: ' + e.args[0]

    try:
        erfa.dat(200, [1, 34, 2], [1, 1, 43], 0.5)
    except erfa.ErfaError as e:
        if 'warning' in e.args[0]:
            assert False, 'Raised the correct type of error, but there were warnings mixed in: ' + e.args[0]


def test_vector_inouts():
    """
    Tests that ERFA functions working with vectors are correctly consumed and spit out
    """

    # values are from test_erfa.c t_ab function
    pnat = [-0.76321968546737951,
            -0.60869453983060384,
            -0.21676408580639883]
    v = [2.1044018893653786e-5,
         -8.9108923304429319e-5,
         -3.8633714797716569e-5]
    s = 0.99980921395708788
    bm1 = 0.99999999506209258

    expected = [-0.7631631094219556269,
                -0.6087553082505590832,
                -0.2167926269368471279]

    res = erfa.ab(pnat, v, s, bm1)
    assert res.shape == (3,)

    np.testing.assert_allclose(res, expected)

    res2 = erfa.ab([pnat]*4, v, s, bm1)
    assert res2.shape == (4, 3)
    np.testing.assert_allclose(res2, [expected]*4)

    # here we stride an array and also do it Fortran-order to make sure
    # it all still works correctly with non-contig arrays
    pnata = np.array(pnat)
    arrin = np.array([pnata, pnata/2, pnata/3, pnata/4, pnata/5]*4, order='F')
    res3 = erfa.ab(arrin[::5], v, s, bm1)
    assert res3.shape == (4, 3)
    np.testing.assert_allclose(res3, [expected]*4)


def test_pv_in():
    jd1 = 2456165.5
    jd2 = 0.401182685

    pv = np.empty((), dtype=erfa.dt_pv)
    pv['p'] = [-6241497.16,
               401346.896,
               -1251136.04]
    pv['v'] = [-29.264597,
               -455.021831,
               0.0266151194]

    astrom = erfa.apcs13(jd1, jd2, pv)
    assert astrom.shape == ()

    # values from t_erfa_c
    np.testing.assert_allclose(astrom['pmt'], 12.65133794027378508)
    np.testing.assert_allclose(astrom['em'], 1.010428384373318379)
    np.testing.assert_allclose(astrom['eb'], [0.9012691529023298391,
                                              -.4173999812023068781,
                                              -.1809906511146821008])
    np.testing.assert_allclose(astrom['bpn'], np.eye(3))

    # first make sure it *fails* if we mess with the input orders
    pvbad = np.empty_like(pv)
    pvbad['p'], pvbad['v'] = pv['v'], pv['p']
    astrombad = erfa.apcs13(jd1, jd2, pvbad)
    assert not np.allclose(astrombad['em'], 1.010428384373318379)

    pvarr = np.array([pv]*3)
    astrom2 = erfa.apcs13(jd1, jd2, pvarr)
    assert astrom2.shape == (3,)
    np.testing.assert_allclose(astrom2['em'], 1.010428384373318379)

    # try striding of the input array to make non-contiguous
    pvmatarr = np.array([pv]*9)[::3]
    astrom3 = erfa.apcs13(jd1, jd2, pvmatarr)
    assert astrom3.shape == (3,)
    np.testing.assert_allclose(astrom3['em'], 1.010428384373318379)


def test_structs():
    """
    Checks producing and consuming of ERFA c structs
    """

    am, eo = erfa.apci13(2456165.5, [0.401182685, 1])
    assert am.shape == (2, )
    assert am.dtype == erfa.dt_eraASTROM
    assert eo.shape == (2, )

    # a few spotchecks from test_erfa.c
    np.testing.assert_allclose(am[0]['pmt'], 12.65133794027378508)
    np.testing.assert_allclose(am[0]['v'], [0.4289638897157027528e-4,
                                            0.8115034002544663526e-4,
                                            0.3517555122593144633e-4])

    ri, di = erfa.atciqz(2.71, 0.174, am[0])
    np.testing.assert_allclose(ri, 2.709994899247599271)
    np.testing.assert_allclose(di, 0.1728740720983623469)


def test_float32_input():
    # Regression test for gh-8615
    xyz = np.array([[1, 0, 0], [0.9, 0.1, 0]])
    out64 = erfa.p2s(xyz)
    out32 = erfa.p2s(xyz.astype('f4'))
    np.testing.assert_allclose(out32, out64, rtol=1.e-5)


class TestAstromNotInplace:
    def setup(self):
        self.mjd_array = np.array(
            [58827.15925499, 58827.15925499, 58827.15925499,
             58827.15925499, 58827.15925499])
        self.mjd_scalar = self.mjd_array[0].item()
        self.utc2mjd = 2400000.5
        paranal_long = -1.228798
        paranal_lat = -0.42982
        paranal_height = 2669.
        self.astrom_scalar, _ = erfa.apco13(
            self.utc2mjd, self.mjd_scalar, 0.0, paranal_long, paranal_lat,
            paranal_height, 0.0, 0.0, 0.0, 0.0, 0.0, 2.5)
        self.astrom_array, _ = erfa.apco13(
            self.utc2mjd, self.mjd_array, 0.0, paranal_long, paranal_lat,
            paranal_height, 0.0, 0.0, 0.0, 0.0, 0.0, 2.5)

    def test_scalar_input(self):
        # Regression test for gh-9799, where astrom0 being a void
        # caused a TypeError, as it was trying to change it in-place.
        assert type(self.astrom_scalar) is np.void
        astrom = erfa.aper13(self.utc2mjd, self.mjd_scalar, self.astrom_scalar)
        assert astrom is not self.astrom_scalar
        assert type(astrom) is np.void

    def test_array_input(self):
        # Trying to fix gh-9799, it became clear that doing things in-place was
        # a bad idea generally (and didn't work), so also for array input we
        # now return a copy.
        assert type(self.astrom_array) is np.ndarray
        astrom = erfa.aper13(self.utc2mjd, self.mjd_array, self.astrom_array)
        assert astrom is not self.astrom_array
        assert type(astrom) is np.ndarray


class TestLeapSecondsBasics:
    def test_get_leap_seconds(self):
        leap_seconds = erfa.leap_seconds.get()
        assert isinstance(leap_seconds, np.ndarray)
        assert leap_seconds.dtype is erfa.dt_eraLEAPSECOND
        # Very basic sanity checks.
        assert np.all((leap_seconds['year'] >= 1960) &
                      (leap_seconds['year'] < 3000))
        assert np.all((leap_seconds['month'] >= 1) &
                      (leap_seconds['month'] <= 12))
        assert np.all(abs(leap_seconds['tai_utc'] < 1000.))

    def test_leap_seconds_expires(self):
        expires = erfa.leap_seconds.expires
        assert isinstance(expires, datetime)
        last_ls = erfa.leap_seconds.get()[-1]
        dt_last = datetime(last_ls['year'], last_ls['month'], 1)
        assert expires > dt_last


class TestLeapSeconds:
    """Test basic methods to control the ERFA leap-second table."""
    def setup(self):
        self.erfa_ls = erfa.leap_seconds.get()
        self.expires = erfa.leap_seconds.expires
        self._expires = erfa.leap_seconds._expires

    def teardown(self):
        erfa.leap_seconds.set(self.erfa_ls)
        erfa.leap_seconds._expires = self._expires

    def test_set_reset_leap_seconds(self):
        erfa.leap_seconds.set()
        leap_seconds = erfa.leap_seconds.get()

        erfa.leap_seconds.set(leap_seconds[:-2])
        new_leap_seconds = erfa.leap_seconds.get()
        assert_array_equal(new_leap_seconds, leap_seconds[:-2])

        erfa.leap_seconds.set()
        reset_leap_seconds = erfa.leap_seconds.get()
        assert_array_equal(reset_leap_seconds, leap_seconds)

    def test_set_leap_seconds(self):
        assert erfa.dat(2018, 1, 1, 0.) == 37.0
        leap_seconds = erfa.leap_seconds.get()
        # Set to a table that misses the 2017 leap second.
        part_leap_seconds = leap_seconds[leap_seconds['year'] < 2017]
        erfa.leap_seconds.set(part_leap_seconds)
        new_leap_seconds = erfa.leap_seconds.get()
        assert_array_equal(new_leap_seconds, part_leap_seconds)
        # Check the 2017 leap second is indeed missing.
        assert erfa.dat(2018, 1, 1, 0.) == 36.0
        # And that this would be expected from the expiration date.
        assert erfa.leap_seconds.expires < datetime(2018, 1, 1)
        assert erfa.leap_seconds.expired
        # Reset and check it is back.
        erfa.leap_seconds.set()
        assert erfa.dat(2018, 1, 1, 0.) == 37.0

    @pytest.mark.parametrize('table,match', [
        ([(2017, 3, 10.)], 'January'),
        ([(2017, 1, 1.),
          (2017, 7, 3.)], 'jump'),
        ([[(2017, 1, 1.)],
          [(2017, 7, 2.)]], 'dimension')])
    def test_validation(self, table, match):
        with pytest.raises(ValueError, match=match):
            erfa.leap_seconds.set(table)
        # Check leap-second table is not corrupted.
        assert_array_equal(erfa.leap_seconds.get(), self.erfa_ls)
        assert erfa.dat(2018, 1, 1, 0.) == 37.0

    def test_update_leap_seconds(self):
        assert erfa.dat(2018, 1, 1, 0.) == 37.0
        leap_seconds = erfa.leap_seconds.get()
        # Get old and new leap seconds
        old_leap_seconds = leap_seconds[leap_seconds['year'] < 2017]
        new_leap_seconds = leap_seconds[leap_seconds['year'] >= 2017]
        # Updating with either of these should do nothing.
        n_update = erfa.leap_seconds.update(new_leap_seconds)
        assert n_update == 0
        assert_array_equal(erfa.leap_seconds.get(), self.erfa_ls)
        n_update = erfa.leap_seconds.update(old_leap_seconds)
        assert n_update == 0
        assert_array_equal(erfa.leap_seconds.get(), self.erfa_ls)

        # But after setting to older part, update with newer should work.
        erfa.leap_seconds.set(old_leap_seconds)
        # Check the 2017 leap second is indeed missing.
        assert erfa.dat(2018, 1, 1, 0.) == 36.0
        # Update with missing leap seconds.
        n_update = erfa.leap_seconds.update(new_leap_seconds)
        assert n_update == len(new_leap_seconds)
        assert erfa.dat(2018, 1, 1, 0.) == 37.0

        # Also a final try with overlapping data.
        erfa.leap_seconds.set(old_leap_seconds)
        n_update = erfa.leap_seconds.update(leap_seconds)
        assert n_update == len(new_leap_seconds)
        assert erfa.dat(2018, 1, 1, 0.) == 37.0

    @pytest.mark.parametrize('expiration', [
        datetime(2345, 1, 1),
        '1 January 2345',
        Time('2345-01-01', scale='tai')])
    def test_with_expiration(self, expiration):
        class ExpiringArray(np.ndarray):
            expires = expiration

        leap_seconds = erfa.leap_seconds.get()
        erfa.leap_seconds.set(leap_seconds.view(ExpiringArray))
        assert erfa.leap_seconds.expires == datetime(2345, 1, 1)

        # Get old and new leap seconds
        old_leap_seconds = leap_seconds[:-10]
        new_leap_seconds = leap_seconds[-10:]

        erfa.leap_seconds.set(old_leap_seconds)
        # Check expiration is reset
        assert erfa.leap_seconds.expires != datetime(2345, 1, 1)
        # Update with missing leap seconds.
        n_update = erfa.leap_seconds.update(
            new_leap_seconds.view(ExpiringArray))
        assert n_update == len(new_leap_seconds)
        assert erfa.leap_seconds.expires == datetime(2345, 1, 1)

    def test_with_expiration_warning(self):
        class ExpiringArray(np.ndarray):
            expires = 'incomprehensible'

        leap_seconds = erfa.leap_seconds.get()
        with pytest.warns(erfa.ErfaWarning,
                          match='non-datetime.*parsing it raised'):
            erfa.leap_seconds.set(leap_seconds.view(ExpiringArray))
