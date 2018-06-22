# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from .. import core as erfa
from ...tests.helper import catch_warnings


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
