# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from .. import core as erfa


def test_erfa_wrapper():
    """
    Runs a set of tests that mostly make sure vectorization is
    working as expected
    """

    jd = np.linspace(2456855.5, 2456855.5+1.0/24.0/60.0, 60*2+1)
    ra = np.linspace(0.0, np.pi*2.0, 5)
    dec = np.linspace(-np.pi/2.0, np.pi/2.0, 4)

    aob, zob, hob, dob, rob, eo = erfa.atco13(0.0,0.0,0.0,0.0,0.0,0.0,jd,0.0,0.0,0.0,np.pi/4.0,0.0,0.0,0.0,1014.0,0.0,0.0,0.5)
    assert aob.shape == (121,)

    aob, zob, hob, dob, rob, eo = erfa.atco13(0.0,0.0,0.0,0.0,0.0,0.0,jd[0],0.0,0.0,0.0,np.pi/4.0,0.0,0.0,0.0,1014.0,0.0,0.0,0.5)
    assert aob.shape == ()

    aob, zob, hob, dob, rob, eo = erfa.atco13(ra[:,None,None],dec[None,:,None],0.0,0.0,0.0,0.0,jd[None,None,:],0.0,0.0,0.0,np.pi/4.0,0.0,0.0,0.0,1014.0,0.0,0.0,0.5)
    (aob.shape) == (5, 4, 121)

    iy, im, id, ihmsf = erfa.d2dtf("UTC", 3, jd, 0.0)
    assert iy.shape == (121,)
    assert ihmsf.shape == (121, 4)
    assert ihmsf.dtype == np.dtype('i4')

    iy, im, id, ihmsf = erfa.d2dtf("UTC", 3, jd[0], 0.0)
    assert iy.shape == ()
    assert ihmsf.shape == (4,)
    assert ihmsf.dtype == np.dtype('i4')


def test_errwarn_reporting(recwarn):
    """
    Test that the ERFA error reporting mechanism works as it should
    """

    # no warning
    erfa.dat(1990, 1, 1, 0.5)

    # check warning is raised for a scalar
    erfa.dat(100, 1, 1, 0.5)
    w = recwarn.pop(erfa.ErfaWarning)
    assert '1 of "dubious year (Note 1)"' in str(w.message)

    # and that the count is right for a vector.
    erfa.dat([100, 200, 1990], 1, 1, 0.5)
    w = recwarn.pop(erfa.ErfaWarning)
    assert '2 of "dubious year (Note 1)"' in str(w.message)

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

    #values are from test_erfa.c t_ab function
    pnat = [-0.76321968546737951,
            -0.60869453983060384,
            -0.21676408580639883]
    v = [ 2.1044018893653786e-5,
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


def test_matrix_in():
    jd1 = 2456165.5
    jd2 = 0.401182685

    pvmat = np.empty((2, 3))
    pvmat[0][0] = -6241497.16
    pvmat[0][1] = 401346.896
    pvmat[0][2] = -1251136.04
    pvmat[1][0] = -29.264597
    pvmat[1][1] = -455.021831
    pvmat[1][2] = 0.0266151194

    astrom = erfa.apcs13(jd1, jd2, pvmat)
    assert astrom.shape == ()

    # values from t_erfa_c
    np.testing.assert_allclose(astrom['pmt'], 12.65133794027378508)
    np.testing.assert_allclose(astrom['em'], 1.010428384373318379)
    np.testing.assert_allclose(astrom['eb'], [.9012691529023298391,
                                             -.4173999812023068781,
                                             -.1809906511146821008])
    np.testing.assert_allclose(astrom['bpn'], np.eye(3))

    #first make sure it *fails* if we mess with the input orders
    pvmatbad = np.roll(pvmat.ravel(), 1).reshape((2, 3))
    astrombad = erfa.apcs13(jd1, jd2, pvmatbad)
    assert not np.allclose(astrombad['em'], 1.010428384373318379)

    pvmatarr = np.array([pvmat]*3)
    astrom2 = erfa.apcs13(jd1, jd2, pvmatarr)
    assert astrom2.shape == (3,)
    np.testing.assert_allclose(astrom2['em'], 1.010428384373318379)

    #try striding of the input array to make non-contiguous
    pvmatarr = np.array([pvmat]*9)[::3]
    astrom3 = erfa.apcs13(jd1, jd2, pvmatarr)
    assert astrom3.shape == (3,)
    np.testing.assert_allclose(astrom3['em'], 1.010428384373318379)

    #try fortran-order
    pvmatarr = np.array([pvmat]*3, order='F')
    astrom4 = erfa.apcs13(jd1, jd2, pvmatarr)
    assert astrom4.shape == (3,)
    np.testing.assert_allclose(astrom4['em'], 1.010428384373318379)


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


def test_structs_2():
    """
    Regression test for this case:
    https://github.com/astropy/astropy/pull/3181#issuecomment-66213070

    This does not test the numerical results, just that the return value has
    the correct type and shape.
    """

    jd = np.linspace(2456855.5, 2456855.5+1.0/24.0/60.0, 100)
    x = np.zeros((1,), erfa.dt_eraASTROM)
    out = erfa.aper(jd, x)

    assert out.shape == (100,)
    assert out.dtype == erfa.dt_eraASTROM
