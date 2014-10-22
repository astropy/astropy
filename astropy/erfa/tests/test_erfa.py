# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from .. import _erfa as erfa


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

    astrom = np.zeros([2], dtype=erfa.dt_eraASTROM)
    theta = np.arange(0, 10.0)
    assert theta.shape == (10,)
    assert astrom.shape == (2,)

    astrom = erfa.aper(theta[:, None], astrom[None, :])
    assert astrom.shape == (10, 2)
    assert astrom.dtype == erfa.dt_eraASTROM


def test_errwarn_reporting(recwarn):
    """
    Test that the ERFA error reporting mechanism works as it should
    """

    erfa.dat(1990, 1, 1, 0.5)

    erfa.dat([100, 200, 1990], 1, 1, 0.5)
    w = recwarn.pop(erfa.ErfaWarning)
    assert '2 of "dubious year (Note 1)"' in w.message[0]

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
