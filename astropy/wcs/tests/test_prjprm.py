# Licensed under a 3-clause BSD style license - see LICENSE.rst
from copy import copy, deepcopy
import pytest

import numpy as np
from astropy import wcs

from . helper import SimModelTAB


def test_prjprm_init():
    # test PyPrjprm_cnew
    assert wcs.WCS().wcs.cel.prj

    # test PyPrjprm_new
    assert wcs.Prjprm()

    with pytest.raises(wcs.InvalidPrjParametersError):
        prj = wcs.Prjprm()
        prj.set()

    # test deletion does not crash
    prj = wcs.Prjprm()
    del prj


def test_prjprm_copy():
    # shallow copy
    prj = wcs.Prjprm()
    prj2 = copy(prj)
    prj3 = copy(prj2)
    prj.pv = [0, 6, 8, 18, 3]
    assert (np.allclose(prj.pv, prj2.pv, atol=1e-12, rtol=0) and
            np.allclose(prj.pv, prj3.pv, atol=1e-12, rtol=0))
    del prj, prj2, prj3

    # deep copy
    prj = wcs.Prjprm()
    prj2 = deepcopy(prj)
    prj.pv = [0, 6, 8, 18, 3]
    assert not np.allclose(prj.pv, prj2.pv, atol=1e-12, rtol=0)
    del prj, prj2


def test_prjprm_flag():
    prj = wcs.Prjprm()
    assert prj._flag == 0


def test_prjprm_code():
    prj = wcs.Prjprm()

    assert prj.code == '   '
    assert prj._flag == 0

    prj.code = 'TAN'
    prj.set()
    assert prj.code == 'TAN'
    assert prj._flag

    prj.code = 'TAN'
    assert prj._flag

    prj.code = None
    assert prj.code == '   '
    assert prj._flag == 0


def test_prjprm_phi0():
    prj = wcs.Prjprm()

    assert prj.phi0 == None
    assert prj._flag == 0

    prj.code = 'TAN'
    prj.phi0 = 2.0
    prj.set()
    assert prj.phi0 == 0

    prj.phi0 = 0.0
    assert prj._flag

    prj.phi0 = 2.0
    assert prj._flag == 0

    prj.phi0 = None
    assert prj.phi0 == None
    assert prj._flag == 0


def test_prjprm_theta0():
    prj = wcs.Prjprm()

    assert prj.theta0 == None
    assert prj._flag == 0

    prj.code = 'TAN'
    prj.phi0 = 2.0
    prj.theta0 = 4.0
    prj.set()
    assert prj.theta0 == 4.0

    prj.theta0 = 4.0
    assert prj._flag

    prj.theta0 = 8.0
    assert prj._flag == 0

    prj.theta0 = None
    assert prj.theta0 == None
    assert prj._flag == 0


def test_prjprm_pv():
    prj = wcs.Prjprm()

    assert prj.pv.size == wcs.PRJ_PVN

    assert prj.theta0 == None
    assert prj._flag == 0

    prj.code = 'ZPN'
    prj.phi0 = 2.0
    prj.theta0 = 4.0
    pv = [float(i) if (i % 2) else i for i in range(wcs.PRJ_PVN)]
    prj.pv = pv
    prj.set()
    assert prj.phi0 == 2.0
    assert prj.theta0 == 4.0
    assert np.allclose(prj.pv, pv, atol=1e-6, rtol=0)

    # test the same with numpy.ndarray (flag not cleared for same values):
    prj.pv = prj.pv
    assert prj._flag != 0
    prj.set()
    assert np.allclose(prj.pv, pv, atol=1e-6, rtol=0)

    # test the same but modify PV values to check that flag is cleared:
    prj.pv = np.array(pv) + 2e-7
    assert prj._flag == 0
    prj.set()
    assert np.allclose(prj.pv, pv, atol=1e-6, rtol=0)

    # check that None in pv list leave value unchanged and values not set
    # by the list are left unchanged:
    prj.code = 'SZP'
    prj.pv = [0.0, 99.0, None]
    assert np.allclose(prj.pv[:4], [0.0, 99.0, 2.0, 3.0], atol=1e-6, rtol=0)

    # check that None resets PV to uninitialized (before prjset()) values:
    prj.pv = None
    assert prj.pv[0] == 0.0
    assert np.all(np.isnan(prj.pv[1:4]))
    assert np.allclose(prj.pv[5:], 0, atol=0, rtol=0)

    # check we can set all PVs to nan:
    nan_pvs = wcs.PRJ_PVN * [np.nan]
    prj.code = 'TAN'
    prj.pv = nan_pvs  # set using a list
    prj.set()
    assert np.all(np.isnan(prj.pv))
    prj.pv = np.array(nan_pvs)  # set using a numpy.ndarray
    prj.set()
    assert np.all(np.isnan(prj.pv))


def test_prjprm_pvrange():
    prj = wcs.Prjprm()
    prj.code = 'ZPN'
    prj.phi0 = 2.0
    prj.theta0 = 4.0
    prj.pv = [0.0, 1.0, 2.0, 3.0]
    prj.set()
    assert prj.pvrange == wcs.PRJ_PVN

    prj.code = 'SZP'
    prj.set()
    assert prj.pvrange == 103


def test_prjprm_bounds(prj_TAB):
    assert prj_TAB.bounds == 7
    prj_TAB.bounds = 0
    assert prj_TAB.bounds == 0


def test_prjprm_category(prj_TAB):
    assert prj_TAB.category == wcs.PRJ_ZENITHAL


def test_prjprm_name(prj_TAB):
    assert prj_TAB.name


def test_prjprm_w(prj_TAB):
    assert np.all(np.isfinite(prj_TAB.w))


def test_prjprm_simplezen(prj_TAB):
    assert prj_TAB.simplezen == 1


def test_prjprm_equiareal(prj_TAB):
    assert prj_TAB.equiareal == 0


def test_prjprm_conformal(prj_TAB):
    assert prj_TAB.conformal == 0


def test_prjprm_global_projection(prj_TAB):
    assert prj_TAB.global_projection == 0


def test_prjprm_divergent(prj_TAB):
    assert prj_TAB.divergent == 1


def test_prjprm_r0(prj_TAB):
    assert prj_TAB.r0 > 0.0


def test_prjprm_x0_y0(prj_TAB):
    assert prj_TAB.x0 == 0.0
    assert prj_TAB.y0 == 0.0


def test_prjprm_n_m(prj_TAB):
    assert prj_TAB.n == 0
    assert prj_TAB.m == 0


def test_prjprm_prj_roundtrips(prj_TAB):
    # some random values:
    x = [-0.002, 0.014, -0.003, 0.015, -0.047, -0.029, -0.042, 0.027, 0.021]
    y = [0.022, -0.017, -0.048, -0.049, -0.043, 0.015, 0.046, 0.031, 0.011]

    xr, yr = prj_TAB.prjs2x(*prj_TAB.prjx2s(x, y))
    assert np.allclose(x, xr, atol=1e-12, rtol=0)
    assert np.allclose(y, yr, atol=1e-12, rtol=0)

    # same test for a Prjprm that was not previously explicitly "set":
    prj = wcs.Prjprm()
    prj.code = 'TAN'
    xr, yr = prj.prjs2x(*prj.prjx2s(x, y))
    assert np.allclose(x, xr, atol=1e-12, rtol=0)
    assert np.allclose(y, yr, atol=1e-12, rtol=0)
