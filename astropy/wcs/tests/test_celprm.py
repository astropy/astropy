# Licensed under a 3-clause BSD style license - see LICENSE.rst
from copy import copy, deepcopy
import pytest

import numpy as np
from astropy import wcs

from . helper import SimModelTAB


_WCS_UNDEFINED = 987654321.0e99


def test_celprm_init():
    # test PyCelprm_cnew
    assert wcs.WCS().wcs.cel

    # test PyCelprm_new
    assert wcs.Celprm()

    with pytest.raises(wcs.InvalidPrjParametersError):
        cel = wcs.Celprm()
        cel.set()

    # test deletion does not crash
    cel = wcs.Celprm()
    del cel


def test_celprm_copy():
    # shallow copy
    cel = wcs.Celprm()
    cel2 = copy(cel)
    cel3 = copy(cel2)
    cel.ref = [6, 8, 18, 3]
    assert (np.allclose(cel.ref, cel2.ref, atol=1e-12, rtol=0) and
            np.allclose(cel.ref, cel3.ref, atol=1e-12, rtol=0))
    del cel, cel2, cel3

    # deep copy
    cel = wcs.Celprm()
    cel2 = deepcopy(cel)
    cel.ref = [6, 8, 18, 3]
    assert not np.allclose(cel.ref, cel2.ref, atol=1e-12, rtol=0)
    del cel, cel2


def test_celprm_offset():
    cel = wcs.Celprm()
    assert not cel.offset
    cel.offset = True
    assert cel.offset


def test_celprm_prj():
    cel = wcs.Celprm()
    assert cel.prj is not None
    cel.prj.code = 'TAN'
    cel.set()
    assert cel._flag


def test_celprm_phi0():
    cel = wcs.Celprm()
    cel.prj.code = 'TAN'

    assert cel.phi0 == None
    assert cel._flag == 0
    cel.set()
    assert cel.phi0 == 0.0

    cel.phi0 = 0.0
    assert cel._flag

    cel.phi0 = 2.0
    assert cel._flag == 0

    cel.phi0 = None
    assert cel.phi0 == None
    assert cel._flag == 0


def test_celprm_theta0():
    cel = wcs.Celprm()
    cel.prj.code = 'TAN'

    assert cel.theta0 == None
    assert cel._flag == 0

    cel.theta0 = 4.0
    cel.set()
    assert cel.theta0 == 4.0

    cel.theta0 = 4.0
    assert cel._flag

    cel.theta0 = 8.0
    assert cel._flag == 0

    cel.theta0 = None
    assert cel.theta0 == None
    assert cel._flag == 0


def test_celprm_ref():
    cel = wcs.Celprm()
    cel.prj.code = 'TAN'
    cel.set()

    assert np.allclose(cel.ref, [0.0, 0.0, 180.0, 0.0], atol=1e-12, rtol=0)

    cel.phi0 = 2.0
    cel.theta0 = 4.0
    cel.ref = [123, 12]
    cel.set()

    assert np.allclose(cel.ref, [123.0, 12.0, 2, 82], atol=1e-12, rtol=0)

    cel.ref = [None, 13, None, None]
    assert np.allclose(cel.ref, [123.0, 13.0, 2, 82], atol=1e-12, rtol=0)


def test_celprm_isolat():
    cel = wcs.Celprm()
    cel.prj.code = 'TAN'
    cel.set()

    assert cel.isolat == 0


def test_celprm_latpreq():
    cel = wcs.Celprm()
    cel.prj.code = 'TAN'
    cel.set()

    assert cel.latpreq == 0


def test_celprm_euler():
    cel = wcs.Celprm()
    cel.prj.code = 'TAN'
    cel.set()

    assert np.allclose(cel.euler, [0.0, 90.0, 180.0, 0.0, 1.0], atol=1e-12, rtol=0)
