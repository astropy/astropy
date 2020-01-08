# Licensed under a 3-clause BSD style license - see LICENSE.rst

from copy import deepcopy

import pytest

import numpy as np

from astropy import wcs

from . helper import SimModelTAB


def test_wcsprm_tab_basic(tab_wcs_2di):
    assert len(tab_wcs_2di.wcs.tab) == 1
    t = tab_wcs_2di.wcs.tab[0]
    assert tab_wcs_2di.wcs.tab[0] is not t


def test_tabprm_coord(tab_wcs_2di_f):
    t = tab_wcs_2di_f.wcs.tab[0]

    c0 = t.coord
    c1 = np.ones_like(c0)
    t.coord = c1
    assert np.allclose(tab_wcs_2di_f.wcs.tab[0].coord, c1)


def test_tabprm_crval_and_deepcopy(tab_wcs_2di_f):
    w = deepcopy(tab_wcs_2di_f)

    t = tab_wcs_2di_f.wcs.tab[0]

    pix = np.array([[2, 3]], dtype=np.float32)

    rd1 = tab_wcs_2di_f.wcs_pix2world(pix, 1)

    c = t.crval.copy()
    d = 0.5 * np.ones_like(c)
    t.crval += d
    assert np.allclose(tab_wcs_2di_f.wcs.tab[0].crval, c + d)

    rd2 = tab_wcs_2di_f.wcs_pix2world(pix - d, 1)
    assert np.allclose(rd1, rd2)

    rd3 = w.wcs_pix2world(pix, 1)
    assert np.allclose(rd1, rd3)


def test_tabprm_delta(tab_wcs_2di):
    t = tab_wcs_2di.wcs.tab[0]
    assert np.allclose([0.0, 0.0], t.delta)


def test_tabprm_K(tab_wcs_2di):
    t = tab_wcs_2di.wcs.tab[0]
    assert np.all(t.K == [4, 2])


def test_tabprm_M(tab_wcs_2di):
    t = tab_wcs_2di.wcs.tab[0]
    assert t.M == 2


def test_tabprm_nc(tab_wcs_2di):
    t = tab_wcs_2di.wcs.tab[0]
    assert t.nc == 8


def test_tabprm_extrema(tab_wcs_2di):
    t = tab_wcs_2di.wcs.tab[0]
    extrema = np.array(
        [[[-0.0026, -0.5], [1.001, -0.5]],
         [[-0.0026, 0.5], [1.001, 0.5]]]
    )
    assert np.allclose(t.extrema, extrema)


def test_tabprm_map(tab_wcs_2di_f):
    t = tab_wcs_2di_f.wcs.tab[0]
    assert np.allclose(t.map, [0, 1])

    t.map[1] = 5
    assert np.all(tab_wcs_2di_f.wcs.tab[0].map == [0, 5])

    t.map = [1, 4]
    assert np.all(tab_wcs_2di_f.wcs.tab[0].map == [1, 4])


def  test_tabprm_sense(tab_wcs_2di):
    t = tab_wcs_2di.wcs.tab[0]
    assert np.all(t.sense == [1, 1])


def  test_tabprm_p0(tab_wcs_2di):
    t = tab_wcs_2di.wcs.tab[0]
    assert np.all(t.p0 == [0, 0])


def test_tabprm_print(tab_wcs_2di_f, capfd):
    tab_wcs_2di_f.wcs.tab[0].print_contents()
    captured = capfd.readouterr()
    s = str(tab_wcs_2di_f.wcs.tab[0])
    out = str(captured.out)
    lout= out.split('\n')
    assert out == s
    assert lout[0] == '       flag: 137'
    assert lout[1] == '          M: 2'


def test_wcstab_copy(tab_wcs_2di_f):
    t = tab_wcs_2di_f.wcs.tab[0]

    c0 = t.coord
    c1 = np.ones_like(c0)
    t.coord = c1
    assert np.allclose(tab_wcs_2di_f.wcs.tab[0].coord, c1)


def test_tabprm_crval(tab_wcs_2di_f):
    w = deepcopy(tab_wcs_2di_f)

    t = tab_wcs_2di_f.wcs.tab[0]

    pix = np.array([[2, 3]], dtype=np.float32)

    rd1 = tab_wcs_2di_f.wcs_pix2world(pix, 1)

    c = t.crval.copy()
    d = 0.5 * np.ones_like(c)
    t.crval += d
    assert np.allclose(tab_wcs_2di_f.wcs.tab[0].crval, c + d)

    rd2 = tab_wcs_2di_f.wcs_pix2world(pix - d, 1)
    assert np.allclose(rd1, rd2)

    rd3 = w.wcs_pix2world(pix, 1)
    assert np.allclose(rd1, rd3)
