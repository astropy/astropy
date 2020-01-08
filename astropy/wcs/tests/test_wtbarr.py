# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

import numpy as np

from astropy import wcs


def test_wtbarr_i(tab_wcs_2di):
    assert tab_wcs_2di.wcs.wtb[0].i == 1


def test_wtbarr_m(tab_wcs_2di):
    assert tab_wcs_2di.wcs.wtb[0].m == 1


def test_wtbarr_m(tab_wcs_2di):
    assert tab_wcs_2di.wcs.wtb[0].kind == 'c'


def test_wtbarr_extnam(tab_wcs_2di):
    assert tab_wcs_2di.wcs.wtb[0].extnam == 'WCS-TABLE'


def test_wtbarr_extver(tab_wcs_2di):
    assert tab_wcs_2di.wcs.wtb[0].extver == 1


def test_wtbarr_extlev(tab_wcs_2di):
    assert tab_wcs_2di.wcs.wtb[0].extlev == 1


def test_wtbarr_ttype(tab_wcs_2di):
    assert tab_wcs_2di.wcs.wtb[0].ttype == 'wavelength'


def test_wtbarr_row(tab_wcs_2di):
    assert tab_wcs_2di.wcs.wtb[0].row == 1


def test_wtbarr_ndim(tab_wcs_2di):
    assert tab_wcs_2di.wcs.wtb[0].ndim == 3


def test_wtbarr_print(tab_wcs_2di, capfd):
    tab_wcs_2di.wcs.wtb[0].print_contents()
    captured = capfd.readouterr()
    s = str(tab_wcs_2di.wcs.wtb[0])
    lines = s.split('\n')
    assert captured.out == s
    assert '     i: 1' == lines[0]
    assert '     m: 1' == lines[1]
    assert '  kind: c' == lines[2]
    assert 'extnam: WCS-TABLE' == lines[3]
    assert 'extver: 1' == lines[4]
    assert 'extlev: 1' == lines[5]
    assert ' ttype: wavelength' == lines[6]
    assert '   row: 1' == lines[7]
    assert '  ndim: 3' == lines[8]
    assert lines[9].startswith('dimlen: ')
    assert '        0:   4' == lines[10]
    assert '        1:   2' == lines[11]
    assert lines[12].startswith('arrayp: ')
