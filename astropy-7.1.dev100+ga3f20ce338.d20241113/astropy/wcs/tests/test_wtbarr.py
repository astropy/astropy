# Licensed under a 3-clause BSD style license - see LICENSE.rst


def test_wtbarr_i(tab_wcs_2di):
    assert tab_wcs_2di.wcs.wtb[0].i == 1


def test_wtbarr_m(tab_wcs_2di):
    assert tab_wcs_2di.wcs.wtb[0].m == 1


def test_wtbarr_kind(tab_wcs_2di):
    assert tab_wcs_2di.wcs.wtb[0].kind == "c"


def test_wtbarr_extnam(tab_wcs_2di):
    assert tab_wcs_2di.wcs.wtb[0].extnam == "WCS-TABLE"


def test_wtbarr_extver(tab_wcs_2di):
    assert tab_wcs_2di.wcs.wtb[0].extver == 1


def test_wtbarr_extlev(tab_wcs_2di):
    assert tab_wcs_2di.wcs.wtb[0].extlev == 1


def test_wtbarr_ttype(tab_wcs_2di):
    assert tab_wcs_2di.wcs.wtb[0].ttype == "wavelength"


def test_wtbarr_row(tab_wcs_2di):
    assert tab_wcs_2di.wcs.wtb[0].row == 1


def test_wtbarr_ndim(tab_wcs_2di):
    assert tab_wcs_2di.wcs.wtb[0].ndim == 3


def test_wtbarr_print(tab_wcs_2di, capfd):
    tab_wcs_2di.wcs.wtb[0].print_contents()
    captured = capfd.readouterr()
    s = str(tab_wcs_2di.wcs.wtb[0])
    lines = s.split("\n")
    assert captured.out == s
    assert lines[0] == "     i: 1"
    assert lines[1] == "     m: 1"
    assert lines[2] == "  kind: c"
    assert lines[3] == "extnam: WCS-TABLE"
    assert lines[4] == "extver: 1"
    assert lines[5] == "extlev: 1"
    assert lines[6] == " ttype: wavelength"
    assert lines[7] == "   row: 1"
    assert lines[8] == "  ndim: 3"
    assert lines[9].startswith("dimlen: ")
    assert lines[10] == "        0:   4"
    assert lines[11] == "        1:   2"
    assert lines[12].startswith("arrayp: ")
