# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from astropy import wcs

from .helper import SimModelTAB


@pytest.fixture(scope="module")
def tab_wcs_2di():
    model = SimModelTAB(nx=150, ny=200)

    # generate FITS HDU list:
    hdulist = model.hdulist

    # create WCS object:
    w = wcs.WCS(hdulist[0].header, hdulist)

    return w


@pytest.fixture(scope="module")
def tab_wcsh_2di():
    model = SimModelTAB(nx=150, ny=200)

    # generate FITS HDU list:
    hdulist = model.hdulist

    # create WCS object:
    w = wcs.WCS(hdulist[0].header, hdulist)

    return w, hdulist


@pytest.fixture(scope="function")
def tab_wcs_2di_f():
    model = SimModelTAB(nx=150, ny=200)

    # generate FITS HDU list:
    hdulist = model.hdulist

    # create WCS object:
    w = wcs.WCS(hdulist[0].header, hdulist)

    return w


@pytest.fixture(scope="function")
def prj_TAB():
    prj = wcs.Prjprm()
    prj.code = "TAN"
    prj.set()
    return prj
