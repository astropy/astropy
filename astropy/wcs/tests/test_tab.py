# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

import numpy as np

from astropy import wcs

from . helper import SimModelTAB


def test_2d_spatial_tab_roundtrip(tab_wcs_2di):
    nx, ny = tab_wcs_2di.pixel_shape
    # generate "random" test coordinates:
    np.random.seed(1)
    xy = 0.51 + [nx + 0.99, ny + 0.99] * np.random.random((100, 2))
    rd = tab_wcs_2di.wcs_pix2world(xy, 1)
    xy_roundtripped = tab_wcs_2di.wcs_world2pix(rd, 1)
    m = np.logical_and(*(np.isfinite(xy_roundtripped).T))
    assert np.allclose(xy[m], xy_roundtripped[m], rtol=0, atol=1e-7)


def test_2d_spatial_tab_vs_model():
    nx = 150
    ny = 200
    model = SimModelTAB(nx=nx, ny=ny)

    # generate FITS HDU list:
    hdulist = model.hdulist

    # create WCS object:
    w = wcs.WCS(hdulist[0].header, hdulist)

    # generate "random" test coordinates:
    np.random.seed(1)
    xy = 0.51 + [nx + 0.99, ny + 0.99] * np.random.random((100, 2))
    rd = w.wcs_pix2world(xy, 1)
    rd_model = model.fwd_eval(xy)
    assert np.allclose(rd, rd_model, rtol=0, atol=1e-7)
