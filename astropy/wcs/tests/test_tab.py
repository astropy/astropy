# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

import pytest
from packaging.version import Version
from astropy import wcs
from astropy.wcs import _wcs  # noqa
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

from . helper import SimModelTAB


_WCSLIB_VER = Version(_wcs.__version__)


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


@pytest.mark.skipif(
    _WCSLIB_VER < Version('7.6'),
    reason="Only in WCSLIB 7.6 a 1D -TAB axis roundtrips unless first axis"
)
def test_mixed_celest_and_1d_tab_roundtrip():
    # Tests WCS roundtripping for the case when there is one -TAB axis and
    # this axis is not the first axis. This tests a bug fixed in WCSLIB 7.6.
    filename = get_pkg_data_filename('data/tab-time-last-axis.fits')
    hdul = fits.open(filename)
    w = wcs.WCS(hdul[0].header, hdul)
    hdul.close()

    pts = np.random.random((10, 3)) * [[2047, 2047, 127]]
    assert np.allclose(pts, w.wcs_world2pix(w.wcs_pix2world(pts, 0), 0))
