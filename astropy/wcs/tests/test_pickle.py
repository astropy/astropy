# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import pickle

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

from astropy import wcs
from astropy.io import fits
from astropy.utils.data import (
    get_pkg_data_contents,
    get_pkg_data_filename,
    get_pkg_data_fileobj,
)
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.utils.misc import NumpyRNGContext
from astropy.wcs.wcs import FITSFixedWarning


def test_basic():
    wcs1 = wcs.WCS()
    s = pickle.dumps(wcs1)
    pickle.loads(s)


def test_dist():
    with get_pkg_data_fileobj(
        os.path.join("data", "dist.fits"), encoding="binary"
    ) as test_file:
        hdulist = fits.open(test_file)
        # The use of ``AXISCORR`` for D2IM correction has been deprecated
        with (
            pytest.warns(AstropyDeprecationWarning),
            pytest.warns(
                wcs.FITSFixedWarning,
                match="The WCS transformation has more axes",
            ),
        ):
            wcs1 = wcs.WCS(hdulist[0].header, hdulist)
        assert wcs1.det2im2 is not None

        s = pickle.dumps(wcs1)
        wcs2 = pickle.loads(s)

        with NumpyRNGContext(123456789):
            x = np.random.rand(2**16, wcs1.wcs.naxis)
            world1 = wcs1.all_pix2world(x, 1)
            world2 = wcs2.all_pix2world(x, 1)

        assert_array_almost_equal(world1, world2)


def test_sip():
    with get_pkg_data_fileobj(
        os.path.join("data", "sip.fits"), encoding="binary"
    ) as test_file:
        hdulist = fits.open(test_file, ignore_missing_end=True)
        with pytest.warns(FITSFixedWarning):
            wcs1 = wcs.WCS(hdulist[0].header)
        assert wcs1.sip is not None
        s = pickle.dumps(wcs1)
        wcs2 = pickle.loads(s)

        with NumpyRNGContext(123456789):
            x = np.random.rand(2**16, wcs1.wcs.naxis)
            world1 = wcs1.all_pix2world(x, 1)
            world2 = wcs2.all_pix2world(x, 1)

        assert_array_almost_equal(world1, world2)


def test_sip2():
    with get_pkg_data_fileobj(
        os.path.join("data", "sip2.fits"), encoding="binary"
    ) as test_file:
        hdulist = fits.open(test_file, ignore_missing_end=True)
        with pytest.warns(FITSFixedWarning):
            wcs1 = wcs.WCS(hdulist[0].header)
        assert wcs1.sip is not None
        s = pickle.dumps(wcs1)
        wcs2 = pickle.loads(s)

        with NumpyRNGContext(123456789):
            x = np.random.rand(2**16, wcs1.wcs.naxis)
            world1 = wcs1.all_pix2world(x, 1)
            world2 = wcs2.all_pix2world(x, 1)

        assert_array_almost_equal(world1, world2)


# Ignore "PV2_2 = 0.209028857410973 invalid keyvalue" warning seen on Windows.
@pytest.mark.filterwarnings(r"ignore:PV2_2")
def test_wcs():
    header = get_pkg_data_contents(
        os.path.join("data", "outside_sky.hdr"), encoding="binary"
    )

    wcs1 = wcs.WCS(header)
    s = pickle.dumps(wcs1)
    wcs2 = pickle.loads(s)

    with NumpyRNGContext(123456789):
        x = np.random.rand(2**16, wcs1.wcs.naxis)
        world1 = wcs1.all_pix2world(x, 1)
        world2 = wcs2.all_pix2world(x, 1)

    assert_array_almost_equal(world1, world2)


class Sub(wcs.WCS):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.foo = 42


def test_subclass():
    wcs1 = Sub()
    wcs1.foo = 45
    s = pickle.dumps(wcs1)
    wcs2 = pickle.loads(s)

    assert isinstance(wcs2, Sub)
    assert wcs1.foo == 45
    assert wcs2.foo == 45
    assert wcs2.wcs is not None


def test_axes_info():
    w = wcs.WCS(naxis=3)
    w.pixel_shape = [100, 200, 300]
    w.pixel_bounds = ((11, 22), (33, 45), (55, 67))
    w.extra = 111

    w2 = pickle.loads(pickle.dumps(w))

    # explicitly test naxis-related info
    assert w.naxis == w2.naxis
    assert w.pixel_shape == w2.pixel_shape
    assert w.pixel_bounds == w2.pixel_bounds

    # test all attributes
    for k, v in w.__dict__.items():
        assert getattr(w2, k) == v


def test_pixlist_wcs_colsel():
    """
    Test selection of a specific pixel list WCS using ``colsel``. See #11412.
    """
    hdr_file = get_pkg_data_filename("data/chandra-pixlist-wcs.hdr")
    hdr = fits.Header.fromtextfile(hdr_file)
    with pytest.warns(wcs.FITSFixedWarning):
        w0 = wcs.WCS(hdr, keysel=["image", "pixel"], colsel=[11, 12])

    with pytest.warns(wcs.FITSFixedWarning):
        w = pickle.loads(pickle.dumps(w0))

    assert w.naxis == 2
    assert list(w.wcs.ctype) == ["RA---TAN", "DEC--TAN"]
    assert np.allclose(w.wcs.crval, [229.38051931869, -58.81108068885])
    assert np.allclose(w.wcs.pc, [[1, 0], [0, 1]])
    assert np.allclose(w.wcs.cdelt, [-0.00013666666666666, 0.00013666666666666])
    assert np.allclose(w.wcs.lonpole, 180.0)


def test_alt_wcskey():
    w = wcs.WCS(key="A")
    w2 = pickle.loads(pickle.dumps(w))

    assert w2.wcs.alt == "A"
