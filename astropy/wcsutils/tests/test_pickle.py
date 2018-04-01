# Licensed under a 3-clause BSD style license - see LICENSE.rst


import os
import pickle

import numpy as np
from numpy.testing import assert_array_almost_equal

from ...utils.data import get_pkg_data_contents, get_pkg_data_fileobj
from ...utils.misc import NumpyRNGContext
from ...io import fits
from ... import wcs


def test_basic():
    wcs1 = wcs.WCS()
    s = pickle.dumps(wcs1)
    wcs2 = pickle.loads(s)


def test_dist():
    with get_pkg_data_fileobj(
            os.path.join("data", "dist.fits"), encoding='binary') as test_file:
        hdulist = fits.open(test_file)
        wcs1 = wcs.WCS(hdulist[0].header, hdulist)
        assert wcs1.det2im2 is not None
        s = pickle.dumps(wcs1)
        wcs2 = pickle.loads(s)

        with NumpyRNGContext(123456789):
            x = np.random.rand(2 ** 16, wcs1.wcs.naxis)
            world1 = wcs1.all_pix2world(x, 1)
            world2 = wcs2.all_pix2world(x, 1)

        assert_array_almost_equal(world1, world2)


def test_sip():
    with get_pkg_data_fileobj(
            os.path.join("data", "sip.fits"), encoding='binary') as test_file:
        hdulist = fits.open(test_file, ignore_missing_end=True)
        wcs1 = wcs.WCS(hdulist[0].header)
        assert wcs1.sip is not None
        s = pickle.dumps(wcs1)
        wcs2 = pickle.loads(s)

        with NumpyRNGContext(123456789):
            x = np.random.rand(2 ** 16, wcs1.wcs.naxis)
            world1 = wcs1.all_pix2world(x, 1)
            world2 = wcs2.all_pix2world(x, 1)

        assert_array_almost_equal(world1, world2)


def test_sip2():
    with get_pkg_data_fileobj(
            os.path.join("data", "sip2.fits"), encoding='binary') as test_file:
        hdulist = fits.open(test_file, ignore_missing_end=True)
        wcs1 = wcs.WCS(hdulist[0].header)
        assert wcs1.sip is not None
        s = pickle.dumps(wcs1)
        wcs2 = pickle.loads(s)

        with NumpyRNGContext(123456789):
            x = np.random.rand(2 ** 16, wcs1.wcs.naxis)
            world1 = wcs1.all_pix2world(x, 1)
            world2 = wcs2.all_pix2world(x, 1)

        assert_array_almost_equal(world1, world2)


def test_wcs():
    header = get_pkg_data_contents(
        os.path.join("data", "outside_sky.hdr"), encoding='binary')

    wcs1 = wcs.WCS(header)
    s = pickle.dumps(wcs1)
    wcs2 = pickle.loads(s)

    with NumpyRNGContext(123456789):
        x = np.random.rand(2 ** 16, wcs1.wcs.naxis)
        world1 = wcs1.all_pix2world(x, 1)
        world2 = wcs2.all_pix2world(x, 1)

    assert_array_almost_equal(world1, world2)


class Sub(wcs.WCS):
    def __init__(self, *args, **kwargs):
        self.foo = 42


def test_subclass():
    wcs = Sub()
    s = pickle.dumps(wcs)
    wcs2 = pickle.loads(s)

    assert isinstance(wcs2, Sub)
    assert wcs.foo == 42
    assert wcs2.foo == 42
    assert wcs2.wcs is not None
