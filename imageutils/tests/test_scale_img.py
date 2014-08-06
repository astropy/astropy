# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
from astropy.tests.helper import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_equal
from .. import scale_img
try:
    import skimage
    HAS_SKIMAGE = True
except ImportError:
    HAS_SKIMAGE = False


DATA = np.array([0, 1., 2.])
DATASCL = 0.5 * DATA


class TestImgCuts(object):
    def test_find_cutlevels(self):
        data = np.arange(101)

        mincut, maxcut = scale_img.find_cutlevels(data, min_percent=20)
        assert_equal([mincut, maxcut], [20, 100])

        mincut, maxcut = scale_img.find_cutlevels(data, max_percent=80)
        assert_equal([mincut, maxcut], [0, 80])

        mincut, maxcut = scale_img.find_cutlevels(data, min_percent=20,
                                                  max_percent=80)
        assert_equal([mincut, maxcut], [20, 80])

        mincut, maxcut = scale_img.find_cutlevels(data, percent=90)
        assert_equal([mincut, maxcut], [5, 95])

        mincut, maxcut = scale_img.find_cutlevels(data, min_percent=20,
                                                  max_percent=80, percent=90)
        assert_equal([mincut, maxcut], [20, 80])


@pytest.mark.skipif('not HAS_SKIMAGE')
class TestImageScaling(object):
    def test_linear(self):
        """Test linear scaling."""
        img = scale_img.scale_image(DATA, scale='linear')
        assert_allclose(img, DATASCL, atol=0, rtol=1.e-5)

    def test_sqrt(self):
        """Test sqrt scaling."""
        img = scale_img.scale_image(DATA, scale='sqrt')
        assert_allclose(img, np.sqrt(DATASCL), atol=0, rtol=1.e-5)

    def test_power(self):
        """Test power scaling."""
        power = 3.0
        img = scale_img.scale_image(DATA, scale='power', power=power)
        assert_allclose(img, DATASCL**power, atol=0, rtol=1.e-5)

    def test_log(self):
        """Test log10 scaling."""
        img = scale_img.scale_image(DATA, scale='log')
        ref = np.log10(DATASCL + 1.0) / np.log10(2.0)
        assert_allclose(img, ref, atol=0, rtol=1.e-5)

    def test_asinh(self):
        """Test arcsinh scaling."""
        img = scale_img.scale_image(DATA, scale='asinh')
        z = 0.658248290464
        ref = np.arcsinh(DATASCL / z) / np.arcsinh(1.0 / z)
        assert_allclose(img, ref, atol=0, rtol=1.e-5)

    def test_asinh_noiselevel(self):
        """Test arcsinh scaling."""
        img = scale_img.scale_image(DATA, scale='asinh', noise_level=1.0)
        z = 0.5
        ref = np.arcsinh(DATASCL / z) / np.arcsinh(1.0 / z)
        assert_allclose(img, ref, atol=0, rtol=1.e-5)

    def test_asinh_noiselevel_zeroz(self):
        """Test arcsinh scaling."""
        img = scale_img.scale_image(DATA, scale='asinh', noise_level=0.0)
        z = 1.e-2
        ref = np.arcsinh(DATASCL / z) / np.arcsinh(1.0 / z)
        assert_allclose(img, ref, atol=0, rtol=1.e-5)
