# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from numpy.testing import assert_allclose

from ..ui import scale_image

DATA = np.array([0, 1., 2.])
DATASCL = 0.5 * DATA


class TestImageScaling(object):

    def test_linear(self):
        """Test linear scaling."""
        img = scale_image(DATA, scale='linear')
        assert_allclose(img, DATASCL, atol=0, rtol=1.e-5)

    def test_sqrt(self):
        """Test sqrt scaling."""
        img = scale_image(DATA, scale='sqrt')
        assert_allclose(img, np.sqrt(DATASCL), atol=0, rtol=1.e-5)

    def test_power(self):
        """Test power scaling."""
        power = 3.0
        img = scale_image(DATA, scale='power', power=power)
        assert_allclose(img, DATASCL**power, atol=0, rtol=1.e-5)

    def test_log(self):
        """Test log10 scaling."""
        img = scale_image(DATA, scale='log')
        ref = np.log10(1000 * DATASCL + 1.0) / np.log10(1001.0)
        assert_allclose(img, ref, atol=0, rtol=1.e-5)

    def test_asinh(self):
        """Test asinh scaling."""
        a = 0.1
        img = scale_image(DATA, scale='asinh', asinh_a=a)
        ref = np.arcsinh(DATASCL / a) / np.arcsinh(1. / a)
        assert_allclose(img, ref, atol=0, rtol=1.e-5)
