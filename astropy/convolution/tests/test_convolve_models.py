# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from ..convolve import convolve, convolve_models
from ...modeling import models
from ...tests.helper import pytest
from numpy.testing import assert_allclose

try:
    import scipy
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True

class TestConvolve1DModels(object):
    @pytest.mark.skipif('not HAS_SCIPY')
    def test_against_scipy(self):
        from scipy.signal import fftconvolve

        g1 = models.Gaussian1D(1, 1, 1)
        g2 = models.Gaussian1D(1, 3, 1)
        x = np.arange(-5, 6)
        model_conv = convolve_models(g1, g2)
        ans = fftconvolve(g1(x), g2(x), 'same')

        assert_allclose(ans, model_conv(x), atol=1e-5)

    def test_sum_of_gaussians(self):
        '''
        Test that convolving N(a, b) with N(c, d) gives N(a + b, c + d),
        where N(.|.) stands for Gaussian probability density function.
        '''

        g1 = models.Gaussian1D(1, 1, 1)
        g2 = models.Gaussian1D(1, 3, 1)
        g3 = convolve_models(g1, g2)
        ans = models.Gaussian1D(1, 4, np.sqrt(2))

        x = np.arange(-5, 6)
        assert_allclose(g3(x) / (2 * np.pi), ans(x) / np.sqrt(2 * np.pi * 2),
                        atol=1e-3)
