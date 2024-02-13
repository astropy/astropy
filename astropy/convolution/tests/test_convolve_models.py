# Licensed under a 3-clause BSD style license - see LICENSE.rst

import math

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_almost_equal

from astropy.convolution.convolve import convolve, convolve_fft, convolve_models
from astropy.modeling import fitting, models
from astropy.utils.compat.optional_deps import HAS_SCIPY
from astropy.utils.misc import NumpyRNGContext


class TestConvolve1DModels:
    @pytest.mark.parametrize("mode", [convolve_fft, convolve])
    @pytest.mark.skipif(not HAS_SCIPY, reason="Requires scipy")
    def test_is_consistency_with_astropy_convolution(self, mode):
        kernel = models.Gaussian1D(1, 0, 1)
        model = models.Gaussian1D(1, 0, 1)
        model_conv = convolve_models(
            model, kernel, mode=mode.__name__, bounding_box=(-5, 5), resolution=1
        )
        x = np.arange(-5, 6)
        ans = mode(model(x), kernel(x))

        assert_allclose(ans, model_conv(x), atol=1e-5)

    @pytest.mark.parametrize("mode", ["convolve_fft", "convolve"])
    @pytest.mark.skipif(not HAS_SCIPY, reason="Requires scipy")
    def test_against_scipy(self, mode):
        from scipy.signal import fftconvolve

        kernel = models.Gaussian1D(1, 0, 1)
        model = models.Gaussian1D(1, 0, 1)
        model_conv = convolve_models(
            model, kernel, mode=mode, bounding_box=(-5, 5), resolution=1
        )
        x = np.arange(-5, 6)
        ans = fftconvolve(kernel(x), model(x), mode="same")

        assert_allclose(ans, model_conv(x) * kernel(x).sum(), atol=1e-5)

    @pytest.mark.parametrize("mode", ["convolve_fft", "convolve"])
    @pytest.mark.skipif(not HAS_SCIPY, reason="Requires scipy")
    def test_against_scipy_with_additional_keywords(self, mode):
        from scipy.signal import fftconvolve

        kernel = models.Gaussian1D(1, 0, 1)
        model = models.Gaussian1D(1, 0, 1)
        model_conv = convolve_models(
            model,
            kernel,
            mode=mode,
            bounding_box=(-5, 5),
            resolution=1,
            normalize_kernel=False,
        )
        x = np.arange(-5, 6)
        ans = fftconvolve(kernel(x), model(x), mode="same")

        assert_allclose(ans, model_conv(x), atol=1e-5)

    @pytest.mark.parametrize("mode", ["convolve_fft", "convolve"])
    @pytest.mark.skipif(not HAS_SCIPY, reason="Requires scipy")
    def test_sum_of_gaussians(self, mode):
        """
        Test that convolving N(a, b) with N(c, d) gives N(a + c, b + d),
        where N(., .) stands for Gaussian probability density function,
        in which a and c are their means and b and d are their variances.
        This test also confirms that the compound model is not reliant
        on the input domain being symmetric about x=0.
        """

        kernel = models.Gaussian1D(1 / math.sqrt(2 * np.pi), 1, 1)
        model = models.Gaussian1D(1 / math.sqrt(2 * np.pi), 3, 1)
        model_conv = convolve_models(
            model,
            kernel,
            mode=mode,
            bounding_box=(-7, 7),
            resolution=1,
            normalize_kernel=False,
        )
        ans = models.Gaussian1D(1 / (2 * math.sqrt(np.pi)), 4, np.sqrt(2))
        x = np.arange(-5, 6)

        assert_allclose(ans(x), model_conv(x), atol=1e-3)

    @pytest.mark.parametrize("mode", ["convolve_fft", "convolve"])
    @pytest.mark.skipif(not HAS_SCIPY, reason="Requires scipy")
    def test_convolve_box_models(self, mode):
        kernel = models.Box1D()
        model = models.Box1D()
        model_conv = convolve_models(
            model,
            kernel,
            mode=mode,
            bounding_box=(-1, 1),
            resolution=0.001,
            normalize_kernel=True,
        )
        x = np.linspace(-1, 1, 99)
        ans = (x + 1) * (x < 0) + (-x + 1) * (x >= 0)

        assert_allclose(ans, model_conv(x), atol=1e-3)

    @pytest.mark.parametrize("mode", ["convolve_fft", "convolve"])
    @pytest.mark.skipif(not HAS_SCIPY, reason="Requires scipy")
    def test_fitting_convolve_models(self, mode):
        """
        test that a convolve model can be fitted
        """
        b1 = models.Box1D()
        g1 = models.Gaussian1D()

        x = np.linspace(-5, 5, 99)
        fake_model = models.Gaussian1D(amplitude=10)
        with NumpyRNGContext(123):
            fake_data = fake_model(x) + np.random.normal(size=len(x))

        init_model = convolve_models(
            b1,
            g1,
            mode=mode,
            bounding_box=(-6 - 10 / 99, 6),
            resolution=10 / 99,
            normalize_kernel=False,
            cache=False,
        )
        fitter = fitting.LevMarLSQFitter()
        fitted_model = fitter(init_model, x, fake_data)

        me = np.mean(fitted_model(x) - fake_data)
        assert_almost_equal(me, 0.0, decimal=2)

    @pytest.mark.parametrize("bounding_box, resolution", [(None, 1), ((-1, 1), None)])
    @pytest.mark.skipif(not HAS_SCIPY, reason="Requires scipy")
    def test_convolve_models_warnings(self, bounding_box, resolution):
        """
        test that convolve models issues warnings when bounding_box or resolution are not specified
        """
        model = models.Gaussian1D()
        kernel = models.Box1D()
        with pytest.warns(UserWarning):
            combined_model = convolve_models(
                model, kernel, bounding_box=bounding_box, resolution=resolution
            )
