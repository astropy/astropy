# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import itertools

import numpy as np
from numpy.testing import assert_allclose

from ...tests.helper import pytest

from ..utils import discretize_model
from ...modeling.functional_models import (
    Gaussian1D, Box1D, MexicanHat1D, Trapezoid1D,
    Gaussian2D, Box2D, MexicanHat2D)
from ...modeling.tests.example_models import models_1D, models_2D
from ...modeling.tests.test_models import create_model

try:
    import scipy
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


modes = ['center', 'linear_interp', 'oversample']
test_models_1D = [Gaussian1D, Box1D, MexicanHat1D]
test_models_2D = [Gaussian2D, Box2D, MexicanHat2D]


@pytest.mark.parametrize(('model_class', 'mode'), list(itertools.product(test_models_1D, modes)))
def test_pixel_sum_1D(model_class, mode):
    """
    Test if the sum of all pixels corresponds nearly to the integral.
    """
    if model_class == Box1D and mode == "center":
        pytest.skip("Non integrating mode. Skip integral test.")
    parameters = models_1D[model_class]
    model = create_model(model_class, parameters)

    values = discretize_model(model, models_1D[model_class]['x_lim'], mode=mode)
    assert_allclose(values.sum(), models_1D[model_class]['integral'], atol=0.0001)


@pytest.mark.parametrize('mode', modes)
def test_gaussian_eval_1D(mode):
    """
    Discretize Gaussian with different modes and check
    if result is at least similar to Gaussian1D.eval().
    """
    model = Gaussian1D(1, 0, 20)
    x = np.arange(-100, 101)
    values = model(x)
    disc_values = discretize_model(model, (-100, 101), mode=mode)
    assert_allclose(values, disc_values, atol=0.001)


@pytest.mark.parametrize(('model_class', 'mode'), list(itertools.product(test_models_2D, modes)))
def test_pixel_sum_2D(model_class, mode):
    """
    Test if the sum of all pixels corresponds nearly to the integral.
    """
    if model_class == Box2D and mode == "center":
        pytest.skip("Non integrating mode. Skip integral test.")

    parameters = models_2D[model_class]
    model = create_model(model_class, parameters)

    values = discretize_model(model, models_2D[model_class]['x_lim'],
                              models_2D[model_class]['y_lim'], mode=mode)
    assert_allclose(values.sum(), models_2D[model_class]['integral'], atol=0.0001)


@pytest.mark.parametrize('mode', modes)
def test_gaussian_eval_2D(mode):
    """
    Discretize Gaussian with different modes and check
    if result is at least similar to Gaussian1D.eval()
    """
    model = Gaussian2D(1, 0, 0, 20, 20)
    x = np.arange(-100, 101)
    y = np.arange(-100, 101)
    x, y = np.meshgrid(x, y)
    values = model(x, y)
    disc_values = discretize_model(model, (-100, 101), (-100, 101), mode=mode)
    assert_allclose(values, disc_values, atol=0.001)


@pytest.mark.skipif('not HAS_SCIPY')
def test_subpixel_gauss_1D():
    """
    Test subpixel accuracy of the oversample mode with gaussian 1D model.
    """
    gauss_1D = Gaussian1D(1, 0, 0.1)
    values = discretize_model(gauss_1D, (-1, 2), mode='integrate', factor=100)
    assert_allclose(values.sum(), np.sqrt(2 * np.pi) * 0.1, atol=0.00001)


@pytest.mark.skipif('not HAS_SCIPY')
def test_subpixel_gauss_2D():
    """
    Test subpixel accuracy of the oversample mode with gaussian 2D model.
    """
    gauss_2D = Gaussian2D(1, 0, 0, 0.1, 0.1)
    values = discretize_model(gauss_2D, (-1, 2), (-1, 2), mode='integrate', factor=100)
    assert_allclose(values.sum(), 2 * np.pi * 0.01, atol=0.00001)
