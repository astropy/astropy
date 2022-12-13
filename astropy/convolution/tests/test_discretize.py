# Licensed under a 3-clause BSD style license - see LICENSE.rst

import itertools

import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy.convolution.utils import discretize_model
from astropy.modeling.functional_models import (
    Box1D,
    Box2D,
    Gaussian1D,
    Gaussian2D,
    RickerWavelet1D,
    RickerWavelet2D,
)
from astropy.modeling.tests.example_models import models_1D, models_2D
from astropy.modeling.tests.test_models import create_model
from astropy.utils.compat.optional_deps import HAS_SCIPY

modes = ["center", "linear_interp", "oversample"]
test_models_1D = [Gaussian1D, Box1D, RickerWavelet1D]
test_models_2D = [Gaussian2D, Box2D, RickerWavelet2D]


@pytest.mark.parametrize(
    ("model_class", "mode"), list(itertools.product(test_models_1D, modes))
)
def test_pixel_sum_1D(model_class, mode):
    """
    Test if the sum of all pixels corresponds nearly to the integral.
    """
    if model_class == Box1D and mode == "center":
        pytest.skip("Non integrating mode. Skip integral test.")
    parameters = models_1D[model_class]
    model = create_model(model_class, parameters)

    values = discretize_model(model, models_1D[model_class]["x_lim"], mode=mode)
    assert_allclose(values.sum(), models_1D[model_class]["integral"], atol=0.0001)


@pytest.mark.parametrize("mode", modes)
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


@pytest.mark.parametrize(
    ("model_class", "mode"), list(itertools.product(test_models_2D, modes))
)
def test_pixel_sum_2D(model_class, mode):
    """
    Test if the sum of all pixels corresponds nearly to the integral.
    """
    if model_class == Box2D and mode == "center":
        pytest.skip("Non integrating mode. Skip integral test.")

    parameters = models_2D[model_class]
    model = create_model(model_class, parameters)

    values = discretize_model(
        model,
        models_2D[model_class]["x_lim"],
        models_2D[model_class]["y_lim"],
        mode=mode,
    )
    assert_allclose(values.sum(), models_2D[model_class]["integral"], atol=0.0001)


@pytest.mark.parametrize(
    ("model_class", "mode"), list(itertools.product(test_models_2D, modes))
)
def test_pixel_sum_compound_2D(model_class, mode):
    """
    Test if the sum of all pixels of a compound model corresponds nearly to the integral.
    """
    if model_class == Box2D and mode == "center":
        pytest.skip("Non integrating mode. Skip integral test.")

    parameters = models_2D[model_class]
    model = create_model(model_class, parameters)

    values = discretize_model(
        model + model,
        models_2D[model_class]["x_lim"],
        models_2D[model_class]["y_lim"],
        mode=mode,
    )

    model_integral = 2 * models_2D[model_class]["integral"]
    assert_allclose(values.sum(), model_integral, atol=0.0001)


@pytest.mark.parametrize("mode", modes)
def test_gaussian_eval_2D(mode):
    """
    Discretize Gaussian with different modes and check
    if result is at least similar to Gaussian2D.eval()
    """
    model = Gaussian2D(0.01, 0, 0, 1, 1)

    x = np.arange(-2, 3)
    y = np.arange(-2, 3)

    x, y = np.meshgrid(x, y)

    values = model(x, y)
    disc_values = discretize_model(model, (-2, 3), (-2, 3), mode=mode)
    assert_allclose(values, disc_values, atol=1e-2)


@pytest.mark.skipif(not HAS_SCIPY, reason="Requires scipy")
@pytest.mark.slow
def test_gaussian_eval_2D_integrate_mode():
    """
    Discretize Gaussian with integrate mode
    """
    model_list = [
        Gaussian2D(0.01, 0, 0, 2, 2),
        Gaussian2D(0.01, 0, 0, 1, 2),
        Gaussian2D(0.01, 0, 0, 2, 1),
    ]

    x = np.arange(-2, 3)
    y = np.arange(-2, 3)

    x, y = np.meshgrid(x, y)

    for model in model_list:
        values = model(x, y)
        disc_values = discretize_model(model, (-2, 3), (-2, 3), mode="integrate")
        assert_allclose(values, disc_values, atol=1e-2)


@pytest.mark.skipif(not HAS_SCIPY, reason="Requires scipy")
def test_subpixel_gauss_1D():
    """
    Test subpixel accuracy of the integrate mode with gaussian 1D model.
    """
    gauss_1D = Gaussian1D(1, 0, 0.1)
    values = discretize_model(gauss_1D, (-1, 2), mode="integrate", factor=100)
    assert_allclose(values.sum(), np.sqrt(2 * np.pi) * 0.1, atol=0.00001)


@pytest.mark.skipif(not HAS_SCIPY, reason="Requires scipy")
def test_subpixel_gauss_2D():
    """
    Test subpixel accuracy of the integrate mode with gaussian 2D model.
    """
    gauss_2D = Gaussian2D(1, 0, 0, 0.1, 0.1)
    values = discretize_model(gauss_2D, (-1, 2), (-1, 2), mode="integrate", factor=100)
    assert_allclose(values.sum(), 2 * np.pi * 0.01, atol=0.00001)


def test_discretize_callable_1d():
    """
    Test discretize when a 1d function is passed.
    """

    def f(x):
        return x**2

    y = discretize_model(f, (-5, 6))
    assert_allclose(y, np.arange(-5, 6) ** 2)


def test_discretize_callable_2d():
    """
    Test discretize when a 2d function is passed.
    """

    def f(x, y):
        return x**2 + y**2

    actual = discretize_model(f, (-5, 6), (-5, 6))
    y, x = np.indices((11, 11)) - 5
    desired = x**2 + y**2
    assert_allclose(actual, desired)


def test_type_exception():
    """
    Test type exception.
    """
    with pytest.raises(TypeError, match=r"Model must be callable\."):
        discretize_model(float(0), (-10, 11))


def test_dim_exception_1d():
    """
    Test dimension exception 1d.
    """

    def f(x):
        return x**2

    with pytest.raises(ValueError, match=r"y range specified, but model is only 1-d\."):
        discretize_model(f, (-10, 11), (-10, 11))


def test_dim_exception_2d():
    """
    Test dimension exception 2d.
    """

    def f(x, y):
        return x**2 + y**2

    with pytest.raises(ValueError, match=r"y range not specified, but model is 2-d"):
        discretize_model(f, (-10, 11))


def test_float_x_range_exception():
    def f(x, y):
        return x**2 + y**2

    with pytest.raises(
        ValueError,
        match=(
            r"The difference between the upper and lower limit of 'x_range' must be a"
            r" whole number\."
        ),
    ):
        discretize_model(f, (-10.002, 11.23))


def test_float_y_range_exception():
    def f(x, y):
        return x**2 + y**2

    with pytest.raises(
        ValueError,
        match=(
            r"The difference between the upper and lower limit of 'y_range' must be a"
            r" whole number\."
        ),
    ):
        discretize_model(f, (-10, 11), (-10.002, 11.23))


def test_discretize_oversample():
    gauss_2D = Gaussian2D(
        amplitude=1.0, x_mean=5.0, y_mean=125.0, x_stddev=0.75, y_stddev=3
    )
    values = discretize_model(
        gauss_2D, x_range=[0, 10], y_range=[100, 135], mode="oversample", factor=10
    )
    vmax = np.max(values)
    vmax_yx = np.unravel_index(values.argmax(), values.shape)
    values_osf1 = discretize_model(
        gauss_2D, x_range=[0, 10], y_range=[100, 135], mode="oversample", factor=1
    )
    values_center = discretize_model(
        gauss_2D, x_range=[0, 10], y_range=[100, 135], mode="center"
    )
    assert values.shape == (35, 10)
    assert_allclose(vmax, 0.927, atol=1e-3)
    assert vmax_yx == (25, 5)
    assert_allclose(values_center, values_osf1)
