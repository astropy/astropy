# Licensed under a 3-clause BSD style license - see LICENSE.rst
# pylint: disable=invalid-name
import numpy as np
import pytest

from astropy.convolution import convolve_models_fft
from astropy.modeling.models import Const1D, Const2D
from astropy.utils.compat.optional_deps import HAS_SCIPY  # noqa: F401


@pytest.mark.skipif('not HAS_SCIPY')
def test_clear_cache():
    m1 = Const1D()
    m2 = Const1D()

    model = convolve_models_fft(m1, m2, (-1, 1), 0.01)
    assert model._kwargs is None
    assert model._convolution is None

    results = model(0)
    assert results.all() == np.array([1.]).all()
    assert model._kwargs is not None
    assert model._convolution is not None

    model.clear_cache()
    assert model._kwargs is None
    assert model._convolution is None


@pytest.mark.skipif('not HAS_SCIPY')
def test_input_shape_1d():
    m1 = Const1D()
    m2 = Const1D()

    model = convolve_models_fft(m1, m2, (-1, 1), 0.01)

    results = model(0)
    assert results.shape == (1,)

    x = np.arange(-1, 1, 0.1)
    results = model(x)
    assert results.shape == x.shape


@pytest.mark.skipif('not HAS_SCIPY')
def test_input_shape_2d():
    m1 = Const2D()
    m2 = Const2D()

    model = convolve_models_fft(m1, m2, ((-1, 1), (-1, 1)), 0.01)

    results = model(0, 0)
    assert results.shape == (1,)

    x = np.arange(-1, 1, 0.1)
    results = model(x, 0)
    assert results.shape == x.shape
    results = model(0, x)
    assert results.shape == x.shape

    grid = np.meshgrid(x, x)
    results = model(*grid)
    assert results.shape == grid[0].shape
    assert results.shape == grid[1].shape


@pytest.mark.skipif('not HAS_SCIPY')
def test__convolution_inputs():
    m1 = Const2D()
    m2 = Const2D()

    model = convolve_models_fft(m1, m2, ((-1, 1), (-1, 1)), 0.01)

    x = np.arange(-1, 1, 0.1)
    y = np.arange(-2, 2, 0.1)
    grid0 = np.meshgrid(x, x)
    grid1 = np.meshgrid(y, y)

    # scalar inputs
    assert (np.array([1]), (1,)) == model._convolution_inputs(1)

    # Multiple inputs
    assert np.all(model._convolution_inputs(*grid0)[0] ==
                  np.reshape([grid0[0], grid0[1]], (2, -1)).T)
    assert model._convolution_inputs(*grid0)[1] == grid0[0].shape
    assert np.all(model._convolution_inputs(*grid1)[0] ==
                  np.reshape([grid1[0], grid1[1]], (2, -1)).T)
    assert model._convolution_inputs(*grid1)[1] == grid1[0].shape

    # Error
    with pytest.raises(ValueError) as err:
        model._convolution_inputs(grid0[0], grid1[1])
    assert str(err.value) == "Values have differing shapes"
