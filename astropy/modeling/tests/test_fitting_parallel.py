# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
from numpy.testing import assert_allclose

from astropy.modeling.models import Gaussian1D
from astropy.modeling.fitting_parallel import parallel_fit_model_nd
from astropy.modeling.fitting import LevMarLSQFitter


def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean) ** 2) / 2 / stddev**2)


def modify_axes(array, shape, swapaxes):
    array = array.reshape(shape)
    if swapaxes is not None:
        array = array.swapaxes(*swapaxes)
    return array


@pytest.mark.parametrize('iterating_shape,swapaxes_data,swapaxes_parameters,fitting_axes', [
    ((120,), None, None, 0),
    ((2, 3, 4, 5), None, None, 0),
    ((2, 3, 4, 5), (0, 1), None, 1),
    ((2, 3, 4, 5), (0, 2), (0, 1), 2),
])
def test_1d_model_fit_axes(iterating_shape, swapaxes_data, swapaxes_parameters, fitting_axes):

    # Test that different iterating and fitting axes work. We start off by
    # creating an array with one fitting and one iterating dimension:

    N = 120
    M = 20

    rng = np.random.default_rng(12345)

    x = np.linspace(-5, 30, M)

    amplitude = rng.uniform(1, 10, N)
    mean = rng.uniform(0, 25, N)
    stddev = rng.uniform(1, 4, N)

    data = gaussian(x[:, None], amplitude, mean, stddev)

    # Modify the axes by reshaping and swaping axes
    data = modify_axes(data, data.shape[0:1] + iterating_shape, swapaxes_data)

    # Set initial parameters to be close to but not exactly equal to true parameters

    model = Gaussian1D(
        amplitude=modify_axes(amplitude * rng.random(N), iterating_shape, swapaxes_parameters),
        mean=modify_axes(mean + rng.random(N), iterating_shape, swapaxes_parameters),
        stddev=modify_axes(stddev + rng.random(N), iterating_shape, swapaxes_parameters)
    )
    fitter = LevMarLSQFitter()

    model_fit = parallel_fit_model_nd(
        data=data,
        model=model,
        fitter=fitter,
        fitting_axes=fitting_axes,
        world={fitting_axes: x},
    )

    # Check that shape and values match

    assert_allclose(model_fit.amplitude.value, modify_axes(amplitude, iterating_shape, swapaxes_parameters))
    assert_allclose(model_fit.mean.value, modify_axes(mean, iterating_shape, swapaxes_parameters))
    assert_allclose(model_fit.stddev.value, modify_axes(stddev, iterating_shape, swapaxes_parameters))
