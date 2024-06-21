# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import re


import pytest
pytest.importorskip("dask")

import numpy as np
from numpy.testing import assert_allclose

from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.fitting_parallel import parallel_fit_model_nd
from astropy.modeling.models import Const1D, Gaussian1D, Linear1D, Planar2D
from astropy.wcs import WCS


def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean) ** 2) / 2 / stddev**2)


def modify_axes(array, shape, swapaxes):
    array = array.reshape(shape)
    if swapaxes is not None:
        array = array.swapaxes(*swapaxes)
    return array


@pytest.mark.parametrize(
    "iterating_shape,swapaxes_data,swapaxes_parameters,fitting_axes",
    [
        ((120,), None, None, 0),
        ((2, 3, 4, 5), None, None, 0),
        ((2, 3, 4, 5), (0, 1), None, 1),
        ((2, 3, 4, 5), (0, 2), (0, 1), 2),
    ],
)
def test_1d_model_fit_axes(
    iterating_shape, swapaxes_data, swapaxes_parameters, fitting_axes
):
    # Test that different iterating and fitting axes work. We start off by
    # creating an array with one fitting and one iterating dimension and we then
    # modify the array dimensionality and order of the dimensions to test
    # different scenarios.

    N = 120
    P = 20

    rng = np.random.default_rng(12345)

    x = np.linspace(-5, 30, P)

    amplitude = rng.uniform(1, 10, N)
    mean = rng.uniform(0, 25, N)
    stddev = rng.uniform(1, 4, N)

    data = gaussian(x[:, None], amplitude, mean, stddev)

    # At this point, the data has shape (P, N)

    # Modify the axes by reshaping and swapping axes
    data = modify_axes(data, data.shape[0:1] + iterating_shape, swapaxes_data)

    # Set initial parameters to be close to but not exactly equal to true parameters

    model = Gaussian1D(
        amplitude=modify_axes(
            amplitude * rng.random(N), iterating_shape, swapaxes_parameters
        ),
        mean=modify_axes(mean + rng.random(N), iterating_shape, swapaxes_parameters),
        stddev=modify_axes(
            stddev + rng.random(N), iterating_shape, swapaxes_parameters
        ),
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

    assert_allclose(
        model_fit.amplitude.value,
        modify_axes(amplitude, iterating_shape, swapaxes_parameters),
    )
    assert_allclose(
        model_fit.mean.value, modify_axes(mean, iterating_shape, swapaxes_parameters)
    )
    assert_allclose(
        model_fit.stddev.value,
        modify_axes(stddev, iterating_shape, swapaxes_parameters),
    )


@pytest.mark.parametrize(
    "iterating_shape,swapaxes_data,swapaxes_parameters,fitting_axes",
    [
        ((120,), None, None, (0, 1)),
        ((2, 3, 4, 5), None, None, (0, 1)),  # make iterating dimensions N-d
        ((2, 3, 4, 5), (0, 1), None, (1, 0)),  # swap data axes to be (y, x)
        ((2, 3, 4, 5), (1, 2), None, (0, 2)),  # start mixing up axes
        ((2, 3, 4, 5), (0, 3), (0, 1), (3, 1)),  # iterating axes out of order
    ],
)
def test_2d_model_fit_axes(
    iterating_shape, swapaxes_data, swapaxes_parameters, fitting_axes
):
    N = 120
    P = 20
    Q = 10

    rng = np.random.default_rng(12345)

    x = np.linspace(-5, 30, P)
    y = np.linspace(5, 20, Q)

    slope_x = rng.uniform(1, 5, N)
    slope_y = rng.uniform(1, 5, N)
    intercept = rng.uniform(-2, 2, N)

    data = slope_x * x[:, None, None] + slope_y * y[None, :, None] + intercept

    # At this point, the data has shape (P, Q, N)

    # Modify the axes by reshaping and swapping axes
    data = modify_axes(data, data.shape[0:2] + iterating_shape, swapaxes_data)

    # Set initial parameters to be close to but not exactly equal to true parameters

    model = Planar2D(
        slope_x=modify_axes(
            slope_x * rng.random(N), iterating_shape, swapaxes_parameters
        ),
        slope_y=modify_axes(
            slope_y + rng.random(N), iterating_shape, swapaxes_parameters
        ),
        intercept=modify_axes(
            intercept + rng.random(N), iterating_shape, swapaxes_parameters
        ),
    )
    fitter = LevMarLSQFitter()

    model_fit = parallel_fit_model_nd(
        data=data,
        model=model,
        fitter=fitter,
        fitting_axes=fitting_axes,
        world={fitting_axes[0]: x, fitting_axes[1]: y},
    )

    # Check that shape and values match

    assert_allclose(
        model_fit.slope_x.value,
        modify_axes(slope_x, iterating_shape, swapaxes_parameters),
    )
    assert_allclose(
        model_fit.slope_y.value,
        modify_axes(slope_y, iterating_shape, swapaxes_parameters),
    )
    assert_allclose(
        model_fit.intercept.value,
        modify_axes(intercept, iterating_shape, swapaxes_parameters),
    )


def test_no_world():
    # This also doubles as a test when there are no iterating dimensions
    data = gaussian(np.arange(20), 2, 10, 1)
    model = Gaussian1D(amplitude=1.5, mean=12, stddev=1.5)
    fitter = LevMarLSQFitter()
    model_fit = parallel_fit_model_nd(
        data=data,
        model=model,
        fitter=fitter,
        fitting_axes=0,
    )
    assert_allclose(model_fit.amplitude.value, 2)
    assert_allclose(model_fit.mean.value, 10)
    assert_allclose(model_fit.stddev.value, 1)


def test_wcs_world_1d():
    # Test specifying world as a WCS, for the 1D model case

    # TODO: decide if WCS dimensionality should match cube or just the fitting
    # dimensions.

    data = gaussian(
        np.arange(20)[:, None],
        np.array([2, 1.8]),
        np.array([10, 11]),
        np.array([1, 1.1]),
    )
    model = Gaussian1D(amplitude=1.5, mean=1.2, stddev=0.15)
    fitter = LevMarLSQFitter()

    wcs = WCS(naxis=2)
    wcs.wcs.ctype = "OFFSET", "WAVE"
    wcs.wcs.crval = 10, 0.1
    wcs.wcs.crpix = 1, 1
    wcs.wcs.cdelt = 10, 0.1

    model_fit = parallel_fit_model_nd(
        data=data, model=model, fitter=fitter, fitting_axes=0, world=wcs
    )
    assert_allclose(model_fit.amplitude.value, [2, 1.8])
    assert_allclose(model_fit.mean.value, [1.1, 1.2])
    assert_allclose(model_fit.stddev.value, [0.1, 0.11])


def test_world_array():
    # Test specifying world as a tuple of arrays with dimensions matching the data

    # TODO: decide if should pass as many arrays as dimensions in cube or
    # fitting dimensions.

    data = gaussian(
        np.arange(21)[:, None],
        np.array([2, 1.8]),
        np.array([10, 10]),
        np.array([1, 1.1]),
    )
    model = Gaussian1D(amplitude=1.5, mean=7, stddev=2)
    fitter = LevMarLSQFitter()

    world1 = np.array([np.linspace(0, 10, 21), np.linspace(5, 15, 21)]).T

    world2 = np.broadcast_to(np.array([1, 2]), data.shape)

    model_fit = parallel_fit_model_nd(
        data=data, model=model, fitter=fitter, fitting_axes=0, world=(world1, world2)
    )
    assert_allclose(model_fit.amplitude.value, [2, 1.8])
    assert_allclose(model_fit.mean.value, [5, 10])
    assert_allclose(model_fit.stddev.value, [0.5, 0.55])


def test_diagnostics(tmp_path):
    data = gaussian(np.arange(20), 2, 10, 1)
    data = np.broadcast_to(data.reshape((20, 1)), (20, 3)).copy()

    data[0, 0] = np.nan

    model = Gaussian1D(amplitude=1.5, mean=12, stddev=1.5)
    fitter = LevMarLSQFitter()

    parallel_fit_model_nd(
        data=data,
        model=model,
        fitter=fitter,
        fitting_axes=0,
        diagnostics="failed",
        diagnostics_path=tmp_path / "diag1",
    )

    assert os.listdir(tmp_path / "diag1") == ["0"]
    assert sorted(os.listdir(tmp_path / "diag1" / "0")) == ["error.log", "fit.png"]

    parallel_fit_model_nd(
        data=data,
        model=model,
        fitter=fitter,
        fitting_axes=0,
        diagnostics="all",
        diagnostics_path=tmp_path / "diag2",
    )

    assert sorted(os.listdir(tmp_path / "diag2")) == ["0", "1", "2"]


@pytest.mark.parametrize(
    "scheduler", ("synchronous", "processes", "threads", "default")
)
def test_dask_scheduler(scheduler):
    N = 120
    P = 20

    rng = np.random.default_rng(12345)

    x = np.linspace(-5, 30, P)

    amplitude = rng.uniform(1, 10, N)
    mean = rng.uniform(0, 25, N)
    stddev = rng.uniform(1, 4, N)

    data = gaussian(x[:, None], amplitude, mean, stddev)

    # At this point, the data has shape (P, N)
    # Set initial parameters to be close to but not exactly equal to true parameters

    model = Gaussian1D(
        amplitude=amplitude * rng.random(N),
        mean=mean + rng.random(N),
        stddev=stddev + rng.random(N),
    )
    fitter = LevMarLSQFitter()

    model_fit = parallel_fit_model_nd(
        data=data,
        model=model,
        fitter=fitter,
        fitting_axes=0,
        world={0: x},
        scheduler=scheduler,
    )

    # Check that shape and values match

    assert_allclose(model_fit.amplitude.value, amplitude)
    assert_allclose(model_fit.mean.value, mean)
    assert_allclose(model_fit.stddev.value, stddev)


def test_compound_model():
    # Compound models have to be treated a little differently so check they
    # work fine.

    data = gaussian(
        np.arange(20)[:, None],
        np.array([2, 1.8]),
        np.array([10, 11]),
        np.array([1, 1.1]),
    )

    data[:, 0] += 2
    data[:, 1] += 3

    model1 = Gaussian1D(amplitude=1.5, mean=1.2, stddev=0.15)
    model2 = Const1D(1)

    model = model1 + model2

    fitter = LevMarLSQFitter()

    wcs = WCS(naxis=2)
    wcs.wcs.ctype = "OFFSET", "WAVE"
    wcs.wcs.crval = 10, 0.1
    wcs.wcs.crpix = 1, 1
    wcs.wcs.cdelt = 10, 0.1

    model_fit = parallel_fit_model_nd(
        data=data, model=model, fitter=fitter, fitting_axes=0, world=wcs
    )
    assert_allclose(model_fit.amplitude_0.value, [2, 1.8])
    assert_allclose(model_fit.mean_0.value, [1.1, 1.2])
    assert_allclose(model_fit.stddev_0.value, [0.1, 0.11])
    assert_allclose(model_fit.amplitude_1.value, [2, 3])

    # Check that constraints work

    model.amplitude_1 = 2
    model.amplitude_1.fixed = True

    model_fit = parallel_fit_model_nd(
        data=data, model=model, fitter=fitter, fitting_axes=0, world=wcs
    )

    assert_allclose(model_fit.amplitude_0.value, [2, 1.63349282])
    assert_allclose(model_fit.mean_0.value, [1.1, 1.145231])
    assert_allclose(model_fit.stddev_0.value, [0.1, 0.73632987])
    assert_allclose(model_fit.amplitude_1.value, [2, 2])


def test_model_dimension_mismatch():
    model = Planar2D()
    data = np.empty((20, 10, 5))
    fitter = LevMarLSQFitter()

    with pytest.raises(
        ValueError,
        match=re.escape("Model is 2-dimensional, but got 1 value(s) in fitting_axes="),
    ):
        parallel_fit_model_nd(
            data=data,
            model=model,
            fitter=fitter,
            fitting_axes=0,
        )

    model = Linear1D()
    data = np.empty((20, 10, 5))
    fitter = LevMarLSQFitter()

    with pytest.raises(
        ValueError,
        match=re.escape("Model is 1-dimensional, but got 2 value(s) in fitting_axes="),
    ):
        parallel_fit_model_nd(
            data=data,
            model=model,
            fitter=fitter,
            fitting_axes=(1, 2),
        )


def test_data_dimension_mismatch():
    model = Planar2D()
    data = np.empty((20, 10, 5))
    fitter = LevMarLSQFitter()

    with pytest.raises(
        ValueError,
        match=re.escape("Fitting index 4 out of range for 3-dimensional data"),
    ):
        parallel_fit_model_nd(
            data=data,
            model=model,
            fitter=fitter,
            fitting_axes=(1, 4),
        )


def test_world_dimension_mismatch():
    model = Planar2D()
    data = np.empty((20, 10, 5))
    fitter = LevMarLSQFitter()

    with pytest.raises(
        ValueError,
        match=re.escape(
            "world[2] has length 6 but data along dimension 2 has length 5"
        ),
    ):
        parallel_fit_model_nd(
            data=data,
            model=model,
            fitter=fitter,
            fitting_axes=(1, 2),
            world={1: np.arange(10), 2: np.arange(6)},
        )
