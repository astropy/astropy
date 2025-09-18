# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import re

import pytest

pytest.importorskip("dask")
import numpy as np
from dask import array as da
from numpy.testing import assert_allclose

from astropy import units as u
from astropy.modeling.fitting import (
    FitInfoArrayContainer,
    LevMarLSQFitter,
    TRFLSQFitter,
    parallel_fit_dask,
)
from astropy.modeling.models import (
    Const1D,
    Gaussian1D,
    Linear1D,
    Planar2D,
)
from astropy.nddata import NDData, StdDevUncertainty
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.compat.optional_deps import HAS_PLT
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

    model_fit = parallel_fit_dask(
        data=data,
        model=model,
        fitter=fitter,
        fitting_axes=fitting_axes,
        world=(x,),
        scheduler="synchronous",
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

    model_fit = parallel_fit_dask(
        data=data,
        model=model,
        fitter=fitter,
        fitting_axes=fitting_axes,
        world=(x, y),
        scheduler="synchronous",
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


class TestWorld:
    def test_no_world(self):
        # This also doubles as a test when there are no iterating dimensions
        data = gaussian(np.arange(20), 2, 10, 1)
        model = Gaussian1D(amplitude=1.5, mean=12, stddev=1.5)
        fitter = LevMarLSQFitter()
        model_fit = parallel_fit_dask(
            data=data,
            model=model,
            fitter=fitter,
            fitting_axes=0,
            scheduler="synchronous",
        )
        assert_allclose(model_fit.amplitude.value, 2)
        assert_allclose(model_fit.mean.value, 10)
        assert_allclose(model_fit.stddev.value, 1)

    def test_wcs_1d(self):
        # Test specifying world as a WCS, for the 1D model case

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

        model_fit = parallel_fit_dask(
            data=data,
            model=model,
            fitter=fitter,
            fitting_axes=0,
            world=wcs,
            scheduler="synchronous",
        )
        assert_allclose(model_fit.amplitude.value, [2, 1.8])
        assert_allclose(model_fit.mean.value, [1.1, 1.2])
        assert_allclose(model_fit.stddev.value, [0.1, 0.11])

    def test_wcs_pixel_dimension_mismatch(self):
        data = np.empty((20, 10))
        model = Gaussian1D(amplitude=1.5, mean=1.2, stddev=0.15)
        fitter = LevMarLSQFitter()

        wcs = WCS(naxis=3)

        with pytest.raises(
            ValueError,
            match=re.escape(
                "The WCS pixel_n_dim (3) does not match the number of "
                "dimensions in the data (2)"
            ),
        ):
            parallel_fit_dask(
                data=data,
                model=model,
                fitter=fitter,
                fitting_axes=0,
                world=wcs,
                scheduler="synchronous",
            )

    def test_wcs_world_dimension_mismatch(self):
        data = np.empty((20, 10))
        model = Gaussian1D(amplitude=1.5, mean=1.2, stddev=0.15)
        fitter = LevMarLSQFitter()

        wcs = WCS(naxis=2)
        wcs.wcs.ctype = "RA---TAN", "DEC--TAN"
        wcs.wcs.crval = 10.0, 20.0
        wcs.wcs.crpix = 1, 1
        wcs.wcs.cdelt = 0.01, 0.01

        with pytest.raises(
            ValueError,
            match=re.escape(
                "The number of WCS world axes corresponding to the fitting axes "
                "(2) does not match the number of fitting axes (1)"
            ),
        ):
            parallel_fit_dask(
                data=data,
                model=model,
                fitter=fitter,
                fitting_axes=0,
                world=wcs,
                scheduler="synchronous",
            )

    def test_array(self):
        # Test specifying world as a tuple of arrays with dimensions matching the data

        data = gaussian(
            np.arange(21)[:, None],
            np.array([2, 1.8]),
            np.array([10, 10]),
            np.array([1, 1.1]),
        )
        model = Gaussian1D(amplitude=1.5, mean=7, stddev=2)
        fitter = LevMarLSQFitter()

        world1 = np.array([np.linspace(0, 10, 21), np.linspace(5, 15, 21)]).T

        model_fit = parallel_fit_dask(
            data=data,
            model=model,
            fitter=fitter,
            fitting_axes=0,
            world=(world1,),
            scheduler="synchronous",
        )
        assert_allclose(model_fit.amplitude.value, [2, 1.8])
        assert_allclose(model_fit.mean.value, [5, 10])
        assert_allclose(model_fit.stddev.value, [0.5, 0.55])

    def test_array_length_mismatch(self):
        data = np.empty((20, 2))
        model = Gaussian1D(amplitude=1.5, mean=7, stddev=2)
        fitter = LevMarLSQFitter()

        world1 = np.array([np.linspace(0, 10, 21), np.linspace(5, 15, 21)]).T

        with pytest.raises(
            ValueError,
            match=re.escape(
                "The number of world arrays (2) must match number of fitting axes (1)"
            ),
        ):
            parallel_fit_dask(
                data=data,
                model=model,
                fitter=fitter,
                fitting_axes=0,
                world=(world1, world1),
            )

    def test_array_shape_mismatch(self):
        data = np.empty((20, 2))
        model = Gaussian1D(amplitude=1.5, mean=7, stddev=2)
        fitter = LevMarLSQFitter()

        world1 = np.array([np.linspace(0, 10, 31), np.linspace(5, 15, 31)]).T

        with pytest.raises(
            ValueError,
            match=re.escape(
                "The arrays in the world tuple should be broadcastable to the "
                "shape of the data (expected (20, 2)), got (31, 2))"
            ),
        ):
            parallel_fit_dask(
                data=data,
                model=model,
                fitter=fitter,
                fitting_axes=0,
                world=(world1,),
            )

    def test_array_dimension_mismatch(self):
        model = Planar2D()
        data = np.empty((20, 10, 5))
        fitter = LevMarLSQFitter()

        with pytest.raises(
            ValueError,
            match=re.escape(
                "world[2] has length 6 but data along dimension 2 has length 5"
            ),
        ):
            parallel_fit_dask(
                data=data,
                model=model,
                fitter=fitter,
                fitting_axes=(1, 2),
                world=(np.arange(10), np.arange(6)),
            )

    def test_invalid_type(self):
        data = np.empty((20, 2))
        model = Gaussian1D(amplitude=1.5, mean=7, stddev=2)
        fitter = LevMarLSQFitter()

        world1 = np.array([np.linspace(0, 10, 21), np.linspace(5, 15, 21)]).T

        with pytest.raises(
            TypeError,
            match=re.escape("world should be None, a WCS object or a tuple of arrays"),
        ):
            parallel_fit_dask(
                data=data,
                model=model,
                fitter=fitter,
                fitting_axes=0,
                world="banana",
            )


def test_fitter_kwargs(tmp_path):
    data = gaussian(np.arange(20), 2, 10, 1)
    data = np.broadcast_to(data.reshape((20, 1)), (20, 3)).copy()

    data[0, 0] = np.nan

    model = Gaussian1D(amplitude=1.5, mean=12, stddev=1.5)
    fitter = LevMarLSQFitter()

    model_fit = parallel_fit_dask(
        data=data,
        model=model,
        fitter=fitter,
        fitting_axes=0,
        fitter_kwargs={"filter_non_finite": True},
        scheduler="synchronous",
    )

    assert_allclose(model_fit.amplitude.value, 2)
    assert_allclose(model_fit.mean.value, 10)
    assert_allclose(model_fit.stddev.value, 1)


class TestDiagnostics:
    def setup_method(self, method):
        self.data = gaussian(np.arange(20), 2, 10, 1)
        self.data = np.broadcast_to(self.data.reshape((20, 1)), (20, 3)).copy()
        self.data_original = self.data.copy()
        self.data[0, 0] = np.nan
        self.model = Gaussian1D(amplitude=1.5, mean=12, stddev=1.5)
        self.fitter = LevMarLSQFitter()

    def test_error(self, tmp_path):
        parallel_fit_dask(
            data=self.data,
            model=self.model,
            fitter=self.fitter,
            fitting_axes=0,
            diagnostics="error",
            diagnostics_path=tmp_path / "diag1",
            scheduler="synchronous",
        )

        assert os.listdir(tmp_path / "diag1") == ["0"]
        assert sorted(os.listdir(tmp_path / "diag1" / "0")) == ["error.log"]

    def test_all(self, tmp_path):
        parallel_fit_dask(
            data=self.data,
            model=self.model,
            fitter=self.fitter,
            fitting_axes=0,
            diagnostics="all",
            diagnostics_path=tmp_path / "diag2",
            scheduler="synchronous",
        )

        assert sorted(os.listdir(tmp_path / "diag2")) == ["0", "1", "2"]

    def test_all_world_wcs(self, tmp_path):
        # Make sure things world also with world=wcs

        parallel_fit_dask(
            data=self.data,
            model=self.model,
            fitter=self.fitter,
            fitting_axes=0,
            world=WCS(naxis=2),
            diagnostics="error",
            diagnostics_path=tmp_path / "diag3",
            scheduler="synchronous",
        )

        assert os.listdir(tmp_path / "diag3") == ["0"]
        assert sorted(os.listdir(tmp_path / "diag3" / "0")) == ["error.log"]

    @pytest.mark.skipif(not HAS_PLT, reason="requires matplotlib.pyplot")
    def test_callable(self, tmp_path):
        # And check that we can pass in a callable

        def custom_callable(path, world, data, weights, model, fitting_kwargs):
            import matplotlib.pyplot as plt

            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(world[0], data, "k.")
            ax.text(0.1, 0.9, "Fit failed!", color="r", transform=ax.transAxes)
            fig.savefig(os.path.join(path, "fit.png"))
            plt.close(fig)

        # Note: here we keep the default scheduler ('processes') to make sure
        # that callables are passed correctly to other processes. This test is
        # not as fast as other ones anyway due to the plotting so this doesn't
        # have a big performance impact.

        parallel_fit_dask(
            data=self.data,
            model=self.model,
            fitter=self.fitter,
            fitting_axes=0,
            world=WCS(naxis=2),
            diagnostics="error",
            diagnostics_path=tmp_path / "diag4",
            diagnostics_callable=custom_callable,
        )

        assert os.listdir(tmp_path / "diag4") == ["0"]
        assert sorted(os.listdir(tmp_path / "diag4" / "0")) == ["error.log", "fit.png"]

    def test_warnings(self, tmp_path):
        # Check that catching warnings works

        parallel_fit_dask(
            data=self.data_original,
            model=self.model,
            fitter=self.fitter,
            fitting_axes=0,
            diagnostics="error+warn",
            diagnostics_path=tmp_path / "diag5",
            fitter_kwargs={"maxiter": 2},
            scheduler="synchronous",
        )

        assert sorted(os.listdir(tmp_path / "diag5")) == ["0", "1", "2"]
        assert sorted(os.listdir(tmp_path / "diag5" / "0")) == ["warn.log"]

    def test_missing_path(self):
        with pytest.raises(ValueError, match="diagnostics_path should be set"):
            parallel_fit_dask(
                data=self.data,
                model=self.model,
                fitter=self.fitter,
                fitting_axes=0,
                diagnostics="error",
            )

    def test_invalid(self):
        with pytest.raises(
            ValueError,
            match=re.escape(
                "diagnostics should be None, 'error', 'error+warn', or 'all'"
            ),
        ):
            parallel_fit_dask(
                data=self.data,
                model=self.model,
                fitter=self.fitter,
                fitting_axes=0,
                diagnostics="spam",
            )


@pytest.mark.parametrize(
    "scheduler", (None, "synchronous", "processes", "threads", "default")
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

    model_fit = parallel_fit_dask(
        data=data,
        model=model,
        fitter=fitter,
        fitting_axes=0,
        world=(x,),
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

    model_fit = parallel_fit_dask(
        data=data,
        model=model,
        fitter=fitter,
        fitting_axes=0,
        world=wcs,
        scheduler="synchronous",
    )
    assert_allclose(model_fit.amplitude_0.value, [2, 1.8])
    assert_allclose(model_fit.mean_0.value, [1.1, 1.2])
    assert_allclose(model_fit.stddev_0.value, [0.1, 0.11])
    assert_allclose(model_fit.amplitude_1.value, [2, 3])

    # Check that constraints work

    model.amplitude_1 = 2
    model.amplitude_1.fixed = True

    model_fit = parallel_fit_dask(
        data=data,
        model=model,
        fitter=fitter,
        fitting_axes=0,
        world=wcs,
        scheduler="synchronous",
    )

    assert_allclose(model_fit.amplitude_0.value, [2, 1.633], atol=0.001)
    assert_allclose(model_fit.mean_0.value, [1.1, 1.145], atol=0.001)
    assert_allclose(model_fit.stddev_0.value, [0.1, 0.736], atol=0.001)
    assert_allclose(model_fit.amplitude_1.value, [2, 2], atol=0.001)


def test_model_dimension_mismatch():
    model = Planar2D()
    data = np.empty((20, 10, 5))
    fitter = LevMarLSQFitter()

    with pytest.raises(
        ValueError,
        match=re.escape("Model is 2-dimensional, but got 1 value(s) in fitting_axes="),
    ):
        parallel_fit_dask(
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
        parallel_fit_dask(
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
        parallel_fit_dask(
            data=data,
            model=model,
            fitter=fitter,
            fitting_axes=(1, 4),
        )


class TestDaskInput:
    # Check that things work correctly when passing dask arrays as input.

    def setup_method(self, method):
        self.base_data = gaussian(np.arange(20), 2, 10, 1)
        self.base_data = np.broadcast_to(
            self.base_data.reshape((20, 1)), (20, 3)
        ).copy()
        self.model = Gaussian1D(amplitude=1.5, mean=12, stddev=1.5)
        self.fitter = LevMarLSQFitter()

    @pytest.mark.parametrize("preserve_native_chunks", (False, True))
    def test_data(self, preserve_native_chunks):
        model_fit = parallel_fit_dask(
            data=da.from_array(self.base_data, chunks=(20, 1)),
            model=self.model,
            fitter=self.fitter,
            fitting_axes=0,
            preserve_native_chunks=preserve_native_chunks,
            scheduler="synchronous",
        )

        assert_allclose(model_fit.amplitude.value, 2)
        assert_allclose(model_fit.mean.value, 10)
        assert_allclose(model_fit.stddev.value, 1)

    @pytest.mark.parametrize("preserve_native_chunks", (False, True))
    def test_data_and_weights(self, preserve_native_chunks):
        data = self.base_data.copy()
        data[0, 0] = 1000.0
        data = da.from_array(data, chunks=(20, 1))

        weights = (data < 100).astype(float)

        model_fit = parallel_fit_dask(
            data=data,
            weights=weights,
            model=self.model,
            fitter=self.fitter,
            fitting_axes=0,
            preserve_native_chunks=preserve_native_chunks,
            scheduler="synchronous",
        )

        assert_allclose(model_fit.amplitude.value, 2)
        assert_allclose(model_fit.mean.value, 10)
        assert_allclose(model_fit.stddev.value, 1)

    def test_preserve_native_chunks_invalid_data_chunks(self):
        with pytest.raises(
            ValueError,
            match=re.escape(
                "When using preserve_native_chunks=True, the chunk size should match the data size along the fitting axe"
            ),
        ):
            parallel_fit_dask(
                data=da.from_array(self.base_data, chunks=(5, 2)),
                model=self.model,
                fitter=self.fitter,
                fitting_axes=0,
                preserve_native_chunks=True,
            )

    def test_preserve_native_chunks_invalid_weight_chunks(self):
        with pytest.raises(
            ValueError,
            match=re.escape(
                "When using preserve_native_chunks=True, the weights should have the same chunk size as the data"
            ),
        ):
            parallel_fit_dask(
                data=da.from_array(self.base_data, chunks=(20, 1)),
                weights=da.random.random(self.base_data.shape).rechunk((18, 2)),
                model=self.model,
                fitter=self.fitter,
                fitting_axes=0,
                preserve_native_chunks=True,
            )

    def test_preserve_native_chunks_invalid_input_data_type(self):
        with pytest.raises(
            TypeError,
            match=re.escape(
                "Can only set preserve_native_chunks=True if input data is a dask array"
            ),
        ):
            parallel_fit_dask(
                data=self.base_data,
                model=self.model,
                fitter=self.fitter,
                fitting_axes=0,
                preserve_native_chunks=True,
            )

    def test_preserve_native_chunks_invalid_input_weights_type(self):
        with pytest.raises(
            TypeError,
            match=re.escape(
                "Can only set preserve_native_chunks=True if input weights is a dask array (if specified)"
            ),
        ):
            parallel_fit_dask(
                data=da.from_array(self.base_data, chunks=(20, 1)),
                weights=np.random.random(self.base_data.shape),
                model=self.model,
                fitter=self.fitter,
                fitting_axes=0,
                preserve_native_chunks=True,
            )


def test_weights():
    # Test specifying world as a tuple of arrays with dimensions matching the data

    data = gaussian(
        np.arange(21)[:, None],
        np.array([2, 1.8]),
        np.array([5, 10]),
        np.array([1, 1.1]),
    )

    # Introduce outliers but then adjust weights to ignore them
    data[10, 0] = 1000.0
    data[20, 1] = 2000.0
    weights = (data < 100.0).astype(float)

    model = Gaussian1D(amplitude=1.5, mean=7, stddev=2)
    fitter = LevMarLSQFitter()

    model_fit = parallel_fit_dask(
        data=data,
        weights=weights,
        model=model,
        fitter=fitter,
        fitting_axes=0,
        scheduler="synchronous",
    )
    assert_allclose(model_fit.amplitude.value, [2, 1.8])
    assert_allclose(model_fit.mean.value, [5, 10])
    assert_allclose(model_fit.stddev.value, [1.0, 1.1])


class TestUnits:
    def test_basic(self):
        # Make sure that fitting with units works

        data = (
            gaussian(
                np.arange(21)[:, None],
                np.array([2, 1.8]),
                np.array([5, 10]),
                np.array([1, 1.1]),
            )
            * u.Jy
        )

        model = Gaussian1D(amplitude=1.5 * u.Jy, mean=7 * u.um, stddev=0.002 * u.mm)
        fitter = LevMarLSQFitter()

        model_fit = parallel_fit_dask(
            data=data,
            model=model,
            fitter=fitter,
            fitting_axes=0,
            world=(1000 * np.arange(21) * u.nm,),
            scheduler="synchronous",
        )

        assert_quantity_allclose(model_fit.amplitude.quantity, [2, 1.8] * u.Jy)
        assert_quantity_allclose(model_fit.mean.quantity, [5, 10] * u.um)
        assert_quantity_allclose(model_fit.stddev.quantity, [1.0, 1.1] * u.um)

    def test_units_no_input_units(self):
        # Make sure that fitting with units works for models without input_units defined

        data = (np.repeat(3, 20)).reshape((20, 1)) * u.Jy

        model = Const1D(1 * u.mJy)
        fitter = LevMarLSQFitter()

        assert not model.input_units

        model_fit = parallel_fit_dask(
            data=data,
            model=model,
            fitter=fitter,
            fitting_axes=0,
            world=(1000 * np.arange(20) * u.nm,),
            scheduler="synchronous",
        )

        assert_quantity_allclose(model_fit.amplitude.quantity, 3 * u.Jy)

    def test_units_with_wcs(self):
        data = gaussian(np.arange(20), 2, 10, 1).reshape((20, 1)) * u.Jy
        model = Gaussian1D(amplitude=1.5 * u.Jy, mean=7 * u.um, stddev=0.002 * u.mm)
        fitter = LevMarLSQFitter()

        wcs = WCS(naxis=2)
        wcs.wcs.ctype = "OFFSET", "WAVE"
        wcs.wcs.crval = 10, 0.1
        wcs.wcs.crpix = 1, 1
        wcs.wcs.cdelt = 10, 0.1
        wcs.wcs.cunit = "deg", "um"

        model_fit = parallel_fit_dask(
            data=data,
            model=model,
            fitter=fitter,
            fitting_axes=0,
            world=wcs,
            scheduler="synchronous",
        )
        assert_allclose(model_fit.amplitude.quantity, 2 * u.Jy)
        assert_allclose(model_fit.mean.quantity, 1.1 * u.um)
        assert_allclose(model_fit.stddev.quantity, 0.1 * u.um)

    def test_units_with_wcs_2d(self):
        # Fit a 2D model with mixed units to a dataset

        N = 3
        P = 20
        Q = 10

        wcs = WCS(naxis=3)
        wcs.wcs.ctype = "OFFSET", "WEIGHTS", "WAVE"
        wcs.wcs.crval = 10, 0.1, 0.2
        wcs.wcs.crpix = 1, 1, 1
        wcs.wcs.cdelt = 10, 0.1, 0.2
        wcs.wcs.cunit = "deg", "kg", "mm"

        rng = np.random.default_rng(12345)

        x = wcs.pixel_to_world(0, 0, np.arange(P))[2]
        y = wcs.pixel_to_world(0, np.arange(Q), 0)[1]

        slope_x = [1, 3, 2] * u.Jy / u.mm
        slope_y = [-1, 2, 3] * u.Jy / u.kg
        intercept = [5, 6, 7] * u.Jy

        data = slope_x * x[:, None, None] + slope_y * y[None, :, None] + intercept

        # At this point, the data has shape (P, Q, N)

        model = Planar2D(
            slope_x=slope_x * rng.uniform(0.9, 1.1, N),
            slope_y=slope_y * rng.uniform(0.9, 1.1, N),
            intercept=intercept * rng.uniform(0.9, 1.1, N),
        )
        fitter = LevMarLSQFitter()

        model_fit = parallel_fit_dask(
            data=data,
            model=model,
            fitter=fitter,
            fitting_axes=(0, 1),
            world=wcs,
            scheduler="synchronous",
        )
        assert_allclose(model_fit.slope_x.quantity, slope_x)
        assert_allclose(model_fit.slope_y.quantity, slope_y)
        assert_allclose(model_fit.intercept.quantity, intercept)


def test_skip_empty_data(tmp_path):
    # Test when one of the datasets being fit is all NaN

    data = gaussian(
        np.arange(21)[:, None],
        np.array([2, 1.8]),
        np.array([5, 10]),
        np.array([1, 1.1]),
    )

    data[:, 1] = np.nan

    model = Gaussian1D(amplitude=1.5, mean=7, stddev=2)
    fitter = TRFLSQFitter()

    model_fit = parallel_fit_dask(
        data=data,
        model=model,
        fitter=fitter,
        fitting_axes=0,
        diagnostics="error+warn",
        diagnostics_path=tmp_path,
        scheduler="synchronous",
    )

    # If we don't properly skip empty data sections, then the fitting would fail
    # and populate the diagnostics directory.
    assert len(sorted(tmp_path.iterdir())) == 0

    assert_allclose(model_fit.amplitude.value, [2, np.nan])
    assert_allclose(model_fit.mean.value, [5, np.nan])
    assert_allclose(model_fit.stddev.value, [1.0, np.nan])


def test_world_wcs_axis_correlation():
    # Regression test for a bug that caused the world coordinates to not be
    # properly extracted from a WCS object if the axis correlation matrix
    # resulted in a difference in order between pixel and world coordinates.

    model = Gaussian1D()
    fitter = TRFLSQFitter()
    data = gaussian(
        np.arange(1, 6)[:, None],
        np.array([5, 5]),
        np.array([3, 2]),
        np.array([1, 1]),
    ).T

    common_kwargs = dict(data=data, model=model, fitter=fitter, scheduler="synchronous")

    # First case: a simple WCS - as the fitting axis is 1 in Numpy order, this
    # means we should use the world coordinates for the first WCS dimension. In
    # this case, the Gaussian means should be 3 and 2.

    wcs1 = WCS(naxis=2)
    wcs1.wcs.cdelt = 1, 2

    model_fit = parallel_fit_dask(fitting_axes=1, world=wcs1, **common_kwargs)
    assert_allclose(model_fit.mean, [3, 2])

    # Second case: as above, but WCS axes swapped. In this case, the means
    # should be 6 and 4.

    wcs1 = WCS(naxis=2)
    wcs1.wcs.cdelt = 2, 1

    model_fit = parallel_fit_dask(fitting_axes=1, world=wcs1, **common_kwargs)
    assert_allclose(model_fit.mean, [6, 4])

    # Third case: as in first case, but this time we set the PC matrix such
    # that the world axes are in a different order to their corresponding pixel
    # axis. In this case, the means should be 6 and 4 because fitting_axes=1
    # should correspond to the second WCS dimension.

    wcs3 = WCS(naxis=2)
    wcs3.wcs.cdelt = 1, 2
    wcs3.wcs.pc = [[0, 1], [1, 0]]

    model_fit = parallel_fit_dask(fitting_axes=1, world=wcs1, **common_kwargs)
    assert_allclose(model_fit.mean, [6, 4])


def test_support_nddata():
    data = gaussian(np.arange(20), 2, 10, 1).reshape((20, 1)) * u.Jy

    # Introduce outliers
    data[10, 0] = 1000.0 * u.Jy
    # Mask the outliers (invalid is True)
    mask = data > 100.0 * u.Jy

    model = Gaussian1D(amplitude=1.5 * u.Jy, mean=7 * u.um, stddev=0.002 * u.mm)
    fitter = LevMarLSQFitter()

    wcs = WCS(naxis=2)
    wcs.wcs.ctype = "OFFSET", "WAVE"
    wcs.wcs.crval = 10, 0.1
    wcs.wcs.crpix = 1, 1
    wcs.wcs.cdelt = 10, 0.1
    wcs.wcs.cunit = "deg", "um"

    nd_data = NDData(
        data=data,
        wcs=wcs,
        mask=mask,
    )

    model_fit = parallel_fit_dask(
        data=nd_data,
        model=model,
        fitter=fitter,
        fitting_axes=0,
        scheduler="synchronous",
    )

    assert_allclose(model_fit.amplitude.quantity, 2 * u.Jy)
    assert_allclose(model_fit.mean.quantity, 1.1 * u.um)
    assert_allclose(model_fit.stddev.quantity, 0.1 * u.um)


def test_support_nddata_uncert():
    data = np.repeat(np.array([1, 2, 4]), 2).reshape((2, -1), order="F").T
    uncert = (
        np.repeat(np.array([1 / 7**0.5, 1 / 2**0.5, 1]), 2)
        .reshape((2, -1), order="F")
        .T
    )

    model = Const1D(0)
    fitter = LevMarLSQFitter()

    wcs = WCS(naxis=2)
    wcs.wcs.ctype = "OFFSET", "WAVE"
    wcs.wcs.crval = 10, 1
    wcs.wcs.crpix = 1, 1
    wcs.wcs.cdelt = 10, 1
    wcs.wcs.cunit = "deg", "m"

    nd_data = NDData(
        data=data,
        wcs=wcs,
        uncertainty=StdDevUncertainty(uncert),
    )

    model_fit = parallel_fit_dask(
        data=nd_data,
        model=model,
        fitter=fitter,
        fitting_axes=0,
        scheduler="synchronous",
    )

    assert_allclose(model_fit.amplitude, 1.5)


class TestFitInfo:
    def setup_method(self, method):
        self.data = gaussian(np.arange(20), 2, 10, 1)
        self.data = np.broadcast_to(self.data.reshape((20, 1)), (20, 3)).copy()
        self.data_original = self.data.copy()
        self.data[0, 0] = np.nan
        self.model = Gaussian1D(amplitude=1.5, mean=12, stddev=1.5)

    def test_default(self, tmp_path):
        fitter = TRFLSQFitter()

        parallel_fit_dask(
            data=self.data,
            model=self.model,
            fitter=fitter,
            fitting_axes=0,
            scheduler="synchronous",
        )

        assert fitter.fit_info is None

    def test_all(self, tmp_path):
        fitter = TRFLSQFitter()

        parallel_fit_dask(
            data=self.data,
            model=self.model,
            fitter=fitter,
            fitting_axes=0,
            scheduler="synchronous",
            fit_info=True,
        )

        assert "message" in fitter.fit_info.properties

        assert_allclose(fitter.fit_info.get_property_as_array("nfev"), [0, 9, 9])

        param_cov_array = fitter.fit_info.get_property_as_array("param_cov")
        assert param_cov_array.shape == (3, 3, 3)
        assert_allclose(param_cov_array[0], 0)
        assert_allclose(param_cov_array[1], param_cov_array[2])
        assert np.any(np.abs(param_cov_array[1]) > 0)

        # Test slicing that returns an array

        assert fitter.fit_info.shape == (3,)
        fit_info_subset = fitter.fit_info[:2]
        assert isinstance(fit_info_subset, FitInfoArrayContainer)
        assert fit_info_subset.shape == (2,)
        assert_allclose(fit_info_subset.get_property_as_array("nfev"), [0, 9])

        # Test slicing that returns a one element array

        fit_info_subset_single = fitter.fit_info[1:2]
        assert isinstance(fit_info_subset_single, FitInfoArrayContainer)
        assert fit_info_subset_single.shape == (1,)
        assert_allclose(fit_info_subset_single.get_property_as_array("nfev"), [9])

        # Test slicing that returns a scalar

        fit_info_indiv = fitter.fit_info[1]
        assert not isinstance(fit_info_indiv, FitInfoArrayContainer)
        assert fit_info_indiv.nfev == 9
        assert fit_info_indiv.message != ""

    def test_subset(self, tmp_path):
        fitter = TRFLSQFitter()

        parallel_fit_dask(
            data=self.data,
            model=self.model,
            fitter=fitter,
            fitting_axes=0,
            scheduler="synchronous",
            fit_info=("message", "nfev", "success"),
        )

        assert fitter.fit_info.properties == ("message", "nfev", "success")

        assert_allclose(fitter.fit_info.get_property_as_array("nfev"), [0, 9, 9])
