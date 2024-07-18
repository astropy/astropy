# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module tests model set evaluation and fitting for some common use cases.
"""

import numpy as np

# pylint: disable=invalid-name
import pytest
from numpy.testing import assert_allclose

from astropy.modeling.core import Model
from astropy.modeling.fitting import LinearLSQFitter
from astropy.modeling.models import (
    Chebyshev1D,
    Chebyshev2D,
    Hermite1D,
    Hermite2D,
    Legendre1D,
    Legendre2D,
    Linear1D,
    Planar2D,
    Polynomial1D,
    Polynomial2D,
)
from astropy.modeling.parameters import Parameter
from astropy.utils import NumpyRNGContext

x = np.arange(4)
xx = np.array([x, x + 10])
xxx = np.arange(24).reshape((3, 4, 2))

_RANDOM_SEED = 0x1337


class TParModel(Model):
    """
    A toy model to test parameters machinery
    """

    # standard_broadasting = False
    n_inputs = 1
    outputs = ("x",)
    coeff = Parameter()
    e = Parameter()

    def __init__(self, coeff, e, **kwargs):
        super().__init__(coeff=coeff, e=e, **kwargs)

    @staticmethod
    def evaluate(x, coeff, e):
        return x * coeff + e


@pytest.mark.parametrize(
    "model_class", [Polynomial1D, Chebyshev1D, Legendre1D, Hermite1D]
)
def test_model1d_axis_1(model_class):
    """
    Test that a model initialized with model_set_axis=1
    can be evaluated with model_set_axis=False.
    """
    n_models = 2
    model_axis = 1

    c0 = [[2, 3]]
    c1 = [[1, 2]]
    t1 = model_class(1, c0=2, c1=1)
    t2 = model_class(1, c0=3, c1=2)

    p1 = model_class(1, c0=c0, c1=c1, n_models=n_models, model_set_axis=model_axis)

    MESSAGE = r"For model_set_axis=1, all inputs must be at least 2-dimensional"
    with pytest.raises(ValueError, match=MESSAGE):
        p1(x)

    y = p1(x, model_set_axis=False)
    assert y.shape[model_axis] == n_models
    assert_allclose(y[:, 0], t1(x))
    assert_allclose(y[:, 1], t2(x))

    y = p1(xx, model_set_axis=False)
    assert y.shape[model_axis] == n_models
    assert_allclose(y[:, 0, :], t1(xx))
    assert_allclose(y[:, 1, :], t2(xx))

    y = p1(xxx, model_set_axis=False)
    assert y.shape[model_axis] == n_models
    assert_allclose(y[:, 0, :, :], t1(xxx))
    assert_allclose(y[:, 1, :, :], t2(xxx))


@pytest.mark.parametrize(
    "model_class", [Polynomial1D, Chebyshev1D, Legendre1D, Hermite1D]
)
def test_model1d_axis_2(model_class):
    """
    Test that a model initialized with model_set_axis=2
    can be evaluated with model_set_axis=False.
    """
    p1 = model_class(
        1, c0=[[[1, 2, 3]]], c1=[[[10, 20, 30]]], n_models=3, model_set_axis=2
    )
    t1 = model_class(1, c0=1, c1=10)
    t2 = model_class(1, c0=2, c1=20)
    t3 = model_class(1, c0=3, c1=30)

    MESSAGE = r"For model_set_axis=2, all inputs must be at least 3-dimensional"
    with pytest.raises(ValueError, match=MESSAGE):
        p1(x)

    with pytest.raises(ValueError, match=MESSAGE):
        p1(xx)

    y = p1(x, model_set_axis=False)
    assert y.shape == (1, 4, 3)
    assert_allclose(y[:, :, 0].ravel(), t1(x))
    assert_allclose(y[:, :, 1].ravel(), t2(x))
    assert_allclose(y[:, :, 2].ravel(), t3(x))


@pytest.mark.parametrize(
    "model_class", [Polynomial1D, Chebyshev1D, Legendre1D, Hermite1D]
)
def test_model1d_axis_0(model_class):
    """
    Test that a model initialized with model_set_axis=0
    can be evaluated with model_set_axis=False.
    """
    p1 = model_class(1, n_models=2, model_set_axis=0)
    p1.c0 = [2, 3]
    p1.c1 = [1, 2]
    t1 = model_class(1, c0=2, c1=1)
    t2 = model_class(1, c0=3, c1=2)

    MESSAGE = r"Input argument 'x' does not have the correct dimensions in .*"
    with pytest.raises(ValueError, match=MESSAGE):
        p1(x)

    y = p1(xx)
    assert len(y) == 2
    assert_allclose(y[0], t1(xx[0]))
    assert_allclose(y[1], t2(xx[1]))

    y = p1(x, model_set_axis=False)
    assert len(y) == 2
    assert_allclose(y[0], t1(x))
    assert_allclose(y[1], t2(x))

    y = p1(xx, model_set_axis=False)
    assert len(y) == 2
    assert_allclose(y[0], t1(xx))
    assert_allclose(y[1], t2(xx))
    y = p1(xxx, model_set_axis=False)
    assert_allclose(y[0], t1(xxx))
    assert_allclose(y[1], t2(xxx))
    assert len(y) == 2


@pytest.mark.parametrize("model_class", [Chebyshev2D, Legendre2D, Hermite2D])
def test_model2d_axis_2(model_class):
    """
    Test that a model initialized with model_set_axis=2
    can be evaluated with model_set_axis=False.
    """
    p2 = model_class(
        1,
        1,
        c0_0=[[[0, 1, 2]]],
        c0_1=[[[3, 4, 5]]],
        c1_0=[[[5, 6, 7]]],
        c1_1=[[[1, 1, 1]]],
        n_models=3,
        model_set_axis=2,
    )
    t1 = model_class(1, 1, c0_0=0, c0_1=3, c1_0=5, c1_1=1)
    t2 = model_class(1, 1, c0_0=1, c0_1=4, c1_0=6, c1_1=1)
    t3 = model_class(1, 1, c0_0=2, c0_1=5, c1_0=7, c1_1=1)

    assert p2.c0_0.shape == (1, 1, 3)
    y = p2(x, x, model_set_axis=False)
    assert y.shape == (1, 4, 3)
    # These are columns along the 2nd axis.
    assert_allclose(y[:, :, 0].ravel(), t1(x, x))
    assert_allclose(y[:, :, 1].ravel(), t2(x, x))
    assert_allclose(y[:, :, 2].ravel(), t3(x, x))


def test_negative_axis():
    p1 = Polynomial1D(1, c0=[1, 2], c1=[3, 4], n_models=2, model_set_axis=-1)
    t1 = Polynomial1D(1, c0=1, c1=3)
    t2 = Polynomial1D(1, c0=2, c1=4)

    MESSAGE = r"Input argument 'x' does not have the correct dimensions in .*"
    with pytest.raises(ValueError, match=MESSAGE):
        p1(x)

    with pytest.raises(ValueError, match=MESSAGE):
        p1(xx)

    xxt = xx.T
    y = p1(xxt)
    assert_allclose(y[:, 0], t1(xxt[:, 0]))
    assert_allclose(y[:, 1], t2(xxt[:, 1]))


def test_shapes():
    p2 = Polynomial1D(1, n_models=3, model_set_axis=2)
    assert p2.c0.shape == (1, 1, 3)
    assert p2.c1.shape == (1, 1, 3)

    p1 = Polynomial1D(1, n_models=2, model_set_axis=1)
    assert p1.c0.shape == (1, 2)
    assert p1.c1.shape == (1, 2)

    p1 = Polynomial1D(1, c0=[1, 2], c1=[3, 4], n_models=2, model_set_axis=-1)
    assert p1.c0.shape == (2,)
    assert p1.c1.shape == (2,)

    e1 = [1, 2]
    e2 = [3, 4]

    a1 = np.array([[10, 20], [30, 40]])
    a2 = np.array([[50, 60], [70, 80]])

    t = TParModel([a1, a2], [e1, e2], n_models=2, model_set_axis=-1)
    assert t.coeff.shape == (2, 2, 2)
    assert t.e.shape == (2, 2)

    t = TParModel([[a1, a2]], [[e1, e2]], n_models=2, model_set_axis=1)
    assert t.coeff.shape == (1, 2, 2, 2)
    assert t.e.shape == (1, 2, 2)

    t = TParModel([a1, a2], [e1, e2], n_models=2, model_set_axis=0)
    assert t.coeff.shape == (2, 2, 2)
    assert t.e.shape == (2, 2)

    t = TParModel([a1, a2], e=[1, 2], n_models=2, model_set_axis=0)
    assert t.coeff.shape == (2, 2, 2)
    assert t.e.shape == (2,)


def test_eval():
    """Tests evaluation of Linear1D and Planar2D with different model_set_axis."""
    model = Linear1D(slope=[1, 2], intercept=[3, 4], n_models=2)
    p = Polynomial1D(1, c0=[3, 4], c1=[1, 2], n_models=2)

    assert_allclose(model(xx), p(xx))
    assert_allclose(model(x, model_set_axis=False), p(x, model_set_axis=False))

    MESSAGE = r"Input argument 'x' does not have the correct dimensions in .*"
    with pytest.raises(ValueError, match=MESSAGE):
        model(x)

    model = Linear1D(slope=[[1, 2]], intercept=[[3, 4]], n_models=2, model_set_axis=1)
    p = Polynomial1D(1, c0=[[3, 4]], c1=[[1, 2]], n_models=2, model_set_axis=1)

    assert_allclose(model(xx.T), p(xx.T))
    assert_allclose(model(x, model_set_axis=False), p(x, model_set_axis=False))
    with pytest.raises(ValueError, match=MESSAGE):
        model(xx)

    model = Planar2D(slope_x=[1, 2], slope_y=[1, 2], intercept=[3, 4], n_models=2)
    y = model(xx, xx)

    assert y.shape == (2, 4)

    MESSAGE = r"Missing input arguments - expected 2, got 1"
    with pytest.raises(ValueError, match=MESSAGE):
        model(x)


# Test fitting


@pytest.mark.parametrize(
    "model_class", [Polynomial1D, Chebyshev1D, Legendre1D, Hermite1D]
)
def test_linearlsqfitter(model_class):
    """
    Issue #7159
    """
    p = model_class(1, n_models=2, model_set_axis=1)

    # Generate data for fitting 2 models and re-stack them along the last axis:
    y = np.array([2 * x + 1, x + 4])
    y = np.rollaxis(y, 0, -1).T

    f = LinearLSQFitter()
    # This seems to fit the model_set correctly:
    fit = f(p, x, y)
    model_y = fit(x, model_set_axis=False)

    m1 = model_class(1, c0=fit.c0[0][0], c1=fit.c1[0][0], domain=fit.domain)
    m2 = model_class(1, c0=fit.c0[0][1], c1=fit.c1[0][1], domain=fit.domain)
    assert_allclose(model_y[:, 0], m1(x))
    assert_allclose(model_y[:, 1], m2(x))

    p = model_class(1, n_models=2, model_set_axis=0)
    fit = f(p, x, y.T)


def test_model_set_axis_outputs():
    fitter = LinearLSQFitter()
    model_set = Polynomial2D(1, n_models=2, model_set_axis=2)
    y2, x2 = np.mgrid[:5, :5]
    # z = np.moveaxis([x2 + y2, 1 - 0.1 * x2 + 0.2 * y2]), 0, 2)
    z = np.rollaxis(np.array([x2 + y2, 1 - 0.1 * x2 + 0.2 * y2]), 0, 3)
    model = fitter(model_set, x2, y2, z)
    res = model(x2, y2, model_set_axis=False)
    assert z.shape == res.shape

    # Test initializing with integer model_set_axis
    # and evaluating with a different model_set_axis
    model_set = Polynomial1D(1, c0=[1, 2], c1=[2, 3], n_models=2, model_set_axis=0)
    y0 = model_set(xx)
    y1 = model_set(xx.T, model_set_axis=1)
    assert_allclose(y0[0], y1[:, 0])
    assert_allclose(y0[1], y1[:, 1])

    model_set = Polynomial1D(1, c0=[[1, 2]], c1=[[2, 3]], n_models=2, model_set_axis=1)
    y0 = model_set(xx.T)
    y1 = model_set(xx, model_set_axis=0)
    assert_allclose(y0[:, 0], y1[0])
    assert_allclose(y0[:, 1], y1[1])

    MESSAGE = r"For model_set_axis=1, all inputs must be at least 2-dimensional"
    with pytest.raises(ValueError, match=MESSAGE):
        model_set(x)


def test_fitting_shapes():
    """Test fitting model sets of Linear1D and Planar2D."""
    fitter = LinearLSQFitter()

    model = Linear1D(slope=[1, 2], intercept=[3, 4], n_models=2)
    y = model(xx)
    fitter(model, x, y)

    model = Linear1D(slope=[[1, 2]], intercept=[[3, 4]], n_models=2, model_set_axis=1)
    fitter(model, x, y.T)

    model = Planar2D(slope_x=[1, 2], slope_y=[1, 2], intercept=[3, 4], n_models=2)
    y = model(xx, xx)
    fitter(model, x, x, y)


def test_compound_model_sets():
    MESSAGE = r"model_set_axis must be False or 0 and consistent for operands"
    with pytest.raises(ValueError, match=MESSAGE):
        (
            Polynomial1D(1, n_models=2, model_set_axis=1)
            | Polynomial1D(1, n_models=2, model_set_axis=0)
        )


def test_linear_fit_model_set_errors():
    init_model = Polynomial1D(degree=2, c0=[1, 1], n_models=2)
    x = np.arange(10)
    y = init_model(x, model_set_axis=False)

    fitter = LinearLSQFitter()

    MESSAGE = r"x and y should have the same shape"
    with pytest.raises(ValueError, match=MESSAGE):
        fitter(init_model, x[:5], y)
    with pytest.raises(ValueError, match=MESSAGE):
        fitter(init_model, x, y[:, :5])


def test_linear_fit_model_set_common_weight():
    """Tests fitting multiple models simultaneously."""

    init_model = Polynomial1D(degree=2, c0=[1, 1], n_models=2)
    x = np.arange(10)
    y_expected = init_model(x, model_set_axis=False)
    assert y_expected.shape == (2, 10)

    # Add a bit of random noise
    with NumpyRNGContext(_RANDOM_SEED):
        y = y_expected + np.random.normal(0, 0.01, size=y_expected.shape)

    fitter = LinearLSQFitter()
    weights = np.ones(10)
    weights[[0, -1]] = 0
    fitted_model = fitter(init_model, x, y, weights=weights)
    assert_allclose(fitted_model(x, model_set_axis=False), y_expected, rtol=1e-1)

    # Check that using null weights raises an error
    # ValueError: On entry to DLASCL parameter number 4 had an illegal value
    with pytest.raises(ValueError, match=r"Found NaNs in the coefficient matrix"):
        with pytest.warns(
            RuntimeWarning, match=r"invalid value encountered in.*divide"
        ):
            fitted_model = fitter(init_model, x, y, weights=np.zeros(10))


def test_linear_fit_model_set_weights():
    """Tests fitting multiple models simultaneously."""

    init_model = Polynomial1D(degree=2, c0=[1, 1], n_models=2)
    x = np.arange(10)
    y_expected = init_model(x, model_set_axis=False)
    assert y_expected.shape == (2, 10)

    # Add a bit of random noise
    with NumpyRNGContext(_RANDOM_SEED):
        y = y_expected + np.random.normal(0, 0.01, size=y_expected.shape)

    weights = np.ones_like(y)
    # Put a null weight for the min and max values
    weights[[0, 1], weights.argmin(axis=1)] = 0
    weights[[0, 1], weights.argmax(axis=1)] = 0
    fitter = LinearLSQFitter()
    fitted_model = fitter(init_model, x, y, weights=weights)
    assert_allclose(fitted_model(x, model_set_axis=False), y_expected, rtol=1e-1)

    # Check that using null weights raises an error
    weights[0] = 0
    with pytest.raises(ValueError, match=r"Found NaNs in the coefficient matrix"):
        with pytest.warns(
            RuntimeWarning, match=r"invalid value encountered in.*divide"
        ):
            fitted_model = fitter(init_model, x, y, weights=weights)

    # Now we mask the values where weight is 0
    with pytest.warns(RuntimeWarning, match=r"invalid value encountered in.*divide"):
        fitted_model = fitter(
            init_model, x, np.ma.array(y, mask=np.isclose(weights, 0)), weights=weights
        )
    # Parameters for the first model are all NaNs
    assert np.all(np.isnan(fitted_model.param_sets[:, 0]))
    assert np.all(np.isnan(fitted_model(x, model_set_axis=False)[0]))
    # Second model is fitted correctly
    assert_allclose(fitted_model(x, model_set_axis=False)[1], y_expected[1], rtol=1e-1)


def test_linear_fit_2d_model_set_errors():
    init_model = Polynomial2D(degree=2, c0_0=[1, 1], n_models=2)
    x = np.arange(10)
    y = np.arange(10)
    z = init_model(x, y, model_set_axis=False)

    fitter = LinearLSQFitter()

    MESSAGE = r"x, y and z should have the same shape"
    with pytest.raises(ValueError, match=MESSAGE):
        fitter(init_model, x[:5], y, z)
    with pytest.raises(ValueError, match=MESSAGE):
        fitter(init_model, x, y, z[:, :5])


def test_linear_fit_2d_model_set_common_weight():
    init_model = Polynomial2D(
        degree=2,
        c1_0=[1, 2],
        c0_1=[-0.5, 1],
        n_models=2,
        fixed={"c1_0": True, "c0_1": True},
    )

    x, y = np.mgrid[0:5, 0:5]
    zz = np.array([1 + x - 0.5 * y + 0.1 * x * x, 2 * x + y - 0.2 * y * y])

    fitter = LinearLSQFitter()
    fitted_model = fitter(init_model, x, y, zz, weights=np.ones((5, 5)))

    assert_allclose(fitted_model(x, y, model_set_axis=False), zz, atol=1e-14)


def test_linear_fit_flat_2d_model_set_common_weight():
    init_model = Polynomial2D(
        degree=2,
        c1_0=[1, 2],
        c0_1=[-0.5, 1],
        n_models=2,
        fixed={"c1_0": True, "c0_1": True},
    )

    x, y = np.mgrid[0:5, 0:5]
    x, y = x.ravel(), y.ravel()
    zz = np.array([1 + x - 0.5 * y + 0.1 * x * x, 2 * x + y - 0.2 * y * y])
    weights = np.ones(25)

    fitter = LinearLSQFitter()
    fitted_model = fitter(init_model, x, y, zz, weights=weights)

    assert_allclose(fitted_model(x, y, model_set_axis=False), zz, atol=1e-14)


def test_linear_fit_2d_model_set_weights():
    init_model = Polynomial2D(
        degree=2,
        c1_0=[1, 2],
        c0_1=[-0.5, 1],
        n_models=2,
        fixed={"c1_0": True, "c0_1": True},
    )

    x, y = np.mgrid[0:5, 0:5]
    zz = np.array([1 + x - 0.5 * y + 0.1 * x * x, 2 * x + y - 0.2 * y * y])

    fitter = LinearLSQFitter()
    weights = [np.ones((5, 5)), np.ones((5, 5))]
    fitted_model = fitter(init_model, x, y, zz, weights=weights)

    assert_allclose(fitted_model(x, y, model_set_axis=False), zz, atol=1e-14)


def test_linear_fit_flat_2d_model_set_weights():
    init_model = Polynomial2D(
        degree=2,
        c1_0=[1, 2],
        c0_1=[-0.5, 1],
        n_models=2,
        fixed={"c1_0": True, "c0_1": True},
    )

    x, y = np.mgrid[0:5, 0:5]
    x, y = x.ravel(), y.ravel()
    zz = np.array([1 + x - 0.5 * y + 0.1 * x * x, 2 * x + y - 0.2 * y * y])
    weights = np.ones((2, 25))

    fitter = LinearLSQFitter()
    fitted_model = fitter(init_model, x, y, zz, weights=weights)

    assert_allclose(fitted_model(x, y, model_set_axis=False), zz, atol=1e-14)


class Test1ModelSet:
    """
    Check that fitting a single model works with a length-1 model set axis.
    It's not clear that this was originally intended usage, but it can be
    convenient, eg. when fitting a range of image rows that may be a single
    row, and some existing scripts might rely on it working.
    Currently this does not work with FittingWithOutlierRemoval.
    """

    def setup_class(self):
        self.x1 = np.arange(0, 10)
        self.y1 = np.array([0.5 + 2.5 * self.x1])
        self.w1 = np.ones((10,))
        self.y1[0, 8] = 100.0
        self.w1[8] = 0.0
        self.y2, self.x2 = np.mgrid[0:10, 0:10]
        self.z2 = np.array([1 - 0.1 * self.x2 + 0.2 * self.y2])
        self.w2 = np.ones((10, 10))
        self.z2[0, 1, 2] = 100.0
        self.w2[1, 2] = 0.0

    def test_linear_1d_common_weights(self):
        model = Polynomial1D(1)
        fitter = LinearLSQFitter()
        model = fitter(model, self.x1, self.y1, weights=self.w1)
        assert_allclose(model.c0, 0.5, atol=1e-12)
        assert_allclose(model.c1, 2.5, atol=1e-12)

    def test_linear_1d_separate_weights(self):
        model = Polynomial1D(1)
        fitter = LinearLSQFitter()
        model = fitter(model, self.x1, self.y1, weights=self.w1[np.newaxis, ...])
        assert_allclose(model.c0, 0.5, atol=1e-12)
        assert_allclose(model.c1, 2.5, atol=1e-12)

    def test_linear_1d_separate_weights_axis_1(self):
        model = Polynomial1D(1, model_set_axis=1)
        fitter = LinearLSQFitter()
        model = fitter(model, self.x1, self.y1.T, weights=self.w1[..., np.newaxis])
        assert_allclose(model.c0, 0.5, atol=1e-12)
        assert_allclose(model.c1, 2.5, atol=1e-12)

    def test_linear_2d_common_weights(self):
        model = Polynomial2D(1)
        fitter = LinearLSQFitter()
        model = fitter(model, self.x2, self.y2, self.z2, weights=self.w2)
        assert_allclose(model.c0_0, 1.0, atol=1e-12)
        assert_allclose(model.c1_0, -0.1, atol=1e-12)
        assert_allclose(model.c0_1, 0.2, atol=1e-12)

    def test_linear_2d_separate_weights(self):
        model = Polynomial2D(1)
        fitter = LinearLSQFitter()
        model = fitter(
            model, self.x2, self.y2, self.z2, weights=self.w2[np.newaxis, ...]
        )
        assert_allclose(model.c0_0, 1.0, atol=1e-12)
        assert_allclose(model.c1_0, -0.1, atol=1e-12)
        assert_allclose(model.c0_1, 0.2, atol=1e-12)

    def test_linear_2d_separate_weights_axis_2(self):
        model = Polynomial2D(1, model_set_axis=2)
        fitter = LinearLSQFitter()
        model = fitter(
            model,
            self.x2,
            self.y2,
            np.rollaxis(self.z2, 0, 3),
            weights=self.w2[..., np.newaxis],
        )
        assert_allclose(model.c0_0, 1.0, atol=1e-12)
        assert_allclose(model.c1_0, -0.1, atol=1e-12)
        assert_allclose(model.c0_1, 0.2, atol=1e-12)
