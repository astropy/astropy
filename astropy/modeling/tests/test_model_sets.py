# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module tests model set evaluation and fitting for some common use cases.
"""
# pylint: disable=invalid-name
import pytest
import numpy as np
from numpy.testing import assert_allclose

from astropy.modeling.models import (Polynomial1D, Polynomial2D, Legendre1D, Legendre2D,
                                     Chebyshev2D, Chebyshev1D, Hermite1D, Hermite2D,
                                     Linear1D, Planar2D)
from astropy.modeling.fitting import LinearLSQFitter
from astropy.modeling.core import Model
from astropy.modeling.parameters import Parameter


x = np.arange(4)
xx = np.array([x, x + 10])
xxx = np.arange(24).reshape((3, 4, 2))


class TParModel(Model):
    """
    A toy model to test parameters machinery
    """
    # standard_broadasting = False
    n_inputs = 1
    outputs = ('x',)
    coeff = Parameter()
    e = Parameter()

    def __init__(self, coeff, e, **kwargs):
        super().__init__(coeff=coeff, e=e, **kwargs)

    @staticmethod
    def evaluate(x, coeff, e):
        return x*coeff + e


@pytest.mark.parametrize('model_class', [Polynomial1D, Chebyshev1D, Legendre1D, Hermite1D])
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

    with pytest.raises(ValueError):
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


@pytest.mark.parametrize('model_class', [Polynomial1D, Chebyshev1D, Legendre1D, Hermite1D])
def test_model1d_axis_2(model_class):
    """
    Test that a model initialized with model_set_axis=2
    can be evaluated with model_set_axis=False.
    """
    p1 = model_class(1, c0=[[[1, 2, 3]]], c1=[[[10, 20, 30]]],
                     n_models=3, model_set_axis=2)
    t1 = model_class(1, c0=1, c1=10)
    t2 = model_class(1, c0=2, c1=20)
    t3 = model_class(1, c0=3, c1=30)

    with pytest.raises(ValueError):
        p1(x)

    with pytest.raises(ValueError):
        p1(xx)

    y = p1(x, model_set_axis=False)
    assert y.shape == (1, 4, 3)
    assert_allclose(y[:, :, 0].flatten(), t1(x))
    assert_allclose(y[:, :, 1].flatten(), t2(x))
    assert_allclose(y[:, :, 2].flatten(), t3(x))


@pytest.mark.parametrize('model_class', [Polynomial1D, Chebyshev1D, Legendre1D, Hermite1D])
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

    with pytest.raises(ValueError):
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


@pytest.mark.parametrize('model_class', [Chebyshev2D, Legendre2D, Hermite2D])
def test_model2d_axis_2(model_class):
    """
    Test that a model initialized with model_set_axis=2
    can be evaluated with model_set_axis=False.
    """
    p2 = model_class(1, 1, c0_0=[[[0, 1, 2]]], c0_1=[[[3, 4, 5]]],
                     c1_0=[[[5, 6, 7]]], c1_1=[[[1,1,1]]], n_models=3, model_set_axis=2)
    t1 = model_class(1, 1, c0_0=0, c0_1=3, c1_0=5, c1_1=1)
    t2 = model_class(1, 1, c0_0=1, c0_1=4, c1_0=6, c1_1=1)
    t3 = model_class(1, 1, c0_0=2, c0_1=5, c1_0=7, c1_1=1)

    assert p2.c0_0.shape == (1, 1, 3)
    y = p2(x, x, model_set_axis=False)
    assert y.shape == (1, 4, 3)
    # These are columns along the 2nd axis.
    assert_allclose(y[:, :, 0].flatten(), t1(x, x))
    assert_allclose(y[:, :, 1].flatten(), t2(x, x))
    assert_allclose(y[:, :, 2].flatten(), t3(x, x))


def test_negative_axis():
    p1 = Polynomial1D(1, c0=[1, 2], c1=[3, 4], n_models=2, model_set_axis=-1)
    t1 = Polynomial1D(1, c0=1, c1=3)
    t2 = Polynomial1D(1, c0=2, c1=4)

    with pytest.raises(ValueError):
        p1(x)

    with pytest.raises(ValueError):
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
    """ Tests evaluation of Linear1D and Planar2D with different model_set_axis."""
    model = Linear1D(slope=[1, 2], intercept=[3, 4], n_models=2)
    p = Polynomial1D(1, c0=[3, 4], c1=[1, 2], n_models=2)

    assert_allclose(model(xx), p(xx))
    assert_allclose(model(x, model_set_axis=False), p(x, model_set_axis=False))

    with pytest.raises(ValueError):
        model(x)

    model = Linear1D(slope=[[1, 2]], intercept=[[3, 4]], n_models=2, model_set_axis=1)
    p = Polynomial1D(1, c0=[[3, 4]], c1=[[1, 2]], n_models=2, model_set_axis=1)

    assert_allclose(model(xx.T), p(xx.T))
    assert_allclose(model(x, model_set_axis=False), p(x, model_set_axis=False))
    with pytest.raises(ValueError):
        model(xx)

    model = Planar2D(slope_x=[1, 2], slope_y=[1, 2], intercept=[3, 4], n_models=2)
    y = model(xx, xx)

    assert y.shape == (2, 4)

    with pytest.raises(ValueError):
        model(x)


# Test fitting


@pytest.mark.parametrize('model_class', [Polynomial1D, Chebyshev1D, Legendre1D, Hermite1D])
def test_linearlsqfitter(model_class):
    """
    Issue #7159
    """
    p = model_class(1, n_models=2, model_set_axis=1)

    # Generate data for fitting 2 models and re-stack them along the last axis:
    y = np.array([2*x+1, x+4])
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
    y2, x2 = np.mgrid[: 5, : 5]
    # z = np.moveaxis([x2 + y2, 1 - 0.1 * x2 + 0.2 * y2]), 0, 2)
    z = np.rollaxis(np.array([x2 + y2, 1 - 0.1 * x2 + 0.2 * y2]), 0, 3)
    model = fitter(model_set, x2, y2, z)
    res = model(x2, y2, model_set_axis=False)
    assert z.shape == res.shape

    # Test initializing with integer model_set_axis
    # and evaluating with a different model_set_axis
    model_set = Polynomial1D(1, c0=[1, 2], c1=[2, 3],
                             n_models=2, model_set_axis=0)
    y0 = model_set(xx)
    y1 = model_set(xx.T, model_set_axis=1)
    assert_allclose(y0[0], y1[:, 0])
    assert_allclose(y0[1], y1[:, 1])

    model_set = Polynomial1D(1, c0=[[1, 2]], c1=[[2, 3]],
                             n_models=2, model_set_axis=1)
    y0 = model_set(xx.T)
    y1 = model_set(xx, model_set_axis=0)
    assert_allclose(y0[:, 0], y1[0])
    assert_allclose(y0[:, 1], y1[1])
    with pytest.raises(ValueError):
        model_set(x)


def test_fitting_shapes():
    """ Test fitting model sets of Linear1D and Planar2D."""
    fitter = LinearLSQFitter()

    model = Linear1D(slope=[1, 2], intercept=[3, 4], n_models=2)
    y = model(xx)
    fit_model = fitter(model, x, y)

    model = Linear1D(slope=[[1, 2]], intercept=[[3, 4]], n_models=2, model_set_axis=1)
    fit_model = fitter(model, x, y.T)

    model = Planar2D(slope_x=[1, 2], slope_y=[1, 2], intercept=[3, 4], n_models=2)
    y = model(xx, xx)
    fit_model = fitter(model, x, x, y)


def test_compound_model_sets():
    with pytest.raises(ValueError):
        Polynomial1D(1, n_models=2, model_set_axis=1) | Polynomial1D(1, n_models=2, model_set_axis=0)
