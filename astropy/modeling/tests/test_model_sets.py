# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module tests model set evaluation for some common use cases.
"""
import pytest
import numpy as np
from numpy.testing.utils import assert_allclose

from ..models import Polynomial1D, Polynomial2D
from ..fitting import LinearLSQFitter
from ..core import Model
from ..parameters import Parameter


x = np.arange(4)
xx = np.array([x, x + 10])
xxx = np.arange(24).reshape((3, 4, 2))


class TParModel(Model):
    """
    A toy model to test parameters machinery
    """
    #standard_broadasting = False
    inputs = ('x',)
    outputs = ('x',)
    coeff = Parameter()
    e = Parameter()

    def __init__(self, coeff, e, **kwargs):
        super().__init__(coeff=coeff, e=e, **kwargs)

    @staticmethod
    def evaluate(x, coeff, e):
        return x*coeff + e


def test_model_axis_1():
    """
    Test that a model initialized with model_set_axis=1
    can be evaluated with model_set_axis=False.
    """

    p1 = Polynomial1D(1, n_models=2, model_set_axis=1)
    p1.c0 = [2, 3]
    p1.c1 = [1, 2]
    t1 = Polynomial1D(1, c0=2, c1=1)
    t2 = Polynomial1D(1, c0=3, c1=2)

    with pytest.raises(ValueError):
        p1(x)

    with pytest.raises(ValueError):
        p1(xx)

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


def test_model_axis_2():
    """
    Test that a model initialized with model_set_axis=2
    can be evaluated with model_set_axis=False.
    """
    p1 = Polynomial1D(1, c0=[[[1, 2,3 ]]], c1=[[[10, 20, 30]]],
                      n_models=3, model_set_axis=2)
    t1 = Polynomial1D(1, c0=1, c1=10)
    t2 = Polynomial1D(1, c0=2, c1=20)
    t3 = Polynomial1D(1, c0=3, c1=30)

    with pytest.raises(ValueError):
        p1(x)

    with pytest.raises(ValueError):
        p1(xx)

    y = p1(x, model_set_axis=False)
    assert y.shape == (4, 1, 3)
    assert_allclose(y[:, :, 0].flatten(), t1(x))
    assert_allclose(y[:, :, 1].flatten(), t2(x))
    assert_allclose(y[:, :, 2].flatten(), t3(x))

    p2 = Polynomial2D(1, c0_0=[[[0,1,2]]], c0_1=[[[3,4,5]]],
                      c1_0=[[[5,6,7]]], n_models=3, model_set_axis=2)
    t1 = Polynomial2D(1, c0_0=0, c0_1=3, c1_0=5)
    t2 = Polynomial2D(1, c0_0=1, c0_1=4, c1_0=6)
    t3 = Polynomial2D(1, c0_0=2, c0_1=5, c1_0=7)

    assert p2.c0_0.shape == ()
    y = p2(x, x, model_set_axis=False)
    assert y.shape == (4, 1, 3)
    # These are columns along the 2nd axis.
    assert_allclose(y[:, :, 0].flatten(), t1(x, x))
    assert_allclose(y[:, :, 1].flatten(), t2(x, x))
    assert_allclose(y[:, :, 2].flatten(), t3(x, x))


def test_axis_0():
    """
    Test that a model initialized with model_set_axis=0
    can be evaluated with model_set_axis=False.
    """
    x = np.arange(4)
    xx = np.array([x, x + 10])
    xxx = np.arange(24).reshape((3, 4, 2))

    p1 = Polynomial1D(1, n_models=2, model_set_axis=0)
    p1.c0 = [2, 3]
    p1.c1 = [1, 2]
    t1 = Polynomial1D(1, c0=2, c1=1)
    t2 = Polynomial1D(1, c0=3, c1=2)

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


def test_negative_axis():
    p1 = Polynomial1D(1, c0=[1, 2], c1=[3, 4], n_models=2, model_set_axis=-1)
    x = np.arange(4)
    xx = np.array([x, x + 10])
    t1 = Polynomial1D(1, c0=1,c1=3)
    t2 = Polynomial1D(1, c0=2,c1=4)

    with pytest.raises(ValueError):
        p1(x)

    with pytest.raises(ValueError):
        p1(xx)

    xxt = xx.T
    y = p1(xxt)
    assert_allclose(y[: ,0], t1(xxt[: ,0]))
    assert_allclose(y[: ,1], t2(xxt[: ,1]))


def test_shapes():
    p2 = Polynomial1D(1, n_models=3, model_set_axis=2)
    assert p2.c0.shape == ()
    assert p2.c1.shape == ()

    p1 = Polynomial1D(1, n_models=2, model_set_axis=1)
    assert p1.c0.shape == ()
    assert p1.c1.shape == ()

    p1 = Polynomial1D(1, c0=[1, 2], c1=[3, 4], n_models=2, model_set_axis=-1)
    assert p1.c0.shape == ()
    assert p1.c1.shape == ()

    e1 = [1, 2]
    e2 = [3, 4]

    a1 = np.array([[10, 20], [30, 40]])
    a2 = np.array([[50, 60], [70, 80]])

    t = TParModel([a1, a2], [e1, e2], n_models=2, model_set_axis=-1)
    assert t.coeff.shape == (2, 2)
    assert t.e.shape == (2,)



    t = TParModel([[a1, a2]], [[e1, e2]], n_models=2, model_set_axis=1)
    assert t.coeff.shape == (2, 2)
    assert t.e.shape == (2,)

    t = TParModel([a1, a2], [e1, e2], n_models=2, model_set_axis=0)
    assert t.coeff.shape == (2, 2)
    assert t.e.shape == (2,)

    t = TParModel([a1, a2], e=[1, 2], n_models=2, model_set_axis=0)
    assert t.coeff.shape == (2, 2)
    assert t.e.shape == ()


def test_linearlsqfitter():
    """
    Issue #7159
    """
    x = np.arange(0,4)
    p = Polynomial1D(1, n_models=2, model_set_axis=1)

    # Generate data for fitting 2 models and re-stack them along the last axis:
    y = np.array([2*x+1, x+4])
    y = np.rollaxis(y, 0, -1).T

    f = LinearLSQFitter()
    # This seems to fit the model_set correctly:
    fit = f(p, x, y)
    model_y = fit(x, model_set_axis=False)

    m1 = Polynomial1D(1, c0=fit.c0[0][0], c1=fit.c1[0][0])
    m2 = Polynomial1D(1, c0=fit.c0[0][1], c1=fit.c1[0][1])
    assert_allclose(model_y[0], m1(x))
    assert_allclose(model_y[1], m2(x))
