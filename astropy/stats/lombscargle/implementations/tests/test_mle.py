import pytest
import numpy as np
from numpy.testing import assert_allclose

from .....extern.six.moves import range
from ..mle import design_matrix, periodic_fit


@pytest.fixture
def t():
    rand = np.random.RandomState(42)
    return 10 * rand.rand(10)


@pytest.mark.parametrize('freq', [1.0, 2])
@pytest.mark.parametrize('dy', [None, 2.0])
@pytest.mark.parametrize('bias', [True, False])
def test_design_matrix(t, freq, dy, bias):
    X = design_matrix(t, freq, dy, bias=bias)
    assert X.shape == (t.shape[0], 2 + bool(bias))
    if bias:
        assert_allclose(X[:, 0], 1. / (dy or 1.0))
    assert_allclose(X[:, -2], np.sin(2 * np.pi * freq * t) / (dy or 1.0))
    assert_allclose(X[:, -1], np.cos(2 * np.pi * freq * t) / (dy or 1.0))


@pytest.mark.parametrize('nterms', range(4))
def test_multiterm_design_matrix(t, nterms):
    dy = 2.0
    freq = 1.5
    X = design_matrix(t, freq, dy=dy, bias=True, nterms=nterms)
    assert X.shape == (t.shape[0], 1 + 2 * nterms)
    assert_allclose(X[:, 0], 1. / dy)
    for i in range(1, nterms + 1):
        assert_allclose(X[:, 2 * i - 1], np.sin(2 * np.pi * i * freq * t) / dy)
        assert_allclose(X[:, 2 * i], np.cos(2 * np.pi * i * freq * t) / dy)


@pytest.mark.parametrize('nterms', range(1, 4))
@pytest.mark.parametrize('freq', [1, 2])
@pytest.mark.parametrize('fit_mean', [True, False])
def test_exact_mle_fit(nterms, freq, fit_mean):
    rand = np.random.RandomState(42)
    t = 10 * rand.rand(30)
    theta = -1 + rand.rand(2 * nterms + 1)
    y = np.zeros(t.shape)
    if fit_mean:
        y = theta[0] * np.ones(t.shape)
    for i in range(1, nterms + 1):
        y += theta[2 * i - 1] * np.sin(2 * np.pi * i * freq * t)
        y += theta[2 * i] * np.cos(2 * np.pi * i * freq * t)

    y_fit = periodic_fit(t, y, dy=1, frequency=freq, t_fit=t, nterms=nterms,
                         center_data=False, fit_mean=fit_mean)
    assert_allclose(y, y_fit)
