import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_equal

from astropy.timeseries.periodograms.lombscargle._testing import (
    assert_not_strictly_equal,
)
from astropy.timeseries.periodograms.lombscargle.implementations.utils import (
    bitceil,
    extirpolate,
    trig_sum,
)


@pytest.mark.parametrize("N", 2 ** np.arange(1, 12))
@pytest.mark.parametrize("offset", [-1, 0, 1])
def test_bitceil(N, offset):
    assert_equal(bitceil(N + offset), int(2 ** np.ceil(np.log2(N + offset))))


@pytest.fixture
def extirpolate_data():
    rng = np.random.default_rng(0)
    x = 100 * rng.random(50)
    y = np.sin(x)
    f = lambda x: np.sin(x / 10)
    return x, y, f


@pytest.mark.parametrize("N", [100, None])
@pytest.mark.parametrize("M", [5])
def test_extirpolate(N, M, extirpolate_data):
    x, y, f = extirpolate_data
    y_hat = extirpolate(x, y, N, M)
    x_hat = np.arange(len(y_hat))
    assert_allclose(np.dot(f(x), y), np.dot(f(x_hat), y_hat), rtol=1.5e-5)


@pytest.fixture
def extirpolate_int_data():
    rng = np.random.default_rng(0)
    x = 100 * rng.random(50)
    x[:25] = x[:25].astype(int)
    y = np.sin(x)
    f = lambda x: np.sin(x / 10)
    return x, y, f


@pytest.mark.parametrize("N", [100, None])
@pytest.mark.parametrize("M", [5])
def test_extirpolate_with_integers(N, M, extirpolate_int_data):
    x, y, f = extirpolate_int_data
    y_hat = extirpolate(x, y, N, M)
    x_hat = np.arange(len(y_hat))
    assert_allclose(np.dot(f(x), y), np.dot(f(x_hat), y_hat), rtol=1.7e-5)


@pytest.fixture
def trig_sum_data():
    rng = np.random.default_rng(0)
    t = 10 * rng.random(50)
    h = np.sin(t)
    return t, h


@pytest.mark.parametrize("f0", [0, 1])
@pytest.mark.parametrize("adjust_t", [True, False])
@pytest.mark.parametrize("freq_factor", [1, 2])
@pytest.mark.parametrize("df", [0.1])
def test_trig_sum(f0, adjust_t, freq_factor, df, trig_sum_data):
    t, h = trig_sum_data

    tfit = t - t.min() if adjust_t else t

    kwargs_ref = {"use_fft": False, "f0": f0, "freq_factor": freq_factor}
    S0, C0 = trig_sum(tfit, h, df, N=1000, **kwargs_ref)

    results = [(S0, C0)]

    fast_methods = [
        (
            {
                "oversampling": 10,
                "algorithm": "fasper",
                "f0": f0,
                "freq_factor": freq_factor,
            },
            1e-2,
        ),
        (
            {"algorithm": "lra", "f0": f0, "freq_factor": freq_factor},
            1e-7,
        ),
    ]

    for alt_kwargs, tol in fast_methods:
        S, C = trig_sum(tfit, h, df, N=1000, **alt_kwargs)

        assert_allclose(S0, S, atol=tol)
        assert_allclose(C0, C, atol=tol)

        for prev_S, prev_C in results:
            assert_not_strictly_equal(prev_S, S)
            assert_not_strictly_equal(prev_C, C)

        results.append((S, C))
