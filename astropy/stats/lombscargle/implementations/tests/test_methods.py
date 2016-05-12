import numpy as np
from numpy.testing import assert_allclose

from .....tests.helper import pytest
from .. import lombscargle, available_methods

ALL_METHODS = available_methods()
BIAS_METHODS = [method for method in ALL_METHODS if method != 'scipy']
FAST_METHODS = [method for method in ALL_METHODS if 'fast' in method]
NTERMS_METHODS = ['auto'] + [method for method in ALL_METHODS
                             if 'chi2' in method]


@pytest.fixture
def data(N=100, period=1, theta=[10, 2, 3], dy=1, rseed=0):
    """Generate some data for testing"""
    rng = np.random.RandomState(rseed)
    t = 20 * period * rng.rand(N)
    omega = 2 * np.pi / period
    y = theta[0] + theta[1] * np.sin(omega * t) + theta[2] * np.cos(omega * t)
    dy = dy * (0.5 + rng.rand(N))
    y += dy * rng.randn(N)

    return t, y, dy


@pytest.mark.parametrize('method', ALL_METHODS)
@pytest.mark.parametrize('center_data', [True, False])
@pytest.mark.parametrize('fit_bias', [True, False])
@pytest.mark.parametrize('with_errors', [True, False])
@pytest.mark.parametrize('normalization', ['standard', 'psd'])
def test_lombscargle_common(method, center_data, fit_bias,
                            with_errors, normalization, data):
    if fit_bias and method not in BIAS_METHODS:
        return
    if with_errors and method == 'scipy':
        return

    t, y, dy = data
    if not with_errors:
        dy = None

    freq = 0.8 + 0.01 * np.arange(40)

    kwds = dict(normalization=normalization,
                center_data=center_data,
                fit_bias=fit_bias)
    expected_output = lombscargle(t, y, dy, frequency=freq, **kwds)

    # don't test fast fft methods here
    if 'fast' in method:
        kwds['method_kwds'] = dict(use_fft=False)

    # check that output matches that of the "auto" method
    output = lombscargle(t, y, dy, frequency=freq, method=method, **kwds)
    assert_allclose(output, expected_output, rtol=1E-7, atol=1E-20)

    # check that output of dy=None matches output of dy=[array of ones]
    output_dy_None = lombscargle(t, y, dy=None,
                                 frequency=freq, method=method, **kwds)
    output_dy_ones = lombscargle(t, y, dy=np.ones_like(t),
                                 frequency=freq, method=method, **kwds)
    assert_allclose(output_dy_None, output_dy_ones)


@pytest.mark.parametrize('method', NTERMS_METHODS)
@pytest.mark.parametrize('center_data', [True, False])
@pytest.mark.parametrize('fit_bias', [True, False])
@pytest.mark.parametrize('with_errors', [True, False])
@pytest.mark.parametrize('nterms', range(5))
@pytest.mark.parametrize('normalization', ['standard', 'psd'])
def test_lombscargle_nterms(method, center_data, fit_bias, with_errors, nterms,
                            normalization, data):
    t, y, dy = data
    if not with_errors:
        dy = None

    freq = 0.8 + 0.01 * np.arange(40)

    kwds = dict(normalization=normalization,
                center_data=center_data,
                fit_bias=fit_bias,
                nterms=nterms)

    if nterms == 0 and not fit_bias:
        with pytest.raises(ValueError) as err:
            lombscargle(t, y, frequency=freq, method=method, **kwds)
        assert 'nterms' in str(err.value) and 'bias' in str(err.value)
    else:
        expected_output = lombscargle(t, y, dy, frequency=freq, **kwds)

        # don't test fast fft methods here
        if 'fast' in method:
            kwds['method_kwds'] = dict(use_fft=False)
        output = lombscargle(t, y, dy, frequency=freq, method=method, **kwds)

        assert_allclose(output, expected_output, rtol=1E-7, atol=1E-20)


@pytest.mark.parametrize('method', FAST_METHODS)
@pytest.mark.parametrize('center_data', [True, False])
@pytest.mark.parametrize('fit_bias', [True, False])
@pytest.mark.parametrize('with_errors', [True, False])
@pytest.mark.parametrize('nterms', range(4))
def test_fast_methods(method, center_data, fit_bias, with_errors,
                      nterms, data):
    # leave out normalization here because we judge via absolute tolerance
    t, y, dy = data
    if not with_errors:
        dy = None

    freq = 0.8 + 0.01 * np.arange(40)

    kwds = dict(method=method,
                center_data=center_data,
                fit_bias=fit_bias,
                nterms=nterms)

    if method == 'fast' and nterms != 1:
        with pytest.raises(ValueError) as err:
            lombscargle(t, y, frequency=freq, **kwds)
        assert 'nterms' in str(err.value)
        return

    if nterms == 0 and not fit_bias:
        with pytest.raises(ValueError) as err:
            lombscargle(t, y, frequency=freq, **kwds)
        assert 'nterms' in str(err.value) and 'bias' in str(err.value)
        return

    output_slow = lombscargle(t, y, dy, frequency=freq,
                              method_kwds=dict(use_fft=False),
                              **kwds)
    output_fast = lombscargle(t, y, dy, frequency=freq,
                              method_kwds=dict(use_fft=True),
                              **kwds)

    assert_allclose(output_slow, output_fast, atol=0.008)
