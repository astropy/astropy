import numpy as np
from numpy.testing import assert_allclose

from ....tests.helper import pytest
from .. import LombScargle
from .. import statistics
from ..statistics import (cdf_single, pdf_single, METHODS,
                          false_alarm_probability,
                          false_alarm_level)
from ..utils import convert_normalization, compute_chi2_ref

METHOD_KWDS = dict(bootstrap={'n_bootstraps': 20, 'random_seed': 42})
NORMALIZATIONS = ['standard', 'psd', 'log', 'model']


@pytest.fixture
def data(N=100, period=1, theta=[10, 2, 3], dy=1, rseed=0):
    """Generate some data for testing"""
    rng = np.random.RandomState(rseed)
    t = 5 * period * rng.rand(N)
    omega = 2 * np.pi / period
    y = theta[0] + theta[1] * np.sin(omega * t) + theta[2] * np.cos(omega * t)
    dy = dy * (0.5 + rng.rand(N))
    y += dy * rng.randn(N)

    return t, y, dy


@pytest.fixture
def null_data(N=1000, dy=1, rseed=0):
    """Generate null hypothesis data"""
    rng = np.random.RandomState(rseed)
    t = 100 * rng.rand(N)
    dy = 0.5 * dy * (1 + rng.rand(N))
    y = dy * rng.randn(N)
    return t, y, dy


@pytest.mark.parametrize('normalization', NORMALIZATIONS)
def test_distribution(null_data, normalization):
    t, y, dy = null_data
    N = len(t)
    ls = LombScargle(t, y, dy)
    freq, power = ls.autopower(normalization=normalization,
                               maximum_frequency=40)
    z = np.linspace(0, power.max(), 1000)

    # Test that pdf and cdf are consistent
    dz = z[1] - z[0]
    z_mid = z[:-1] + 0.5 * dz
    pdf = pdf_single(z_mid, N, normalization=normalization)
    cdf = cdf_single(z, N, normalization=normalization)
    assert_allclose(pdf, np.diff(cdf) / dz, rtol=1E-5, atol=1E-8)

    # Test that observed power is distributed according to the theoretical pdf
    hist, bins = np.histogram(power, 30, normed=True)
    midpoints = 0.5 * (bins[1:] + bins[:-1])
    pdf = pdf_single(midpoints, N, normalization=normalization)
    assert_allclose(hist, pdf, rtol=0.05, atol=0.05 * pdf[0])


@pytest.mark.parametrize('N', [10, 100, 1000])
@pytest.mark.parametrize('normalization', NORMALIZATIONS)
def test_inverse_single(N, normalization):
    fap = np.linspace(0, 1, 100)

    z = statistics.inv_fap_single(fap, N, normalization)
    fap_out = statistics.fap_single(z, N, normalization)
    assert_allclose(fap, fap_out)


@pytest.mark.parametrize('normalization', NORMALIZATIONS)
def test_inverse_bootstrap(null_data, normalization, fmax=5):
    t, y, dy = null_data

    fap = np.linspace(0, 1, 10)
    method = 'bootstrap'
    method_kwds = METHOD_KWDS['bootstrap']

    z = false_alarm_level(fap, fmax, t, y, dy, normalization,
                          method=method, method_kwds=method_kwds)
    fap_out = false_alarm_probability(z, fmax, t, y, dy, normalization,
                                      method=method, method_kwds=method_kwds)

    # atol = 1 / n_bootstraps
    assert_allclose(fap, fap_out, atol=0.05)


@pytest.mark.parametrize('method', set(METHODS) - {'bootstrap'})
@pytest.mark.parametrize('normalization', NORMALIZATIONS)
@pytest.mark.parametrize('N', [10, 100, 1000])
def test_inverses(method, normalization, N, T=5, fmax=5):
    t, y, dy = data(N, rseed=543)
    method_kwds = METHOD_KWDS.get(method, None)

    fap = np.logspace(-10, 0, 10)

    z = false_alarm_level(fap, fmax, t, y, dy, normalization,
                          method=method, method_kwds=method_kwds)
    fap_out = false_alarm_probability(z, fmax, t, y, dy, normalization,
                                      method=method, method_kwds=method_kwds)

    assert_allclose(fap, fap_out)


@pytest.mark.parametrize('method', METHODS)
@pytest.mark.parametrize('normalization', NORMALIZATIONS)
def test_false_alarm_smoketest(method, normalization, data):
    kwds = METHOD_KWDS.get(method, None)
    t, y, dy = data
    fmax = 5

    freq, power = LombScargle(t, y, dy).autopower(normalization=normalization,
                                                  maximum_frequency=fmax)
    Z = np.linspace(power.min(), power.max(), 30)
    print(normalization)
    print(Z.max())

    fap = false_alarm_probability(Z, fmax, t, y, dy,
                                  normalization=normalization,
                                  method=method,
                                  method_kwds=METHOD_KWDS.get(method, None))
    assert len(fap) == len(Z)
    if method != 'davies':
        assert np.all(fap <= 1)
        assert np.all(fap[:-1] >= fap[1:])  # monotonically decreasing


@pytest.mark.parametrize('method', METHODS)
@pytest.mark.parametrize('normalization', set(NORMALIZATIONS) - {'psd'})
def test_false_alarm_equivalence(method, normalization, data):
    # Note: the PSD normalization is not equivalent to the others, in that it
    # depends on the absolute errors rather than relative errors. Because the
    # scaling contributes to the distribution, it cannot be converted directly
    # from any of the three normalized versions.

    kwds = METHOD_KWDS.get(method, None)
    t, y, dy = data
    fmax = 5

    freq, power = LombScargle(t, y, dy).autopower(normalization=normalization,
                                                  maximum_frequency=fmax)
    Z = np.linspace(power.min(), power.max(), 30)
    fap = false_alarm_probability(Z, fmax, t, y, dy,
                                  normalization=normalization,
                                  method=method,
                                  method_kwds=METHOD_KWDS.get(method, None))

    # Compute the equivalent Z values in the standard normalization
    # and check that the FAP is consistent
    Z_std = convert_normalization(Z, len(t),
                                  from_normalization=normalization,
                                  to_normalization='standard',
                                  chi2_ref=compute_chi2_ref(y, dy))
    fap_std = false_alarm_probability(Z_std, fmax, t, y, dy,
                                      normalization='standard',
                                      method=method,
                                      method_kwds=METHOD_KWDS.get(method, None))

    assert_allclose(fap, fap_std, rtol=0.1)
