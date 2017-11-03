import numpy as np
import pytest
from numpy.testing import assert_allclose

try:
    import scipy
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True

from .. import LombScargle
from .._statistics import (cdf_single, pdf_single, fap_single, inv_fap_single,
                           METHODS)
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
@pytest.mark.parametrize('with_errors', [True, False])
def test_distribution(null_data, normalization, with_errors, fmax=40):
    t, y, dy = null_data
    if not with_errors:
        dy = None

    N = len(t)
    ls = LombScargle(t, y, dy, normalization=normalization)
    freq, power = ls.autopower(maximum_frequency=fmax)
    z = np.linspace(0, power.max(), 1000)

    # Test that pdf and cdf are consistent
    dz = z[1] - z[0]
    z_mid = z[:-1] + 0.5 * dz
    pdf = ls.distribution(z_mid)
    cdf = ls.distribution(z, cumulative=True)
    assert_allclose(pdf, np.diff(cdf) / dz, rtol=1E-5, atol=1E-8)

    # psd normalization without specified errors produces bad results
    if not (normalization == 'psd' and not with_errors):
        # Test that observed power is distributed according to the theoretical pdf
        hist, bins = np.histogram(power, 30, normed=True)
        midpoints = 0.5 * (bins[1:] + bins[:-1])
        pdf = ls.distribution(midpoints)
        assert_allclose(hist, pdf, rtol=0.05, atol=0.05 * pdf[0])


@pytest.mark.parametrize('N', [10, 100, 1000])
@pytest.mark.parametrize('normalization', NORMALIZATIONS)
def test_inverse_single(N, normalization):
    fap = np.linspace(0, 1, 100)

    z = inv_fap_single(fap, N, normalization)
    fap_out = fap_single(z, N, normalization)
    assert_allclose(fap, fap_out)


@pytest.mark.parametrize('normalization', NORMALIZATIONS)
@pytest.mark.parametrize('use_errs', [True, False])
def test_inverse_bootstrap(null_data, normalization, use_errs, fmax=5):
    t, y, dy = null_data
    if not use_errs:
        dy = None

    fap = np.linspace(0, 1, 10)
    method = 'bootstrap'
    method_kwds = METHOD_KWDS['bootstrap']

    ls = LombScargle(t, y, dy, normalization=normalization)

    z = ls.false_alarm_level(fap, maximum_frequency=fmax,
                             method=method, method_kwds=method_kwds)
    fap_out = ls.false_alarm_probability(z, maximum_frequency=fmax,
                                         method=method,
                                         method_kwds=method_kwds)

    # atol = 1 / n_bootstraps
    assert_allclose(fap, fap_out, atol=0.05)


@pytest.mark.parametrize('method', sorted(set(METHODS) - {'bootstrap'}))
@pytest.mark.parametrize('normalization', NORMALIZATIONS)
@pytest.mark.parametrize('use_errs', [True, False])
@pytest.mark.parametrize('N', [10, 100, 1000])
def test_inverses(method, normalization, use_errs, N, T=5, fmax=5):
    if not HAS_SCIPY and method in ['baluev', 'davies']:
        pytest.skip("SciPy required")

    t, y, dy = data(N, rseed=543)
    if not use_errs:
        dy = None
    method_kwds = METHOD_KWDS.get(method, None)

    fap = np.logspace(-10, 0, 10)

    ls = LombScargle(t, y, dy, normalization=normalization)
    z = ls.false_alarm_level(fap, maximum_frequency=fmax,
                             method=method,
                             method_kwds=method_kwds)
    fap_out = ls.false_alarm_probability(z, maximum_frequency=fmax,
                                         method=method,
                                         method_kwds=method_kwds)
    assert_allclose(fap, fap_out)


@pytest.mark.parametrize('method', sorted(METHODS))
@pytest.mark.parametrize('normalization', NORMALIZATIONS)
def test_false_alarm_smoketest(method, normalization, data):
    if not HAS_SCIPY and method in ['baluev', 'davies']:
        pytest.skip("SciPy required")

    kwds = METHOD_KWDS.get(method, None)
    t, y, dy = data
    fmax = 5

    ls = LombScargle(t, y, dy, normalization=normalization)
    freq, power = ls.autopower(maximum_frequency=fmax)
    Z = np.linspace(power.min(), power.max(), 30)

    fap = ls.false_alarm_probability(Z, maximum_frequency=fmax,
                                     method=method, method_kwds=kwds)

    assert len(fap) == len(Z)
    if method != 'davies':
        assert np.all(fap <= 1)
        assert np.all(fap[:-1] >= fap[1:])  # monotonically decreasing


@pytest.mark.parametrize('method', sorted(METHODS))
@pytest.mark.parametrize('use_errs', [True, False])
@pytest.mark.parametrize('normalization', sorted(set(NORMALIZATIONS) - {'psd'}))
def test_false_alarm_equivalence(method, normalization, use_errs, data):
    # Note: the PSD normalization is not equivalent to the others, in that it
    # depends on the absolute errors rather than relative errors. Because the
    # scaling contributes to the distribution, it cannot be converted directly
    # from any of the three normalized versions.
    if not HAS_SCIPY and method in ['baluev', 'davies']:
        pytest.skip("SciPy required")

    kwds = METHOD_KWDS.get(method, None)
    t, y, dy = data
    if not use_errs:
        dy = None
    fmax = 5

    ls = LombScargle(t, y, dy, normalization=normalization)
    freq, power = ls.autopower(maximum_frequency=fmax)
    Z = np.linspace(power.min(), power.max(), 30)
    fap = ls.false_alarm_probability(Z, maximum_frequency=fmax,
                                     method=method, method_kwds=kwds)

    # Compute the equivalent Z values in the standard normalization
    # and check that the FAP is consistent
    Z_std = convert_normalization(Z, len(t),
                                  from_normalization=normalization,
                                  to_normalization='standard',
                                  chi2_ref=compute_chi2_ref(y, dy))
    ls = LombScargle(t, y, dy, normalization='standard')
    fap_std = ls.false_alarm_probability(Z_std, maximum_frequency=fmax,
                                         method=method, method_kwds=kwds)

    assert_allclose(fap, fap_std, rtol=0.1)
