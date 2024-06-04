import numpy as np
import pytest
from numpy.testing import assert_allclose

import astropy.units as u
from astropy.timeseries.periodograms.lombscargle import LombScargle
from astropy.timeseries.periodograms.lombscargle._statistics import (
    METHODS,
    fap_single,
    inv_fap_single,
)
from astropy.timeseries.periodograms.lombscargle.utils import (
    compute_chi2_ref,
    convert_normalization,
)
from astropy.utils.compat.optional_deps import HAS_SCIPY

METHOD_KWDS = dict(bootstrap={"n_bootstraps": 20, "random_seed": 42})
NORMALIZATIONS = ["standard", "psd", "log", "model"]


def make_data(N=100, period=1, theta=[10, 2, 3], dy=1, rseed=0, units=False):
    """Generate some data for testing"""
    rng = np.random.default_rng(rseed)
    t = 5 * period * rng.random(N)
    omega = 2 * np.pi / period
    y = theta[0] + theta[1] * np.sin(omega * t) + theta[2] * np.cos(omega * t)
    dy = dy * (0.5 + rng.random(N))
    y += dy * rng.standard_normal(N)

    fmax = 5

    if units:
        return t * u.day, y * u.mag, dy * u.mag, fmax / u.day
    else:
        return t, y, dy, fmax


def null_data(N=1000, dy=1, rseed=0, units=False):
    """Generate null hypothesis data"""
    rng = np.random.default_rng(rseed)
    t = 100 * rng.random(N)
    dy = 0.5 * dy * (1 + rng.random(N))
    y = dy * rng.standard_normal(N)
    fmax = 40

    if units:
        return t * u.day, y * u.mag, dy * u.mag, fmax / u.day
    else:
        return t, y, dy, fmax


@pytest.mark.parametrize("normalization", NORMALIZATIONS)
@pytest.mark.parametrize("with_errors", [True, False])
@pytest.mark.parametrize("units", [False, True])
def test_distribution(normalization, with_errors, units):
    t, y, dy, fmax = null_data(units=units)

    if not with_errors:
        dy = None

    ls = LombScargle(t, y, dy, normalization=normalization)
    freq, power = ls.autopower(maximum_frequency=fmax)
    z = np.linspace(0, power.max(), 1000)

    # Test that pdf and cdf are consistent
    dz = z[1] - z[0]
    z_mid = z[:-1] + 0.5 * dz
    pdf = ls.distribution(z_mid)
    cdf = ls.distribution(z, cumulative=True)
    if isinstance(dz, u.Quantity):
        dz = dz.value
    assert_allclose(pdf, np.diff(cdf) / dz, rtol=1e-5, atol=1e-8)

    # psd normalization without specified errors produces bad results
    if not (normalization == "psd" and not with_errors):
        # Test that observed power is distributed according to the theoretical pdf
        hist, bins = np.histogram(power, 30, density=True)
        midpoints = 0.5 * (bins[1:] + bins[:-1])
        pdf = ls.distribution(midpoints)
        assert_allclose(hist, pdf, rtol=0.05, atol=0.05 * pdf[0])


@pytest.mark.parametrize("N", [10, 100, 1000])
@pytest.mark.parametrize("normalization", NORMALIZATIONS)
def test_inverse_single(N, normalization):
    fap = np.linspace(0, 1, 11)

    z = inv_fap_single(fap, N, normalization)
    fap_out = fap_single(z, N, normalization)
    assert_allclose(fap, fap_out)


@pytest.mark.parametrize("normalization", NORMALIZATIONS)
@pytest.mark.parametrize("use_errs", [True, False])
@pytest.mark.parametrize("units", [False, True])
def test_inverse_bootstrap(normalization, use_errs, units):
    t, y, dy, fmax = null_data(units=units)
    if not use_errs:
        dy = None

    fap = np.linspace(0, 1, 11)
    method = "bootstrap"
    method_kwds = METHOD_KWDS["bootstrap"]

    ls = LombScargle(t, y, dy, normalization=normalization)

    z = ls.false_alarm_level(
        fap, maximum_frequency=fmax, method=method, method_kwds=method_kwds
    )
    fap_out = ls.false_alarm_probability(
        z, maximum_frequency=fmax, method=method, method_kwds=method_kwds
    )

    # atol = 1 / n_bootstraps
    assert_allclose(fap, fap_out, atol=0.05)


@pytest.mark.parametrize("method", sorted(set(METHODS) - {"bootstrap"}))
@pytest.mark.parametrize("normalization", NORMALIZATIONS)
@pytest.mark.parametrize("use_errs", [True, False])
@pytest.mark.parametrize("N", [10, 100, 1000])
@pytest.mark.parametrize("units", [False, True])
def test_inverses(method, normalization, use_errs, N, units, T=5):
    if not HAS_SCIPY and method in ["baluev", "davies"]:
        pytest.skip("SciPy required")

    t, y, dy, fmax = make_data(N, rseed=543, units=units)
    if not use_errs:
        dy = None
    method_kwds = METHOD_KWDS.get(method)

    fap = np.logspace(-10, 0, 11)

    ls = LombScargle(t, y, dy, normalization=normalization)
    z = ls.false_alarm_level(
        fap, maximum_frequency=fmax, method=method, method_kwds=method_kwds
    )
    fap_out = ls.false_alarm_probability(
        z, maximum_frequency=fmax, method=method, method_kwds=method_kwds
    )
    assert_allclose(fap, fap_out)


@pytest.mark.parametrize("method", sorted(METHODS))
@pytest.mark.parametrize("normalization", NORMALIZATIONS)
@pytest.mark.parametrize("units", [False, True])
def test_false_alarm_smoketest(method, normalization, units):
    if not HAS_SCIPY and method in ["baluev", "davies"]:
        pytest.skip("SciPy required")

    kwds = METHOD_KWDS.get(method)
    t, y, dy, fmax = make_data(rseed=42, units=units)

    ls = LombScargle(t, y, dy, normalization=normalization)
    freq, power = ls.autopower(maximum_frequency=fmax)
    Z = np.linspace(power.min(), power.max(), 30)

    fap = ls.false_alarm_probability(
        Z, maximum_frequency=fmax, method=method, method_kwds=kwds
    )

    assert len(fap) == len(Z)
    if method != "davies":
        assert np.all(fap <= 1)
        assert np.all(fap[:-1] >= fap[1:])  # monotonically decreasing


@pytest.mark.parametrize("method", sorted(METHODS))
@pytest.mark.parametrize("use_errs", [True, False])
@pytest.mark.parametrize("normalization", sorted(set(NORMALIZATIONS) - {"psd"}))
@pytest.mark.parametrize("units", [False, True])
def test_false_alarm_equivalence(method, normalization, use_errs, units):
    # Note: the PSD normalization is not equivalent to the others, in that it
    # depends on the absolute errors rather than relative errors. Because the
    # scaling contributes to the distribution, it cannot be converted directly
    # from any of the three normalized versions.
    if not HAS_SCIPY and method in ["baluev", "davies"]:
        pytest.skip("SciPy required")

    kwds = METHOD_KWDS.get(method)
    t, y, dy, fmax = make_data(units=units)
    if not use_errs:
        dy = None

    ls = LombScargle(t, y, dy, normalization=normalization)
    freq, power = ls.autopower(maximum_frequency=fmax)
    Z = np.linspace(power.min(), power.max(), 30)
    fap = ls.false_alarm_probability(
        Z, maximum_frequency=fmax, method=method, method_kwds=kwds
    )

    # Compute the equivalent Z values in the standard normalization
    # and check that the FAP is consistent
    Z_std = convert_normalization(
        Z,
        len(t),
        from_normalization=normalization,
        to_normalization="standard",
        chi2_ref=compute_chi2_ref(y, dy),
    )
    ls = LombScargle(t, y, dy, normalization="standard")
    fap_std = ls.false_alarm_probability(
        Z_std, maximum_frequency=fmax, method=method, method_kwds=kwds
    )

    assert_allclose(fap, fap_std, rtol=0.1)
