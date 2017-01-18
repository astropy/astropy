import numpy as np
from numpy.testing import assert_allclose

from .... import units
from ....tests.helper import pytest
from .. import LombScargle, statistics


NORMALIZATIONS = ['standard', 'psd', 'log', 'model']


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
    pdf = statistics.pdf(z_mid, N, normalization=normalization)
    cdf = statistics.cdf(z, N, normalization=normalization)
    assert_allclose(pdf, np.diff(cdf) / dz, rtol=1E-5, atol=1E-8)

    # Test that observed power is distributed according to the theoretical pdf
    hist, bins = np.histogram(power, 30, normed=True)
    midpoints = 0.5 * (bins[1:] + bins[:-1])
    pdf = statistics.pdf(midpoints, N, normalization=normalization)
    assert_allclose(hist, pdf, rtol=0.05, atol=0.05 * pdf[0])


def test_distribution_units(null_data):
    t, y, dy = null_data
    N = len(t)

    t_days = t * units.day
    y_mag = y * units.mag
    dy_mag = dy * units.mag

    # unnormalized: this should fail
    ls = LombScargle(t_days, y_mag)
    freq, power = ls.autopower(normalization='psd',
                               maximum_frequency=1 / units.day)
    assert power.unit == units.mag ** 2

    with pytest.raises(ValueError) as err:
        pdf = statistics.pdf(power, N, normalization='psd')
    assert str(err.value).startswith('The distribution can be computed')

    with pytest.raises(ValueError) as err:
        cdf = statistics.cdf(power, N, normalization='psd')
    assert str(err.value).startswith('The distribution can be computed')

    # normalized: case with units should match case without units
    ls = LombScargle(t_days, y_mag, dy_mag)
    freq, power = ls.autopower(normalization='psd',
                               maximum_frequency=1 / units.day)
    assert power.unit == units.dimensionless_unscaled

    pdf_with_units = statistics.pdf(power, N, normalization='psd')
    cdf_with_units = statistics.cdf(power, N, normalization='psd')

    ls = LombScargle(t, y, dy)
    freq, power = ls.autopower(normalization='psd',
                               maximum_frequency=1)
    pdf_no_units = statistics.pdf(power, N, normalization='psd')
    cdf_no_units = statistics.cdf(power, N, normalization='psd')

    assert_allclose(pdf_with_units, pdf_no_units)
    assert_allclose(cdf_with_units, cdf_no_units)
