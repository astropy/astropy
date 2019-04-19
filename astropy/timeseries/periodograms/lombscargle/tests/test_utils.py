import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy.timeseries.periodograms.lombscargle.utils import convert_normalization, compute_chi2_ref
from astropy.timeseries.periodograms.lombscargle.core import LombScargle


NORMALIZATIONS = ['standard', 'model', 'log', 'psd']


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


@pytest.mark.parametrize('norm_in', NORMALIZATIONS)
@pytest.mark.parametrize('norm_out', NORMALIZATIONS)
def test_convert_normalization(norm_in, norm_out, data):
    t, y, dy = data

    _, power_in = LombScargle(t, y, dy).autopower(maximum_frequency=5,
                                                  normalization=norm_in)
    _, power_out = LombScargle(t, y, dy).autopower(maximum_frequency=5,
                                                   normalization=norm_out)
    power_in_converted = convert_normalization(power_in, N=len(t),
                                               from_normalization=norm_in,
                                               to_normalization=norm_out,
                                               chi2_ref = compute_chi2_ref(y, dy))
    assert_allclose(power_in_converted, power_out)
