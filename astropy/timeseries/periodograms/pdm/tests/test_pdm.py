import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy.timeseries.periodograms.pdm import PhaseDispersionMinimization


@pytest.fixture
def data(n_samples=1000, period=1.0, offset=1.0, amplitude=2.0, dy=0.05, rseed=0):
    rng = np.random.default_rng(rseed)

    t = 20 * period * rng.random(n_samples)
    y = offset + amplitude * np.sin(2 * np.pi * t / period)

    y += dy * rng.standard_normal(size=n_samples)

    return t, y


def test_pdm(data):
    t, y = data

    pdm = PhaseDispersionMinimization(t, y)
    periods = np.linspace(0.5, 1.5, 201)
    theta = pdm.run(periods)

    min_period = periods[np.argmin(theta)]

    assert_allclose(min_period, 1.0, atol=0.01)
