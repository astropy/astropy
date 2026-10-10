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


def test_pdm_recovers_period(data):
    t, y = data

    pdm = PhaseDispersionMinimization(t, y)
    periods = np.linspace(0.5, 1.5, 201)
    theta = pdm.run(periods)

    min_period = periods[np.argmin(theta)]

    assert_allclose(min_period, 1.0, atol=0.01)


def test_pdm_theta():
    t = np.array(
        [
            0.0,
            0.1,
            1.0,
            1.1,
            0.5,
            0.6,
            1.5,
            1.6,
        ]
    )

    y = np.array(
        [
            1.0,
            2.0,
            1.0,
            2.0,
            10.0,
            12.0,
            10.0,
            12.0,
        ]
    )

    pdm = PhaseDispersionMinimization(t, y)

    theta = pdm.run(
        np.array([1.0]),
        bins=2,
    )

    assert_allclose(theta[0], 5 / 159)


def test_pdm_perfect_period():

    # for test purposes we put the data in the center of the bins
    phase = (np.arange(10) + 0.5) / 10

    t = np.concatenate(
        [
            phase,
            phase + 1,
            phase + 2,
            phase + 3,
        ]
    )

    y_cycle = np.sin(2 * np.pi * phase)
    y = np.tile(y_cycle, 4)

    pdm = PhaseDispersionMinimization(t, y)

    theta = pdm.run(
        np.array([1.0]),
        bins=10,
    )

    # since every bin has a single value the variance is zero and therefore theta is zero
    assert_allclose(theta[0], 0.0, atol=1e-14)
