from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import pytest

from astropy.timeseries.periodograms.lombscargle.implementations.cython_impl import (
    lombscargle_cython,
)

type _F64Array = npt.NDArray[np.float64]


@dataclass(kw_only=True, slots=True, frozen=True)
class SyntheticSignal:
    t: _F64Array
    y: _F64Array
    dy: _F64Array
    true_frequency: float
    frequency_grid: _F64Array


@pytest.fixture
def perfect_sine_wave() -> SyntheticSignal:
    """Generates a synthetic sine wave signal with a known period to verify power recovery."""
    t = np.linspace(0, 10, 101, dtype=np.float64)
    true_frequency = 0.5

    y = np.sin(2 * np.pi * true_frequency * t)
    dy = np.ones_like(t)

    # We tune the frequency grid to exactly 91 steps so the peak frequency (0.5)
    # lands perfectly on a discrete float coordinate. This avoids false-negative
    # Cython interpolation errors when comparing the recovered peak.
    frequency_grid = np.linspace(0.1, 1.0, 91, dtype=np.float64)

    return SyntheticSignal(
        t=t, y=y, dy=dy, true_frequency=true_frequency, frequency_grid=frequency_grid
    )


def test_cython_lombscargle_exact_recovery(perfect_sine_wave: SyntheticSignal):
    """Validates the Cython Lomb-Scargle engine mathematically recovers a known signal."""
    power = lombscargle_cython(
        perfect_sine_wave.t,
        perfect_sine_wave.y,
        perfect_sine_wave.dy,
        perfect_sine_wave.frequency_grid,
        normalization="standard",
        fit_mean=False,
        center_data=False,
        assume_regular_frequency=True,
    )

    assert type(power) is np.ndarray
    assert power.dtype == np.float64

    # Ensure the C-engine successfully allocated the memory and returned valid numbers
    # across the entire grid without overflowing into NaNs or infinities.
    assert np.all(np.isfinite(power))
    assert power.shape == perfect_sine_wave.frequency_grid.shape

    peak_idx = np.argmax(power)
    recovered_frequency = perfect_sine_wave.frequency_grid[peak_idx]

    assert np.isclose(recovered_frequency, perfect_sine_wave.true_frequency, rtol=1e-4)


@pytest.mark.parametrize("norm_type", ["standard", "model", "log", "psd"])
def test_cython_lombscargle_normalizations(
    perfect_sine_wave: SyntheticSignal, norm_type: str
):
    """Validates that the C-engine successfully calculates all normalizations."""
    power = lombscargle_cython(
        perfect_sine_wave.t,
        perfect_sine_wave.y,
        perfect_sine_wave.dy,
        perfect_sine_wave.frequency_grid,
        normalization=norm_type,
    )

    assert type(power) is np.ndarray
    assert np.all(np.isfinite(power))


def test_cython_lombscargle_exceptions(perfect_sine_wave: SyntheticSignal):
    """Triggers C-level validation exceptions to ensure dimensions are strictly enforced."""
    t_2d = np.ones((10, 2), dtype=np.float64)
    freq_2d = np.ones((10, 2), dtype=np.float64)

    with pytest.raises(ValueError, match="t, y, dy should be one dimensional"):
        lombscargle_cython(t_2d, t_2d, t_2d, perfect_sine_wave.frequency_grid)

    with pytest.raises(ValueError, match="frequency should be one-dimensional"):
        lombscargle_cython(
            perfect_sine_wave.t, perfect_sine_wave.y, perfect_sine_wave.dy, freq_2d
        )

    with pytest.raises(ValueError, match="normalization='invalid' not recognized"):
        lombscargle_cython(
            perfect_sine_wave.t,
            perfect_sine_wave.y,
            perfect_sine_wave.dy,
            perfect_sine_wave.frequency_grid,
            normalization="invalid",
        )
