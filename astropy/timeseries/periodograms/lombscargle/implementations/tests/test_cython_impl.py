from dataclasses import dataclass
from typing import TypeAlias

import numpy as np
import numpy.typing as npt
import pytest

from astropy.timeseries.periodograms.lombscargle.implementations.cython_impl import (
    lombscargle_cython,
)

_F64Array: TypeAlias = npt.NDArray[np.float64]


@dataclass(kw_only=True, slots=True, frozen=True)
class SyntheticSignal:
    t: _F64Array
    y: _F64Array
    dy: _F64Array
    true_frequency: float
    frequency_grid: _F64Array


@pytest.fixture
def perfect_sine_wave() -> SyntheticSignal:
    """Generates a mathematically perfect synthetic signal to verify period recovery."""
    t = np.linspace(0, 10, 101, dtype=np.float64)
    true_frequency = 0.5

    y = np.sin(2 * np.pi * true_frequency * t)
    dy = np.ones_like(t)

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
    assert np.any(np.isfinite(power))
    assert power.shape == perfect_sine_wave.frequency_grid.shape

    peak_idx = np.argmax(power)
    recovered_frequency = perfect_sine_wave.frequency_grid[peak_idx]

    assert np.isclose(recovered_frequency, perfect_sine_wave.true_frequency, rtol=1e-4)


@dataclass(kw_only=True, slots=True, frozen=True)
class NormMatrix:
    norm_type: str


@pytest.mark.parametrize(
    "matrix",
    [
        pytest.param(NormMatrix(norm_type="standard"), id="norm_standard"),
        pytest.param(NormMatrix(norm_type="model"), id="norm_model"),
        pytest.param(NormMatrix(norm_type="log"), id="norm_log"),
        pytest.param(NormMatrix(norm_type="psd"), id="norm_psd"),
    ],
)
def test_cython_lombscargle_normalizations(
    perfect_sine_wave: SyntheticSignal, matrix: NormMatrix
):
    """Validates that the C-engine successfully calculates all normalizations."""
    power = lombscargle_cython(
        perfect_sine_wave.t,
        perfect_sine_wave.y,
        perfect_sine_wave.dy,
        perfect_sine_wave.frequency_grid,
        normalization=matrix.norm_type,
    )

    assert type(power) is np.ndarray
    assert np.any(np.isfinite(power))


def test_cython_lombscargle_trapdoors(perfect_sine_wave: SyntheticSignal):
    """Triggers C-level validation trapdoors to ensure dimensions are strictly enforced."""
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
