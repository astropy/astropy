"""Tests for the low-level Box Least Squares (BLS) C-extension implementation."""

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import pytest

from astropy.timeseries.periodograms.bls._impl import bls_impl

_F64Array = npt.NDArray[np.float64]


@dataclass(kw_only=True, slots=True, frozen=True)
class SyntheticTransit:
    """Dataclass for a synthetic transit."""

    t: _F64Array
    y: _F64Array
    ivar: _F64Array
    true_period: float
    true_duration: float
    true_depth: float


@pytest.fixture
def perfect_transit() -> SyntheticTransit:
    """Generate a synthetic transit signal."""
    true_period = 2.0
    true_duration = 0.1
    true_depth = 0.5

    t = np.linspace(0, 10, 1000, dtype=np.float64)
    y = np.ones_like(t, dtype=np.float64)
    ivar = np.ones_like(t, dtype=np.float64) * 1e6

    phase = t % true_period
    in_transit = (phase > (true_period / 2 - true_duration / 2)) & (
        phase < (true_period / 2 + true_duration / 2)
    )
    y[in_transit] -= true_depth

    return SyntheticTransit(
        t=t,
        y=y,
        ivar=ivar,
        true_period=true_period,
        true_duration=true_duration,
        true_depth=true_depth,
    )


def test_bls_exact_types_and_shapes(perfect_transit: SyntheticTransit) -> None:
    """Check memory allocation, types, and shapes."""
    periods = np.array([1.9, 2.0, 2.1], dtype=np.float64)
    durations = np.array([0.05, 0.1, 0.15], dtype=np.float64)

    results = bls_impl(
        perfect_transit.t,
        perfect_transit.y,
        perfect_transit.ivar,
        periods,
        durations,
        10,  # oversample
        0,  # obj_flag
    )

    assert type(results) is tuple
    assert len(results) == 7

    for array_out in results:
        assert type(array_out) is np.ndarray
        assert array_out.dtype == np.float64
        assert array_out.shape == periods.shape


def test_bls_mathematical_recovery(perfect_transit: SyntheticTransit) -> None:
    """Check if the C-engine recovers the planted signal."""
    periods = np.linspace(1.5, 2.5, 100, dtype=np.float64)
    durations = np.linspace(0.05, 0.2, 10, dtype=np.float64)

    results = bls_impl(
        perfect_transit.t,
        perfect_transit.y,
        perfect_transit.ivar,
        periods,
        durations,
        10,
        0,
    )

    out_objective, out_depth, _, out_duration, _, _, _ = results
    best_idx = np.argmax(out_objective)

    np.testing.assert_allclose(
        periods[best_idx], perfect_transit.true_period, rtol=1e-2
    )
    np.testing.assert_allclose(
        out_duration[best_idx], perfect_transit.true_duration, rtol=1e-1
    )
    np.testing.assert_allclose(
        out_depth[best_idx], perfect_transit.true_depth, rtol=5e-2
    )


@pytest.mark.parametrize(
    ("bad_periods", "bad_durations", "expected_err_msg"),
    [
        pytest.param(
            np.array([0.0], dtype=np.float64),
            np.array([0.1], dtype=np.float64),
            "invalid inputs",
            id="zero_period_flag_1",
        ),
        pytest.param(
            np.array([1.0], dtype=np.float64),
            np.array([1.5], dtype=np.float64),
            "invalid inputs",
            id="duration_exceeds_period_flag_2",
        ),
    ],
)
def test_bls_c_level_exceptions(
    perfect_transit: SyntheticTransit,
    bad_periods: _F64Array,
    bad_durations: _F64Array,
    expected_err_msg: str,
) -> None:
    """Test C-level early-exit flags."""
    with pytest.raises(ValueError) as exc_info:
        bls_impl(
            perfect_transit.t,
            perfect_transit.y,
            perfect_transit.ivar,
            bad_periods,
            bad_durations,
            10,
            0,
        )

    assert expected_err_msg in str(exc_info.value).lower()
