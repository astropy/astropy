# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for the _sigma_clip_fast C extension in astropy.stats.

Note: The underlying C implementation assumes input data arrays are sorted.
"""

import warnings

import numpy as np
import pytest

from astropy.stats._fast_sigma_clip import _sigma_clip_fast


def test_basic_returns_finite_bounds():
    data = np.array([[1.0, 2.0, 3.0]], dtype=np.float64)
    mask = np.array([False, False, False], dtype=np.bool_)
    bound_low, bound_high = _sigma_clip_fast(data, mask, False, False, 5, 3.0, 3.0)
    assert np.isfinite(bound_low.item())
    assert np.isfinite(bound_high.item())
    assert bound_low.item() < data.min()
    assert bound_high.item() > data.max()


def test_median_vs_mean_differ():
    data = np.array([[1.0, 2.0, 3.0, 100.0]], dtype=np.float64)
    mask = np.array([False, False, False, False], dtype=np.bool_)
    bound_low_mean, _ = _sigma_clip_fast(data, mask, False, False, 5, 3.0, 3.0)
    bound_low_median, _ = _sigma_clip_fast(data, mask, True, False, 5, 3.0, 3.0)
    assert bound_low_mean.item() != bound_low_median.item()


def test_pre_masked_outlier_excluded():
    data = np.array([[1.0, 2.0, 3.0, 100.0]], dtype=np.float64)
    mask = np.array([False, False, False, True], dtype=np.bool_)
    _, bound_high = _sigma_clip_fast(data, mask, False, False, 5, 3.0, 3.0)
    assert bound_high.item() < data.max()


def test_mad_std_vs_std_differ():
    data = np.array([[1.0, 2.0, 3.0, 100.0]], dtype=np.float64)
    mask = np.array([False, False, False, False], dtype=np.bool_)
    bound_low_std, _ = _sigma_clip_fast(data, mask, False, False, 5, 3.0, 3.0)
    # RuntimeWarning may be emitted as an implementation detail when use_mad_std=True
    # converges to a state where MAD=0. We suppress it here as it is not a guaranteed
    # contract of the function.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        bound_low_mad, _ = _sigma_clip_fast(data, mask, False, True, 5, 3.0, 3.0)
    assert bound_low_std.item() != bound_low_mad.item()


def test_max_iter_zero():
    data = np.array([[1.0, 2.0, 3.0, 100.0]], dtype=np.float64)
    mask = np.array([False, False, False, False], dtype=np.bool_)
    bound_low, bound_high = _sigma_clip_fast(data, mask, False, False, 0, 3.0, 3.0)
    assert np.isfinite(bound_low.item())
    assert np.isfinite(bound_high.item())


def test_single_element():
    data = np.array([[5.0]], dtype=np.float64)
    mask = np.array([False], dtype=np.bool_)
    bound_low, bound_high = _sigma_clip_fast(data, mask, False, False, 5, 3.0, 3.0)
    assert bound_low.item() == 5.0
    assert bound_high.item() == 5.0


def test_all_masked_returns_nan():
    data = np.array([[1.0, 2.0, 3.0]], dtype=np.float64)
    mask = np.array([True, True, True], dtype=np.bool_)
    bound_low, bound_high = _sigma_clip_fast(data, mask, False, False, 5, 3.0, 3.0)
    assert np.isnan(bound_low.item())
    assert np.isnan(bound_high.item())


@pytest.mark.parametrize(
    "sigma_lower, sigma_upper, check",
    [
        pytest.param(
            1.0,
            3.0,
            lambda lo_sym, hi_sym, lo, hi: (
                lo.item() != lo_sym.item() and hi.item() == hi_sym.item()
            ),
            id="lower_changes_upper_unchanged_asymmetric_sigma",
        ),
        pytest.param(
            5.0,
            1.0,
            lambda lo_sym, hi_sym, lo, hi: hi.item() < hi_sym.item(),
            id="upper_tighter_when_sigma_lower_gt_sigma_upper",
        ),
        pytest.param(
            3.0,
            0.0,
            lambda _lo_sym, _hi_sym, lo, hi: lo.item() == hi.item(),
            id="sigma_upper_zero_bounds_converge",
        ),
        pytest.param(
            0.0,
            3.0,
            lambda _lo_sym, _hi_sym, lo, hi: lo.item() == hi.item(),
            id="sigma_lower_zero_bounds_converge",
        ),
    ],
)
def test_asymmetric_sigma(sigma_lower, sigma_upper, check):
    """Test sigma clipping with various asymmetric sigma_lower/sigma_upper combinations."""
    data = np.array([[1.0, 2.0, 3.0, 100.0]], dtype=np.float64)
    mask = np.array([False, False, False, False], dtype=np.bool_)
    bound_low_sym, bound_high_sym = _sigma_clip_fast(
        data, mask, False, False, 5, 3.0, 3.0
    )
    bound_low, bound_high = _sigma_clip_fast(
        data, mask, False, False, 5, sigma_lower, sigma_upper
    )
    assert check(bound_low_sym, bound_high_sym, bound_low, bound_high)
