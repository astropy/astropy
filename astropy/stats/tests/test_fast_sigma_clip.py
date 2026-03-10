# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest

from astropy.stats._fast_sigma_clip import _sigma_clip_fast


def test_basic_returns_finite_bounds():
    data = np.array([[1.0, 2.0, 3.0]], dtype=np.float64)
    mask = np.array([False, False, False], dtype=np.bool_)
    bound_low, bound_high = _sigma_clip_fast(data, mask, False, False, 5, 3.0, 3.0)
    assert np.isfinite(bound_low[0])
    assert np.isfinite(bound_high[0])
    assert bound_low[0] < 1.0
    assert bound_high[0] > 3.0


def test_median_vs_mean_differ():
    data = np.array([[1.0, 2.0, 3.0, 100.0]], dtype=np.float64)
    mask = np.array([False, False, False, False], dtype=np.bool_)
    bound_low_mean, _ = _sigma_clip_fast(data, mask, False, False, 5, 3.0, 3.0)
    bound_low_median, _ = _sigma_clip_fast(data, mask, True, False, 5, 3.0, 3.0)
    assert bound_low_mean[0] != bound_low_median[0]


def test_pre_masked_outlier_excluded():
    data = np.array([[1.0, 2.0, 3.0, 100.0]], dtype=np.float64)
    mask = np.array([False, False, False, True], dtype=np.bool_)
    _, bound_high = _sigma_clip_fast(data, mask, False, False, 5, 3.0, 3.0)
    assert bound_high[0] < 100.0


def test_mad_std_vs_std_differ():
    data = np.array([[1.0, 2.0, 3.0, 100.0]], dtype=np.float64)
    mask = np.array([False, False, False, False], dtype=np.bool_)
    bound_low_std, _ = _sigma_clip_fast(data, mask, False, False, 5, 3.0, 3.0)
    with pytest.warns(RuntimeWarning, match="invalid value"):
        bound_low_mad, _ = _sigma_clip_fast(data, mask, False, True, 5, 3.0, 3.0)
    assert bound_low_std[0] != bound_low_mad[0]


def test_max_iter_zero():
    data = np.array([[1.0, 2.0, 3.0, 100.0]], dtype=np.float64)
    mask = np.array([False, False, False, False], dtype=np.bool_)
    bound_low, bound_high = _sigma_clip_fast(data, mask, False, False, 0, 3.0, 3.0)
    assert np.isfinite(bound_low[0])
    assert np.isfinite(bound_high[0])


def test_single_element():
    data = np.array([[5.0]], dtype=np.float64)
    mask = np.array([False], dtype=np.bool_)
    bound_low, bound_high = _sigma_clip_fast(data, mask, False, False, 5, 3.0, 3.0)
    assert bound_low[0] == 5.0
    assert bound_high[0] == 5.0


def test_all_masked_returns_nan():
    data = np.array([[1.0, 2.0, 3.0]], dtype=np.float64)
    mask = np.array([True, True, True], dtype=np.bool_)
    bound_low, bound_high = _sigma_clip_fast(data, mask, False, False, 5, 3.0, 3.0)
    assert np.isnan(bound_low[0])
    assert np.isnan(bound_high[0])
