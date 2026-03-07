import numpy as np
import pytest
from scipy.stats import ks_2samp as scipy_ks_2samp

from astropy.stats._stats import ks_2samp


def test_identical_arrays_returns_zero():
    a = np.sort(np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64))
    result = ks_2samp(a, a)
    assert result == pytest.approx(0.0)


def test_completely_different_returns_one():
    a = np.sort(np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64))
    b = np.sort(np.array([10.0, 20.0, 30.0, 40.0], dtype=np.float64))
    result = ks_2samp(a, b)
    assert result == pytest.approx(1.0)


def test_partial_overlap_matches_scipy():
    a = np.sort(np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float64))
    b = np.sort(np.array([3.0, 4.0, 5.0, 6.0, 7.0], dtype=np.float64))
    result = ks_2samp(a, b)
    scipy_result = scipy_ks_2samp(a, b).statistic
    assert abs(result - scipy_result) < 1e-10


def test_different_sizes_matches_scipy():
    a = np.sort(np.array([1.0, 2.0, 3.0], dtype=np.float64))
    b = np.sort(np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], dtype=np.float64))
    result = ks_2samp(a, b)
    scipy_result = scipy_ks_2samp(a, b).statistic
    assert abs(result - scipy_result) < 1e-10


def test_integer_dtype():
    a = np.sort(np.array([1, 2, 3, 4], dtype=np.int32))
    b = np.sort(np.array([5, 6, 7, 8], dtype=np.int32))
    result = ks_2samp(a, b)
    assert result == pytest.approx(1.0)


def test_uint8_dtype():
    a = np.sort(np.array([1, 2, 3], dtype=np.uint8))
    b = np.sort(np.array([4, 5, 6], dtype=np.uint8))
    result = ks_2samp(a, b)
    assert result == pytest.approx(1.0)


def test_single_element_each():
    a = np.sort(np.array([1.0], dtype=np.float64))
    b = np.sort(np.array([2.0], dtype=np.float64))
    result = ks_2samp(a, b)
    assert result == pytest.approx(1.0)


def test_result_between_zero_and_one():
    a = np.sort(np.array([1.0, 2.0, 2.5, 3.0, 4.0], dtype=np.float64))
    b = np.sort(np.array([0.5, 1.5, 3.5, 5.0], dtype=np.float64))
    result = ks_2samp(a, b)
    assert 0.0 <= result <= 1.0
