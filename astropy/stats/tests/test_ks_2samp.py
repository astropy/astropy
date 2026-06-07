# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy.stats._stats import ks_2samp


def test_identical_array():
    a = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    assert_allclose(ks_2samp(a, a), 0.0, atol=1e-15)


def test_completely_different():
    a = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    b = np.array([10.0, 20.0, 30.0, 40.0], dtype=np.float64)
    assert_allclose(ks_2samp(a, b), 1.0, atol=1e-15)


def test_partial_overlap():
    # KS statistic computed analytically: max CDF difference = 2/5
    a = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float64)
    b = np.array([3.0, 4.0, 5.0, 6.0, 7.0], dtype=np.float64)
    assert_allclose(ks_2samp(a, b), 0.4, atol=1e-15)


def test_different_sizes():
    # KS statistic computed analytically: max CDF difference = 1/2
    a = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    b = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], dtype=np.float64)
    assert_allclose(ks_2samp(a, b), 0.5, atol=1e-15)


@pytest.mark.parametrize("dtype", [np.int32, np.uint8, np.float64])
def test_dtype(dtype):
    a = np.array([1, 2, 3, 4], dtype=dtype)
    b = np.array([5, 6, 7, 8], dtype=dtype)
    assert_allclose(ks_2samp(a, b), 1.0, atol=1e-15)


def test_mixed_overlap():
    # KS statistic computed analytically: max CDF difference = 3/5
    a = np.array([1.0, 2.0, 2.5, 3.0, 4.0], dtype=np.float64)
    b = np.array([0.5, 1.5, 3.5, 5.0], dtype=np.float64)
    assert_allclose(ks_2samp(a, b), 0.6, atol=1e-15)
