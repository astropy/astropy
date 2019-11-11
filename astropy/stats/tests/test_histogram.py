# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings

import pytest
import numpy as np
from numpy.testing import assert_allclose

from astropy.stats import (histogram, calculate_bin_edges, scott_bin_width,
                           freedman_bin_width, knuth_bin_width)
from astropy.utils.exceptions import AstropyUserWarning

try:
    import scipy  # pylint: disable=W0611 # noqa
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True


def test_scott_bin_width(N=10000, rseed=0):
    rng = np.random.RandomState(rseed)
    X = rng.randn(N)

    delta = scott_bin_width(X)
    assert_allclose(delta, 3.5 * np.std(X) / N ** (1 / 3))

    delta, bins = scott_bin_width(X, return_bins=True)
    assert_allclose(delta, 3.5 * np.std(X) / N ** (1 / 3))

    with pytest.raises(ValueError):
        scott_bin_width(rng.rand(2, 10))


def test_freedman_bin_width(N=10000, rseed=0):
    rng = np.random.RandomState(rseed)
    X = rng.randn(N)

    v25, v75 = np.percentile(X, [25, 75])

    delta = freedman_bin_width(X)
    assert_allclose(delta, 2 * (v75 - v25) / N ** (1 / 3))

    delta, bins = freedman_bin_width(X, return_bins=True)
    assert_allclose(delta, 2 * (v75 - v25) / N ** (1 / 3))

    with pytest.raises(ValueError):
        freedman_bin_width(rng.rand(2, 10))

    # data with too small IQR
    test_x = [1, 2, 3] + [4] * 100 + [5, 6, 7]
    with pytest.raises(ValueError) as e:
        with pytest.warns(RuntimeWarning, match=r'divide by zero encountered'):
            freedman_bin_width(test_x, return_bins=True)
        assert 'Please use another bin method' in str(e.value)

    # data with small IQR but not too small
    test_x = np.asarray([1, 2, 3] * 100 + [4] + [5, 6, 7], dtype=np.float32)
    test_x *= 1.5e-6
    delta, bins = freedman_bin_width(test_x, return_bins=True)
    assert_allclose(delta, 8.923325554510689e-07)


@pytest.mark.skipif('not HAS_SCIPY')
def test_knuth_bin_width(N=10000, rseed=0):
    rng = np.random.RandomState(rseed)
    X = rng.randn(N)

    dx, bins = knuth_bin_width(X, return_bins=True)
    assert_allclose(len(bins), 59)

    dx2 = knuth_bin_width(X)
    assert dx == dx2

    with pytest.raises(ValueError):
        knuth_bin_width(rng.rand(2, 10))


@pytest.mark.skipif('not HAS_SCIPY')
def test_knuth_histogram(N=1000, rseed=0):
    rng = np.random.RandomState(rseed)
    x = rng.randn(N)
    counts, bins = histogram(x, 'knuth')
    assert (counts.sum() == len(x))
    assert (len(counts) == len(bins) - 1)


_bin_types_to_test = [30, 'scott', 'freedman', 'blocks']

if HAS_SCIPY:
    _bin_types_to_test += ['knuth']


@pytest.mark.parametrize('bin_type',
                         _bin_types_to_test + [np.linspace(-5, 5, 31)])
def test_histogram(bin_type, N=1000, rseed=0):
    rng = np.random.RandomState(rseed)
    x = rng.randn(N)
    counts, bins = histogram(x, bin_type)
    assert (counts.sum() == len(x))
    assert (len(counts) == len(bins) - 1)


# Don't include a list of bins as a bin_type here because the effect
# of range is different in that case
@pytest.mark.parametrize('bin_type', _bin_types_to_test)
def test_histogram_range(bin_type, N=1000, rseed=0):
    # Regression test for #8010
    rng = np.random.RandomState(rseed)
    x = rng.randn(N)
    range = (0.1, 0.8)

    bins = calculate_bin_edges(x, bin_type, range=range)
    assert bins.max() == range[1]
    assert bins.min() == range[0]


def test_histogram_range_with_bins_list(N=1000, rseed=0):
    # The expected result when the input bins is a list is
    # the same list on output.
    rng = np.random.RandomState(rseed)
    x = rng.randn(N)
    range = (0.1, 0.8)

    input_bins = np.linspace(-5, 5, 31)
    bins = calculate_bin_edges(x, input_bins, range=range)
    assert all(bins == input_bins)


@pytest.mark.skipif('not HAS_SCIPY')
def test_histogram_output_knuth():
    rng = np.random.RandomState(0)
    X = rng.randn(100)

    counts, bins = histogram(X, bins='knuth')
    assert_allclose(counts, [1, 6, 9, 14, 21, 22, 12, 8, 7])
    assert_allclose(bins, [-2.55298982, -2.01712932, -1.48126883, -0.94540834,
                           -0.40954784, 0.12631265, 0.66217314, 1.19803364,
                           1.73389413, 2.26975462])


def test_histogram_output():
    rng = np.random.RandomState(0)
    X = rng.randn(100)

    counts, bins = histogram(X, bins=10)
    assert_allclose(counts, [1, 5, 7, 13, 17, 18, 16, 11, 7, 5])
    assert_allclose(bins, [-2.55298982, -2.07071537, -1.58844093, -1.10616648,
                           -0.62389204, -0.1416176, 0.34065685, 0.82293129,
                           1.30520574, 1.78748018, 2.26975462])

    counts, bins = histogram(X, bins='scott')
    assert_allclose(counts, [2, 13, 23, 34, 16, 10, 2])
    assert_allclose(bins, [-2.55298982, -1.79299405, -1.03299829, -0.27300252,
                           0.48699324, 1.24698901, 2.00698477, 2.76698054])

    counts, bins = histogram(X, bins='freedman')
    assert_allclose(counts, [2, 7, 13, 20, 26, 14, 11, 5, 2])
    assert_allclose(bins, [-2.55298982, -1.95796338, -1.36293694, -0.7679105,
                           -0.17288406, 0.42214237, 1.01716881, 1.61219525,
                           2.20722169, 2.80224813])

    counts, bins = histogram(X, bins='blocks')
    assert_allclose(counts, [10, 61, 29])
    assert_allclose(bins, [-2.55298982, -1.24381059, 0.46422235, 2.26975462])


def test_histogram_badargs(N=1000, rseed=0):
    rng = np.random.RandomState(rseed)
    x = rng.randn(N)

    # weights is not supported
    for bins in ['scott', 'freedman', 'blocks']:
        with pytest.raises(NotImplementedError):
            histogram(x, bins, weights=x)

    # bad bins arg gives ValueError
    with pytest.raises(ValueError):
        histogram(x, bins='bad_argument')
