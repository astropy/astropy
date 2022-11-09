# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy.stats import (
    calculate_bin_edges,
    freedman_bin_width,
    histogram,
    knuth_bin_width,
    scott_bin_width,
)
from astropy.utils.compat.optional_deps import HAS_SCIPY


def test_scott_bin_width(N=10000, rseed=0):
    rng = np.random.default_rng(rseed)
    X = rng.standard_normal(N)

    delta = scott_bin_width(X)
    assert_allclose(delta, 3.5 * np.std(X) / N ** (1 / 3))

    delta, bins = scott_bin_width(X, return_bins=True)
    assert_allclose(delta, 3.5 * np.std(X) / N ** (1 / 3))

    with pytest.raises(ValueError):
        scott_bin_width(rng.random((2, 10)))


def test_freedman_bin_width(N=10000, rseed=0):
    rng = np.random.default_rng(rseed)
    X = rng.standard_normal(N)

    v25, v75 = np.percentile(X, [25, 75])

    delta = freedman_bin_width(X)
    assert_allclose(delta, 2 * (v75 - v25) / N ** (1 / 3))

    delta, bins = freedman_bin_width(X, return_bins=True)
    assert_allclose(delta, 2 * (v75 - v25) / N ** (1 / 3))

    with pytest.raises(ValueError):
        freedman_bin_width(rng.random((2, 10)))

    # data with too small IQR
    test_x = [1, 2, 3] + [4] * 100 + [5, 6, 7]
    with pytest.raises(ValueError, match=r"Please use another bin method"):
        with pytest.warns(RuntimeWarning, match=r"divide by zero encountered"):
            freedman_bin_width(test_x, return_bins=True)

    # data with small IQR but not too small
    test_x = np.asarray([1, 2, 3] * 100 + [4] + [5, 6, 7], dtype=np.float32)
    test_x *= 1.5e-6
    delta, bins = freedman_bin_width(test_x, return_bins=True)
    assert_allclose(delta, 8.923325554510689e-07)


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
def test_knuth_bin_width(N=10000, rseed=0):
    rng = np.random.default_rng(rseed)
    X = rng.standard_normal(N)

    dx, bins = knuth_bin_width(X, return_bins=True)
    assert_allclose(len(bins), 58)

    dx2 = knuth_bin_width(X)
    assert dx == dx2

    with pytest.raises(ValueError):
        knuth_bin_width(rng.random((2, 10)))


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
def test_knuth_histogram(N=1000, rseed=0):
    rng = np.random.default_rng(rseed)
    x = rng.standard_normal(N)
    counts, bins = histogram(x, "knuth")
    assert counts.sum() == len(x)
    assert len(counts) == len(bins) - 1


_bin_types_to_test = [30, "scott", "freedman", "blocks"]

if HAS_SCIPY:
    _bin_types_to_test += ["knuth"]


@pytest.mark.parametrize("bin_type", _bin_types_to_test + [np.linspace(-5, 5, 31)])
def test_histogram(bin_type, N=1000, rseed=0):
    rng = np.random.default_rng(rseed)
    x = rng.standard_normal(N)
    counts, bins = histogram(x, bin_type)
    assert counts.sum() == len(x)
    assert len(counts) == len(bins) - 1


# Don't include a list of bins as a bin_type here because the effect
# of range is different in that case
@pytest.mark.parametrize("bin_type", _bin_types_to_test)
def test_histogram_range(bin_type, N=1000, rseed=0):
    # Regression test for #8010
    rng = np.random.default_rng(rseed)
    x = rng.standard_normal(N)
    range = (0.1, 0.8)

    bins = calculate_bin_edges(x, bin_type, range=range)
    assert bins.max() == range[1]
    assert bins.min() == range[0]


def test_histogram_range_with_bins_list(N=1000, rseed=0):
    # The expected result when the input bins is a list is
    # the same list on output.
    rng = np.random.default_rng(rseed)
    x = rng.standard_normal(N)
    range = (0.1, 0.8)

    input_bins = np.linspace(-5, 5, 31)
    bins = calculate_bin_edges(x, input_bins, range=range)
    assert all(bins == input_bins)


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
def test_histogram_output_knuth():
    rng = np.random.default_rng(0)
    X = rng.standard_normal(100)

    counts, bins = histogram(X, bins="knuth")
    assert_allclose(counts, [2, 1, 13, 19, 15, 18, 14, 10, 8])

    # fmt: off
    assert_allclose(bins, [-2.32503077, -1.84420596, -1.36338114, -0.88255632, -0.4017315,
                           0.07909331, 0.55991813, 1.04074295, 1.52156777, 2.00239258])
    # fmt: on


def test_histogram_output():
    rng = np.random.default_rng(0)
    X = rng.standard_normal(100)

    counts, bins = histogram(X, bins=10)
    assert_allclose(counts, [2, 0, 12, 14, 14, 17, 16, 8, 9, 8])

    # fmt: off
    assert_allclose(bins, [-2.32503077, -1.89228844, -1.4595461, -1.02680377, -0.59406143,
                           -0.1613191, 0.27142324, 0.70416558, 1.13690791, 1.56965025,
                           2.00239258])
    # fmt: on

    counts, bins = histogram(X, bins="scott")
    assert_allclose(counts, [2, 14, 27, 25, 16, 16])

    # fmt: off
    assert_allclose(bins, [-2.32503077, -1.59953424, -0.87403771, -0.14854117, 0.57695536,
                           1.3024519, 2.02794843])
    # fmt: on

    counts, bins = histogram(X, bins="freedman")
    assert_allclose(counts, [2, 11, 16, 18, 22, 14, 13, 4])

    # fmt: off
    assert_allclose(bins, [-2.32503077, -1.74087192, -1.15671306, -0.5725542,  0.01160465,
                           0.59576351, 1.17992237, 1.76408122, 2.34824008], rtol=2e-7)
    # fmt: on

    counts, bins = histogram(X, bins="blocks")
    assert_allclose(counts, [3, 97])
    assert_allclose(bins, [-2.32503077, -1.37136996, 2.00239258])


def test_histogram_badargs(N=1000, rseed=0):
    rng = np.random.default_rng(rseed)
    x = rng.standard_normal(N)

    # weights is not supported
    for bins in ["scott", "freedman", "blocks"]:
        with pytest.raises(NotImplementedError):
            histogram(x, bins, weights=x)

    # bad bins arg gives ValueError
    with pytest.raises(ValueError):
        histogram(x, bins="bad_argument")
