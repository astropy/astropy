# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from numpy.testing import assert_allclose, assert_, assert_raises

from ...tests.helper import pytest
from .. import histogram, scotts_bin_width, freedman_bin_width, knuth_bin_width

try:
    import scipy  # pylint: disable=W0611
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True


def test_scotts_bin_width(N=10000, rseed=0):
    np.random.seed(rseed)
    X = np.random.normal(size=N)

    delta = scotts_bin_width(X)
    assert_allclose(delta,  3.5 * np.std(X) / N ** (1 / 3))

    delta, bins = scotts_bin_width(X, return_bins=True)
    assert_allclose(delta,  3.5 * np.std(X) / N ** (1 / 3))

    assert_raises(ValueError, scotts_bin_width, X.reshape(2, -1))


def test_freedman_bin_width(N=10000, rseed=0):
    np.random.seed(rseed)
    X = np.random.normal(size=N)
    v25, v75 = np.percentile(X, [25, 75])

    delta = freedman_bin_width(X)
    assert_allclose(delta, 2 * (v75 - v25) / N ** (1 / 3))

    delta, bins = freedman_bin_width(X, return_bins=True)
    assert_allclose(delta, 2 * (v75 - v25) / N ** (1 / 3))

    assert_raises(ValueError, freedman_bin_width, X.reshape(2, -1))


@pytest.mark.skipif('not HAS_SCIPY')
def test_knuth_bin_width(N=10000, rseed=0):
    np.random.seed(0)
    X = np.random.normal(size=N)

    dx, bins = knuth_bin_width(X, return_bins=True)
    assert_allclose(len(bins), 59)

    dx2 = knuth_bin_width(X)
    assert dx == dx2

    assert_raises(ValueError, knuth_bin_width, X.reshape(2, -1))


@pytest.mark.skipif('not HAS_SCIPY')
def test_knuth_histogram(N=1000, rseed=0):
    np.random.seed(rseed)
    x = np.random.normal(0, 1, N)
    counts, bins = histogram(x, 'knuth')
    assert_(counts.sum() == len(x))
    assert_(len(counts) == len(bins) - 1)


def test_histogram(N=1000, rseed=0):
    np.random.seed(0)
    x = np.random.normal(0, 1, N)

    for bins in [30, np.linspace(-5, 5, 31),
                 'scotts', 'freedman', 'blocks', 'adaptive']:
        counts, bins = histogram(x, bins)
        assert_(counts.sum() == len(x))
        assert_(len(counts) == len(bins) - 1)


def test_histogram_range(N=1000, rseed=0):
    np.random.seed(rseed)
    x = np.random.normal(0, 1, N)
    rng = (0.1, 0.8)
    
    for bins in ['scotts', 'freedman', 'blocks', 'adaptive']:
        print(bins)
        counts, bins = histogram(x, bins, range=rng)
        

def test_histogram_badargs(N=1000, rseed=0):
    np.random.seed(rseed)
    x = np.random.normal(0, 1, N)

    # weights is not supported
    for bins in ['scotts', 'freedman', 'blocks', 'adaptive']:
        assert_raises(NotImplementedError, histogram, x, bins, weights=x)

    # bad bins arg gives ValueError
    assert_raises(ValueError, histogram, x, 'blahblah')

    
