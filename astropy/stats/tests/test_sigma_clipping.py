# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from numpy.random import randn
from numpy.testing import assert_equal

try:
    from scipy import stats  # used in testing
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True

from ...tests.helper import pytest

from ..sigma_clipping import sigma_clip, sigma_clipped_stats
from ...utils.misc import NumpyRNGContext


def test_sigma_clip():
    #need to seed the numpy RNG to make sure we don't get some amazingly flukey
    #random number that breaks one of the tests

    with NumpyRNGContext(12345):
        # Amazing, I've got the same combination on my luggage!
        randvar = randn(10000)

        filtered_data = sigma_clip(randvar, sigma=1, iters=2)

        assert sum(filtered_data.mask) > 0
        assert sum(~filtered_data.mask) < randvar.size

        # this is actually a silly thing to do, because it uses the
        # standard deviation as the variance, but it tests to make sure
        # these arguments are actually doing something
        filtered_data2 = sigma_clip(randvar, sigma=1, iters=2, stdfunc=np.var)
        assert not np.all(filtered_data.mask == filtered_data2.mask)

        filtered_data3 = sigma_clip(randvar, sigma=1, iters=2,
                                    cenfunc=np.mean)
        assert not np.all(filtered_data.mask == filtered_data3.mask)

        # make sure the iters=None method works at all.
        filtered_data = sigma_clip(randvar, sigma=3, iters=None)

        # test copying
        assert filtered_data.data[0] == randvar[0]
        filtered_data.data[0] += 1.
        assert filtered_data.data[0] != randvar[0]

        filtered_data = sigma_clip(randvar, sigma=3, iters=None, copy=False)
        assert filtered_data.data[0] == randvar[0]
        filtered_data.data[0] += 1.
        assert filtered_data.data[0] == randvar[0]

        # test axis
        data = np.arange(5) + np.random.normal(0., 0.05, (5, 5)) + \
            np.diag(np.ones(5))
        filtered_data = sigma_clip(data, axis=0, sigma=2.3)
        assert filtered_data.count() == 20
        filtered_data = sigma_clip(data, axis=1, sigma=2.3)
        assert filtered_data.count() == 25


@pytest.mark.skipif('not HAS_SCIPY')
def test_compare_to_scipy_sigmaclip():
    #need to seed the numpy RNG to make sure we don't get some amazingly flukey
    #random number that breaks one of the tests

    with NumpyRNGContext(12345):

        randvar = randn(10000)

        astropyres = sigma_clip(randvar, sigma=3, iters=None, cenfunc=np.mean)
        scipyres = stats.sigmaclip(randvar, 3, 3)[0]

        assert astropyres.count() == len(scipyres)
        assert_equal(astropyres[~astropyres.mask].data, scipyres)


def test_sigma_clipped_stats():
    """Test list data with input mask or mask_value (#3268)."""
    # test list data with mask
    data = [0, 1]
    mask = np.array([True, False])
    result = sigma_clipped_stats(data, mask=mask)
    assert result[0] == 1.
    assert result[1] == 1.
    assert result[2] == 0.

    # test list data with mask_value
    result2 = sigma_clipped_stats(data, mask_value=0.)
    assert result2[0] == 1.
    assert result2[1] == 1.
    assert result2[2] == 0.
