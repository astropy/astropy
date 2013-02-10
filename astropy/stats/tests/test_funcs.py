# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division

import numpy as np
from numpy.testing import assert_equal, assert_almost_equal
from ...tests.helper import pytest

from .. import funcs
from ...utils.misc import NumpyRNGContext

try:
    from scipy import stats  # used in testing
    from scipy.integrate import quad  # used in testing
    from scipy.special import betainc, betaincinv  # used in funcs
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True


def test_sigma_clip():
    from numpy.random import randn

    #need to seed the numpy RNG to make sure we don't get some amazingly flukey
    #random number that breaks one of the tests

    with NumpyRNGContext(12345):  # Amazing, I've got the same combination on my luggage!

        randvar = randn(10000)

        data, mask = funcs.sigma_clip(randvar, 1, 2)
        maskedarr = funcs.sigma_clip(randvar, 1, 2, maout=True)

        assert sum(mask) > 0
        assert data.size < randvar.size
        assert np.all(mask == ~maskedarr.mask)

        #this is actually a silly thing to do, because it uses the standard
        #deviation as the variance, but it tests to make sure these arguments
        #are actually doing something
        data2, mask2 = funcs.sigma_clip(randvar, 1, 2, varfunc=np.std)
        assert not np.all(data == data2)
        assert not np.all(mask == mask2)

        data3, mask3 = funcs.sigma_clip(randvar, 1, 2, cenfunc=np.mean)
        assert not np.all(data == data3)
        assert not np.all(mask == mask3)

        #now just make sure the iters=None method works at all.
        maskedarr = funcs.sigma_clip(randvar, 3, None, maout=True)


@pytest.mark.skipif('not HAS_SCIPY')
def test_compare_to_scipy_sigmaclip():
    from numpy.random import randn

    #need to seed the numpy RNG to make sure we don't get some amazingly flukey
    #random number that breaks one of the tests

    with NumpyRNGContext(12345):  # Amazing, I've got the same combination on my luggage!

        randvar = randn(10000)

        astropyres = funcs.sigma_clip(randvar, 3, None, np.mean)[0]
        scipyres = stats.sigmaclip(randvar, 3, 3)[0]

        assert_equal(astropyres, scipyres)


@pytest.mark.skipif('not HAS_SCIPY')
def test_binom_conf_interval():

    # Test for corner cases k = 0, k = n, conf = 0., conf = 1.
    n = 5
    k = [0, 4, 5]
    for conf in [0., 0.5, 1.]:
        interval = funcs.binom_conf_interval(k, n, conf=conf)

    # Test that confidence interval contains the probability we request.
    conf = 0.5
    for k in [0, 1, 2, 3, 4, 5]:
        lo, hi = funcs.binom_conf_interval(k, n, conf=conf)

        # The posterior PDF is proportional to binom.pmf. Compare the 
        # integratral of binom.pmf over [0, 1] to the integral over [lo, hi]
        p_tot = quad(lambda e: stats.binom.pmf(k, n, e), 0., 1.)[0]
        p_inside = quad(lambda e: stats.binom.pmf(k, n, e), lo, hi)[0]
        assert_almost_equal(p_inside / p_tot, conf)

        # check correct percentage below the peak
        p_peak = k / n
        if k == 0:
            assert_almost_equal(lo, 0.)
        else:
            p_tot_below = quad(lambda e: stats.binom.pmf(k, n, e),
                               0., p_peak)[0]
            p_inside_below = quad(lambda e: stats.binom.pmf(k, n, e),
                                  lo, p_peak)[0]
            assert_almost_equal(p_inside_below / p_tot_below, conf)

        if k == n:
            assert_almost_equal(hi, 1.)


@pytest.mark.skipif('not HAS_SCIPY')
def test_effhist():

    # Check that it works.
    nbins = 20
    x = 10. * np.random.rand(100)
    success = np.zeros(len(x), dtype=np.bool)
    success[0:50] = True
    bin_ctr, bin_hw, p, perr = funcs.effhist(x, success, bins=nbins)
    assert bin_ctr.shape == (nbins,)
    assert bin_hw.shape == (nbins,)
    assert p.shape == (nbins,)
    assert perr.shape == (2, nbins)
