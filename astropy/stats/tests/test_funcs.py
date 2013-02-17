# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division

import numpy as np
from numpy.testing import assert_equal, assert_almost_equal, assert_allclose
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

    # Corner cases k = 0, k = n, conf = 0., conf = 1.
    n = 5
    k = [0, 4, 5]
    for conf in [0., 0.5, 1.]:
        res = funcs.binom_conf_interval(k, n, conf=conf, interval='wilson')
        assert ((res >= 0.) & (res <= 1.)).all()
        res = funcs.binom_conf_interval(k, n, conf=conf, interval='jeffreys')
        assert ((res >= 0.) & (res <= 1.)).all()

    # Test against table in Brown et al. (2001). (See function docstring.)
    k = [0, 1, 2, 3, 4]
    n = 7
    conf = 0.95
    result = funcs.binom_conf_interval(k, n, conf=conf,
                                       interval='jeffreys')
    table = np.array([[   0., 0.016, 0.065, 0.139, 0.234],
                      [0.292, 0.501, 0.648, 0.766, 0.861]])
    assert_allclose(result, table, atol=1.e-3, rtol=0.)


@pytest.mark.skipif('not HAS_SCIPY')
def test_binned_efficiency():

    # Check that it works.
    nbins = 20
    x = np.linspace(0., 10., 100)  # Guarantee an `x` in every bin.
    success = np.ones(len(x), dtype=np.bool)
    bin_ctr, bin_hw, p, perr = funcs.binned_efficiency(x, success, bins=nbins)

    # Check shape of outputs
    assert bin_ctr.shape == (nbins,)
    assert bin_hw.shape == (nbins,)
    assert p.shape == (nbins,)
    assert perr.shape == (2, nbins)

    # Check that p is 1 in all bins. (success = True for all `x`)
    assert (p == 1.).all()

    # Check that p is 0 in all bins if success = False for all `x`
    success[:] = False
    bin_ctr, bin_hw, p, perr = funcs.binned_efficiency(x, success, bins=nbins)
    assert (p == 0.).all()
