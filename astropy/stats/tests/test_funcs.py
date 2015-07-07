# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from numpy.random import randn, normal
from numpy.testing import assert_equal
from numpy.testing.utils import assert_allclose

try:
    import scipy
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True

from ...tests.helper import pytest

from .. import funcs
from ...utils.misc import NumpyRNGContext


def test_median_absolute_deviation():
    #need to seed the numpy RNG to make sure we don't get some amazingly flukey
    #random number that breaks one of the tests

    with NumpyRNGContext(12345):

        #test that it runs
        randvar = randn(10000)
        mad = funcs.median_absolute_deviation(randvar)

        #test whether an array is returned if an axis is used
        randvar = randvar.reshape((10, 1000))
        mad = funcs.median_absolute_deviation(randvar, axis=1)
        assert len(mad) == 10
        assert mad.size < randvar.size
        mad = funcs.median_absolute_deviation(randvar, axis=0)
        assert len(mad) == 1000
        assert mad.size < randvar.size
        # Test some actual values in a 3 dimensional array
        x = np.arange(3*4*5)
        a = np.array([sum(x[:i+1]) for i in range(len(x))]).reshape(3, 4, 5)
        mad = funcs.median_absolute_deviation(a)
        assert mad == 389.5
        mad = funcs.median_absolute_deviation(a, axis=0)
        assert_allclose(mad, [[ 210.,  230.,  250.,  270.,  290.],
                              [ 310.,  330.,  350.,  370.,  390.],
                              [ 410.,  430.,  450.,  470.,  490.],
                              [ 510.,  530.,  550.,  570.,  590.]])
        mad = funcs.median_absolute_deviation(a, axis=1)
        assert_allclose(mad, [[ 27.5,   32.5,   37.5,   42.5,   47.5],
                              [ 127.5,  132.5,  137.5,  142.5,  147.5],
                              [ 227.5,  232.5,  237.5,  242.5,  247.5]])
        mad = funcs.median_absolute_deviation(a, axis=2)
        assert_allclose(mad, [[  3.,   8.,  13.,  18.],
                              [ 23.,  28.,  33.,  38.],
                              [ 43.,  48.,  53.,  58.]])


def test_biweight_location():
    #need to seed the numpy RNG to make sure we don't get some amazingly flukey
    #random number that breaks one of the tests

    with NumpyRNGContext(12345):

        #test that it runs
        randvar = randn(10000)
        cbl = funcs.biweight_location(randvar)

        assert abs(cbl-0) < 1e-2


def test_biweight_location_small():

    cbl = funcs.biweight_location([1, 3, 5, 500, 2])
    assert abs(cbl-2.745) < 1e-3


def test_biweight_midvariance():
    #need to seed the numpy RNG to make sure we don't get some amazingly flukey
    #random number that breaks one of the tests

    with NumpyRNGContext(12345):

        #test that it runs
        randvar = randn(10000)
        scl = funcs.biweight_midvariance(randvar)

        assert abs(scl-1) < 1e-2


def test_biweight_midvariance_small():
    scl = funcs.biweight_midvariance([1, 3, 5, 500, 2])
    assert abs(scl-1.529) < 1e-3


@pytest.mark.skipif('not HAS_SCIPY')
def test_binom_conf_interval():

    # Test Wilson and Jeffreys interval for corner cases:
    # Corner cases: k = 0, k = n, conf = 0., conf = 1.
    n = 5
    k = [0, 4, 5]
    for conf in [0., 0.5, 1.]:
        res = funcs.binom_conf_interval(k, n, conf=conf, interval='wilson')
        assert ((res >= 0.) & (res <= 1.)).all()
        res = funcs.binom_conf_interval(k, n, conf=conf, interval='jeffreys')
        assert ((res >= 0.) & (res <= 1.)).all()

    # Test Jeffreys interval accuracy against table in Brown et al. (2001).
    # (See `binom_conf_interval` docstring for reference.)
    k = [0, 1, 2, 3, 4]
    n = 7
    conf = 0.95
    result = funcs.binom_conf_interval(k, n, conf=conf, interval='jeffreys')
    table = np.array([[0.000, 0.016, 0.065, 0.139, 0.234],
                      [0.292, 0.501, 0.648, 0.766, 0.861]])
    assert_allclose(result, table, atol=1.e-3, rtol=0.)

    # Test scalar version
    result = np.array([funcs.binom_conf_interval(kval, n, conf=conf,
                                                 interval='jeffreys')
                       for kval in k]).transpose()
    assert_allclose(result, table, atol=1.e-3, rtol=0.)

    # Test flat
    result = funcs.binom_conf_interval(k, n, conf=conf, interval='flat')
    table = np.array([[0., 0.03185, 0.08523, 0.15701, 0.24486],
                      [0.36941, 0.52650, 0.65085, 0.75513, 0.84298]])
    assert_allclose(result, table, atol=1.e-3, rtol=0.)

    # Test scalar version
    result = np.array([funcs.binom_conf_interval(kval, n, conf=conf,
                                                 interval='flat')
                       for kval in k]).transpose()
    assert_allclose(result, table, atol=1.e-3, rtol=0.)

    # Test Wald interval
    result = funcs.binom_conf_interval(0, 5, interval='wald')
    assert_allclose(result, 0.)  # conf interval is [0, 0] when k = 0
    result = funcs.binom_conf_interval(5, 5, interval='wald')
    assert_allclose(result, 1.)  # conf interval is [1, 1] when k = n
    result = funcs.binom_conf_interval(500, 1000, conf=0.68269,
                                       interval='wald')
    assert_allclose(result[0], 0.5 - 0.5 / np.sqrt(1000.))
    assert_allclose(result[1], 0.5 + 0.5 / np.sqrt(1000.))

    # Test shapes
    k = 3
    n = 7
    for interval in ['wald', 'wilson', 'jeffreys', 'flat']:
        result = funcs.binom_conf_interval(k, n, interval=interval)
        assert result.shape == (2,)

    k = np.array(k)
    for interval in ['wald', 'wilson', 'jeffreys', 'flat']:
        result = funcs.binom_conf_interval(k, n, interval=interval)
        assert result.shape == (2,)

    n = np.array(n)
    for interval in ['wald', 'wilson', 'jeffreys', 'flat']:
        result = funcs.binom_conf_interval(k, n, interval=interval)
        assert result.shape == (2,)

    k = np.array([1, 3, 5])
    for interval in ['wald', 'wilson', 'jeffreys', 'flat']:
        result = funcs.binom_conf_interval(k, n, interval=interval)
        assert result.shape == (2, 3)

    n = np.array([5, 5, 5])
    for interval in ['wald', 'wilson', 'jeffreys', 'flat']:
        result = funcs.binom_conf_interval(k, n, interval=interval)
        assert result.shape == (2, 3)

@pytest.mark.skipif('not HAS_SCIPY')
def test_binned_binom_proportion():

    # Check that it works.
    nbins = 20
    x = np.linspace(0., 10., 100)  # Guarantee an `x` in every bin.
    success = np.ones(len(x), dtype=np.bool)
    bin_ctr, bin_hw, p, perr = funcs.binned_binom_proportion(x, success,
                                                             bins=nbins)

    # Check shape of outputs
    assert bin_ctr.shape == (nbins,)
    assert bin_hw.shape == (nbins,)
    assert p.shape == (nbins,)
    assert perr.shape == (2, nbins)

    # Check that p is 1 in all bins, since success = True for all `x`.
    assert (p == 1.).all()

    # Check that p is 0 in all bins if success = False for all `x`.
    success[:] = False
    bin_ctr, bin_hw, p, perr = funcs.binned_binom_proportion(x, success,
                                                             bins=nbins)
    assert (p == 0.).all()


def test_signal_to_noise_oir_ccd():

    result = funcs.signal_to_noise_oir_ccd(1, 25, 0, 0, 0, 1)
    assert 5.0 == result
    #check to make sure gain works
    result = funcs.signal_to_noise_oir_ccd(1, 5, 0, 0, 0, 1, 5)
    assert 5.0 == result

    #now add in sky, dark current, and read noise
    #make sure the snr goes down
    result = funcs.signal_to_noise_oir_ccd(1, 25, 1, 0, 0, 1)
    assert result < 5.0
    result = funcs.signal_to_noise_oir_ccd(1, 25, 0, 1, 0, 1)
    assert result < 5.0
    result = funcs.signal_to_noise_oir_ccd(1, 25, 0, 0, 1, 1)
    assert result < 5.0

    #make sure snr increases with time
    result = funcs.signal_to_noise_oir_ccd(2, 25, 0, 0, 0, 1)
    assert result > 5.0


def test_bootstrap():
    bootarr = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 0])
    #test general bootstrapping
    answer = np.array([[7, 4, 8, 5, 7, 0, 3, 7, 8, 5],
                       [4, 8, 8, 3, 6, 5, 2, 8, 6, 2]])
    with NumpyRNGContext(42):
        assert_equal(answer, funcs.bootstrap(bootarr, 2))

    #test with a bootfunction
    with NumpyRNGContext(42):
        bootresult = np.mean(funcs.bootstrap(bootarr, 10000, bootfunc=np.mean))
        assert_allclose(np.mean(bootarr), bootresult, atol=0.01)


def test_mad_std():
    with NumpyRNGContext(12345):
        data = normal(5, 2, size=(100, 100))
        assert_allclose(funcs.mad_std(data), 2.0, rtol=0.05)


def test_gaussian_fwhm_to_sigma():
    fwhm = (2.0 * np.sqrt(2.0 * np.log(2.0)))
    assert_allclose(funcs.gaussian_fwhm_to_sigma * fwhm, 1.0, rtol=1.0e-6)


def test_gaussian_sigma_to_fwhm():
    sigma = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    assert_allclose(funcs.gaussian_sigma_to_fwhm * sigma, 1.0, rtol=1.0e-6)


def test_gaussian_sigma_to_fwhm_to_sigma():
    assert_allclose(funcs.gaussian_fwhm_to_sigma *
                    funcs.gaussian_sigma_to_fwhm, 1.0)

def test_akaike_information_criterion():
    #testing with arrays
    aicc = funcs.akaike_information_criterion(np.array([-10,-5]), np.array([2,3]), 20, True)
    assert_allclose(aicc, np.array([ 24.70588235, 17.5]))
    
    #testing with single float
    aicc = funcs.akaike_information_criterion(-10,2, 20, True)
    assert_allclose(aicc, 24.70588235)
    
    #test regular aic
    aicc = funcs.akaike_information_criterion(np.array([-10,-5]), np.array([2,3]), 20, False)
    assert_allclose(aicc, np.array([24, 16]))
    aicc = funcs.akaike_information_criterion(-10, 2, 20, False)
    assert_allclose(aicc, 24)
    
    #make sure AssertionError in AICc if num_samples - num_params - 1 <= 0
    with pytest.raises(AssertionError):
        aicc = funcs.akaike_information_criterion(-10, 30, 31, True)
    with pytest.raises(AssertionError):
        aicc = funcs.akaike_information_criterion(-10, 40, 31, True)
    with pytest.raises(AssertionError):
        aicc = funcs.akaike_information_criterion(np.array([-10,-5]), np.array([2,30]), 20, True)

def test_bayesian_information_criterion():
    #testing with arrays
    bic = funcs.bayesian_information_criterion(np.array([-10,-5]), np.array([2,3]), 20)
    assert_allclose(bic, np.array([ 25.99146455, 18.98719682]))
    
    #testing with single float
    bic = funcs.bayesian_information_criterion(-10,2, 20)
    assert_allclose(bic, 25.99146455)
