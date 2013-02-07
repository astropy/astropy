# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test `astropy.utils.timer`.

.. note::

    The tests might fail if function being timed
    deviates from expected run time by more than
    `ACCURACY_DECIMAL` decimals.

"""
# STDLIB
import time

# THIRD-PARTY
import numpy as np

# LOCAL
from ..timer import timeit, SimpleRunTimePredictor
from ...tests.helper import pytest


ACCURACY_DECIMAL = 3  # For np.testing.assert_almost_equal()


def func_to_time(x):
    """This is sleeps for x seconds for timing tests."""
    time.sleep(x)
    return 'Slept for {0} second(s)'.format(x)


@timeit(3)
def timed_func(*args, **kwargs):
    """Run `func_to_time` 3 times and time it."""
    return func_to_time(*args, **kwargs)


def test_timeit(object):
    """Test `astropy.utils.timer.timeit`."""
    t_ave, out_str = timed_func(1)
    assert out_str == 'Slept for {0} second(s)'
    np.testing.assert_almost_equal(t_ave, 1, ACCURACY_DECIMAL)


class TestSimpleRunTimePredictor(object):
    """Test `astropy.utils.timer.SimpleRunTimePredictor`."""
    def setup_class(self):
        p = SimpleRunTimePredictor(func_to_time)

    def test_expected_errors(self):
        try:
            p.do_fit()
        except AssertionError as e:
            assert e.message == 'Requires 3 points but as 0'

        try:
            p.predict_time(100)
        except AssertionError as e:
            assert e.message == 'No fitted data for prediction'

    def test_baseline(self):
        p.time_func([0.1, 0.2, 0.5, -1, 1.5])
        p.time_func(1.0)

        assert p._funcname == 'func_to_time'
        assert p._cache_bad == [-1]
        assert p.results == {0.1: 'Slept for 0.1 second(s)',
                             0.2: 'Slept for 0.2 second(s)',
                             0.5: 'Slept for 0.5 second(s)',
                             1.5: 'Slept for 1.5 second(s)',
                             1.0: 'Slept for 1.0 second(s)'}

    def test_fitting(self):
        a = p.do_fit()
        assert p._power == 1
        np.testing.assert_almost_equal(a, (1, 0), ACCURACY_DECIMAL)

    def test_prediction(self):
        t = p.predict_time(100)
        assert round(t) == 100

        # Repeated call to access cached run time
        t2 = p.predict_time(100)
        assert t == t2
