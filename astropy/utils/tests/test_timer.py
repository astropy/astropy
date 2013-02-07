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


ACCURACY_DECIMAL = 3  # For np.testing.assert_almost_equal()


def func_to_time(x):
    """This is sleeps for x seconds for timing tests."""
    time.sleep(x)
    return 'Slept for {0} second(s)'.format(x)


class TestSimpleRunTimePredictor(object):
    """Test `astropy.utils.timer.SimpleRunTimePredictor`."""
    def setup_class(self):
        self.p = SimpleRunTimePredictor(func_to_time)

    def test_expected_errors(self):
        try:
            self.p.do_fit()
        except AssertionError as e:
            assert e.message == 'Requires 3 points but has 0'

        try:
            self.p.predict_time(100)
        except AssertionError as e:
            assert e.message == 'No fitted data for prediction'

    def test_baseline(self):
        self.p.time_func([0.1, 0.2, 0.5, -1, 1.5])
        self.p.time_func(1.0)

        assert self.p._funcname == 'func_to_time'
        assert self.p._cache_bad == [-1]
        assert self.p.results == {0.1: 'Slept for 0.1 second(s)',
                                  0.2: 'Slept for 0.2 second(s)',
                                  0.5: 'Slept for 0.5 second(s)',
                                  1.5: 'Slept for 1.5 second(s)',
                                  1.0: 'Slept for 1.0 second(s)'}

    def test_fitting(self):
        a = self.p.do_fit()
        assert self.p._power == 1
        np.testing.assert_almost_equal(a, (1, 0), ACCURACY_DECIMAL)

    def test_prediction(self):
        t = self.p.predict_time(100)
        assert round(t) == 100

        # Repeated call to access cached run time
        t2 = self.p.predict_time(100)
        assert t == t2
