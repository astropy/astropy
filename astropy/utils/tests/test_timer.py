# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test `astropy.utils.timer`.

.. note::

    The tests only compare fitted results rounded to
    nearest integer. More accurate comparisons might
    fail on some machines.

"""
# STDLIB
import time

# THIRD-PARTY
import numpy as np

# LOCAL
from ..timer import RunTimePredictor


def func_to_time(x):
    """This is sleeps for x seconds for timing tests."""
    time.sleep(x)
    return 'Slept for {0} second(s)'.format(x)


class TestRunTimePredictor(object):
    """Test `astropy.utils.timer.RunTimePredictor`."""
    def setup_class(self):
        self.p = RunTimePredictor(func_to_time)

    def test_expected_errors(self):
        try:
            self.p.do_fit()
        except AssertionError as e:
            assert str(e) == 'Requires 3 points but has 0'

        try:
            self.p.predict_time(100)
        except AssertionError as e:
            assert str(e) == 'No fitted data for prediction'

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
        np.testing.assert_array_equal(np.round(np.array(a)), (1, 0))

    def test_prediction(self):
        t = self.p.predict_time(100)
        assert round(t) == 100

        # Repeated call to access cached run time
        t2 = self.p.predict_time(100)
        assert t == t2
