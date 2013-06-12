# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test `astropy.utils.timer`.

.. note::

    The tests only compare rough estimates as
    performance is machine-dependent.

"""
# STDLIB
import time

# LOCAL
from ..timer import RunTimePredictor


def func_to_time(x):
    """This is sleeps for x seconds for timing tests."""
    time.sleep(x)
    return 'Slept for {0} second(s)'.format(x)


def test_timer():
    p = RunTimePredictor(func_to_time)
    try:
        p.do_fit()
    except AssertionError as e:
        assert str(e) == 'Requires 3 points but has 0'

    try:
        p.predict_time(100)
    except AssertionError as e:
        assert str(e) == 'No fitted data for prediction'

    p.time_func([0.1, 0.2, 0.5, 'a', 1.5])
    p.time_func(1.0)

    assert p._funcname == 'func_to_time'
    assert p._cache_bad == ['a']
    assert p.results == {0.1: 'Slept for 0.1 second(s)',
                         0.2: 'Slept for 0.2 second(s)',
                         0.5: 'Slept for 0.5 second(s)',
                         1.5: 'Slept for 1.5 second(s)',
                         1.0: 'Slept for 1.0 second(s)'}

    a = p.do_fit()
    assert p._power == 1

    # Perfect slope is 1, with 10% uncertainty
    assert 0.9 <= a[0] <= 1.1

    # Perfect intercept is 0, with 1-sec uncertainty
    assert -1 <= a[1] <= 1

    # Perfect answer is 100, with 10% uncertainty
    t = p.predict_time(100)
    assert 90 <= t <= 110

    # Repeated call to access cached run time
    t2 = p.predict_time(100)
    assert t == t2
