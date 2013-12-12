# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test `astropy.utils.timer`.

.. note::

    The tests only compare rough estimates as
    performance is machine-dependent.

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ...extern import six

# STDLIB
import time

# THIRD-PARTY
import numpy as np

# LOCAL
from ..timer import RunTimePredictor
from ...modeling.fitting import ModelsError
from ...tests.helper import pytest


def func_to_time(x):
    """This sleeps for y seconds for use with timing tests.

    .. math::

        y = 5 * x - 10

    """
    y = 5.0 * np.asarray(x) - 10
    time.sleep(y)
    return y


def test_timer():
    """Test function timer."""
    p = RunTimePredictor(func_to_time)

    # --- These must run before data points are introduced. ---

    with pytest.raises(AssertionError):
        p.do_fit()

    with pytest.raises(AssertionError):
        p.predict_time(100)

    # --- These must run next to set up data points. ---

    p.time_func([2.02, 2.04, 2.1, 'a', 2.3])
    p.time_func(2.2)  # Test OrderedDict

    assert p._funcname == 'func_to_time'
    assert p._cache_bad == ['a']

    k = list(six.iterkeys(p.results))
    v = list(six.itervalues(p.results))
    np.testing.assert_array_equal(k, [2.02, 2.04, 2.1, 2.3, 2.2])
    np.testing.assert_allclose(v, [0.1, 0.2, 0.5, 1.5, 1.0])

    # --- These should only run once baseline is established. ---

    with pytest.raises(ModelsError):
        a = p.do_fit(model='foo')

    with pytest.raises(ModelsError):
        a = p.do_fit(fitter='foo')

    a = p.do_fit()

    assert p._power == 1

    # Perfect slope is 5, with 10% uncertainty
    assert 4.5 <= a[1] <= 5.5

    # Perfect intercept is -10, with 1-sec uncertainty
    assert -11 <= a[0] <= -9

    # --- These should only run once fitting is completed. ---

    # Perfect answer is 490, with 10% uncertainty
    t = p.predict_time(100)
    assert 441 <= t <= 539

    # Repeated call to access cached run time
    t2 = p.predict_time(100)
    assert t == t2
