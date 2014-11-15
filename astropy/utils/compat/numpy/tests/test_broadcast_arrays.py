# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test broadcast_arrays replacement on Quantity class.
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np

from astropy.tests.helper import pytest
import astropy.units as u

from ..lib.stride_tricks import broadcast_arrays, GE1P10


def test_import():
    """Check that what is imported from code is what we are testing."""
    from ... import numpy as anp
    assert anp.broadcast_arrays is broadcast_arrays


def test_test_function():
    """Test the test function

    The possibly patched version of broadcast_arrays should always be OK
    The numpy version may be, in which case we just use it, or it may not,
    it which case we use the patched version.
    """
    assert GE1P10(broadcast_arrays) is True
    if GE1P10():
        assert broadcast_arrays is np.broadcast_arrays
    else:
        assert broadcast_arrays is not np.broadcast_arrays


def test_broadcast_quantity():
    q1 = u.Quantity([1., 2., 3., 4.], u.m)
    q2 = u.Quantity([5., 6.], u.deg).reshape(-1, 1)
    ba1, ba2 = broadcast_arrays(q1, q2)
    assert type(ba1) is np.ndarray
    assert type(ba2) is np.ndarray
    bq1, bq2 = broadcast_arrays(q1, q2, subok=True)
    assert type(bq1) is u.Quantity
    assert type(bq2) is u.Quantity
    assert bq1.shape == bq2.shape == (2, 4)
    with pytest.raises(ValueError):
        broadcast_arrays(q1, q1[:2], subok=True)
