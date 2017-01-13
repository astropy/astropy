# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test broadcast_arrays replacement on Quantity class.
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import pytest
import numpy as np

import astropy.units as u

from ..lib.stride_tricks import broadcast_arrays, broadcast_to, GE1P10


def test_import():
    """Check that what is imported from code is what we are testing."""
    from ... import numpy as anp
    assert anp.broadcast_arrays is broadcast_arrays
    assert anp.broadcast_to is broadcast_to


def test_test_function():
    """Test the test function

    The possibly patched version of broadcast_arrays should always be OK
    The numpy version may be, in which case we just use it, or it may not,
    it which case we use the patched version.
    """
    from ... import numpy as anp
    assert GE1P10(module=anp) is True
    if GE1P10(module=np):
        assert broadcast_arrays is np.broadcast_arrays
        assert broadcast_to is np.broadcast_to
    else:
        assert broadcast_arrays is not np.broadcast_arrays
        assert not hasattr(np, 'broadcast_to')


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


def test_broadcast_to():
    q1 = u.Quantity([1., 2., 3., 4.], u.m)
    bq1 = broadcast_to(q1, (2, 4), subok=True)
    assert type(bq1) is type(q1)
    assert bq1.shape == (2, 4)
    assert bq1.strides[0] == 0
