# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test matmul replacement.
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np

import pytest

from ..core.multiarray import matmul, GE1P10


def test_import():
    """Check that what is imported from code is what we are testing."""
    from ... import numpy as anp
    assert anp.matmul is matmul


def test_test_function():
    """Test the test function

    The possibly patched version of broadcast_arrays should always be OK
    The numpy version may be, in which case we just use it, or it may not,
    it which case we use the patched version.
    """
    from ... import numpy as anp
    assert GE1P10(module=anp) is True
    if GE1P10(module=np):
        assert matmul is np.matmul
    else:
        assert not hasattr(np, 'matmul')


def test_matmul():
    a = np.arange(18).reshape(2, 3, 3)
    # with another matrix
    b1 = np.identity(3)
    assert np.all(matmul(a, b1) == a)
    b2 = 1. - b1
    assert np.all(matmul(a[0], b2) == np.array([[3, 2, 1],
                                                [9, 8, 7],
                                                [15, 14, 13]]))
    b3 = np.ones((4, 1, 3, 2))
    out = np.zeros((4, 2, 3, 2))
    res = matmul(a, b3, out=out)
    assert res is out
    with pytest.raises(ValueError):  # wrong shape
        matmul(b3, a)
    out2 = np.zeros((4, 1, 3, 2))
    with pytest.raises(ValueError):
        matmul(a, b3, out=out2)

    # with a vector
    b4 = np.ones((3,))
    assert np.all(matmul(a, b4) == a.sum(-1))
    out = np.zeros((a.shape[0], a.shape[2]))
    res = matmul(b4, a, out=out)
    assert res is out
    assert np.all(out == a.sum(-2))

    with pytest.raises(ValueError):
        matmul(a, 1.)
    with pytest.raises(ValueError):
        matmul(1., a)
    with pytest.raises(ValueError):
        matmul(a, 1.)
