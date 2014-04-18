# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test broadcast_arrays replacement on Quantity class.
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np

from .. import broadcast_arrays, PR4622
from ....tests.helper import pytest
from .... import units as u


def test_PR():
    assert PR4622(broadcast_arrays) is True


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
