# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from ..masked import masked_arrays_equal


def test_masked_arrays_equal():

    assert masked_arrays_equal(np.ma.array([1, 2, 3]),
                               np.ma.array([1, 3, 3])).tolist() == [True, False, True]

    assert masked_arrays_equal(np.ma.array([1, 2, 3], mask=[1, 0, 0]),
                               np.ma.array([1, 3, 3])).tolist() == [False, False, True]

    assert masked_arrays_equal(np.ma.array([1, 2, 3], mask=[0, 1, 0]),
                               np.ma.array([1, 3, 3], mask=[0, 1, 0])).tolist() == [True, True, True]

    assert masked_arrays_equal(np.ma.array([1, 2, 3], mask=[0, 1, 0], dtype=[('a', '?')]),
                               np.ma.array([1, 3, 3], mask=[0, 1, 0], dtype=[('a', '?')])).tolist() == [True, True, True]
