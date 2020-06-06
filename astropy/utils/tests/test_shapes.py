# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np

from astropy.utils.shapes import check_broadcast, unbroadcast


def test_check_broadcast():
    assert check_broadcast((10, 1), (3,)) == (10, 3)
    assert check_broadcast((10, 1), (3,), (4, 1, 1, 3)) == (4, 1, 10, 3)
    with pytest.raises(ValueError):
        check_broadcast((10, 2), (3,))

    with pytest.raises(ValueError):
        check_broadcast((10, 1), (3,), (4, 1, 2, 3))


def test_unbroadcast():

    x = np.array([1, 2, 3])
    y = np.broadcast_to(x, (2, 4, 3))
    z = unbroadcast(y)
    assert z.shape == (3,)
    np.testing.assert_equal(z, x)

    x = np.ones((3, 5))
    y = np.broadcast_to(x, (5, 3, 5))
    z = unbroadcast(y)
    assert z.shape == (3, 5)
