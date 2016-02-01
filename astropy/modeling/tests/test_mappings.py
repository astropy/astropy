# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
from ...tests.helper import pytest
from ..models import Shift, Rotation2D, Identity, Mapping
import numpy as np
from numpy.testing.utils import (assert_allclose, assert_array_equal)


def test_swap_axes():
    x = np.zeros((2, 3))
    y = np.ones((2, 3))

    mapping = Mapping((1, 0))
    assert(mapping(1, 2) == (2.0, 1.0))
    assert(mapping.inverse(2, 1) == (1, 2))
    assert_array_equal(mapping(x, y), (y, x))
    assert_array_equal(mapping.inverse(y, x), (x, y))


def test_duplicate_axes():
    mapping = Mapping((0, 1, 0, 1))
    assert(mapping(1, 2) == (1.0, 2., 1., 2))
    assert(mapping.inverse(1, 2, 1, 2) == (1, 2))
    assert(mapping.inverse.n_inputs == 4)
    assert(mapping.inverse.n_outputs == 2)


def test_drop_axes_1():
    mapping = Mapping((0,), n_inputs=2)
    assert(mapping(1, 2) == (1.))

def test_drop_axes_2():
    mapping = Mapping((1, ))
    assert(mapping(1, 2) == (2.))
    with pytest.raises(NotImplementedError):
        mapping.inverse


def test_drop_axes_3():
    mapping = Mapping((1,), n_inputs=2)
    assert(mapping.n_inputs == 2)
    rotation = Rotation2D(60)
    model = rotation | mapping
    assert_allclose(model(1, 2), 1.86602540378)


def test_identity():
    x = np.zeros((2, 3))
    y = np.ones((2, 3))

    ident1 = Identity(1)
    shift = Shift(1)
    rotation = Rotation2D(angle=60)
    model = ident1 & shift | rotation
    assert_allclose(model(1, 2), (-2.098076211353316, 2.3660254037844393))
    res_x, res_y = model(x, y)
    assert_allclose((res_x, res_y),
                       (np.array([[-1.73205081, -1.73205081, -1.73205081],
                                  [-1.73205081, -1.73205081, -1.73205081]]),
                        np.array([[ 1.,  1.,  1.],
                                  [ 1.,  1.,  1.]])))
    assert_allclose(model.inverse(res_x, res_y), (x, y), atol=1.e-10)
