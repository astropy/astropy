# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
from numpy.testing.utils import assert_allclose, assert_array_equal

from ... import units as u
from ..matrix_utilities import rotation_matrix, angle_axis


def test_rotation_matrix():
    assert_array_equal(rotation_matrix(0*u.deg, 'x'), np.eye(3))

    assert_allclose(rotation_matrix(90*u.deg, 'y'), [[ 0, 0,-1],
                                                     [ 0, 1, 0],
                                                     [ 1, 0, 0]], atol=1e-12)

    assert_allclose(rotation_matrix(-90*u.deg, 'z'), [[ 0,-1, 0],
                                                      [ 1, 0, 0],
                                                      [ 0, 0, 1]], atol=1e-12)

    assert_allclose(rotation_matrix(45*u.deg, 'x'),
                    rotation_matrix(45*u.deg, [1, 0, 0]))
    assert_allclose(rotation_matrix(125*u.deg, 'y'),
                    rotation_matrix(125*u.deg, [0, 1, 0]))
    assert_allclose(rotation_matrix(-30*u.deg, 'z'),
                    rotation_matrix(-30*u.deg, [0, 0, 1]))

    assert_allclose(np.dot(rotation_matrix(180*u.deg, [1, 1, 0]), [1, 0, 0]),
                    [0, 1, 0], atol=1e-12)

    #make sure it also works for very small angles
    assert_allclose(rotation_matrix(0.000001*u.deg, 'x'),
                    rotation_matrix(0.000001*u.deg, [1, 0, 0]))


def test_angle_axis():
    m1 = rotation_matrix(35*u.deg, 'x')
    an1, ax1 = angle_axis(m1)

    assert an1 - 35*u.deg < 1e-10*u.deg
    assert_allclose(ax1, [1, 0, 0])


    m2 = rotation_matrix(-89*u.deg, [1, 1, 0])
    an2, ax2 = angle_axis(m2)

    assert an2 - 89*u.deg < 1e-10*u.deg
    assert_allclose(ax2, [-2**-0.5, -2**-0.5, 0])
