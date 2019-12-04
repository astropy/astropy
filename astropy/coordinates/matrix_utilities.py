# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains utililies used for constructing rotation matrices.
"""
from functools import reduce
import numpy as np

from astropy import units as u
from .angles import Angle


def matrix_product(*matrices):
    """Matrix multiply all arguments together.

    Arguments should have dimension 2 or larger. Larger dimensional objects
    are interpreted as stacks of matrices residing in the last two dimensions.

    This function mostly exists for readability: using `~numpy.matmul`
    directly, one would have ``matmul(matmul(m1, m2), m3)``, etc. For even
    better readability, one might consider using `~numpy.matrix` for the
    arguments (so that one could write ``m1 * m2 * m3``), but then it is not
    possible to handle stacks of matrices. Once only python >=3.5 is supported,
    this function can be replaced by ``m1 @ m2 @ m3``.
    """
    return reduce(np.matmul, matrices)


def matrix_transpose(matrix):
    """Transpose a matrix or stack of matrices by swapping the last two axes.

    This function mostly exists for readability; seeing ``.swapaxes(-2, -1)``
    it is not that obvious that one does a transpose.  Note that one cannot
    use `~numpy.ndarray.T`, as this transposes all axes and thus does not
    work for stacks of matrices.
    """
    return matrix.swapaxes(-2, -1)


def rotation_matrix(angle, axis='z', unit=None):
    """
    Generate matrices for rotation by some angle around some axis.

    Parameters
    ----------
    angle : convertible to `~astropy.coordinates.Angle`
        The amount of rotation the matrices should represent.  Can be an array.
    axis : str or array_like
        Either ``'x'``, ``'y'``, ``'z'``, or a (x,y,z) specifying the axis to
        rotate about. If ``'x'``, ``'y'``, or ``'z'``, the rotation sense is
        counterclockwise looking down the + axis (e.g. positive rotations obey
        left-hand-rule).  If given as an array, the last dimension should be 3;
        it will be broadcast against ``angle``.
    unit : UnitBase, optional
        If ``angle`` does not have associated units, they are in this
        unit.  If neither are provided, it is assumed to be degrees.

    Returns
    -------
    rmat : `numpy.matrix`
        A unitary rotation matrix.
    """
    if unit is None:
        unit = u.degree

    angle = Angle(angle, unit=unit)

    s = np.sin(angle)
    c = np.cos(angle)

    # use optimized implementations for x/y/z
    try:
        i = 'xyz'.index(axis)
    except TypeError:
        axis = np.asarray(axis)
        axis = axis / np.sqrt((axis * axis).sum(axis=-1, keepdims=True))
        R = (axis[..., np.newaxis] * axis[..., np.newaxis, :] *
             (1. - c)[..., np.newaxis, np.newaxis])

        for i in range(0, 3):
            R[..., i, i] += c
            a1 = (i + 1) % 3
            a2 = (i + 2) % 3
            R[..., a1, a2] += axis[..., i] * s
            R[..., a2, a1] -= axis[..., i] * s

    else:
        a1 = (i + 1) % 3
        a2 = (i + 2) % 3
        R = np.zeros(angle.shape + (3, 3))
        R[..., i, i] = 1.
        R[..., a1, a1] = c
        R[..., a1, a2] = s
        R[..., a2, a1] = -s
        R[..., a2, a2] = c

    return R


def angle_axis(matrix):
    """
    Angle of rotation and rotation axis for a given rotation matrix.

    Parameters
    ----------
    matrix : array_like
        A 3 x 3 unitary rotation matrix (or stack of matrices).

    Returns
    -------
    angle : `~astropy.coordinates.Angle`
        The angle of rotation.
    axis : array
        The (normalized) axis of rotation (with last dimension 3).
    """
    m = np.asanyarray(matrix)
    if m.shape[-2:] != (3, 3):
        raise ValueError('matrix is not 3x3')

    axis = np.zeros(m.shape[:-1])
    axis[..., 0] = m[..., 2, 1] - m[..., 1, 2]
    axis[..., 1] = m[..., 0, 2] - m[..., 2, 0]
    axis[..., 2] = m[..., 1, 0] - m[..., 0, 1]
    r = np.sqrt((axis * axis).sum(-1, keepdims=True))
    angle = np.arctan2(r[..., 0],
                       m[..., 0, 0] + m[..., 1, 1] + m[..., 2, 2] - 1.)
    return Angle(angle, u.radian), -axis / r
