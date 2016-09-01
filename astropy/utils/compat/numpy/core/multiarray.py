# coding: utf-8
# Licensed like numpy; see licenses/NUMPY_LICENSE.rst
"""
Replacement for matmul in 'numpy.core.multiarray'.

Notes
-----
The pure python version here allows matrix multiplications for numpy <= 1.10
"""
from __future__ import division, absolute_import, print_function

import numpy as np

__all__ = ['matmul', 'GE1P10']


def GE1P10(module=np):
    return hasattr(module, 'matmul')


if GE1P10():
    from numpy import matmul

else:
    def matmul(a, b, out=None):
        """Matrix product of two arrays.

        The behavior depends on the arguments in the following way.

        - If both arguments are 2-D they are multiplied like conventional
          matrices.
        - If either argument is N-D, N > 2, it is treated as a stack of
          matrices residing in the last two indexes and broadcast accordingly.
        - If the first argument is 1-D, it is promoted to a matrix by
          prepending a 1 to its dimensions. After matrix multiplication
          the prepended 1 is removed.
        - If the second argument is 1-D, it is promoted to a matrix by
          appending a 1 to its dimensions. After matrix multiplication
          the appended 1 is removed.

        Multiplication by a scalar is not allowed, use ``*`` instead. Note that
        multiplying a stack of matrices with a vector will result in a stack of
        vectors, but matmul will not recognize it as such.

        ``matmul`` differs from ``dot`` in two important ways.

        - Multiplication by scalars is not allowed.
        - Stacks of matrices are broadcast together as if the matrices
          were elements.

        Parameters
        ----------
        a : array_like
            First argument.
        b : array_like
            Second argument.
        out : ndarray, optional
            Output argument. This must have the exact kind that would be returned
            if it was not used. In particular, it must have the right type, must be
            C-contiguous, and its dtype must be the dtype that would be returned
            for `dot(a,b)`. This is a performance feature. Therefore, if these
            conditions are not met, an exception is raised, instead of attempting

        Notes
        -----
        This routine mimicks ``matmul`` using ``einsum``.  See
        http://docs.scipy.org/doc/numpy/reference/generated/numpy.matmul.html
        """
        a = np.asanyarray(a)
        b = np.asanyarray(b)

        if out is None:
            kwargs = {}
        else:
            kwargs = {'out': out}

        if a.ndim >= 2:
            if b.ndim >= 2:
                return np.einsum('...ij,...jk->...ik', a, b, **kwargs)

            if b.ndim == 1:
                return np.einsum('...ij,...j->...i', a, b, **kwargs)

        elif a.ndim == 1 and b.ndim >= 2:
            return np.einsum('...i,...ik->...k', a, b, **kwargs)

        raise ValueError("Scalar operands are not allowed, use '*' instead.")
