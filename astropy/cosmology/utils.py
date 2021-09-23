# Licensed under a 3-clause BSD style license - see LICENSE.rst

from math import inf

import numpy as np

from astropy.utils import isiterable

__all__ = []  # nothing is publicly scoped


def vectorize_if_needed(f, *x, **vkw):
    """Helper function to vectorize scalar functions on array inputs.

    Parameters
    ----------
    f : callable
        'f' must accept positional arguments and no mandatory keyword
        arguments.
    *x
        Arguments into ``f``.
    **vkw
        Keyword arguments into :class:`numpy.vectorize`.

    Examples
    --------
    >>> func = lambda x: x ** 2
    >>> vectorize_if_needed(func, 2)
    4
    >>> vectorize_if_needed(func, [2, 3])
    array([4, 9])
    """
    return np.vectorize(f, **vkw)(*x) if any(map(isiterable, x)) else f(*x)


def inf_like(x):
    """Return the shape of x with value infinity and dtype='float'.

    Preserves 'shape' for both array and scalar inputs.
    But always returns a float array, even if x is of integer type.

    Parameters
    ----------
    x : scalar or array-like
        Must work with functions `numpy.isscalar` and `numpy.full_like` (if `x`
        is not a scalar`

    Returns
    -------
    `math.inf` or ndarray[float] thereof
        Returns a scalar `~math.inf` if `x` is a scalar, an array of floats
        otherwise.

    Examples
    --------
    >>> inf_like(0.)  # float scalar
    inf
    >>> inf_like(1)  # integer scalar should give float output
    inf
    >>> inf_like([0., 1., 2., 3.])  # float list
    array([inf, inf, inf, inf])
    >>> inf_like([0, 1, 2, 3])  # integer list should give float output
    array([inf, inf, inf, inf])
    """
    return inf if np.isscalar(x) else np.full_like(x, inf, dtype=float)
