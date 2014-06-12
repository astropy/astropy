# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module provides utility functions for the models package
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from itertools import izip_longest

import numpy as np

from ..extern.six.moves import xrange

__all__ = ['check_broadcast', 'poly_map_domain', 'comb']



def check_broadcast(shape_a, shape_b, *shapes):
    """
    Determines whether two or more Numpy arrays can be broadcast with each
    other based on their shape tuple alone.

    Parameters
    ----------
    shape_a : tuple
        The shape tuple of the first array.
    shape_b : tuple
        The shape tuple of the second array.
    *shapes : tuple
        Any subsequent shapes to include in the comparison.

    Returns
    -------
    can_broadcast : `tuple`, `None`
        If all the supplied shapes can broadcast with each other, returns
        `None`, otherwise returns a two-tuple of the indices of the first two
        shapes that were determined not to broadcast with each other.
    """

    all_shapes = ((reversed(shape_a), reversed(shape_b)) +
                  tuple(reversed(shape) for shape in shapes))

    for dims in izip_longest(*all_shapes, fillvalue=1):
        max_dim = None
        max_dim_idx = None
        for idx, dim in enumerate(dims):
            if dim == 1:
                continue

            if max_dim is None:
                # The first dimension of size greater than 1
                max_dim = dim
                max_dim_idx = idx
            elif dim != max_dim:
                return (max_dim_idx, idx)

    return None


def poly_map_domain(oldx, domain, window):
    """
    Map domain into window by shifting and scaling.

    Parameters
    ----------
    oldx : array
          original coordinates
    domain : list or tuple of length 2
          function domain
    window : list or tuple of length 2
          range into which to map the domain
    """
    domain = np.array(domain, dtype=np.float64)
    window = np.array(window, dtype=np.float64)
    scl = (window[1] - window[0]) / (domain[1] - domain[0])
    off = (window[0] * domain[1] - window[1] * domain[0]) / (domain[1] - domain[0])
    return off + scl * oldx


def comb(N, k):
    """
    The number of combinations of N things taken k at a time.

    Parameters
    ----------
    N : int, array
        Number of things.
    k : int, array
        Number of elements taken.

    """
    if (k > N) or (N < 0) or (k < 0):
        return 0
    val = 1
    for j in xrange(min(k, N - k)):
        val = (val * (N - j)) / (j + 1)
    return val


def array_repr_oneline(array):
    """
    Represents a multi-dimensional Numpy array flattened onto a single line.
    """

    r = np.array2string(array, separator=',', suppress_small=True)
    return ' '.join(l.strip() for l in r.splitlines())
