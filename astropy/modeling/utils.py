# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module provides utility functions for the models package
"""
from __future__ import division
import numpy as np

__all__ = ['poly_map_domain', 'comb', 'InputParameterError', 'ModelDefinitionError']


class ModelDefinitionError(Exception):
    """
    Called when models are defined in a wrong way
    """
    def __init__(self, message):
        self._message = message

    def __str__(self):
        return self._message


class InputParameterError(Exception):

    """
    Called when there's a problem with input parameters.
    """
    def __init__(self, message):
        self._message = message

    def __str__(self):
        return self._message


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
