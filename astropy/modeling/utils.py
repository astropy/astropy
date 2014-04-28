# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module provides utility functions for the models package
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import re

import numpy as np

from ..extern.six.moves import xrange

__all__ = ['poly_map_domain', 'comb']


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


def format_formula(templ, **parameters):
    """
    Format a model formula with the given LaTeX symbols for its parameters.

    The template strings for formulae use a special syntax that uses pipe
    (``|``) characters around a parameter name to indicate where a given
    parameter's LaTeX symbol should be substituted.

    This format is used instead of the standard Python string template
    formatting due to the prevalence of curly braces (``{}``) in that format as
    well as in LaTeX.
    """

    sub_re = r'\|(?P<param>{0})\|'.format('|'.join(parameters.keys()))

    def sub_repl(m):
        return parameters[m.group('param')]

    return re.sub(sub_re, sub_repl, templ)
