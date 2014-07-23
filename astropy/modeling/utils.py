# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module provides utility functions for the models package
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import inspect
import itertools
import keyword
import re
import textwrap

import numpy as np

from ..extern.six import iteritems
from ..extern.six.moves import xrange, zip_longest
from ..utils import find_current_module

__all__ = ['check_broadcast', 'poly_map_domain', 'comb']


class IncompatibleShapeError(ValueError):
    def __init__(self, shape_a, shape_a_idx, shape_b, shape_b_idx):
        super(IncompatibleShapeError, self).__init__(
                shape_a, shape_a_idx, shape_b, shape_b_idx)


def check_broadcast(*shapes):
    """
    Determines whether two or more Numpy arrays can be broadcast with each
    other based on their shape tuple alone.

    Parameters
    ----------
    *shapes : tuple
        All shapes to include in the comparison.  If only one shape is given it
        is passed through unmodified.  If no shapes are given returns an empty
        `tuple`.

    Returns
    -------
    broadcast : `tuple`
        If all shapes are mutually broadcastable, returns a tuple of the full
        broadcast shape.
    """

    if len(shapes) == 0:
        return ()
    elif len(shapes) == 1:
        return shapes[0]

    reversed_shapes = (reversed(shape) for shape in shapes)

    full_shape = []

    for dims in zip_longest(*reversed_shapes, fillvalue=1):
        max_dim = 1
        max_dim_idx = None
        for idx, dim in enumerate(dims):
            if dim == 1:
                continue

            if max_dim == 1:
                # The first dimension of size greater than 1
                max_dim = dim
                max_dim_idx = idx
            elif dim != max_dim:
                raise IncompatibleShapeError(
                    shapes[max_dim_idx], max_dim_idx, shapes[idx], idx)

        full_shape.append(max_dim)

    return tuple(full_shape[::-1])


_ARGNAME_RE = re.compile(r'^[A-Za-z][A-Za-z_]*')
"""
Regular expression used my make_func which limits the allowed argument
names for the created function.  Only valid Python variable names in
the ASCII range and not beginning with '_' are allowed, currently.
"""


def make_func_with_sig(func, *args, **kwargs):
    """
    Make a new function from an existing function but with the desired
    signature.

    The desired signature must of course be compatible with the arguments
    actually accepted by the input function.

    The ``*args`` are strings that should be the names of the positional
    arguments.  ``**kwargs`` can map names of keyword arguments to their
    default values.

    Alternatively, ``*args`` may be a list of zero or more strings *followed*
    by 2-tuples of ``(keyword, value)`` pairs representing the keyword
    arguments and their default values.  This ensures the order of keyword
    arguments in the signature.

    Note, the names may only be valid Python variable names.
    """

    pos_args = []
    key_args = []

    # Check that all the argument names are valid
    for item in itertools.chain(args, iteritems(kwargs)):
        if isinstance(item, tuple):
            argname = item[0]
            key_args.append(item)
        else:
            argname = item
            pos_args.append(item)

        if keyword.iskeyword(argname) or not _ARGNAME_RE.match(argname):
            raise SyntaxError('invalid argument name: {0}'.format(argname))

    def_signature = [', '.join(pos_args)]
    call_signature = def_signature[:]

    name = func.__name__

    global_vars = {'__{0}__func'.format(name): func}
    local_vars = {}
    # Make local variables to handle setting the default args
    for idx, item in enumerate(key_args):
        key, value = item
        default_var = '_kwargs{0}'.format(idx)
        local_vars[default_var] = value
        def_signature.append(', {0}={1}'.format(key, default_var))
        call_signature.append(', {0}={0}'.format(key))

    def_signature = ''.join(def_signature).lstrip(', ')
    call_signature = ''.join(call_signature).lstrip(', ')
    # The lstrip is in case there were *no* positional arguments (a rare case)
    # in any context this will actually be used...
    template = textwrap.dedent("""\
    def {name}({sig1}):
        return __{name}__func({sig2})
    """.format(name=name, sig1=def_signature, sig2=call_signature))

    mod = find_current_module(2)
    if mod:
        filename = mod.__file__
        modname = mod.__name__
    else:
        filename = '<string>'
        modname = '__main__'

    code = compile(template, filename, 'single')

    eval(code, global_vars, local_vars)

    new_func = local_vars[name]
    new_func.__module__ = modname
    new_func.__doc__ = func.__doc__

    return new_func


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
