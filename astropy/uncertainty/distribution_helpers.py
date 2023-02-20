# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Helpers for letting numpy functions interact with Distributions."""

import numpy as np

from astropy.units.quantity_helper.function_helpers import FunctionAssigner

# This module should not really be imported, but we define __all__
# such that sphinx can typeset the functions with docstrings.
# The latter are added to __all__ at the end.
__all__ = [
    "DISTRIBUTION_SAFE_FUNCTIONS",
    "FUNCTION_HELPERS",
    "DISPATCHED_FUNCTIONS",
    "UNSUPPORTED_FUNCTIONS",
]


DISTRIBUTION_SAFE_FUNCTIONS = set()
"""Set of functions that work fine on Distributions already.

Most of these internally use `numpy.ufunc` or other functions that
are already covered.
"""

FUNCTION_HELPERS = {}
"""Functions with implementations usable with conversion to regular arrays.

Should return args, kwargs, out, with the first two passed on for a super() call
and out the actual distribution in which output is meant to be stored.

It should raise `NotImplementedError` if one of the arguments is
a Distribution when it should not be or vice versa.
"""

DISPATCHED_FUNCTIONS = {}
"""Dict of functions that provide the numpy function's functionality.

These should return ``result, out``, with result the regular array,
and ``out`` a possible Distribution to put the output in.

It should raise `NotImplementedError` if one of the arguments is
a Distribution when it should not be or vice versa.
"""

UNSUPPORTED_FUNCTIONS = set()
"""Set of numpy functions that are not supported for Distributions.

For most, distributions simply makes no sense, but for others it may have
been lack of time.  Issues or PRs for support for functions are welcome.
"""


function_helper = FunctionAssigner(FUNCTION_HELPERS)

dispatched_function = FunctionAssigner(DISPATCHED_FUNCTIONS)


def is_distribution(x):
    from astropy.uncertainty import Distribution

    return isinstance(x, Distribution)


def get_n_samples(*arrays):
    for array in arrays:
        if is_distribution(array):
            return array.n_samples

    raise RuntimeError("no Distribution found! Please raise an issue.")


@function_helper
def concatenate(arrays, axis=0, out=None, dtype=None, casting="same_kind"):
    n_samples = get_n_samples(*arrays, out)
    converted = tuple(
        array.distribution
        if is_distribution(array)
        else (
            np.broadcast_to(
                array[..., np.newaxis], array.shape + (n_samples,), subok=True
            )
            if getattr(array, "shape", ())
            else array
        )
        for array in arrays
    )
    if axis < 0:
        axis = axis - 1  # not in-place, just in case.
    kwargs = dict(axis=axis, dtype=dtype, casting=casting)
    if is_distribution(out):
        kwargs["out"] = out.distribution
    return (converted,), kwargs, out


# Add any dispatched or helper function that has a docstring to
# __all__, so they will be typeset by sphinx. The logic is that for
# those presumably the use of the mask is not entirely obvious.
__all__ += sorted(
    helper.__name__
    for helper in (set(FUNCTION_HELPERS.values()) | set(DISPATCHED_FUNCTIONS.values()))
    if helper.__doc__
)
