# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Helpers for letting numpy functions interact with distributions.

The module supplies helper routines for numpy functions that propagate
distributions appropriately., for use in the ``__array_function__``
implementation of `~astropy.uncertainty.core.Distribution`.  They are not
very useful on their own, but the ones with docstrings are included in
the documentation so that there is a place to find out how the distributions
are interpreted.

"""
import numpy as np

from astropy.units.quantity_helper.function_helpers import FunctionAssigner

# This module should not really be imported, but we define __all__
# such that sphinx can typeset the functions with docstrings.
# The latter are added to __all__ at the end.
__all__ = [
    "DISTRIBUTION_SAFE_FUNCTIONS",
    "DISPATCHED_FUNCTIONS",
    "UNSUPPORTED_FUNCTIONS",
]


DISTRIBUTION_SAFE_FUNCTIONS = set()
"""Set of functions that work fine on Distribution classes already.

Most of these internally use `numpy.ufunc` or other functions that
are already covered.
"""

DISPATCHED_FUNCTIONS = {}
"""Dict of functions that provide the numpy function's functionality.

These are for more complicated versions where the numpy function itself
cannot easily be used.  It should return the result of the function.

It should raise `NotImplementedError` if one of the arguments is a
distribution when it should not be or vice versa.
"""


FUNCTION_HELPERS = {}
"""Dict of functions for which Distribution can be used after some conversions.

The `dict` is keyed by the numpy function and the values are functions
that take the input arguments of the numpy function and organize these
for passing the distribution data to the numpy function.
"""


UNSUPPORTED_FUNCTIONS = set()
"""Set of numpy functions that are not supported for distributions.

For most, distributions simply make no sense, but for others it may have
been lack of time.  Issues or PRs for support for functions are welcome.
"""


function_helper = FunctionAssigner(FUNCTION_HELPERS)
dispatched_function = FunctionAssigner(DISPATCHED_FUNCTIONS)


@function_helper
def empty_like(prototype, dtype=None, *args, **kwargs):
    dtype = prototype._get_distribution_dtype(
        prototype.dtype if dtype is None else dtype, prototype.n_samples
    )
    return (prototype, dtype) + args, kwargs, None


@function_helper
def broadcast_arrays(*args, subok=True):
    """Broadcast arrays to a common shape.

    Like `numpy.broadcast_arrays`, applied to both distributions and other data.
    Note that ``subok`` is taken to mean whether or not subclasses of
    the distribution are allowed, i.e., for ``subok=False``,
    `~astropy.uncertainty.NdarrayDistribution` instances will be returned.
    """
    if not subok:
        args = tuple(
            arg.view(np.ndarray) if isinstance(arg, np.ndarray) else np.array(arg)
            for arg in args
        )
    return args, {"subok": True}, None
