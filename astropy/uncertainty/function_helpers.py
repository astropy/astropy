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
for passing the distribution data to the numpy function, by returning
``args, kwargs, out``. Here, the former two are passed on, while ``out``
is used to indicate whether there was an output argument.  If ``out`` is
set to `True`, then no further processing should be done; otherwise, it
it is assumed that the function operates on unwrapped distributions and
that the results need to be rewrapped as |Distribution|.

The function should raise `NotImplementedError` if one of the arguments is a
distribution when it should not be or vice versa.

"""


UNSUPPORTED_FUNCTIONS = set()
"""Set of numpy functions that are not supported for distributions.

For most, distributions simply make no sense, but for others it may have
been lack of time.  Issues or PRs for support for functions are welcome.
"""


function_helper = FunctionAssigner(FUNCTION_HELPERS)
dispatched_function = FunctionAssigner(DISPATCHED_FUNCTIONS)


def is_distribution(x):
    from astropy.uncertainty import Distribution

    return isinstance(x, Distribution)


def get_n_samples(*arrays):
    """Get n_samples from the first Distribution amount arrays.

    The logic of getting ``n_samples`` from the first |Distribution|
    is that the code will raise an appropriate exception later if
    distributions do not have the same ``n_samples``.
    """
    # TODO: add verification if another function needs it.
    for array in arrays:
        if is_distribution(array):
            return array.n_samples

    raise RuntimeError("no Distribution found! Please raise an issue.")


@function_helper
def empty_like(prototype, dtype=None, *args, **kwargs):
    dtype = prototype._get_distribution_dtype(
        prototype.dtype if dtype is None else dtype, prototype.n_samples
    )
    return (prototype, dtype) + args, kwargs, None


@function_helper
def broadcast_arrays(*args, subok=False):
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
    return args, {"subok": True}, True


@function_helper
def concatenate(arrays, axis=0, out=None, dtype=None, casting="same_kind"):
    """Concatenate arrays.

    Like `numpy.concatenate`, but any array that is not already a |Distribution|
    is turned into one with identical samples.
    """
    n_samples = get_n_samples(*arrays, out)
    converted = tuple(
        array.distribution
        if is_distribution(array)
        else (
            np.broadcast_to(
                array[..., np.newaxis], array.shape + (n_samples,), subok=True
            )
            if getattr(array, "shape", False)
            else array
        )
        for array in arrays
    )
    if axis < 0:
        axis = axis - 1  # not in-place, just in case.
    kwargs = dict(axis=axis, dtype=dtype, casting=casting)
    if out is not None:
        if is_distribution(out):
            kwargs["out"] = out.distribution
        else:
            raise NotImplementedError
    return (converted,), kwargs, out


# Add any dispatched or helper function that has a docstring to __all__, so
# they will be typeset by sphinx. The logic is that for those presumably the
# way distributions are dealt with is not entirely obvious.
__all__ += sorted(  # noqa: PLE0605
    helper.__name__
    for helper in (set(FUNCTION_HELPERS.values()) | set(DISPATCHED_FUNCTIONS.values()))
    if helper.__doc__
)
