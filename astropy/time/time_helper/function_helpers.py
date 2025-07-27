"""
Helpers for overriding numpy functions in
`~astropy.time.Time.__array_function__`.
"""

from typing import TYPE_CHECKING, TypeVar, Union

import numpy as np

from astropy.units.quantity_helper.function_helpers import FunctionAssigner
from astropy.utils.compat import NUMPY_LT_2_0

if TYPE_CHECKING:
    from astropy.time import Time, TimeDelta

if NUMPY_LT_2_0:
    from numpy.core.multiarray import normalize_axis_index
else:
    from numpy.lib.array_utils import normalize_axis_index

# TODO: Fill this in with functions that don't make sense for times
UNSUPPORTED_FUNCTIONS = {}
# Functions that return the final result of the numpy function
CUSTOM_FUNCTIONS = {}

custom_function = FunctionAssigner(CUSTOM_FUNCTIONS)


@custom_function
def linspace(tstart, tstop, *args, **kwargs):
    from astropy.time import Time

    if isinstance(tstart, Time):
        if not isinstance(tstop, Time):
            return NotImplemented

    if kwargs.get("retstep"):
        offsets, step = np.linspace(
            np.zeros(tstart.shape), np.ones(tstop.shape), *args, **kwargs
        )
        tdelta = tstop - tstart
        return tstart + tdelta * offsets, tdelta * step
    else:
        offsets = np.linspace(
            np.zeros(tstart.shape), np.ones(tstop.shape), *args, **kwargs
        )
        return tstart + (tstop - tstart) * offsets


TimeLike = TypeVar("TimeLike", bound=Union["Time", "TimeDelta"])


@custom_function
def zeros_like(
    a: TimeLike,
    dtype=None,
    order="K",
    subok=True,
    shape=None,
    **kwargs,
) -> TimeLike:
    """Create a new "zero" Time or TimeDelta object with the properties of a.

    The other parameters allow to override what is used to create ``jd1`` and ``jd2``
    (but one cannot override ``dtype``).

    Parameters
    ----------
    a : Time or TimeDelta
        The object to mimic.
    dtype : None
        This parameter is not used, as the dtype is determined by the class of `a`.
    order : {'C', 'F', 'A', 'K'}, optional
        Memory layout order for the output array. Default is 'K', which means the output
        will have the same order as ``a``.
    subok : bool, optional
        If True, the output will be a subclass of `a`'s class. If False, it will always
        return a base class instance (e.g., ``Time`` or ``TimeDelta``).
    shape : tuple, optional
        The shape of the output array. If None, it will match the shape of `a`.
    kwargs : dict, optional
        Additional keyword arguments passed to the array creation functions
        ``np.zeros_like`` and ``np.full_like``.

    Returns
    -------
    Time or TimeDelta
        A new instance of the same class as ``a``, with ``jd1`` and ``jd2`` initialized
        to J2000.0 for ``Time`` or zero for ``TimeDelta`` and the specified shape and
        order.
    """
    from astropy.time import Time

    if dtype is not None:
        raise ValueError(f"Cannot set dtype for {a.__class__.__name__} creation.")
    if shape is None:
        shape = a.shape

    # Get class, where the only other choice is TimeDelta due to the way custom_function
    is_time = isinstance(a, Time)

    # Create jd1 and jd2 arrays, with the same shape as `a`, but filled with the
    # appropriate value for the class. For Time, this is the arbitrary time J2000.0 that
    # will work with ERFA.
    fill = 2451544.5 if is_time else 0.0
    jd1 = np.full_like(a.jd1, fill, shape=shape, order=order, subok=subok, **kwargs)
    jd2 = np.zeros_like(a.jd2, shape=shape, order=order, subok=subok, **kwargs)

    # We first create a JD format instance, and then convert the
    # format and add back other attributes. This will break if location was
    # an array that cannot be broadcast to the new shape, but that seems
    # reasonable; it is unclear what the user wants in that case.
    tkw = {"location": a.location} if is_time else {}
    out = a.__class__(jd1, jd2, format="jd", scale=a.scale, copy=False, **tkw)
    for attr in ("format", "precision", "in_subfmt", "out_subfmt"):
        setattr(out, attr, getattr(a, attr))

    return out


def _combine_helper(func, arrays, axis, out, dtype):
    # Apply on arrays of bool with the same shape, to get the final shape,
    # and to avoid having test that shapes match, etc.
    empties = [np.empty(shape=np.shape(array), dtype=bool) for array in arrays]
    shape = func(empties, axis=axis).shape
    axis = normalize_axis_index(axis, len(shape))
    if out is None:
        out = np.zeros_like(arrays[0], shape=shape, dtype=dtype)
    return axis, out


@custom_function
def concatenate(arrays, axis=0, out=None, dtype=None, casting="same_kind"):
    axis, out = _combine_helper(np.concatenate, arrays, axis, out, dtype)

    offset = 0
    for array in arrays:
        n_el = array.shape[axis]
        out[(slice(None),) * axis + (slice(offset, offset + n_el),)] = array
        offset += n_el

    return out


@custom_function
def stack(arrays, axis=0, out=None, *, dtype=None, casting="same_kind"):
    axis, out = _combine_helper(np.stack, arrays, axis, out, dtype)
    for i, array in enumerate(arrays):
        out[(slice(None),) * axis + (i,)] = array

    return out
