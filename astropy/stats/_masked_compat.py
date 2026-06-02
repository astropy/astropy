# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Private helpers that let ``astropy.stats`` functions written against
`numpy.ma` transparently accept astropy ``Masked`` inputs
(``MaskedNDArray`` / ``MaskedQuantity``).
"""

import functools
import inspect

import numpy as np

from astropy.units import Quantity
from astropy.utils.masked import Masked

__all__ = ["support_masked"]


def _to_numpy_masked(data):
    """
    Convert an astropy ``Masked`` instance to ``(np.ma.MaskedArray, unit)``.

    ``unit`` is ``None`` unless ``data`` was a ``MaskedQuantity``.  Only plain
    numeric data (optionally a `~astropy.units.Quantity`) is supported; other
    flavours (structured dtypes, ``Time``, etc.) raise `NotImplementedError`
    rather than silently round-tripping through `numpy.ma`.
    """
    unmasked = data.unmasked
    if isinstance(unmasked, Quantity):
        unit = unmasked.unit
        values = np.asarray(unmasked.value)
    else:
        unit = None
        values = np.asarray(unmasked)

    # ``kind`` is one of "biufc" only for bool/int/uint/float/complex; strings,
    # datetimes, objects and structured dtypes (kind "V") are all excluded.
    if values.dtype.kind not in "biufc":
        raise NotImplementedError(
            "astropy.stats only supports Masked inputs backed by plain numeric "
            f"data, not dtype {values.dtype!r}."
        )

    mask = np.broadcast_to(np.asarray(data.mask), values.shape)
    return np.ma.MaskedArray(values, mask=mask.copy(), copy=True), unit


def _attach_unit(value, unit):
    return value if unit is None else value << unit


def _as_masked(result, unit):
    """Wrap a single array result back into a ``Masked`` instance."""
    if isinstance(result, np.ma.MaskedArray):
        values, mask = np.asarray(result.data), np.ma.getmaskarray(result)
    else:
        # ``masked=False`` style outputs: clipped/invalid entries are NaN
        # (when an axis is given) or have been dropped (axis=None); flag any
        # remaining non-finite values as masked.
        values = np.asarray(result)
        mask = ~np.isfinite(values)
    return Masked(_attach_unit(values, unit), mask=mask)


def _repack(result, unit):
    """
    Re-wrap a function result as ``Masked``.

    A bare array result becomes a ``Masked`` array.  A tuple is treated as
    ``(primary_array, *extras)`` -- the primary array is re-masked while the
    extras (e.g. the lower/upper bounds from ``sigma_clip``) only have the unit
    re-attached.
    """
    if isinstance(result, tuple):
        primary, *extras = result
        return (_as_masked(primary, unit), *(_attach_unit(e, unit) for e in extras))
    return _as_masked(result, unit)


def support_masked(data_arg="data"):
    """
    Decorate a ``numpy.ma``-based function to also accept astropy ``Masked``.

    When the ``data_arg`` argument is an astropy ``Masked`` instance it is
    converted to a `numpy.ma.MaskedArray` (with any unit stripped), the wrapped
    function runs unchanged, and the result is converted back so that ``Masked``
    input yields ``Masked`` output (and ``MaskedQuantity`` input yields
    ``MaskedQuantity`` output).  Non-``Masked`` inputs are passed through
    untouched, so existing behaviour is unaffected.

    Parameters
    ----------
    data_arg : str, optional
        Name of the parameter holding the data array.
    """

    def decorator(func):
        sig = inspect.signature(func)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            bound = sig.bind(*args, **kwargs)
            value = bound.arguments.get(data_arg)
            if not isinstance(value, Masked):
                return func(*args, **kwargs)
            converted, unit = _to_numpy_masked(value)
            bound.arguments[data_arg] = converted
            return _repack(func(*bound.args, **bound.kwargs), unit)

        return wrapper

    return decorator
