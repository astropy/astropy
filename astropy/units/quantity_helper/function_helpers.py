# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Helpers for overriding numpy functions for Quantity."""

import numpy as np

from astropy.units.core import (
    UnitsError, UnitConversionError, UnitTypeError,
    dimensionless_unscaled, get_current_unit_registry)
from .helpers import _d, get_converter


UNSUPPORTED_FUNCTIONS = set()
FUNCTION_HELPERS = {}
DISPATCHED_FUNCTIONS = {}


def function_helper(f):
    FUNCTION_HELPERS[getattr(np, f.__name__)] = f
    return f


def dispatched_function(f):
    DISPATCHED_FUNCTIONS[getattr(np, f.__name__)] = f
    return f


def invariant_a_helper(a, *args, **kwargs):
    return (a.view(np.ndarray),) + args, kwargs, a.unit, None


FUNCTION_HELPERS[np.copy] = invariant_a_helper
FUNCTION_HELPERS[np.asfarray] = invariant_a_helper
FUNCTION_HELPERS[np.zeros_like] = invariant_a_helper
FUNCTION_HELPERS[np.ones_like] = invariant_a_helper
FUNCTION_HELPERS[np.real_if_close] = invariant_a_helper
FUNCTION_HELPERS[np.sort_complex] = invariant_a_helper
FUNCTION_HELPERS[np.resize] = invariant_a_helper


def invariant_m_helper(m, *args, **kwargs):
    return (m.view(np.ndarray),) + args, kwargs, m.unit, None


FUNCTION_HELPERS[np.tril] = invariant_m_helper
FUNCTION_HELPERS[np.triu] = invariant_m_helper


@function_helper
def empty_like(prototype, *args, **kwargs):
    return (prototype.view(np.ndarray),) + args, kwargs, prototype.unit, None


@function_helper
def sinc(x):
    from astropy.units.si import radian
    try:
        x = x.to_value(radian)
    except UnitsError:
        raise UnitTypeError("Can only apply 'sinc' function to "
                            "quantities with angle units")
    return (x,), {}, dimensionless_unscaled, None


@dispatched_function
def unwrap(p, discont=None, axis=-1):
    from astropy.units.si import radian
    if discont is None:
        discont = np.pi << radian

    try:
        p = p << radian
        discont = discont.to_value(radian)
    except UnitsError:
        raise UnitTypeError("Can only apply 'unwrap' function to "
                            "quantities with angle units")

    return p._wrap_function(np.unwrap.__wrapped__, discont, axis=axis)


@function_helper
def argpartition(a, *args, **kwargs):
    return (a.view(np.ndarray),) + args, kwargs, None, None


@function_helper
def full_like(a, fill_value, *args, **kwargs):
    unit = a.unit if kwargs.get('subok', True) else None
    return (a.view(np.ndarray),
            a._to_own_unit(fill_value)) + args, kwargs, unit, None


@function_helper
def putmask(a, mask, values):
    from astropy.units import Quantity
    if isinstance(a, Quantity):
        return (a.view(np.ndarray), mask,
                a._to_own_unit(values)), {}, a.unit, None
    elif isinstance(values, Quantity):
        return (a, mask,
                values.to_value(dimensionless_unscaled)), {}, None, None
    else:
        return NotImplemented


@function_helper
def place(arr, mask, vals):
    from astropy.units import Quantity
    if isinstance(arr, Quantity):
        return (arr.view(np.ndarray), mask,
                arr._to_own_unit(vals)), {}, arr.unit, None
    elif isinstance(vals, Quantity):
        return (arr, mask,
                vals.to_value(dimensionless_unscaled)), {}, None, None
    else:
        return NotImplemented


@function_helper
def copyto(dst, src, *args, **kwargs):
    from astropy.units import Quantity
    if isinstance(dst, Quantity):
        return ((dst.view(np.ndarray), dst._to_own_unit(src)) + args,
                kwargs, None, None)
    elif isinstance(src, Quantity):
        return ((dst,  src.to_value(dimensionless_unscaled)) + args,
                kwargs, None, None)
    else:
        return NotImplemented


@function_helper
def nan_to_num(x, copy=True, nan=0.0, posinf=None, neginf=None):
    nan = x._to_own_unit(nan)
    if posinf is not None:
        posinf = x._to_own_unit(posinf)
    if neginf is not None:
        neginf = x._to_own_unit(neginf)
    return ((x.view(np.ndarray),),
            dict(copy=True, nan=nan, posinf=posinf, neginf=neginf),
            x.unit, None)


def _iterable_helper(arrays, out=None, **kwargs):
    # Get unit from input if possible.
    from astropy.units import Quantity

    for a in arrays:
        if hasattr(a, 'unit'):
            unit = a.unit
            arrays = tuple(Quantity(_a, subok=True, copy=False).to_value(unit)
                           for _a in arrays)
            break
    else:
        # No input has even a unit, so the only way we can have gotten
        # here is for output to be a Quantity.  Assume all arrays are
        # are standard ndarray and set unit to dimensionless.
        unit = dimensionless_unscaled

    if out is not None:
        if isinstance(out, Quantity):
            kwargs['out'] = out.view(np.ndarray)
        else:
            # TODO: for an ndarray output, we could in principle
            # try converting all Quantity to dimensionless.
            return NotImplemented

    return arrays, unit, out, kwargs


@function_helper
def concatenate(arrays, axis=0, out=None):
    # TODO: make this smarter by creating an appropriately shaped
    # empty output array and just filling it.
    arrays, unit, out, kwargs = _iterable_helper(arrays, out,
                                                 axis=axis)
    return (arrays,), kwargs, unit, out


@function_helper
def choose(a, choices, out=None, **kwargs):
    choices, unit, out, kwargs = _iterable_helper(choices, out, **kwargs)
    return (a, choices,), kwargs, unit, out


@function_helper
def select(condlist, choicelist, default=0):
    choicelist, unit, _, _ = _iterable_helper(choicelist, None)
    if default != 0:
        default = (1 * unit)._to_own_unit(default)
    return (condlist, choicelist, default), {}, unit, None


@function_helper
def append(arr, values, *args, **kwargs):
    from astropy.units import Quantity
    if isinstance(arr, Quantity):
        return (arr.view(np.ndarray),
                arr._to_own_unit(values)) + args, kwargs, arr.unit, None
    else:  # values must be Quantity
        unit = getattr(arr, 'unit', dimensionless_unscaled)
        return (arr, values.to_value(unit)) + args, kwargs, unit, None


@function_helper
def insert(arr, obj, values, *args, **kwargs):
    from astropy.units import Quantity, dimensionless_unscaled

    if isinstance(obj, Quantity):
        return NotImplemented

    if isinstance(arr, Quantity):
        return (arr.view(np.ndarray), obj,
                arr._to_own_unit(values)) + args, kwargs, arr.unit, None
    else:  # values must be Quantity
        unit = getattr(arr, 'unit', dimensionless_unscaled)
        return (arr, obj, values.to_value(unit)) + args, kwargs, unit, None


@function_helper
def pad(array, pad_width, mode='constant', **kwargs):
    # pad dispatches only on array, so that must be a Quantity.
    for key in 'constant_values', 'end_values':
        value = kwargs.pop(key, None)
        if value is None:
            continue
        if not isinstance(value, tuple):
            value = (value,)

        new_value = []
        for v in value:
            new_value.append(
                tuple(array._to_own_unit(_v) for _v in v)
                if isinstance(v, tuple) else array._to_own_unit(v))
        kwargs[key] = new_value

    return (array, pad_width, mode), kwargs, array.unit, None


@function_helper
def where(condition, *args):
    from astropy.units import Quantity
    if isinstance(condition, Quantity) or len(args) != 2:
        return NotImplemented

    one, two = args
    if isinstance(one, Quantity):
        return ((condition, one.value, one._to_own_unit(args[1])), {},
                one.unit, None)
    else:
        unit = getattr(one, 'unit', dimensionless_unscaled)
        return (condition, one, two.to_value(unit)), {}, unit, None


def _qp_helper(a, q, *args, q_unit=dimensionless_unscaled, **kwargs):
    if len(args) > 2:
        out = args[1]
        args = args[:1] + args[2:]
    else:
        out = kwargs.pop('out', None)

    from astropy.units import Quantity
    if isinstance(q, Quantity):
        q = q.to_value(q_unit)

    if isinstance(a, Quantity):
        unit = a.unit
        a = a.value
    else:
        unit = getattr(a, 'unit', dimensionless_unscaled)

    if out is not None:
        if isinstance(out, Quantity):
            kwargs['out'] = out.view(np.ndarray)
        else:
            # TODO: for an ndarray output, we could in principle
            # try converting all Quantity to dimensionless.
            return NotImplemented

    return (a, q) + args, kwargs, unit, out


@function_helper
def quantile(a, q, *args, **kwargs):
    return _qp_helper(a, q, *args, **kwargs)


@function_helper
def percentile(a, q, *args, **kwargs):
    from astropy.units import percent
    return _qp_helper(a, q, *args, q_unit=percent, **kwargs)


@function_helper
def count_nonzero(a, *args, **kwargs):
    return (a.value,) + args, kwargs, None, None
