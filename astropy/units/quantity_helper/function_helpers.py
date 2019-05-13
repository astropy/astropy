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
    from astropy.units.quantity import Quantity
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
    from astropy.units.quantity import Quantity
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
    from astropy.units.quantity import Quantity
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


@function_helper
def concatenate(arrays, axis=0, out=None):
    # TODO: make this smarter by creating an appropriately shaped
    # empty output array and just filling it.
    from astropy.units import Quantity
    kwargs = {'axis': axis}

    # Get unit from input if possible.
    for a in arrays:
        if hasattr(a, 'unit'):
            unit = a.unit
            break
    else:
        # No input has even a unit, so the only way we can have gotten
        # here is for output to be a Quantity.  Assume all arrays are
        # are standard ndarray and set unit to dimensionless.
        kwargs['out'] = out.view(np.ndarray)
        return (arrays,), kwargs, dimensionless_unscaled, out

    if out is not None:
        if isinstance(out, Quantity):
            kwargs['out'] = out.view(np.ndarray)
        else:
            # TODO: for an ndarray output, we could in principle
            # try converting all Quantity to dimensionless.
            return NotImplemented

    arrays = tuple(Quantity(a, subok=True, copy=False).to_value(unit)
                   for a in arrays)
    return (arrays,), kwargs, unit, out
