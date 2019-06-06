# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Helpers for overriding numpy functions for Quantity."""

from functools import reduce
import operator

import numpy as np

from astropy.units.core import (
    UnitsError, UnitConversionError, UnitTypeError,
    dimensionless_unscaled, get_current_unit_registry)
from .helpers import _d, get_converter


UNSUPPORTED_FUNCTIONS = {np.packbits, np.unpackbits, np.unravel_index,
                         np.ravel_multi_index, np.ix_, np.cov,
                         np.corrcoef}
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
        raise NotImplementedError


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
        raise NotImplementedError


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
        raise NotImplementedError


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


def _as_quantity(a):
    from astropy.units import Quantity

    try:
        return Quantity(a, copy=False, subok=True)
    except Exception:
        # If we cannot convert to Quantity, we should just bail.
        raise NotImplementedError


def _as_quantities(*args):
    from astropy.units import Quantity

    try:
        return tuple(Quantity(a, copy=False, subok=True)
                     for a in args)
    except Exception:
        # If we cannot convert to Quantity, we should just bail.
        raise NotImplementedError


def _quantities2arrays(*args):
    qs = _as_quantities(*args)
    unit = qs[0].unit
    # Allow any units error to be raised.
    arrays = tuple(q.to_value(unit) for q in qs)
    return arrays, unit


def _iterable_helper(*args, out=None, **kwargs):
    from astropy.units import Quantity

    if out is not None:
        if isinstance(out, Quantity):
            kwargs['out'] = out.view(np.ndarray)
        else:
            # TODO: for an ndarray output, we could in principle
            # try converting all Quantity to dimensionless.
            raise NotImplementedError

    if not args:
        return args, kwargs, None if out is None else out.unit, out

    arrays, unit = _quantities2arrays(*args)
    return arrays, kwargs, unit, out


@function_helper
def concatenate(arrays, axis=0, out=None):
    # TODO: make this smarter by creating an appropriately shaped
    # empty output array and just filling it.
    arrays, kwargs, unit, out = _iterable_helper(*arrays, out=out, axis=axis)
    return (arrays,), kwargs, unit, out


@function_helper
def choose(a, choices, out=None, **kwargs):
    choices, kwargs, unit, out = _iterable_helper(*choices, out=out, **kwargs)
    return (a, choices,), kwargs, unit, out


@function_helper
def select(condlist, choicelist, default=0):
    choicelist, kwargs, unit, out = _iterable_helper(*choicelist)
    if default != 0:
        default = (1 * unit)._to_own_unit(default)
    return (condlist, choicelist, default), kwargs, unit, out


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
        raise NotImplementedError

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

    return (array.view(np.ndarray), pad_width, mode), kwargs, array.unit, None


@function_helper
def where(condition, *args):
    from astropy.units import Quantity
    if isinstance(condition, Quantity) or len(args) != 2:
        raise NotImplementedError

    one, two = args
    if isinstance(one, Quantity):
        return ((condition, one.value, one._to_own_unit(args[1])), {},
                one.unit, None)
    else:
        unit = getattr(one, 'unit', dimensionless_unscaled)
        return (condition, one, two.to_value(unit)), {}, unit, None


@function_helper
def quantile(a, q, *args, q_unit=dimensionless_unscaled, **kwargs):
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
            raise NotImplementedError

    return (a, q) + args, kwargs, unit, out


@function_helper
def percentile(a, q, *args, **kwargs):
    from astropy.units import percent
    return quantile(a, q, *args, q_unit=percent, **kwargs)


FUNCTION_HELPERS[np.nanquantile] = quantile
FUNCTION_HELPERS[np.nanpercentile] = percentile


@function_helper
def count_nonzero(a, *args, **kwargs):
    return (a.value,) + args, kwargs, None, None


@function_helper
def array_equal(a1, a2):
    args, unit = _quantities2arrays(a1, a2)
    return args, {}, None, None


@function_helper
def array_equiv(a1, a2):
    args, unit = _quantities2arrays(a1, a2)
    return args, {}, None, None


def _dot_like(a, b, out=None):
    from astropy.units import Quantity

    a, b = _as_quantities(a, b)
    unit = a.unit * b.unit
    if out is not None:
        if not isinstance(out, Quantity):
            raise NotImplementedError
        return tuple(x.view(np.ndarray) for x in (a, b, out)), {}, unit, out
    else:
        return (a.view(np.ndarray), b.view(np.ndarray)), {}, unit, None


FUNCTION_HELPERS[np.dot] = _dot_like
FUNCTION_HELPERS[np.outer] = _dot_like


def _cross_like(a, b, *args, **kwargs):
    a, b = _as_quantities(a, b)
    unit = a.unit * b.unit
    return (a.view(np.ndarray), b.view(np.ndarray)) + args, kwargs, unit, None


FUNCTION_HELPERS[np.cross] = _cross_like
FUNCTION_HELPERS[np.inner] = _cross_like
FUNCTION_HELPERS[np.vdot] = _cross_like
FUNCTION_HELPERS[np.tensordot] = _cross_like
FUNCTION_HELPERS[np.kron] = _cross_like
FUNCTION_HELPERS[np.correlate] = _cross_like
FUNCTION_HELPERS[np.convolve] = _cross_like


@function_helper
def einsum(subscripts, *operands, out=None, **kwargs):
    from astropy.units import Quantity

    if not isinstance(subscripts, str):
        raise ValueError('only "subscripts" string mode supported for einsum.')

    if out is not None:
        if not isinstance(out, Quantity):
            raise NotImplementedError

        else:
            kwargs['out'] = out.view(np.ndarray)

    qs = _as_quantities(*operands)
    unit = reduce(operator.mul, (q.unit for q in qs),
                  dimensionless_unscaled)
    arrays = tuple(q.view(np.ndarray) for q in qs)
    return (subscripts,) + arrays, kwargs, unit, out


@function_helper
def bincount(x, weights=None, minlength=0):
    from astropy.units import Quantity
    if isinstance(x, Quantity):
        return None
    return (x, weights.value, minlength), {}, weights.unit, None


@function_helper
def digitize(x, bins, *args, **kwargs):
    arrays, unit = _quantities2arrays(x, bins)
    return arrays + args, kwargs, None, None


def _check_bins(bins, unit):
    check = _as_quantity(bins)
    if check.ndim > 0:
        return check.to_value(unit)
    else:
        return bins


@function_helper
def histogram(a, bins=10, range=None, normed=None, weights=None,
              density=None):
    if weights is not None:
        weights = _as_quantity(weights)
        unit = weights.unit
        weights = weights.value
    else:
        unit = None

    a = _as_quantity(a)
    if not isinstance(bins, str):
        bins = _check_bins(bins, a.unit)

    if density or normed:
        unit = (unit or 1) / a.unit

    return ((a.value, bins, range, normed, weights, density), {},
            (unit, a.unit), None)


@function_helper
def histogram_bin_edges(a, bins=10, range=None, weights=None):
    # weights is currently unused
    a = _as_quantity(a)
    if not isinstance(bins, str):
        bins = _check_bins(bins, a.unit)

    return (a.value, bins, range, weights), {}, a.unit, None


@function_helper
def histogram2d(x, y, bins=10, range=None, normed=None, weights=None,
                density=None):

    if weights is not None:
        weights = _as_quantity(weights)
        unit = weights.unit
        weights = weights.value
    else:
        unit = None

    x, y = _as_quantities(x, y)
    if not isinstance(bins, str):
        try:
            n = len(bins)
        except TypeError:
            pass
        else:
            if n == 2:
                bins = tuple(_check_bins(b, unit) for (b, unit) in
                             zip(bins, (x.unit, y.unit)))
            elif n == 1:
                return NotImplementedError
            else:
                bins = _check_bins(bins, x.unit)
                y = y.to(x.unit)

    if density or normed:
        unit = (unit or 1) / x.unit / y.unit

    return ((x.value, y.value, bins, range, normed, weights, density), {},
            (unit, x.unit, y.unit), None)


@function_helper
def diff(a, n=1, axis=-1, prepend=np._NoValue, append=np._NoValue):
    a = _as_quantity(a)
    if prepend is not np._NoValue:
        prepend = _as_quantity(prepend).to_value(a.unit)
    if append is not np._NoValue:
        append = _as_quantity(append).to_value(a.unit)
    return (a.value, n, axis, prepend, append), {}, a.unit, None


@function_helper
def gradient(f, *varargs, **kwargs):
    f = _as_quantity(f)
    axis = kwargs.get('axis', None)
    if axis is None:
        n_axis = f.ndim
    elif isinstance(axis, tuple):
        n_axis = len(axis)
    else:
        n_axis = 1

    if varargs:
        varargs = _as_quantities(*varargs)
        if len(varargs) == 1 and n_axis > 1:
            varargs = varargs * n_axis

    if varargs:
        units = [f.unit / q.unit for q in varargs]
        varargs = tuple(q.value for q in varargs)
    else:
        units = [f.unit] * n_axis

    if len(units) == 1:
        units = units[0]

    return (f.value,) + varargs, kwargs, units, None


@function_helper
def logspace(start, stop, *args, **kwargs):
    from astropy.units import LogQuantity, dex
    if (not isinstance(start, LogQuantity) or
            not isinstance(stop, LogQuantity)):
        raise NotImplementedError

    # Get unit from end point as for linspace.
    stop = stop.to(dex(stop.unit.physical_unit))
    start = start.to(stop.unit)
    unit = stop.unit.physical_unit
    return (start.value, stop.value) + args, kwargs, unit, None


@function_helper
def geomspace(start, stop, *args, **kwargs):
    # Get unit from end point as for linspace.
    (stop, start), unit = _quantities2arrays(stop, start)
    return (start, stop) + args, kwargs, unit, None


@function_helper
def interp(x, xp, fp, *args, **kwargs):
    from astropy.units import Quantity

    (x, xp), _ = _quantities2arrays(x, xp)
    if isinstance(fp, Quantity):
        unit = fp.unit
        fp = fp.value
    else:
        unit = None

    return (x, xp, fp) + args, kwargs, unit, None
