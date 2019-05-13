# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Helpers for overriding numpy functions for Quantity."""

import numpy as np

from astropy.units.core import (
    UnitsError, UnitConversionError, UnitTypeError,
    dimensionless_unscaled, get_current_unit_registry)
from .helpers import _d, get_converter


INVARIANT_FUNCTIONS = {
        np.copy, np.asfarray, np.empty_like, np.zeros_like,
        np.real_if_close, np.tril, np.triu,
        np.sort_complex, np.broadcast_to}
UNSUPPORTED_FUNCTIONS = set()
FUNCTION_HELPERS = {}
DISPATCHED_FUNCTIONS = {}


def function_helper(f):
    FUNCTION_HELPERS[getattr(np, f.__name__.replace('_helper', ''))] = f
    return f


def dispatched_function(f):
    DISPATCHED_FUNCTIONS[getattr(np, f.__name__)] = f
    return f


@dispatched_function
def broadcast_arrays(*args, subok=True):
    return np.broadcast_arrays.__wrapped__(*args, subok=subok)


@function_helper
def sinc_helper(x):
    from astropy.units.si import radian
    try:
        x = x << radian
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
def argpartition(a, kth, **kwargs):
    return (a.value, kth), kwargs, None, None


@function_helper
def full_like(a, fill_value, **kwargs):
    unit = a.unit if kwargs.get('subok', True) else None
    return (a.value, a._to_own_unit(fill_value)), kwargs, unit, None


@function_helper
def putmask(a, mask, values):
    from astropy.units.quantity import Quantity
    if isinstance(a, Quantity):
        return (a.value, mask, a._to_own_unit(values)), {}, a.unit, None
        a = a.value
    elif isinstance(values, Quantity):
        return (a, mask, values.to_value(dimensionless_unscaled)), {}, None, None
    else:
        return NotImplemented
