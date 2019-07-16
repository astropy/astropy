# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license. See LICENSE.rst except
# for parts explicitly labelled as being (largely) copies of numpy
# implementations; for those, see licenses/NUMPY_LICENSE.rst.
"""Helpers for overriding numpy functions.

We override numpy functions in `~astropy.units.Quantity.__array_function__`.
In this module, the numpy functions are split in four groups, each of
which has an associated `set` or `dict`:

1. SUBCLASS_SAFE_FUNCTIONS (set), if the numpy implementation
   supports Quantity; we pass on to ndarray.__array_function__.
2. FUNCTION_HELPERS (dict), if the numpy implementation is usable
   after converting quantities to arrays with suitable units,
   and possibly setting units on the result.
3. DISPATCHED_FUNCTIONS (dict), if the function makes sense but
   requires a Quantity-specific implementation
4. UNSUPPORTED_FUNCTIONS (set), if the function does not make sense.

For the FUNCTION_HELPERS `dict`, the value is a function that does the
unit conversion.  It should take the same arguments as the numpy
function would (though one can use ``*args`` and ``**kwargs``) and
return a tuple of ``args, kwargs, unit, out``, where ``args`` and
``kwargs`` will be will be passed on to the numpy implementation,
``unit`` is a possible unit of the result (`None` if it should not be
converted to Quantity), and ``out`` is a possible output Quantity passed
in, which will be filled in-place.

For the DISPATCHED_FUNCTIONS `dict`, the value is a function that
implements the numpy functionality for Quantity input. It should
return a tuple of ``result, unit, out``, where ``result`` is generally
a plain array with the result, and ``unit`` and ``out`` are as above.
If unit is `None`, result gets returned directly, so one can also
return a Quantity directly using ``quantity_result, None, None``.

"""

import functools
import operator

import numpy as np

from astropy.units.core import (
    UnitsError, UnitTypeError, dimensionless_unscaled)
from astropy.utils.compat import NUMPY_LT_1_17, NUMPY_LT_1_15
from astropy.utils import isiterable


if NUMPY_LT_1_17:  # pragma: no cover
    # pre 1.16, overrides.py did not exist; in 1.16, it existed, but
    # __array_function__ overrides were turned off by default.
    ARRAY_FUNCTION_ENABLED = (hasattr(np.core, 'overrides') and
                              np.core.overrides.ENABLE_ARRAY_FUNCTION)
else:
    # In 1.17, overrides are enabled by default, but it is still possible to
    # turn them off using an environment variable.  We use getattr since it
    # is planned to remove that possibility in later numpy versions.
    ARRAY_FUNCTION_ENABLED = getattr(np.core.overrides,
                                     'ENABLE_ARRAY_FUNCTION', True)
SUBCLASS_SAFE_FUNCTIONS = set()
"""Functions with implementations supporting subclasses like Quantity."""
FUNCTION_HELPERS = {}
"""Functions with implementations usable with proper unit conversion."""
DISPATCHED_FUNCTIONS = {}
"""Functions for which we provide our own implementation."""
UNSUPPORTED_FUNCTIONS = set()
"""Functions that cannot sensibly be used with quantities."""

SUBCLASS_SAFE_FUNCTIONS |= {
    np.alen, np.shape, np.size, np.ndim,
    np.reshape, np.ravel, np.moveaxis, np.rollaxis, np.swapaxes,
    np.transpose, np.atleast_1d, np.atleast_2d, np.atleast_3d,
    np.expand_dims, np.squeeze, np.broadcast_to, np.broadcast_arrays,
    np.flip, np.fliplr, np.flipud, np.rot90,
    np.argmin, np.argmax, np.argsort, np.lexsort, np.searchsorted,
    np.nonzero, np.argwhere, np.flatnonzero,
    np.diag_indices_from, np.triu_indices_from, np.tril_indices_from,
    np.real, np.imag, np.diag, np.diagonal, np.diagflat,
    np.empty_like, np.zeros_like,
    np.compress, np.extract, np.delete, np.trim_zeros, np.roll, np.take,
    np.put, np.fill_diagonal, np.tile, np.repeat,
    np.split, np.array_split, np.hsplit, np.vsplit, np.dsplit,
    np.stack, np.column_stack, np.hstack, np.vstack, np.dstack,
    np.amax, np.amin, np.ptp, np.sum, np.cumsum,
    np.prod, np.product, np.cumprod, np.cumproduct,
    np.round, np.around,
    np.fix, np.angle, np.i0, np.clip,
    np.isposinf, np.isneginf, np.isreal, np.iscomplex,
    np.average, np.mean, np.std, np.var, np.median, np.trace,
    np.nanmax, np.nanmin, np.nanargmin, np.nanargmax, np.nanmean,
    np.nanmedian, np.nansum, np.nancumsum, np.nanstd, np.nanvar,
    np.nanprod, np.nancumprod,
    np.einsum_path, np.trapz, np.linspace,
    np.sort, np.msort, np.partition, np.meshgrid,
    np.common_type, np.result_type, np.can_cast, np.min_scalar_type,
    np.iscomplexobj, np.isrealobj,
    np.shares_memory, np.may_share_memory,
    np.apply_along_axis}

if not NUMPY_LT_1_15:
    SUBCLASS_SAFE_FUNCTIONS |= {np.take_along_axis, np.put_along_axis}


# Implemented as methods on Quantity:
# np.ediff1d is from setops, but we support it anyway; the others
# currently return NotImplementedError.
# TODO: move latter to UNSUPPORTED? Would raise TypeError instead.
SUBCLASS_SAFE_FUNCTIONS |= {
    np.ediff1d,
    np.all, np.any, np.sometrue, np.alltrue}

# Subclass safe, but possibly better if overridden (e.g., with different
# default arguments for isclose, allclose).
# TODO: decide on desired behaviour.
SUBCLASS_SAFE_FUNCTIONS |= {
    np.isclose, np.allclose,
    np.array2string, np.array_repr, np.array_str}

# Nonsensical for quantities.
UNSUPPORTED_FUNCTIONS |= {
    np.packbits, np.unpackbits, np.unravel_index,
    np.ravel_multi_index, np.ix_, np.cov, np.corrcoef,
    np.busday_count, np.busday_offset, np.datetime_as_string,
    np.is_busday}

# The following are not just unsupported, but so unlikely to be thought
# to be supported that we ignore them in testing.  (Kept in a separate
# variable so that we can check consistency in the test routine -
# test_quantity_non_ufuncs.py)
IGNORED_FUNCTIONS = {
    # Deprecated
    np.rank, np.asscalar,
    # I/O - useless for Quantity, since no way to store the unit.
    np.save, np.savez, np.savetxt, np.savez_compressed,
    # Polynomials
    np.poly, np.polyadd, np.polyder, np.polydiv, np.polyfit, np.polyint,
    np.polymul, np.polysub, np.polyval, np.roots, np.vander,
    # financial
    np.fv, np.ipmt, np.irr, np.mirr, np.nper, np.npv, np.pmt, np.ppmt,
    np.pv, np.rate}
UNSUPPORTED_FUNCTIONS |= IGNORED_FUNCTIONS


class FunctionAssigner:
    def __init__(self, assignments):
        self.assignments = assignments

    def __call__(self, f=None, helps=None):
        """Add a helper to a numpy function.

        Normally used as a decorator.

        If ``helps`` is given, it should be the numpy function helped (or an
        iterable of numpy functions helped).

        If ``helps`` is not given, it is assumed the function helped is the
        numpy function with the same name as the decorated function.
        """
        if f is not None:
            if helps is None:
                helps = getattr(np, f.__name__)
            if not isiterable(helps):
                helps = (helps,)
            for h in helps:
                self.assignments[h] = f
            return f
        elif helps is not None:
            return functools.partial(self.__call__, helps=helps)
        else:  # pragma: no cover
            raise ValueError("function_helper requires at least one argument.")


function_helper = FunctionAssigner(FUNCTION_HELPERS)

dispatched_function = FunctionAssigner(DISPATCHED_FUNCTIONS)


@function_helper(helps={np.copy, np.asfarray, np.real_if_close,
                        np.sort_complex, np.resize})
def invariant_a_helper(a, *args, **kwargs):
    return (a.view(np.ndarray),) + args, kwargs, a.unit, None


@function_helper(helps={np.tril, np.triu})
def invariant_m_helper(m, *args, **kwargs):
    return (m.view(np.ndarray),) + args, kwargs, m.unit, None


# Note that ones_like does *not* work by default (unlike zeros_like) since if
# one creates an empty array with a unit, one cannot just fill it with unity.
# Indeed, in this respect, it is a bit of an odd function for Quantity. On the
# other hand, it matches the idea that a unit is the same as the quantity with
# that unit and value of 1. Also, it used to work without __array_function__.
@function_helper
def ones_like(a, *args, **kwargs):
    subok = args[2] if len(args) > 2 else kwargs.pop('subok', True)
    unit = a.unit if subok else None
    return (a.view(np.ndarray),) + args, kwargs, unit, None


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

    p, discont = _as_quantities(p, discont)
    result = np.unwrap.__wrapped__(p.to_value(radian),
                                   discont.to_value(radian), axis=axis)
    result = radian.to(p.unit, result)
    return result, p.unit, None


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
    """Convert argument to a Quantity (or raise NotImplementedError)."""
    from astropy.units import Quantity

    try:
        return Quantity(a, copy=False, subok=True)
    except Exception:
        # If we cannot convert to Quantity, we should just bail.
        raise NotImplementedError


def _as_quantities(*args):
    """Convert arguments to Quantity (or raise NotImplentedError)."""
    from astropy.units import Quantity

    try:
        return tuple(Quantity(a, copy=False, subok=True)
                     for a in args)
    except Exception:
        # If we cannot convert to Quantity, we should just bail.
        raise NotImplementedError


def _quantities2arrays(*args):
    """Convert to Quantities with the unit of the first argument."""
    qs = _as_quantities(*args)
    unit = qs[0].unit
    # Allow any units error to be raised.
    arrays = tuple(q.to_value(unit) for q in qs)
    return arrays, unit


def _iterable_helper(*args, out=None, **kwargs):
    """Convert arguments to Quantity, and treat possible 'out'."""
    from astropy.units import Quantity

    if out is not None:
        if isinstance(out, Quantity):
            kwargs['out'] = out.view(np.ndarray)
        else:
            # TODO: for an ndarray output, we could in principle
            # try converting all Quantity to dimensionless.
            raise NotImplementedError

    arrays, unit = _quantities2arrays(*args)
    return arrays, kwargs, unit, out


@function_helper
def concatenate(arrays, axis=0, out=None):
    # TODO: make this smarter by creating an appropriately shaped
    # empty output array and just filling it.
    arrays, kwargs, unit, out = _iterable_helper(*arrays, out=out, axis=axis)
    return (arrays,), kwargs, unit, out


@dispatched_function
def block(arrays):
    # We need to override block since the numpy implementation can take two
    # different paths, one for concatenation, one for creating a large empty
    # result array in which parts are set.  Each assumes array input and
    # cannot be used directly.  Since it would be very costly to inspect all
    # arrays and then turn them back into a nested list, we just copy here the
    # second implementation, np.core.shape_base._block_slicing, since it is
    # shortest and easiest.
    (arrays, list_ndim, result_ndim,
     final_size) = np.core.shape_base._block_setup(arrays)
    shape, slices, arrays = np.core.shape_base._block_info_recursion(
        arrays, list_ndim, result_ndim)
    # Here, one line of difference!
    arrays, unit = _quantities2arrays(*arrays)
    # Back to _block_slicing
    dtype = np.result_type(*[arr.dtype for arr in arrays])
    F_order = all(arr.flags['F_CONTIGUOUS'] for arr in arrays)
    C_order = all(arr.flags['C_CONTIGUOUS'] for arr in arrays)
    order = 'F' if F_order and not C_order else 'C'
    result = np.empty(shape=shape, dtype=dtype, order=order)
    for the_slice, arr in zip(slices, arrays):
        result[(Ellipsis,) + the_slice] = arr
    return result, unit, None


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


@dispatched_function
def piecewise(x, condlist, funclist, *args, **kw):
    from astropy.units import Quantity

    # Copied implementation from numpy.lib.function_base.piecewise,
    # taking care of units of function outputs.
    n2 = len(funclist)
    # undocumented: single condition is promoted to a list of one condition
    if np.isscalar(condlist) or (
            not isinstance(condlist[0], (list, np.ndarray)) and x.ndim != 0):
        condlist = [condlist]

    if any(isinstance(c, Quantity) for c in condlist):
        raise NotImplementedError

    condlist = np.array(condlist, dtype=bool)
    n = len(condlist)

    if n == n2 - 1:  # compute the "otherwise" condition.
        condelse = ~np.any(condlist, axis=0, keepdims=True)
        condlist = np.concatenate([condlist, condelse], axis=0)
        n += 1
    elif n != n2:
        raise ValueError(
            "with {} condition(s), either {} or {} functions are expected"
            .format(n, n, n+1)
        )

    y = np.zeros(x.shape, x.dtype)
    where = []
    what = []
    for k in range(n):
        item = funclist[k]
        if not callable(item):
            where.append(condlist[k])
            what.append(item)
        else:
            vals = x[condlist[k]]
            if vals.size > 0:
                where.append(condlist[k])
                what.append(item(vals, *args, **kw))

    what, unit = _quantities2arrays(*what)
    for item, value in zip(where, what):
        y[item] = value

    return y, unit, None


@function_helper
def append(arr, values, *args, **kwargs):
    arrays, unit = _quantities2arrays(arr, values)
    return arrays + args, kwargs, unit, None


@function_helper
def insert(arr, obj, values, *args, **kwargs):
    from astropy.units import Quantity

    if isinstance(obj, Quantity):
        raise NotImplementedError

    (arr, values), unit = _quantities2arrays(arr, values)
    return (arr, obj, values) + args, kwargs, unit, None


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

    args, unit = _quantities2arrays(*args)
    return (condition,) + args, {}, unit, None


# Quantile was only introduced in numpy 1.15.
@function_helper(helps=({np.quantile, np.nanquantile}
                        if not NUMPY_LT_1_15 else ()))
def quantile(a, q, *args, _q_unit=dimensionless_unscaled, **kwargs):
    if len(args) >= 2:
        out = args[1]
        args = args[:1] + args[2:]
    else:
        out = kwargs.pop('out', None)

    from astropy.units import Quantity
    if isinstance(q, Quantity):
        q = q.to_value(_q_unit)

    (a,), kwargs, unit, out = _iterable_helper(a, out=out, **kwargs)

    return (a, q) + args, kwargs, unit, out


@function_helper(helps={np.percentile, np.nanpercentile})
def percentile(a, q, *args, **kwargs):
    from astropy.units import percent
    return quantile(a, q, *args, _q_unit=percent, **kwargs)


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


@function_helper(helps={np.dot, np.outer})
def dot_like(a, b, out=None):
    from astropy.units import Quantity

    a, b = _as_quantities(a, b)
    unit = a.unit * b.unit
    if out is not None:
        if not isinstance(out, Quantity):
            raise NotImplementedError
        return tuple(x.view(np.ndarray) for x in (a, b, out)), {}, unit, out
    else:
        return (a.view(np.ndarray), b.view(np.ndarray)), {}, unit, None


@function_helper(helps={np.cross, np.inner, np.vdot, np.tensordot, np.kron,
                        np.correlate, np.convolve})
def cross_like(a, b, *args, **kwargs):
    a, b = _as_quantities(a, b)
    unit = a.unit * b.unit
    return (a.view(np.ndarray), b.view(np.ndarray)) + args, kwargs, unit, None


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
    unit = functools.reduce(operator.mul, (q.unit for q in qs),
                            dimensionless_unscaled)
    arrays = tuple(q.view(np.ndarray) for q in qs)
    return (subscripts,) + arrays, kwargs, unit, out


@function_helper
def bincount(x, weights=None, minlength=0):
    from astropy.units import Quantity
    if isinstance(x, Quantity):
        raise NotImplementedError
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


# histogram_bin_edges was only introduced in numpy 1.15.
@function_helper(helps=np.histogram_bin_edges if not NUMPY_LT_1_15 else ())
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
                bins = [_check_bins(b, unit)
                        for (b, unit) in zip(bins, (x.unit, y.unit))]
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
def histogramdd(sample, bins=10, range=None, normed=None, weights=None,
                density=None):
    if weights is not None:
        weights = _as_quantity(weights)
        unit = weights.unit
        weights = weights.value
    else:
        unit = None

    try:
        # Sample is an ND-array.
        _, D = sample.shape
    except (AttributeError, ValueError):
        # Sample is a sequence of 1D arrays.
        sample = _as_quantities(*sample)
        sample_units = [s.unit for s in sample]
        sample = [s.value for s in sample]
        D = len(sample)
    else:
        sample = _as_quantity(sample)
        sample_units = [sample.unit] * D

    try:
        M = len(bins)
    except TypeError:
        # bins should be an integer
        from astropy.units import Quantity

        if isinstance(bins, Quantity):
            raise NotImplementedError
    else:
        if M != D:
            raise ValueError(
                'The dimension of bins must be equal to the dimension of the '
                ' sample x.')
        bins = [_check_bins(b, unit)
                for (b, unit) in zip(bins, sample_units)]

    if density or normed:
        unit = functools.reduce(operator.truediv, sample_units, (unit or 1))

    return ((sample, bins, range, normed, weights, density), {},
            (unit, sample_units), None)


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


@function_helper
def unique(ar, return_index=False, return_inverse=False,
           return_counts=False, axis=None):
    unit = ar.unit
    n_index = sum(bool(i) for i in
                  (return_index, return_inverse, return_counts))
    if n_index:
        unit = [unit] + n_index * [None]

    return (ar.value, return_index, return_inverse, return_counts,
            axis), {}, unit, None


@function_helper
def intersect1d(ar1, ar2, assume_unique=False, return_indices=False):
    (ar1, ar2), unit = _quantities2arrays(ar1, ar2)
    if return_indices:
        unit = [unit, None, None]
    return (ar1, ar2, assume_unique, return_indices), {}, unit, None


@function_helper(helps=(np.setxor1d, np.union1d, np.setdiff1d))
def twosetop(ar1, ar2, *args, **kwargs):
    (ar1, ar2), unit = _quantities2arrays(ar1, ar2)
    return (ar1, ar2) + args, kwargs, unit, None


@function_helper(helps=(np.isin, np.in1d))
def setcheckop(ar1, ar2, *args, **kwargs):
    # This tests whether ar1 is in ar2, so we should change the unit of
    # a1 to that of a2.
    (ar2, ar1), unit = _quantities2arrays(ar2, ar1)
    return (ar1, ar2) + args, kwargs, None, None


@dispatched_function
def apply_over_axes(func, a, axes):
    # Copied straight from numpy/lib/shape_base, just to omit its
    # val = asarray(a); if only it had been asanyarray, or just not there
    # since a is assumed to an an array in the next line...
    # Which is what we do here - we can only get here if it is a Quantity.
    val = a
    N = a.ndim
    if np.array(axes).ndim == 0:
        axes = (axes,)
    for axis in axes:
        if axis < 0:
            axis = N + axis
        args = (val, axis)
        res = func(*args)
        if res.ndim == val.ndim:
            val = res
        else:
            res = np.expand_dims(res, axis)
            if res.ndim == val.ndim:
                val = res
            else:
                raise ValueError("function is not returning "
                                 "an array of the correct shape")
    # Returning unit is None to signal nothing should happen to
    # the output.
    return val, None, None
