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
import warnings

import numpy as np
from numpy.lib import recfunctions as rfn

from astropy.coordinates import funcs
from astropy.units.errors import UnitConversionError, UnitsError, UnitTypeError
from astropy.utils.compat import (
    COPY_IF_NEEDED,
    NUMPY_LT_2_0,
    NUMPY_LT_2_1,
    NUMPY_LT_2_2,
    NUMPY_LT_2_4,
)
from astropy.utils.exceptions import AstropyUserWarning


from .helpers import (
    _as_quantities,
    _as_quantity,
    _get_np_func_name,
    _quantities2arrays,
)
def _as_quantities(*args):
    """Convert all arguments to quantities."""
    from astropy.units import Quantity
    return tuple(Quantity(arg, copy=COPY_IF_NEEDED) for arg in args)

def _as_quantity(a):
    """Convert an argument to a quantity."""
    from astropy.units import Quantity
    return Quantity(a, copy=COPY_IF_NEEDED)

def _quantities2arrays(*args, unit_from_first=False):
    """Convert all quantities to arrays, ensuring they have the same unit."""
    from astropy.units import Quantity
    if unit_from_first:
        unit = getattr(args[0], "unit", None)
        args = tuple(arg.to_value(unit) if isinstance(arg, Quantity) else arg for arg in args)
        return args, unit
    
    # Standard logic: find a common unit
    arg_units = [getattr(arg, "unit", None) for arg in args if isinstance(arg, Quantity)]
    if not arg_units:
        return args, None
    unit = arg_units[0]
    args = tuple(arg.to_value(unit) if isinstance(arg, Quantity) else arg for arg in args)
    return args, unit

def _interpret_tol(tol, unit):
    """
    Convert a tolerance to a specific unit and return its value.
    Used in linalg functions like matrix_rank and pinv.
    """
    # Using local import of Quantity inside to avoid potential 
    # circularity if this file is imported by units/__init__.py
    from astropy.units import Quantity
    return Quantity(tol, unit).value


if NUMPY_LT_2_0:
    import numpy.core as np_core
else:
    import numpy._core as np_core

SUBCLASS_SAFE_FUNCTIONS = set()
"""Functions with implementations supporting subclasses like Quantity."""
FUNCTION_HELPERS = {}
"""Functions with implementations usable with proper unit conversion."""
DISPATCHED_FUNCTIONS = {}
"""Functions for which we provide our own implementation."""

if NUMPY_LT_2_2:
    # in numpy 2.2 these are auto detected by numpy itself
    # xref https://github.com/numpy/numpy/issues/27451
    SUPPORTED_NEP35_FUNCTIONS = {
        "arange",
        "empty", "ones", "zeros", "full",
        "array", "asarray", "asanyarray", "ascontiguousarray", "asfortranarray",
        "frombuffer", "fromfile", "fromfunction", "fromiter", "fromstring",
        "require", "identity", "eye", "tri", "genfromtxt", "loadtxt",
    }  # fmt: skip
    """Functions that support a 'like' keyword argument and dispatch on it (NEP 35)"""
else:
    # When our minimum becomes numpy>=2.2, this can be removed, here and in the tests
    SUPPORTED_NEP35_FUNCTIONS = set()

"""Functions that support a 'like' keyword argument and dispatch on it (NEP 35)"""
UNSUPPORTED_FUNCTIONS = set()
"""Functions that cannot sensibly be used with quantities."""
SUBCLASS_SAFE_FUNCTIONS |= {
    "shape", "size", "ndim",
    "reshape", "ravel", "moveaxis", "rollaxis", "swapaxes",
    "transpose", "atleast_1d", "atleast_2d", "atleast_3d",
    "expand_dims", "squeeze", "broadcast_to", "broadcast_arrays",
    "flip", "fliplr", "flipud", "rot90",
    "argmin", "argmax", "argsort", "lexsort", "searchsorted",
    "nonzero", "argwhere", "flatnonzero",
    "diag_indices_from", "triu_indices_from", "tril_indices_from",
    "real", "imag", "diagonal", "diagflat", "empty_like",
    "compress", "extract", "delete", "trim_zeros", "roll", "take",
    "put", "fill_diagonal", "tile", "repeat",
    "split", "array_split", "hsplit", "vsplit", "dsplit",
    "stack", "column_stack", "hstack", "vstack", "dstack",
    "max", "min", "amax", "amin", "ptp", "sum", "cumsum",
    "prod", "cumprod",
    "round", "around",
    "fix", "angle", "i0", "clip",
    "isposinf", "isneginf", "isreal", "iscomplex",
    "average", "mean", "std", "var", "trace",
    "nanmax", "nanmin", "nanargmin", "nanargmax", "nanmean",
    "nansum", "nancumsum", "nanprod", "nancumprod",
    "einsum_path", "linspace",
    "sort", "partition", "meshgrid",
    "common_type", "result_type", "can_cast", "min_scalar_type",
    "iscomplexobj", "isrealobj",
    "shares_memory", "may_share_memory",
    "apply_along_axis", "take_along_axis", "put_along_axis",
    "linalg.cond", "linalg.multi_dot",
}  # fmt: skip

SUBCLASS_SAFE_FUNCTIONS |= {"median"}

if NUMPY_LT_2_0:
    # functions (re)moved in numpy 2.0
    SUBCLASS_SAFE_FUNCTIONS |= {
        "msort",
        "round_",  # noqa: NPY003, NPY201
        "trapz",  # noqa: NPY201
        "product",  # noqa: NPY003, NPY201
        "cumproduct",  # noqa: NPY003, NPY201
    }
if not NUMPY_LT_2_0:
    # Array-API compatible versions (matrix axes always at end).
    SUBCLASS_SAFE_FUNCTIONS |= {
        "matrix_transpose", "linalg.matrix_transpose",
        "linalg.diagonal", "linalg.trace",
        "linalg.matrix_norm", "linalg.vector_norm", "linalg.vecdot",
    }  # fmt: skip

    # these work out of the box (and are tested), because they
    # delegate to other, already wrapped functions from the np namespace
    SUBCLASS_SAFE_FUNCTIONS |= {
        "linalg.cross", "linalg.svdvals", "linalg.tensordot", "linalg.matmul",
        "unique_all", "unique_counts", "unique_inverse", "unique_values",
        "astype",
    }  # fmt: skip

    # trapz was renamed to trapezoid
    SUBCLASS_SAFE_FUNCTIONS |= {"trapezoid"}
if not NUMPY_LT_2_1:
    SUBCLASS_SAFE_FUNCTIONS |= {"unstack", "cumulative_prod", "cumulative_sum"}

# Implemented as methods on Quantity:
# np.ediff1d is from setops, but we support it anyway; the others
# currently return NotImplementedError.
# TODO: move latter to UNSUPPORTED? Would raise TypeError instead.
UNSUPPORTED_FUNCTIONS |= {
    "packbits", "unpackbits", "unravel_index",
    "ravel_multi_index", "ix_", "cov", "corrcoef",
    "busday_count", "busday_offset", "datetime_as_string",
    "is_busday", "all", "any",
}  # fmt: skip

if NUMPY_LT_2_0:
    UNSUPPORTED_FUNCTIONS |= {  # removed in numpy 2.0
        "sometrue", "alltrue",  # noqa: NPY003, NPY201
    }  # fmt: skip

# Could be supported if we had a natural logarithm unit.
UNSUPPORTED_FUNCTIONS |= {"linalg.slogdet"}
TBD_FUNCTIONS = {
    "lib.recfunctions.drop_fields", "lib.recfunctions.rename_fields", 
    "lib.recfunctions.append_fields", "lib.recfunctions.join_by",
    "lib.recfunctions.apply_along_fields", "lib.recfunctions.assign_fields_by_name",
    "lib.recfunctions.find_duplicates", "lib.recfunctions.recursive_fill_fields", 
    "lib.recfunctions.require_fields", "lib.recfunctions.repack_fields", 
    "lib.recfunctions.stack_arrays",
}  # fmt: skip
UNSUPPORTED_FUNCTIONS |= TBD_FUNCTIONS
IGNORED_FUNCTIONS = {
    # I/O - useless for Quantity, since no way to store the unit.
    "save", "savez", "savetxt", "savez_compressed",
    # Polynomials
    "poly", "polyadd", "polyder", "polydiv", "polyfit", "polyint",
    "polymul", "polysub", "polyval", "roots", "vander",
    # functions taking record arrays (which are deprecated)
    "lib.recfunctions.rec_append_fields", "lib.recfunctions.rec_drop_fields", "lib.recfunctions.rec_join",
}  # fmt: skip
UNSUPPORTED_FUNCTIONS |= IGNORED_FUNCTIONS


class FunctionAssigner:
    def __init__(self, assignments):
        self.assignments = assignments

    def __call__(self, helps):
        """Add a helper to a numpy function.

        Normally used as a decorator.

        If ``helps`` is given, it should be the numpy function helped (or an
        iterable of numpy functions helped).

        If ``helps`` is not given, it is assumed the function helped is the
        numpy function with the same name as the decorated function.
        """
        def decorator(f):
            # Normalize 'helps' to a collection
            targets = helps if isinstance(helps, (set, list, tuple)) else (helps,)
            
            for h in targets:
                if isinstance(h, str):
                    # Direct string registration (e.g., "linalg.inv")
                    name = h
                else:
                    # Get name from function object (e.g., np.sin -> "sin")
                    name = h.__name__
                    
                    # Minimal logic to prefix submodules if using the function object
                    mod = getattr(h, "__module__", "")
                    if "linalg" in mod:
                        name = f"linalg.{name}"
                    elif "fft" in mod:
                        name = f"fft.{name}"
                    elif "recfunctions" in mod:
                        name = f"lib.recfunctions.{name}"
                
                self.assignments[name] = f
            return f
        return decorator

    function_helper = FunctionAssigner(FUNCTION_HELPERS)

    dispatched_function = FunctionAssigner(DISPATCHED_FUNCTIONS)


    @function_helper(
        helps={
            "copy", "real_if_close", "sort_complex", "resize",
            "fft.fft", "fft.ifft", "fft.rfft", "fft.irfft",
            "fft.fft2", "fft.ifft2", "fft.rfft2", "fft.irfft2",
            "fft.fftn", "fft.ifftn", "fft.rfftn", "fft.irfftn",
            "fft.hfft", "fft.ihfft",
            "linalg.eigvals", "linalg.eigvalsh",
        } | ({"asfarray"} if NUMPY_LT_2_0 else set())  
    )  # fmt: skip
    def invariant_a_helper(a, *args, **kwargs):
        return (a.view(np.ndarray),) + args, kwargs, a.unit, None


    @function_helper(helps={np.tril, np.triu})
    def invariant_m_helper(m, *args, **kwargs):
        return (m.view(np.ndarray),) + args, kwargs, m.unit, None


    @function_helper(helps={np.fft.fftshift, np.fft.ifftshift})
    def invariant_x_helper(x, *args, **kwargs):
        return (x.view(np.ndarray),) + args, kwargs, x.unit, None


    # Note that ones_like does *not* work by default since if one creates an empty
    # array with a unit, one cannot just fill it with unity.  Indeed, in this
    # respect, it is a bit of an odd function for Quantity. On the other hand, it
    # matches the idea that a unit is the same as the quantity with that unit and
    # value of 1. Also, it used to work without __array_function__.
    # zeros_like does work by default for regular quantities, because numpy first
    # creates an empty array with the unit and then fills it with 0 (which can have
    # any unit), but for structured dtype this fails (0 cannot have an arbitrary
    # structured unit), so we include it here too.
    @function_helper(helps={np.ones_like, np.zeros_like})
    def like_helper(a, *args, **kwargs):
        subok = args[2] if len(args) > 2 else kwargs.pop("subok", True)
        unit = a.unit if subok else None
        return (a.view(np.ndarray),) + args, kwargs, unit, None


    def _quantity_out_as_array(out):
        from astropy.units import Quantity

        if isinstance(out, Quantity):
            return out.view(np.ndarray)
        else:
            # TODO: for an ndarray output, one could in principle
            # try converting the input to dimensionless.
            raise NotImplementedError


    # nanvar is safe for Quantity and was previously in SUBCLASS_FUNCTIONS, but it
    # is not safe for Angle, since the resulting unit is inconsistent with being
    # an Angle. By using FUNCTION_HELPERS, the unit gets passed through
    # _result_as_quantity, which will correctly drop to Quantity.
    # A side effect would be that np.nanstd then also produces Quantity; this
    # is avoided by it being helped below.
    @function_helper(helps=np.nanvar)
    def nanvar(a, axis=None, dtype=None, out=None, ddof=0, keepdims=np._NoValue, **kwargs):
        a = _as_quantity(a)
        out_array = None if out is None else _quantity_out_as_array(out)
        return (
            (a.view(np.ndarray), axis, dtype, out_array, ddof, keepdims),
            kwargs,
            a.unit**2,
            out,
        )


    @function_helper(helps=np.nanstd)
    def nanstd(a, axis=None, dtype=None, out=None, ddof=0, keepdims=np._NoValue, **kwargs):
        a = _as_quantity(a)
        out_array = None if out is None else _quantity_out_as_array(out)
        return (
            (a.view(np.ndarray), axis, dtype, out_array, ddof, keepdims),
            kwargs,
            a.unit,
            out,
        )


    @function_helper(helps=np.sinc)
    def sinc(x):
        from astropy.units.si import radian

        try:
            x = x.to_value(radian)
        except UnitsError:
            raise UnitTypeError(
                "Can only apply 'sinc' function to quantities with angle units"
            )
        return (x,), {}, dimensionless_unscaled, None


    @dispatched_function(helps=np.unwrap)
    def unwrap(p, discont=None, axis=-1, *, period=2 * np.pi):
        from astropy.units.si import radian

        if discont is None:
            discont = np.pi << radian

        if period == 2 * np.pi:
            period <<= radian

        p, discont, period = _as_quantities(p, discont, period)
        result = np.unwrap.__wrapped__(
            p.to_value(radian),
            discont.to_value(radian),
            axis=axis,
            period=period.to_value(radian),
        )
        result = radian.to(p.unit, result)
        return result, p.unit, None


    @function_helper(helps=np.argpartition)
    def argpartition(a, *args, **kwargs):
        return (a.view(np.ndarray),) + args, kwargs, None, None


    @function_helper(helps=np.full_like)
    def full_like(a, fill_value, *args, **kwargs):
        unit = a.unit if kwargs.get("subok", True) else None
        return (a.view(np.ndarray), a._to_own_unit(fill_value)) + args, kwargs, unit, None


    def putmask_impl(a, /, mask, values):
        from astropy.units import Quantity

        if isinstance(a, Quantity):
            return (a.view(np.ndarray), mask, a._to_own_unit(values)), {}, a.unit, None
        elif isinstance(values, Quantity):
            return (a, mask, values.to_value(dimensionless_unscaled)), {}, None, None
        else:
            raise NotImplementedError


    if not NUMPY_LT_2_4:

        @function_helper(helps=np.putmask)
        def putmask(a, /, mask, values):
            return putmask_impl(a, mask=mask, values=values)
    else:

        @function_helper(helps=np.putmask)
        def putmask(a, mask, values):
            return putmask_impl(a, mask=mask, values=values)


    @function_helper(helps=np.place)
    def place(arr, mask, vals):
        from astropy.units import Quantity

        if isinstance(arr, Quantity):
            return (arr.view(np.ndarray), mask, arr._to_own_unit(vals)), {}, arr.unit, None
        elif isinstance(vals, Quantity):
            return (arr, mask, vals.to_value(dimensionless_unscaled)), {}, None, None
        else:
            raise NotImplementedError


    @function_helper(helps=np.copyto)
    def copyto(dst, src, *args, **kwargs):
        from astropy.units import Quantity

        if isinstance(dst, Quantity):
            return (dst.view(np.ndarray), dst._to_own_unit(src)) + args, kwargs, None, None
        elif isinstance(src, Quantity):
            return (dst, src.to_value(dimensionless_unscaled)) + args, kwargs, None, None
        else:
            raise NotImplementedError


    @function_helper(helps=np.nan_to_num)
    def nan_to_num(x, copy=True, nan=0.0, posinf=None, neginf=None):
        nan = x._to_own_unit(nan)
        if posinf is not None:
            posinf = x._to_own_unit(posinf)
        if neginf is not None:
            neginf = x._to_own_unit(neginf)
        return (
            (x.view(np.ndarray),),
            dict(copy=True, nan=nan, posinf=posinf, neginf=neginf),
            x.unit,
            None,
        )


def _iterable_helper(*args, out=None, **kwargs):
    """Convert arguments to Quantity, and treat possible 'out'."""
    if out is not None:
        kwargs["out"] = _quantity_out_as_array(out)  # raises if not Quantity.

    arrays, unit = _quantities2arrays(*args)
    return arrays, kwargs, unit, out


if NUMPY_LT_2_4:

    @function_helper(helps=np.concatenate)
    def concatenate(arrays, axis=0, out=None, **kwargs):
        # TODO: make this smarter by creating an appropriately shaped
        # empty output array and just filling it.
        arrays, kwargs, unit, out = _iterable_helper(
            *arrays, out=out, axis=axis, **kwargs
        )
        return (arrays,), kwargs, unit, out
else:

    @function_helper(helps=np.concatenate)
    def concatenate(arrays, /, axis=0, out=None, **kwargs):
        # TODO: make this smarter by creating an appropriately shaped
        # empty output array and just filling it.
        arrays, kwargs, unit, out = _iterable_helper(
            *arrays, out=out, axis=axis, **kwargs
        )
        return (arrays,), kwargs, unit, out


def _block(arrays, max_depth, result_ndim, depth=0):
    # Block by concatenation, copied from np._core.shape_base,
    # but ensuring that we call regular concatenate.
    if depth < max_depth:
        arrs = [_block(arr, max_depth, result_ndim, depth + 1) for arr in arrays]
        # The one difference with the numpy code.
        return np.concatenate(arrs, axis=-(max_depth - depth))
    else:
        return np_core.shape_base._atleast_nd(arrays, result_ndim)


UNIT_FROM_LIKE_ARG = object()


if not NUMPY_LT_2_0: 

    @function_helper(helps=np.arange)
    def arange(
        start_or_stop,
        /,
        stop=None,
        step=1,
        *,
        dtype=None,
        device=None,
    ):
        return arange_impl(
            start_or_stop, stop=stop, step=step, dtype=dtype, device=device
        )

else:

    @function_helper(helps=np.arange)
    def arange(start_or_stop, /, stop=None, step=1, *, dtype=None):
        return arange_impl(start_or_stop, stop=stop, step=step, dtype=dtype)


def unwrap_arange_args(*, start_or_stop, stop_, step_):
    # handle the perilous task of disentangling original arguments
    # This isn't trivial because start_or_stop may actually bind to two
    # different inputs, as the name suggests.
    # We (ab)use structural pattern matching here to bind output variables
    # (start, stop, step), so no additional logic is actually needed after
    # a match is found.
    match (start_or_stop, stop_, step_):
        case (stop, None as start, step):
            pass
        case (start, stop, step):
            pass

    # purely defensive programming
    assert stop is not None, "Please report this."
    return start, stop, step


def wrap_arange_args(*, start, stop, step, expected_out_unit):
    # do the reverse operation than unwrap_arange_args
    # this is needed because start_or_stop *must* be passed as positional

    # purely defensive programming
    assert stop is not None, "Please report this."

    match start, stop:
        case (None, _):
            qty_args = (stop,)
        case _:
            qty_args = (start, stop)

    step_val = step.to_value(expected_out_unit) if hasattr(step, "unit") else step

    kwargs = {} if step == 1 else {"step": step_val}

    # reverse positional arguments so `stop` always comes first
    # this is done to ensure that the arrays are first converted to the
    # expected unit, which we guarantee should be stop's
    args_rev, out_unit = _quantities2arrays(*qty_args[::-1])
    if expected_out_unit is not UNIT_FROM_LIKE_ARG:
        assert out_unit == expected_out_unit
    if hasattr(stop, "unit"):
        assert out_unit == stop.unit

    # reverse args again to restore initial order
    args = args_rev[::-1]

    if "step" in kwargs:
        kwargs["step"] = step_val
    return args, kwargs


def arange_impl(start_or_stop, /, *, stop, step, dtype, device=None):
    # Because this wrapper requires exceptional amounts of additional logic
    # to unwrap/wrap its complicated signature, we'll sprinkle a few
    # sanity checks in the form of `assert` statements, which should help making
    # heads or tails of what's happening in the event of an unexpected exception.
    #
    # also note that we intentionally choose to match numpy.arange's signature
    # at typecheck time, as opposed to its actual runtime signature, which is
    # even richer (as of numpy 2.4). For instance, this means we don't support
    # `start` being passed as keyword, or `dtype` being passed as positional.
    # This is done to improve the overall stability and maintainability of this
    # complicated wrapper function.
    start, stop, step = unwrap_arange_args(
        start_or_stop=start_or_stop, stop_=stop, step_=step
    )
    out_unit = getattr(stop, "unit", UNIT_FROM_LIKE_ARG)

    if out_unit is UNIT_FROM_LIKE_ARG and (
        hasattr(start, "unit") or hasattr(step, "unit")
    ):
        raise TypeError(
            "stop without a unit cannot be combined with start or step with a unit."
        )

    args, kwargs = wrap_arange_args(
        start=start, stop=stop, step=step, expected_out_unit=out_unit
    )

    kwargs["dtype"] = dtype
    if not NUMPY_LT_2_0:
        kwargs["device"] = device

    return args, kwargs, out_unit, None


if NUMPY_LT_2_0:

    @function_helper((helps={np.empty, np.ones, np.zeros})
    def creation_helper(shape, dtype=None, order="C"):
        return (shape, dtype, order), {}, UNIT_FROM_LIKE_ARG, None
else:

    @function_helper(helps={np.empty, np.ones, np.zeros})
    def creation_helper(shape, dtype=None, order="C", *, device=None):
        return (shape, dtype, order), {"device": device}, UNIT_FROM_LIKE_ARG, None


if NUMPY_LT_2_0:

    @function_helper(helps=np.full)
    def full(shape, fill_value, dtype=None, order="C"):
        return full_impl(shape, fill_value, dtype, order)
else:

    @function_helper(helps=np.full)
    def full(shape, fill_value, dtype=None, order="C", *, device=None):
        return full_impl(shape, fill_value, dtype, order, device=device)


def full_impl(shape, fill_value, *args, **kwargs):
    out_unit = getattr(fill_value, "unit", UNIT_FROM_LIKE_ARG)
    if out_unit is not UNIT_FROM_LIKE_ARG:
        fill_value = _as_quantity(fill_value).value
    return (shape, fill_value) + args, kwargs, out_unit, None


@function_helper(helps=np.require)
def require(a, dtype=None, requirements=None):
    out_unit = getattr(a, "unit", UNIT_FROM_LIKE_ARG)
    if out_unit is not UNIT_FROM_LIKE_ARG:
        a = _as_quantity(a).value
    return (a, dtype, requirements), {}, out_unit, None


if not NUMPY_LT_2_4:

    @function_helper(helps=np.array)
    def array(
        object, dtype=None, *, copy=True, order="K", subok=False, ndmin=0, ndmax=0
    ):
        return array_impl(
            object,
            dtype=dtype,
            copy=copy,
            order=order,
            subok=subok,
            ndmin=ndmin,
            ndmax=ndmax,
        )
else:

    @function_helper
    def array(object, dtype=None, *, copy=True, order="K", subok=False, ndmin=0):
        return array_impl(
            object, dtype=dtype, copy=copy, order=order, subok=subok, ndmin=ndmin
        )


def array_impl(object, *, dtype, copy, order, subok, ndmin, ndmax=0):
    out_unit = getattr(object, "unit", UNIT_FROM_LIKE_ARG)
    if out_unit is not UNIT_FROM_LIKE_ARG:
        object = _as_quantity(object).value
    kwargs = {"copy": copy, "order": order, "subok": subok, "ndmin": ndmin}
    if not NUMPY_LT_2_4:
        kwargs |= {"ndmax": ndmax}
    return (object, dtype), kwargs, out_unit, None



if NUMPY_LT_2_0:
    asarray_impl_1_helps = {np.asarray, np.asanyarray}
    asarray_impl_2_helps = {}
elif NUMPY_LT_2_1:
    asarray_impl_1_helps = {np.asanyarray}
    asarray_impl_2_helps = {np.asarray}
else:
    asarray_impl_1_helps = {}
    asarray_impl_2_helps = {np.asarray, np.asanyarray}

@function_helper(helps=asarray_impl_1_helps)
def asarray_impl_1(a, dtype=None, order=None):
    out_unit = getattr(a, "unit", UNIT_FROM_LIKE_ARG)
    if out_unit is not UNIT_FROM_LIKE_ARG:
        a = _as_quantity(a).value
    return (a, dtype, order), {}, out_unit, None


@function_helper(helps=asarray_impl_2_helps)
def asarray_impl_2(a, dtype=None, order=None, *, device=None, copy=None):
    out_unit = getattr(a, "unit", UNIT_FROM_LIKE_ARG)
    if out_unit is not UNIT_FROM_LIKE_ARG:
        a = _as_quantity(a).value
    return (a, dtype, order), {"device": device, "copy": copy}, out_unit, None


@function_helper(helps={"ascontiguousarray", "asfortranarray"})
def aslayoutarray_helper(a, dtype=None):
    out_unit = getattr(a, "unit", UNIT_FROM_LIKE_ARG)
    if out_unit is not UNIT_FROM_LIKE_ARG:
        a = _as_quantity(a).value
    return (a, dtype), {}, out_unit, None


@function_helper(helps=np.fromfunction)
def fromfunction(function, shape, *, dtype=float, **kwargs):
    zero_arg = np.zeros(len(shape), dtype)
    try:
        out_unit = function(*zero_arg).unit
    except Exception:
        out_unit = UNIT_FROM_LIKE_ARG
    return (function, shape), {"dtype": dtype, **kwargs}, out_unit, None


@function_helper(helps={
        np.frombuffer, np.fromfile, np.fromiter, np.fromstring,
        np.identity, np.eye, np.tri,
        np.genfromtxt, np.loadtxt,
    }
)  # fmt: skip
def generic_like_array_function_helper(*args, **kwargs):
    return args, kwargs, UNIT_FROM_LIKE_ARG, None


@dispatched_function(helps=np.block)
def block(arrays):
    # We need to override block since the numpy implementation can take two
    # different paths, one for concatenation, one for creating a large empty
    # result array in which parts are set.  Each assumes array input and
    # cannot be used directly.  Since it would be very costly to inspect all
    # arrays and then turn them back into a nested list, we just copy here the
    # first implementation, np.core.shape_base._block, which is the easiest to
    # adjust while making sure that both units and class are properly kept.
    (arrays, list_ndim, result_ndim, final_size) = np_core.shape_base._block_setup(
        arrays
    )
    result = _block(arrays, list_ndim, result_ndim)
    if list_ndim == 0:
        result = result.copy()
    return result, None, None


@function_helper(helps=np.choose)
def choose(a, choices, out=None, mode="raise"):
    choices, kwargs, unit, out = _iterable_helper(*choices, out=out, mode=mode)
    return (a, choices), kwargs, unit, out


@function_helper(helps=np.select)
def select(condlist, choicelist, default=0):
    choicelist, kwargs, unit, out = _iterable_helper(*choicelist)
    if default != 0:
        default = (1 * unit)._to_own_unit(default)
    return (condlist, choicelist, default), kwargs, unit, out


@dispatched_function(helps=np.piecewise)
def piecewise(x, condlist, funclist, *args, **kw):
    from astropy.units import Quantity

    # Copied implementation from numpy.lib._function_base_impl.piecewise,
    # taking care of units of function outputs.
    n2 = len(funclist)
    # undocumented: single condition is promoted to a list of one condition
    if np.isscalar(condlist) or (
        not isinstance(condlist[0], (list, np.ndarray)) and x.ndim != 0
    ):
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
            f"with {n} condition(s), either {n} or {n + 1} functions are expected"
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


@function_helper(helps=np.append)
def append(arr, values, *args, **kwargs):
    arrays, unit = _quantities2arrays(arr, values, unit_from_first=True)
    return arrays + args, kwargs, unit, None


@function_helper(helps=np.insert)
def insert(arr, obj, values, *args, **kwargs):
    from astropy.units import Quantity

    if isinstance(obj, Quantity):
        raise NotImplementedError

    (arr, values), unit = _quantities2arrays(arr, values, unit_from_first=True)
    return (arr, obj, values) + args, kwargs, unit, None


@function_helper(helps=np.pad)
def pad(array, pad_width, mode="constant", **kwargs):
    # pad dispatches only on array, so that must be a Quantity.
    for key in "constant_values", "end_values":
        value = kwargs.pop(key, None)
        if value is None:
            continue
        if not isinstance(value, tuple):
            value = (value,)

        new_value = []
        for v in value:
            new_value.append(
                tuple(array._to_own_unit(_v) for _v in v)
                if isinstance(v, tuple)
                else array._to_own_unit(v)
            )
        kwargs[key] = new_value

    return (array.view(np.ndarray), pad_width, mode), kwargs, array.unit, None


@function_helper(helps=np.where)
def where(condition, /, *args):
    from astropy.units import Quantity

    if isinstance(condition, Quantity) or len(args) != 2:
        raise NotImplementedError

    args, unit = _quantities2arrays(*args)
    return (condition,) + args, {}, unit, None


@function_helper(helps={np.quantile, np.nanquantile})
def quantile(a, q, *args, _q_unit=None, **kwargs):
    if _q_unit is None:
        from astropy.units import dimensionless_unscaled
        _q_unit = dimensionless_unscaled
    if len(args) >= 2:
        out = args[1]
        args = args[:1] + args[2:]
    else:
        out = kwargs.pop("out", None)

    from astropy.units import Quantity

    if isinstance(q, Quantity):
        q = q.to_value(_q_unit)

    (a,), kwargs, unit, out = _iterable_helper(a, out=out, **kwargs)

    return (a, q) + args, kwargs, unit, out


@function_helper(helps={np.percentile, np.nanpercentile})
def percentile(a, q, *args, **kwargs):
    from astropy.units import percent

    return quantile(a, q, *args, _q_unit=percent, **kwargs)


@function_helper(helps=np.nanmedian)
def nanmedian(a, axis=None, out=None, overwrite_input=False, keepdims=np._NoValue):
    return _iterable_helper(
        a, axis=axis, out=out, overwrite_input=overwrite_input, keepdims=keepdims
    )


@function_helper(helps=np.count_nonzero)    
def count_nonzero(a, *args, **kwargs):
    return (a.value,) + args, kwargs, None, None


@function_helper(helps={np.isclose, np.allclose})
def close(a, b, rtol=1e-05, atol=1e-08, *args, **kwargs):
    from astropy.units import Quantity

    (a, b), unit = _quantities2arrays(a, b, unit_from_first=True)
    # Allow number without a unit as having the unit.
    atol = Quantity(atol, unit).value

    return (a, b, rtol, atol) + args, kwargs, None, None


@dispatched_function(helps=np.array_equal)
def array_equal(a1, a2, equal_nan=False):
    try:
        args, unit = _quantities2arrays(a1, a2)
    except UnitConversionError:
        return False, None, None
    return np.array_equal(*args, equal_nan=equal_nan), None, None


@dispatched_function(helps=np.array_equiv)
def array_equiv(a1, a2):
    try:
        args, unit = _quantities2arrays(a1, a2)
    except UnitConversionError:
        return False, None, None
    return np.array_equiv(*args), None, None


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


@function_helper(
    helps={
        np.cross,
        np.kron,
        np.tensordot,
    }
)
def cross_like_a_b(a, b, *args, **kwargs):
    a, b = _as_quantities(a, b)
    unit = a.unit * b.unit
    return (a.view(np.ndarray), b.view(np.ndarray)) + args, kwargs, unit, None


@function_helper(helps={np.inner, np.vdot})
def cross_like_a_b_posonly(a, b, /):
    a, b = _as_quantities(a, b)
    unit = a.unit * b.unit
    return (a.view(np.ndarray), b.view(np.ndarray)), {}, unit, None


@function_helper(helps={np.correlate, np.convolve})
def cross_like_a_v(a, v, *args, **kwargs):
    a, v = _as_quantities(a, v)
    unit = a.unit * v.unit
    return (a.view(np.ndarray), v.view(np.ndarray)) + args, kwargs, unit, None


@function_helper(helps=np.einsum)
def einsum(*operands, out=None, **kwargs):
    subscripts, *operands = operands

    if not isinstance(subscripts, str):
        raise ValueError('only "subscripts" string mode supported for einsum.')

    if out is not None:
        kwargs["out"] = _quantity_out_as_array(out)

    qs = _as_quantities(*operands)
    unit = functools.reduce(operator.mul, (q.unit for q in qs), dimensionless_unscaled)
    arrays = tuple(q.view(np.ndarray) for q in qs)
    return (subscripts,) + arrays, kwargs, unit, out


@function_helper(helps=np.bincount)
def bincount(x, /, weights=None, minlength=0):
    from astropy.units import Quantity

    if isinstance(x, Quantity):
        raise NotImplementedError
    return (x, weights.value, minlength), {}, weights.unit, None


@function_helper(helps=np.digitize)
def digitize(x, bins, *args, **kwargs):
    arrays, unit = _quantities2arrays(x, bins, unit_from_first=True)
    return arrays + args, kwargs, None, None


def _check_bins(bins, unit):
    from astropy.units import Quantity

    check = _as_quantity(bins)
    if check.ndim > 0:
        return check.to_value(unit)
    elif isinstance(bins, Quantity):
        # bins should be an integer (or at least definitely not a Quantity).
        raise NotImplementedError
    else:
        return bins


def _check_range(range, unit):
    range = _as_quantity(range)
    range = range.to_value(unit)
    return range


@function_helper(helps=np.histogram_bin_edges)
def histogram_bin_edges(a, bins=10, range=None, weights=None):
    # weights is currently unused
    a = _as_quantity(a)
    if not isinstance(bins, str):
        bins = _check_bins(bins, a.unit)

    if range is not None:
        range = _check_range(range, a.unit)

    return (a.value, bins, range, weights), {}, a.unit, None


@function_helper(helps=np.histogram)
def histogram(a, bins=10, range=None, density=None, weights=None):
    if weights is not None:
        weights = _as_quantity(weights)
        unit = weights.unit
        weights = weights.value
    else:
        unit = None

    a = _as_quantity(a)
    if not isinstance(bins, str):
        bins = _check_bins(bins, a.unit)

    if range is not None:
        range = _check_range(range, a.unit)

    if density:
        unit = (unit or 1) / a.unit

    return (
        (a.value, bins, range),
        {"weights": weights, "density": density},
        (unit, a.unit),
        None,
    )


@function_helper(helps=np.histogram2d)
def histogram2d(x, y, bins=10, range=None, density=None, weights=None):
    from astropy.units import Quantity

    if weights is not None:
        weights = _as_quantity(weights)
        unit = weights.unit
        weights = weights.value
    else:
        unit = None

    x, y = _as_quantities(x, y)
    try:
        n = len(bins)
    except TypeError:
        # bins should be an integer (or at least definitely not a Quantity).
        if isinstance(bins, Quantity):
            raise NotImplementedError

    else:
        if n == 1:
            raise NotImplementedError
        elif n == 2 and not isinstance(bins, Quantity):
            bins = [_check_bins(b, unit) for (b, unit) in zip(bins, (x.unit, y.unit))]
        else:
            bins = _check_bins(bins, x.unit)
            y = y.to(x.unit)

    if range is not None:
        range = tuple(
            _check_range(r, unit) for (r, unit) in zip(range, (x.unit, y.unit))
        )

    if density:
        unit = (unit or 1) / x.unit / y.unit

    return (
        (x.value, y.value, bins, range),
        {"weights": weights, "density": density},
        (unit, x.unit, y.unit),
        None,
    )


@function_helper(helps=np.histogramdd)
def histogramdd(sample, bins=10, range=None, density=None, weights=None):
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
                "The dimension of bins must be equal to the dimension of the  sample x."
            )
        bins = [_check_bins(b, unit) for (b, unit) in zip(bins, sample_units)]

    if range is not None:
        range = tuple(_check_range(r, unit) for (r, unit) in zip(range, sample_units))

    if density:
        unit = functools.reduce(operator.truediv, sample_units, (unit or 1))

    return (
        (sample, bins, range),
        {"weights": weights, "density": density},
        (unit, sample_units),
        None,
    )


@function_helper(helps=np.diff)
def diff(a, n=1, axis=-1, prepend=np._NoValue, append=np._NoValue):
    a = _as_quantity(a)
    if prepend is not np._NoValue:
        prepend = _as_quantity(prepend).to_value(a.unit)
    if append is not np._NoValue:
        append = _as_quantity(append).to_value(a.unit)
    return (a.value, n, axis, prepend, append), {}, a.unit, None


@function_helper(helps=np.gradient)
def gradient(f, *varargs, **kwargs):
    f = _as_quantity(f)
    axis = kwargs.get("axis")
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


@function_helper(helps=np.logspace)
def logspace(start, stop, *args, **kwargs):
    from astropy.units import LogQuantity, dex

    if not isinstance(start, LogQuantity) or not isinstance(stop, LogQuantity):
        raise NotImplementedError

    # Get unit from end point as for linspace.
    stop = stop.to(dex(stop.unit.physical_unit))
    start = start.to(stop.unit)
    unit = stop.unit.physical_unit
    return (start.value, stop.value) + args, kwargs, unit, None


@function_helper(helps=np.geomspace)
def geomspace(start, stop, *args, **kwargs):
    # Get unit from end point as for linspace.
    (stop, start), unit = _quantities2arrays(stop, start)
    return (start, stop) + args, kwargs, unit, None


@function_helper(helps=np.interp)
def interp(x, xp, fp, *args, **kwargs):
    from astropy.units import Quantity

    (x, xp), _ = _quantities2arrays(x, xp)
    if isinstance(fp, Quantity):
        unit = fp.unit
        fp = fp.value
    else:
        unit = None

    return (x, xp, fp) + args, kwargs, unit, None


@function_helper(helps=np.unique)
def unique(
    ar,
    return_index=False,
    return_inverse=False,
    return_counts=False,
    axis=None,
    *,
    equal_nan=True,
    **kwargs,
):
    # having **kwargs allows to support sorted (for not NUMPY_LT_2_3) without
    # introducing it pre-maturely in older supported numpy versions
    unit = ar.unit
    n_index = sum(bool(i) for i in (return_index, return_inverse, return_counts))
    if n_index:
        unit = [unit] + n_index * [None]
# Handle version compatibility for equal_nan
    return_kwargs = kwargs.copy()
    return_kwargs["equal_nan"] = equal_nan

    return (
        (ar.value, return_index, return_inverse, return_counts, axis),
        return_kwargs,
        unit,
        None,
    )


@function_helper(helps=np.intersect1d)
def intersect1d(ar1, ar2, assume_unique=False, return_indices=False):
    (ar1, ar2), unit = _quantities2arrays(ar1, ar2)
    if return_indices:
        unit = [unit, None, None]
    return (ar1, ar2, assume_unique, return_indices), {}, unit, None


@function_helper(helps={np.setxor1d, np.union1d, np.setdiff1d})
def twosetop(ar1, ar2, *args, **kwargs):
    (ar1, ar2), unit = _quantities2arrays(ar1, ar2)
    return (ar1, ar2) + args, kwargs, unit, None


@function_helper(helps=np.isin)
def isin(element, test_elements, *args, **kwargs):
    # The unit of 'element' is the reference.
    (ar1, ar2), unit = _quantities2arrays(element, test_elements)
    return (ar1, ar2) + args, kwargs, None, None


if NUMPY_LT_2_4:
    # np.in1d deprecated in not NUMPY_LT_2_0, removed in not NUMPY_LT_24
    @function_helper(helps=np.in1d)
    def in1d(ar1, ar2, *args, **kwargs):
        # This tests whether ar1 is in ar2, so we should change the unit of
        # ar1 to that of ar2.
        (ar1, ar2), unit = _quantities2arrays(ar1, ar2)
        return (ar1, ar2) + args, kwargs, None, None


@dispatched_function(helps=np.apply_over_axes)
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
                raise ValueError(
                    "function is not returning an array of the correct shape"
                )
    # Returning unit is None to signal nothing should happen to
    # the output.
    return val, None, None


@dispatched_function(helps=np.array_repr)
def array_repr(arr, *args, **kwargs):
    # TODO: The addition of "unit='...'" doesn't worry about line
    # length.  Could copy & adapt _array_repr_implementation from
    # numpy.core.arrayprint.py
    cls_name = arr.__class__.__name__
    fake_name = "_" * len(cls_name)
    fake_cls = type(fake_name, (np.ndarray,), {})
    no_unit = np.array_repr(arr.view(fake_cls), *args, **kwargs).replace(
        fake_name, cls_name
    )
    unit_part = f"unit='{arr.unit}'"
    pre, dtype, post = no_unit.rpartition("dtype")
    if dtype:
        return f"{pre}{unit_part}, {dtype}{post}", None, None
    else:
        return f"{no_unit[:-1]}, {unit_part})", None, None


@dispatched_function(helps=np.array_str)
def array_str(a, *args, **kwargs):
    # TODO: The addition of the unit doesn't worry about line length.
    # Could copy & adapt _array_repr_implementation from
    # numpy.core.arrayprint.py
    no_unit = np.array_str(a.value, *args, **kwargs)
    return no_unit + a._unitstr, None, None


@function_helper(helps=np.array2string)
def array2string(a, *args, **kwargs):
    # array2string breaks on quantities as it tries to turn individual
    # items into float, which works only for dimensionless.  Since the
    # defaults would not keep any unit anyway, this is rather pointless -
    # we're better off just passing on the array view.  However, one can
    # also work around this by passing on a formatter (as is done in Angle).
    # So, we do nothing if the formatter argument is present and has the
    # relevant formatter for our dtype.
    if NUMPY_LT_2_4:
        # In NumPy < 2.0, index 6 is 'style' and index 7 is 'formatter'.
        # In NumPy 2.0 to 2.3, index 6 is 'formatter'.
        f_idx = 7 if NUMPY_LT_2_0 else 6
        formatter = args[f_idx] if len(args) > f_idx else kwargs.get("formatter")
    else:
        # In 2.4+, assuming it's strictly keyword-only or following 2.0 logic
        formatter = kwargs.get("formatter")

    if formatter is None:
        a = a.view(np.ndarray)
    else:
        # See whether it covers our dtype.
        if NUMPY_LT_2_0:
            from numpy.core.arrayprint import _get_format_function, _make_options_dict
        else:
            from numpy._core.arrayprint import _get_format_function, _make_options_dict

        with np.printoptions(formatter=formatter) as options:
            options = _make_options_dict(**options)
            try:
                ff = _get_format_function(a.view(np.ndarray), **options)
            except Exception:
                # Shouldn't happen, but possibly we're just not being smart
                # enough, so let's pass things on as is.
                pass
            else:
                # If the selected format function is that of numpy, we know
                # things will fail if we pass in the Quantity, so use .value.
                if "numpy" in ff.__module__:
                    a = a.view(np.ndarray)

    return (a,) + args, kwargs, None, None


@function_helper(helps=np.diag)
def diag(v, *args, **kwargs):
    # Function works for *getting* the diagonal, but not *setting*.
    # So, override always.
    return (v.value,) + args, kwargs, v.unit, None


@function_helper(helps=np.linalg.svd)
def svd(a, full_matrices=True, compute_uv=True, hermitian=False):
    unit = a.unit
    if compute_uv:
        unit = (None, unit, None)

    return ((a.view(np.ndarray), full_matrices, compute_uv, hermitian), {}, unit, None)


def _interpret_tol(tol, unit):
    from astropy.units import Quantity

    return Quantity(tol, unit).value


@function_helper(helps=np.linalg.matrix_rank)
def matrix_rank(A, tol=None, *args, **kwargs):
    from astropy.units import dimensionless_unscaled
    if tol is not None:
        tol = _interpret_tol(tol, A.unit)

    return (A.view(np.ndarray), tol) + args, kwargs, None, None


@function_helper(helps={np.linalg.inv, np.linalg.tensorinv})
def inv(a, *args, **kwargs):
    return (a.view(np.ndarray),) + args, kwargs, 1 / a.unit, None


if NUMPY_LT_2_0:

    @function_helper(helps=np.linalg.pinv)
    def pinv(a, rcond=1e-15, *args, **kwargs):
        from astropy.units import dimensionless_unscaled
        rcond = _interpret_tol(rcond, dimensionless_unscaled)

        return (a.view(np.ndarray), rcond) + args, kwargs, 1 / a.unit, None

else:

    @function_helper(helps=np.linalg.pinv)
    def pinv(a, rcond=None, hermitian=False, *, rtol=np._NoValue):
        from astropy.units import dimensionless_unscaled
        if rcond is not None:
            rcond = _interpret_tol(rcond, dimensionless_unscaled)
        if rtol is not np._NoValue and rtol is not None:
            rtol = _interpret_tol(rtol, dimensionless_unscaled)

        return (
            (a.view(np.ndarray),),
            dict(rcond=rcond, hermitian=hermitian, rtol=rtol),
            1 / a.unit,
            None,
        )


@function_helper(helps=np.linalg.det)
def det(a):
    return (a.view(np.ndarray),), {}, a.unit ** a.shape[-1], None


@function_helper(helps={"linalg.solve", "linalg.tensorsolve"})
def solve(a, b, *args, **kwargs):
    a, b = _as_quantities(a, b)

    return (
        (a.view(np.ndarray), b.view(np.ndarray)) + args,
        kwargs,
        b.unit / a.unit,
        None,
    )


@function_helper(helps=np.linalg.lstsq)
def lstsq(a, b, rcond="warn" if NUMPY_LT_2_0 else None):
    from astropy.units import dimensionless_unscaled
    a, b = _as_quantities(a, b)

    if rcond not in (None, "warn", -1):
        rcond = _interpret_tol(rcond, dimensionless_unscaled)

    return (
        (a.view(np.ndarray), b.view(np.ndarray), rcond),
        {},
        (b.unit / a.unit, b.unit**2, None, a.unit),
        None,
    )


@function_helper(helps=np.linalg.norm)
def norm(x, ord=None, *args, **kwargs):
    if ord == 0:
        from astropy.units import dimensionless_unscaled

        unit = dimensionless_unscaled
    else:
        unit = x.unit
    return (x.view(np.ndarray), ord) + args, kwargs, unit, None


@function_helper(helps=np.linalg.matrix_power)
def matrix_power(a, n):
    return (a.value, n), {}, a.unit**n, None


if NUMPY_LT_2_0:

    @function_helper(helps="linalg.cholesky")
    def cholesky(a):
        return (a.value,), {}, a.unit**0.5, None

else:

    @function_helper(helps=np.linalg.cholesky)
    def cholesky(a, /, *, upper=False):
        return (a.value,), {"upper": upper}, a.unit**0.5, None


@function_helper(helps=np.linalg.qr)
def qr(a, mode="reduced"):
    if mode.startswith("e"):
        units = None
    elif mode == "r":
        units = a.unit
    else:
        from astropy.units import dimensionless_unscaled

        units = (dimensionless_unscaled, a.unit)

    return (a.value, mode), {}, units, None


@function_helper(helps={np.linalg.eig, np.linalg.eigh})
def eig(a, *args, **kwargs):
    from astropy.units import dimensionless_unscaled

    return (a.value,) + args, kwargs, (a.unit, dimensionless_unscaled), None


if not NUMPY_LT_2_0:
    # these functions were added in numpy 2.0

    @function_helper(helps=np.linalg.outer)
    def outer(x1, x2, /):
        # maybe this one can be marked as subclass-safe in the near future ?
        # see https://github.com/numpy/numpy/pull/25101#discussion_r1419879122
        x1, x2 = _as_quantities(x1, x2)
        return (x1.view(np.ndarray), x2.view(np.ndarray)), {}, x1.unit * x2.unit, None


# ======================= np.lib.recfunctions =======================


@function_helper(helps="lib.recfunctions.structured_to_unstructured")
def structured_to_unstructured(arr, *args, **kwargs):
    """
    Convert a structured quantity to an unstructured one.
    This only works if all the units are compatible.

    """
    from astropy.units import StructuredUnit
 
    target_unit = arr.unit.values()[0]

    def replace_unit(x):
        if isinstance(x, StructuredUnit):
            return x._recursively_apply(replace_unit)
        else:
            return target_unit

    to_unit = arr.unit._recursively_apply(replace_unit)
    return (arr.to_value(to_unit),) + args, kwargs, target_unit, None

def _build_structured_unit(dtype, unit):
    """Build structured unit from dtype.

    Parameters
    ----------
    dtype : `numpy.dtype`
    unit : `astropy.units.Unit`

    Returns
    -------
    `astropy.units.Unit` or tuple
    """
    if dtype.fields is None:
        return unit

    return tuple(_build_structured_unit(v[0], unit) for v in dtype.fields.values())


@function_helper(helps="lib.recfunctions.unstructured_to_structured")
def unstructured_to_structured(arr, dtype=None, *args, **kwargs):
    from astropy.units import StructuredUnit

    target_unit = StructuredUnit(_build_structured_unit(dtype, arr.unit))

    return (arr.to_value(arr.unit), dtype) + args, kwargs, target_unit, None


def _izip_units_flat(iterable):
    """Returns an iterator of collapsing any nested unit structure.

    Parameters
    ----------
    iterable : Iterable[StructuredUnit | Unit] or StructuredUnit
        A structured unit or iterable thereof.

    Yields
    ------
    unit
    """
    from astropy.units import StructuredUnit

    # Make Structured unit (pass-through if it is already).
    units = StructuredUnit(iterable)

    # Yield from structured unit.
    for v in units.values():
        if isinstance(v, StructuredUnit):
            yield from _izip_units_flat(v)
        else:
            yield v


@function_helper(helps="lib.recfunctions.merge_arrays")
def merge_arrays(
    seqarrays,
    fill_value=-1,
    flatten=False,
    usemask=False,
    asrecarray=False,
):
    """Merge structured Quantities field by field.

    Like :func:`numpy.lib.recfunctions.merge_arrays`. Note that ``usemask`` and
    ``asrecarray`` are not supported at this time and will raise a ValueError if
    not `False`.
    """
    from astropy.units import Quantity, StructuredUnit

    if asrecarray:
        # TODO? implement if Quantity ever supports rec.array
        raise ValueError("asrecarray=True is not supported.")
    if usemask:
        # TODO: use MaskedQuantity for this case
        raise ValueError("usemask=True is not supported.")

    # Do we have a single Quantity as input?
    if isinstance(seqarrays, Quantity):
        seqarrays = (seqarrays,)

    # Note: this also converts ndarray -> Quantity[dimensionless]
    seqarrays = _as_quantities(*seqarrays)
    arrays = tuple(q.value for q in seqarrays)
    units = tuple(q.unit for q in seqarrays)

    if flatten:
        unit = StructuredUnit(tuple(_izip_units_flat(units)))
    elif len(arrays) == 1:
        unit = StructuredUnit(units[0])
    else:
        unit = StructuredUnit(units)

    return (
        (arrays,),
        dict(
            fill_value=fill_value,
            flatten=flatten,
            usemask=usemask,
            asrecarray=asrecarray,
        ),
        unit,
        None,
    )
