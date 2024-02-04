# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Helpers for letting numpy functions interact with Masked arrays.

The module supplies helper routines for numpy functions that propagate
masks appropriately., for use in the ``__array_function__``
implementation of `~astropy.utils.masked.MaskedNDArray`.  They are not
very useful on their own, but the ones with docstrings are included in
the documentation so that there is a place to find out how the mask is
interpreted.

"""
import numpy as np

from astropy.units.quantity_helper.function_helpers import FunctionAssigner
from astropy.utils.compat import NUMPY_LT_1_23, NUMPY_LT_1_24, NUMPY_LT_2_0

if NUMPY_LT_2_0:
    import numpy.core as np_core
else:
    import numpy._core as np_core

# This module should not really be imported, but we define __all__
# such that sphinx can typeset the functions with docstrings.
# The latter are added to __all__ at the end.
__all__ = [
    "MASKED_SAFE_FUNCTIONS",
    "APPLY_TO_BOTH_FUNCTIONS",
    "DISPATCHED_FUNCTIONS",
    "UNSUPPORTED_FUNCTIONS",
]


MASKED_SAFE_FUNCTIONS = set()
"""Set of functions that work fine on Masked classes already.

Most of these internally use `numpy.ufunc` or other functions that
are already covered.
"""

APPLY_TO_BOTH_FUNCTIONS = {}
"""Dict of functions that should apply to both data and mask.

The `dict` is keyed by the numpy function and the values are functions
that take the input arguments of the numpy function and organize these
for passing the data and mask to the numpy function.

Returns
-------
data_args : tuple
    Arguments to pass on to the numpy function for the unmasked data.
mask_args : tuple
    Arguments to pass on to the numpy function for the masked data.
kwargs : dict
    Keyword arguments to pass on for both unmasked data and mask.
out : `~astropy.utils.masked.Masked` instance or None
    Optional instance in which to store the output.

Raises
------
NotImplementedError
   When an arguments is masked when it should not be or vice versa.
"""

DISPATCHED_FUNCTIONS = {}
"""Dict of functions that provide the numpy function's functionality.

These are for more complicated versions where the numpy function itself
cannot easily be used.  It should return either the result of the
function, or a tuple consisting of the unmasked result, the mask for the
result and a possible output instance.

It should raise `NotImplementedError` if one of the arguments is masked
when it should not be or vice versa.
"""

UNSUPPORTED_FUNCTIONS = set()
"""Set of numpy functions that are not supported for masked arrays.

For most, masked input simply makes no sense, but for others it may have
been lack of time.  Issues or PRs for support for functions are welcome.
"""

# Almost all from np.core.fromnumeric defer to methods so are OK.
MASKED_SAFE_FUNCTIONS |= {
    getattr(np, name)
    for name in np_core.fromnumeric.__all__
    if name not in {"choose", "put", "resize", "searchsorted", "where", "alen", "ptp"}
}
MASKED_SAFE_FUNCTIONS |= {
    # built-in from multiarray
    np.may_share_memory, np.can_cast, np.min_scalar_type, np.result_type,
    np.shares_memory,
    # np.core.arrayprint
    np.array_repr,
    # np.core.function_base
    np.linspace, np.logspace, np.geomspace,
    # np.core.numeric
    np.isclose, np.allclose, np.flatnonzero, np.argwhere,
    # np.core.shape_base
    np.atleast_1d, np.atleast_2d, np.atleast_3d, np.stack, np.hstack, np.vstack,
    # np.lib._function_base_impl
    np.average, np.diff, np.extract, np.meshgrid, np.gradient,
    # np.lib.index_tricks
    np.diag_indices_from, np.triu_indices_from, np.tril_indices_from,
    np.fill_diagonal,
    # np.lib.shape_base
    np.column_stack, np.dstack,
    np.array_split, np.split, np.hsplit, np.vsplit, np.dsplit,
    np.expand_dims, np.apply_along_axis, np.kron, np.tile,
    np.take_along_axis, np.put_along_axis,
    # np.lib.type_check (all but nan_to_num)
    np.iscomplexobj, np.isrealobj, np.imag, np.isreal, np.real,
    np.real_if_close, np.common_type,
    # np.lib.ufunclike
    np.fix, np.isneginf, np.isposinf,
    # np.lib._function_base_impl
    np.angle, np.i0,
}  # fmt: skip


if NUMPY_LT_2_0:
    # Safe in < 2.0, because it deferred to the method. Overridden in >= 2.0.
    MASKED_SAFE_FUNCTIONS |= {np.ptp}
    # Removed in numpy 2.0.  Just an alias to vstack.
    MASKED_SAFE_FUNCTIONS |= {np.row_stack}
    # renamed in numpy 2.0
    MASKED_SAFE_FUNCTIONS |= {np.trapz}
else:
    # new in numpy 2.0
    MASKED_SAFE_FUNCTIONS |= {np.astype, np.trapezoid}

IGNORED_FUNCTIONS = {
    # I/O - useless for Masked, since no way to store the mask.
    np.save, np.savez, np.savetxt, np.savez_compressed,
    # Polynomials
    np.poly, np.polyadd, np.polyder, np.polydiv, np.polyfit, np.polyint,
    np.polymul, np.polysub, np.polyval, np.roots, np.vander,
}  # fmt: skip
IGNORED_FUNCTIONS |= {
    np.pad, np.searchsorted, np.digitize,
    np.is_busday, np.busday_count, np.busday_offset,
    # numpy.lib._function_base_impl
    np.cov, np.corrcoef, np.trim_zeros,
    # numpy.core.numeric
    np.correlate, np.convolve,
    # numpy.lib.histograms
    np.histogram, np.histogram2d, np.histogramdd, np.histogram_bin_edges,
    # TODO!!
    np.dot, np.vdot, np.inner, np.tensordot, np.cross,
    np.einsum, np.einsum_path,
}  # fmt: skip

# Really should do these...
if NUMPY_LT_2_0:
    from numpy.lib import arraysetops
else:
    # Public set operations have been moved to the top-level namespace in numpy 2.0
    # (numpy/numpy#24507), raising an AttributeError when accessed through np.lib.arraysetops.
    from numpy.lib import _arraysetops_impl as arraysetops

IGNORED_FUNCTIONS |= {getattr(np, setopsname) for setopsname in arraysetops.__all__}

if NUMPY_LT_1_23:
    IGNORED_FUNCTIONS |= {
        # Deprecated, removed in numpy 1.23
        np.asscalar,
        np.alen,
    }

# Explicitly unsupported functions
UNSUPPORTED_FUNCTIONS |= {
    np.unravel_index,
    np.ravel_multi_index,
    np.ix_,
}

# No support for the functions also not supported by Quantity
# (io, polynomial, etc.).
UNSUPPORTED_FUNCTIONS |= IGNORED_FUNCTIONS


apply_to_both = FunctionAssigner(APPLY_TO_BOTH_FUNCTIONS)
dispatched_function = FunctionAssigner(DISPATCHED_FUNCTIONS)


def _get_data_and_masks(*args):
    """Separate out arguments into tuples of data and masks.

    An all-False mask is created if an argument does not have a mask.
    """
    from .core import Masked

    data, masks = Masked._get_data_and_masks(*args)
    masks = tuple(
        m if m is not None else np.zeros(np.shape(d), bool) for d, m in zip(data, masks)
    )
    return data, masks


# Following are simple ufunc-like functions which should just copy the mask.
@dispatched_function
def datetime_as_string(arr, *args, **kwargs):
    return (np.datetime_as_string(arr.unmasked, *args, **kwargs), arr.mask.copy(), None)


@dispatched_function
def sinc(x):
    return np.sinc(x.unmasked), x.mask.copy(), None


@dispatched_function
def iscomplex(x):
    return np.iscomplex(x.unmasked), x.mask.copy(), None


@dispatched_function
def unwrap(p, *args, **kwargs):
    return np.unwrap(p.unmasked, *args, **kwargs), p.mask.copy(), None


@dispatched_function
def nan_to_num(x, copy=True, nan=0.0, posinf=None, neginf=None):
    data = np.nan_to_num(x.unmasked, copy=copy, nan=nan, posinf=posinf, neginf=neginf)
    if copy:
        return (data, x.mask.copy(), None)
    else:
        return (x, None, None)


# Following are simple functions related to shapes, where the same function
# should be applied to the data and the mask.  They cannot all share the
# same helper, because the first arguments have different names.
@apply_to_both(
    helps=(
        {np.copy, np.resize, np.moveaxis, np.rollaxis, np.roll}
        | ({np.asfarray} if NUMPY_LT_2_0 else set())
    )
)
def masked_a_helper(a, *args, **kwargs):
    data, mask = _get_data_and_masks(a)
    return data + args, mask + args, kwargs, None


@apply_to_both(helps={np.flip, np.flipud, np.fliplr, np.rot90, np.triu, np.tril})
def masked_m_helper(m, *args, **kwargs):
    data, mask = _get_data_and_masks(m)
    return data + args, mask + args, kwargs, None


@apply_to_both(helps={np.diag, np.diagflat})
def masked_v_helper(v, *args, **kwargs):
    data, mask = _get_data_and_masks(v)
    return data + args, mask + args, kwargs, None


@apply_to_both(helps={np.delete})
def masked_arr_helper(array, *args, **kwargs):
    data, mask = _get_data_and_masks(array)
    return data + args, mask + args, kwargs, None


@apply_to_both
def broadcast_to(array, shape, subok=False):
    """Broadcast array to the given shape.

    Like `numpy.broadcast_to`, and applied to both unmasked data and mask.
    Note that ``subok`` is taken to mean whether or not subclasses of
    the unmasked data and mask are allowed, i.e., for ``subok=False``,
    a `~astropy.utils.masked.MaskedNDArray` will be returned.
    """
    data, mask = _get_data_and_masks(array)
    return data, mask, dict(shape=shape, subok=subok), None


@dispatched_function
def outer(a, b, out=None):
    return np.multiply.outer(np.ravel(a), np.ravel(b), out=out), None, None


@dispatched_function
def empty_like(prototype, dtype=None, order="K", subok=True, shape=None):
    """Return a new array with the same shape and type as a given array.

    Like `numpy.empty_like`, but will add an empty mask.
    """
    unmasked = np.empty_like(
        prototype.unmasked, dtype=dtype, order=order, subok=subok, shape=shape
    )
    if dtype is not None:
        dtype = (
            np.ma.make_mask_descr(unmasked.dtype)
            if unmasked.dtype.names
            else np.dtype("?")
        )
    mask = np.empty_like(
        prototype.mask, dtype=dtype, order=order, subok=subok, shape=shape
    )

    return unmasked, mask, None


@dispatched_function
def zeros_like(a, dtype=None, order="K", subok=True, shape=None):
    """Return an array of zeros with the same shape and type as a given array.

    Like `numpy.zeros_like`, but will add an all-false mask.
    """
    unmasked = np.zeros_like(
        a.unmasked, dtype=dtype, order=order, subok=subok, shape=shape
    )
    return unmasked, False, None


@dispatched_function
def ones_like(a, dtype=None, order="K", subok=True, shape=None):
    """Return an array of ones with the same shape and type as a given array.

    Like `numpy.ones_like`, but will add an all-false mask.
    """
    unmasked = np.ones_like(
        a.unmasked, dtype=dtype, order=order, subok=subok, shape=shape
    )
    return unmasked, False, None


@dispatched_function
def full_like(a, fill_value, dtype=None, order="K", subok=True, shape=None):
    """Return a full array with the same shape and type as a given array.

    Like `numpy.full_like`, but with a mask that is also set.
    If ``fill_value`` is `numpy.ma.masked`, the data will be left unset
    (i.e., as created by `numpy.empty_like`).
    """
    result = np.empty_like(a, dtype=dtype, order=order, subok=subok, shape=shape)
    result[...] = fill_value
    return result, None, None


@dispatched_function
def put(a, ind, v, mode="raise"):
    """Replaces specified elements of an array with given values.

    Like `numpy.put`, but for masked array ``a`` and possibly masked
    value ``v``.  Masked indices ``ind`` are not supported.
    """
    from astropy.utils.masked import Masked

    if isinstance(ind, Masked) or not isinstance(a, Masked):
        raise NotImplementedError

    v_data, v_mask = a._get_data_and_mask(v)
    if v_data is not None:
        np.put(a.unmasked, ind, v_data, mode=mode)
    # v_mask of None will be correctly interpreted as False.
    np.put(a.mask, ind, v_mask, mode=mode)


@dispatched_function
def putmask(a, mask, values):
    """Changes elements of an array based on conditional and input values.

    Like `numpy.putmask`, but for masked array ``a`` and possibly masked
    ``values``.  Masked ``mask`` is not supported.
    """
    from astropy.utils.masked import Masked

    if isinstance(mask, Masked) or not isinstance(a, Masked):
        raise NotImplementedError

    values_data, values_mask = a._get_data_and_mask(values)
    if values_data is not None:
        np.putmask(a.unmasked, mask, values_data)
    np.putmask(a.mask, mask, values_mask)


@dispatched_function
def place(arr, mask, vals):
    """Change elements of an array based on conditional and input values.

    Like `numpy.place`, but for masked array ``a`` and possibly masked
    ``values``.  Masked ``mask`` is not supported.
    """
    from astropy.utils.masked import Masked

    if isinstance(mask, Masked) or not isinstance(arr, Masked):
        raise NotImplementedError

    vals_data, vals_mask = arr._get_data_and_mask(vals)
    if vals_data is not None:
        np.place(arr.unmasked, mask, vals_data)
    np.place(arr.mask, mask, vals_mask)


@dispatched_function
def copyto(dst, src, casting="same_kind", where=True):
    """Copies values from one array to another, broadcasting as necessary.

    Like `numpy.copyto`, but for masked destination ``dst`` and possibly
    masked source ``src``.
    """
    from astropy.utils.masked import Masked

    if not isinstance(dst, Masked) or isinstance(where, Masked):
        raise NotImplementedError

    src_data, src_mask = dst._get_data_and_mask(src)

    if src_data is not None:
        np.copyto(dst.unmasked, src_data, casting=casting, where=where)
    if src_mask is not None:
        np.copyto(dst.mask, src_mask, where=where)


@dispatched_function
def packbits(a, *args, **kwargs):
    result = np.packbits(a.unmasked, *args, **kwargs)
    mask = np.packbits(a.mask, *args, **kwargs).astype(bool)
    return result, mask, None


@dispatched_function
def unpackbits(a, *args, **kwargs):
    result = np.unpackbits(a.unmasked, *args, **kwargs)
    mask = np.zeros(a.shape, dtype="u1")
    mask[a.mask] = 255
    mask = np.unpackbits(mask, *args, **kwargs).astype(bool)
    return result, mask, None


@dispatched_function
def bincount(x, weights=None, minlength=0):
    """Count number of occurrences of each value in array of non-negative ints.

    Like `numpy.bincount`, but masked entries in ``x`` will be skipped.
    Any masked entries in ``weights`` will lead the corresponding bin to
    be masked.
    """
    from astropy.utils.masked import Masked

    if weights is not None:
        weights = np.asanyarray(weights)
    if isinstance(x, Masked) and x.ndim <= 1:
        # let other dimensions lead to errors.
        if weights is not None and weights.ndim == x.ndim:
            weights = weights[~x.mask]
        x = x.unmasked[~x.mask]
    mask = None
    if weights is not None:
        weights, w_mask = Masked._get_data_and_mask(weights)
        if w_mask is not None:
            mask = np.bincount(x, w_mask.astype(int), minlength=minlength).astype(bool)
    result = np.bincount(x, weights, minlength=0)
    return result, mask, None


if NUMPY_LT_2_0:

    @dispatched_function
    def msort(a):
        result = a.copy()
        result.sort(axis=0)
        return result, None, None

else:
    # Used to work via ptp method, but now need to override, otherwise
    # plain reduction is used, which gives different mask.
    @dispatched_function
    def ptp(a, axis=None, out=None, keepdims=False):
        result = a.max(axis=axis, out=out, keepdims=keepdims)
        result -= a.min(axis=axis, keepdims=keepdims)
        return result, None, None


@dispatched_function
def sort_complex(a):
    # Just a copy of np.lib._function_base_impl.sort_complex, to avoid the asarray.
    b = a.copy()
    b.sort()
    if not issubclass(b.dtype.type, np.complexfloating):  # pragma: no cover
        if b.dtype.char in "bhBH":
            result = b.astype("F")
        elif b.dtype.char == "g":
            result = b.astype("G")
        else:
            result = b.astype("D")
    else:
        result = b

    return result, None, None


@dispatched_function
def concatenate(arrays, axis=0, out=None, dtype=None, casting="same_kind"):
    data, masks = _get_data_and_masks(*arrays)
    if out is None:
        return (
            np.concatenate(data, axis=axis, dtype=dtype, casting=casting),
            np.concatenate(masks, axis=axis),
            None,
        )
    else:
        from astropy.utils.masked import Masked

        if not isinstance(out, Masked):
            raise NotImplementedError
        np.concatenate(masks, out=out.mask, axis=axis)
        np.concatenate(data, out=out.unmasked, axis=axis, dtype=dtype, casting=casting)
        return out, None, None


@apply_to_both
def append(arr, values, axis=None):
    data, masks = _get_data_and_masks(arr, values)
    return data, masks, dict(axis=axis), None


@dispatched_function
def block(arrays):
    # We need to override block since the numpy implementation can take two
    # different paths, one for concatenation, one for creating a large empty
    # result array in which parts are set.  Each assumes array input and
    # cannot be used directly.  Since it would be very costly to inspect all
    # arrays and then turn them back into a nested list, we just copy here the
    # second implementation, np.core.shape_base._block_slicing, since it is
    # shortest and easiest.
    from astropy.utils.masked import Masked

    arrays, list_ndim, result_ndim, final_size = np_core.shape_base._block_setup(arrays)
    shape, slices, arrays = np_core.shape_base._block_info_recursion(
        arrays, list_ndim, result_ndim
    )
    dtype = np.result_type(*[arr.dtype for arr in arrays])
    F_order = all(arr.flags["F_CONTIGUOUS"] for arr in arrays)
    C_order = all(arr.flags["C_CONTIGUOUS"] for arr in arrays)
    order = "F" if F_order and not C_order else "C"
    result = Masked(np.empty(shape=shape, dtype=dtype, order=order))
    for the_slice, arr in zip(slices, arrays):
        result[(Ellipsis,) + the_slice] = arr
    return result, None, None


@dispatched_function
def broadcast_arrays(*args, subok=True):
    """Broadcast arrays to a common shape.

    Like `numpy.broadcast_arrays`, applied to both unmasked data and masks.
    Note that ``subok`` is taken to mean whether or not subclasses of
    the unmasked data and masks are allowed, i.e., for ``subok=False``,
    `~astropy.utils.masked.MaskedNDArray` instances will be returned.
    """
    from .core import Masked

    are_masked = [isinstance(arg, Masked) for arg in args]
    data = [
        (arg.unmasked if is_masked else arg) for arg, is_masked in zip(args, are_masked)
    ]
    results = np.broadcast_arrays(*data, subok=subok)

    return_type = list if NUMPY_LT_2_0 else tuple
    shape = results[0].shape if isinstance(results, return_type) else results.shape
    masks = [
        (np.broadcast_to(arg.mask, shape, subok=subok) if is_masked else None)
        for arg, is_masked in zip(args, are_masked)
    ]
    results = return_type(
        (Masked(result, mask) if mask is not None else result)
        for (result, mask) in zip(results, masks)
    )
    return (results if len(results) > 1 else results[0]), None, None


@apply_to_both
def insert(arr, obj, values, axis=None):
    """Insert values along the given axis before the given indices.

    Like `numpy.insert` but for possibly masked ``arr`` and ``values``.
    Masked ``obj`` is not supported.
    """
    from astropy.utils.masked import Masked

    if isinstance(obj, Masked) or not isinstance(arr, Masked):
        raise NotImplementedError

    (arr_data, val_data), (arr_mask, val_mask) = _get_data_and_masks(arr, values)
    return ((arr_data, obj, val_data, axis), (arr_mask, obj, val_mask, axis), {}, None)


@dispatched_function
def count_nonzero(a, axis=None, *, keepdims=False):
    """Counts the number of non-zero values in the array ``a``.

    Like `numpy.count_nonzero`, with masked values counted as 0 or `False`.
    """
    filled = a.filled(np.zeros((), a.dtype))
    return np.count_nonzero(filled, axis, keepdims=keepdims), None, None


def _masked_median_1d(a, overwrite_input):
    # TODO: need an in-place mask-sorting option.
    unmasked = a.unmasked[~a.mask]
    if unmasked.size:
        return a.from_unmasked(np.median(unmasked, overwrite_input=overwrite_input))
    else:
        return a.from_unmasked(np.zeros_like(a.unmasked, shape=(1,))[0], mask=True)


def _masked_median(a, axis=None, out=None, overwrite_input=False):
    # As for np.nanmedian, but without a fast option as yet.
    if axis is None or a.ndim == 1:
        part = a.ravel()
        result = _masked_median_1d(part, overwrite_input)
    else:
        result = np.apply_along_axis(_masked_median_1d, axis, a, overwrite_input)
    if out is not None:
        out[...] = result
    return result


@dispatched_function
def median(a, axis=None, out=None, **kwargs):
    from astropy.utils.masked import Masked

    if out is not None and not isinstance(out, Masked):
        raise NotImplementedError

    a = Masked(a)

    if NUMPY_LT_1_24:
        keepdims = kwargs.pop("keepdims", False)
        r, k = np.lib.function_base._ureduce(
            a, func=_masked_median, axis=axis, out=out, **kwargs
        )
        result = (r.reshape(k) if keepdims else r) if out is None else out

    elif NUMPY_LT_2_0:
        result = np.lib.function_base._ureduce(
            a, func=_masked_median, axis=axis, out=out, **kwargs
        )

    else:
        result = np.lib._function_base_impl._ureduce(
            a, func=_masked_median, axis=axis, out=out, **kwargs
        )
    return result, None, None


def _masked_quantile_1d(a, q, **kwargs):
    """
    Private function for rank 1 arrays. Compute quantile ignoring NaNs.
    See nanpercentile for parameter usage.
    """
    unmasked = a.unmasked[~a.mask]
    if unmasked.size:
        if NUMPY_LT_2_0:
            result = np.lib.function_base._quantile_unchecked(unmasked, q, **kwargs)
        else:
            result = np.lib._function_base_impl._quantile_unchecked(
                unmasked, q, **kwargs
            )
        return a.from_unmasked(result)
    else:
        return a.from_unmasked(np.zeros_like(a.unmasked, shape=q.shape), True)


def _masked_quantile(a, q, axis=None, out=None, **kwargs):
    # As for np.nanmedian, but without a fast option as yet.
    if axis is None or a.ndim == 1:
        part = a.ravel()
        result = _masked_quantile_1d(part, q, **kwargs)
    else:
        result = np.apply_along_axis(_masked_quantile_1d, axis, a, q, **kwargs)
        # apply_along_axis fills in collapsed axis with results.
        # Move that axis to the beginning to match percentile's
        # convention.
        if q.ndim != 0:
            result = np.moveaxis(result, axis, 0)

    if out is not None:
        out[...] = result
    return result


@dispatched_function
def quantile(a, q, axis=None, out=None, **kwargs):
    from astropy.utils.masked import Masked

    if isinstance(q, Masked) or out is not None and not isinstance(out, Masked):
        raise NotImplementedError

    a = Masked(a)
    q = np.asanyarray(q)
    if (NUMPY_LT_2_0 and not np.lib.function_base._quantile_is_valid(q)) or (
        not NUMPY_LT_2_0 and not np.lib._function_base_impl._quantile_is_valid(q)
    ):
        raise ValueError("Quantiles must be in the range [0, 1]")

    if NUMPY_LT_1_24:
        keepdims = kwargs.pop("keepdims", False)
        r, k = np.lib.function_base._ureduce(
            a, func=_masked_quantile, q=q, axis=axis, out=out, **kwargs
        )
        result = (r.reshape(q.shape + k) if keepdims else r) if out is None else out

    elif NUMPY_LT_2_0:
        result = np.lib.function_base._ureduce(
            a, func=_masked_quantile, q=q, axis=axis, out=out, **kwargs
        )
    else:
        result = np.lib._function_base_impl._ureduce(
            a, func=_masked_quantile, q=q, axis=axis, out=out, **kwargs
        )

    return result, None, None


@dispatched_function
def percentile(a, q, *args, **kwargs):
    q = np.true_divide(q, 100)
    return quantile(a, q, *args, **kwargs)


@dispatched_function
def array_equal(a1, a2, equal_nan=False):
    (a1d, a2d), (a1m, a2m) = _get_data_and_masks(a1, a2)
    if a1d.shape != a2d.shape:
        return False, None, None

    equal = a1d == a2d
    if equal_nan:
        equal |= np.isnan(a1d) & np.isnan(a2d)
    return bool((equal | a1m | a2m).all()), None, None


@dispatched_function
def array_equiv(a1, a2):
    return bool((a1 == a2).all()), None, None


@dispatched_function
def where(condition, *args):
    from astropy.utils.masked import Masked

    if not args:
        return condition.nonzero(), None, None

    condition, c_mask = Masked._get_data_and_mask(condition)

    data, masks = _get_data_and_masks(*args)
    unmasked = np.where(condition, *data)
    mask = np.where(condition, *masks)
    if c_mask is not None:
        mask |= c_mask
    return Masked(unmasked, mask=mask), None, None


@dispatched_function
def choose(a, choices, out=None, mode="raise"):
    """Construct an array from an index array and a set of arrays to choose from.

    Like `numpy.choose`.  Masked indices in ``a`` will lead to masked output
    values and underlying data values are ignored if out of bounds (for
    ``mode='raise'``).  Any values masked in ``choices`` will be propagated
    if chosen.

    """
    from astropy.utils.masked import Masked

    a_data, a_mask = Masked._get_data_and_mask(a)
    if a_mask is not None and mode == "raise":
        # Avoid raising on masked indices.
        a_data = a.filled(fill_value=0)

    kwargs = {"mode": mode}
    if out is not None:
        if not isinstance(out, Masked):
            raise NotImplementedError
        kwargs["out"] = out.unmasked

    data, masks = _get_data_and_masks(*choices)
    data_chosen = np.choose(a_data, data, **kwargs)
    if out is not None:
        kwargs["out"] = out.mask

    mask_chosen = np.choose(a_data, masks, **kwargs)
    if a_mask is not None:
        mask_chosen |= a_mask

    return Masked(data_chosen, mask_chosen) if out is None else out, None, None


@apply_to_both
def select(condlist, choicelist, default=0):
    """Return an array drawn from elements in choicelist, depending on conditions.

    Like `numpy.select`, with masks in ``choicelist`` are propagated.
    Any masks in ``condlist`` are ignored.

    """
    from astropy.utils.masked import Masked

    condlist = [c.unmasked if isinstance(c, Masked) else c for c in condlist]

    data_list, mask_list = _get_data_and_masks(*choicelist)
    default = Masked(default) if default is not np.ma.masked else Masked(0, mask=True)
    return (
        (condlist, data_list, default.unmasked),
        (condlist, mask_list, default.mask),
        {},
        None,
    )


@dispatched_function
def piecewise(x, condlist, funclist, *args, **kw):
    """Evaluate a piecewise-defined function.

    Like `numpy.piecewise` but for masked input array ``x``.
    Any masks in ``condlist`` are ignored.

    """
    # Copied implementation from numpy.lib._function_base_impl.piecewise,
    # just to ensure output is Masked.
    n2 = len(funclist)
    # undocumented: single condition is promoted to a list of one condition
    if np.isscalar(condlist) or (
        not isinstance(condlist[0], (list, np.ndarray)) and x.ndim != 0
    ):  # pragma: no cover
        condlist = [condlist]

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

    # The one real change...
    y = np.zeros_like(x)
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

    for item, value in zip(where, what):
        y[item] = value

    return y, None, None


@dispatched_function
def interp(x, xp, fp, *args, **kwargs):
    """One-dimensional linear interpolation.

    Like `numpy.interp`, but any masked points in ``xp`` and ``fp``
    are ignored.  Any masked values in ``x`` will still be evaluated,
    but masked on output.
    """
    from astropy.utils.masked import Masked

    xd, xm = Masked._get_data_and_mask(x)
    if isinstance(xp, Masked) or isinstance(fp, Masked):
        (xp, fp), (xpm, fpm) = _get_data_and_masks(xp, fp)
        if xp.ndim == fp.ndim == 1:
            # Avoid making arrays 1-D; will just raise below.
            m = xpm | fpm
            xp = xp[~m]
            fp = fp[~m]

    result = np.interp(xd, xp, fp, *args, **kwargs)
    return (result if xm is None else Masked(result, xm.copy())), None, None


@dispatched_function
def lexsort(keys, axis=-1):
    """Perform an indirect stable sort using a sequence of keys.

    Like `numpy.lexsort` but for possibly masked ``keys``.  Masked
    values are sorted towards the end for each key.
    """
    # Sort masks to the end.
    from .core import Masked

    new_keys = []
    for key in keys:
        if isinstance(key, Masked):
            # If there are other keys below, want to be sure that
            # for masked values, those other keys set the order.
            new_key = key.unmasked
            if new_keys and key.mask.any():
                new_key = new_key.copy()
                new_key[key.mask] = new_key.flat[0]
            new_keys.extend([new_key, key.mask])
        else:
            new_keys.append(key)

    return np.lexsort(new_keys, axis=axis), None, None


@dispatched_function
def apply_over_axes(func, a, axes):
    # Copied straight from numpy/lib/shape_base, just to omit its
    # val = asarray(a); if only it had been asanyarray, or just not there
    # since a is assumed to an an array in the next line...
    # Which is what we do here - we can only get here if it is Masked.
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

    return val, None, None


class MaskedFormat:
    """Formatter for masked array scalars.

    For use in `numpy.array2string`, wrapping the regular formatters such
    that if a value is masked, its formatted string is replaced.

    Typically initialized using the ``from_data`` class method.
    """

    def __init__(self, format_function):
        self.format_function = format_function
        # Special case for structured void and subarray: we need to make all the
        # format functions for the items masked as well.
        # TODO: maybe is a separate class is more logical?
        ffs = getattr(format_function, "format_functions", None)
        if ffs:
            # StructuredVoidFormat: multiple format functions to be changed.
            self.format_function.format_functions = [MaskedFormat(ff) for ff in ffs]

        ff = getattr(format_function, "format_function", None)
        if ff:
            # SubarrayFormat: change format function for the elements.
            self.format_function.format_function = MaskedFormat(ff)

    def __call__(self, x):
        if x.dtype.names:
            # The replacement of x with a list is needed because the function
            # inside StructuredVoidFormat iterates over x, which works for an
            # np.void but not an array scalar.
            return self.format_function([x[field] for field in x.dtype.names])

        if x.shape:
            # For a subarray pass on the data directly, since the
            # items will be iterated on inside the function.
            return self.format_function(x)

        # Single element: first just typeset it normally, replace with masked
        # string if needed.
        string = self.format_function(x.unmasked[()])
        if x.mask:
            # Strikethrough would be neat, but terminal needs a different
            # formatting than, say, jupyter notebook.
            # return "\x1B[9m"+string+"\x1B[29m"
            # return ''.join(s+'\u0336' for s in string)
            n = min(3, max(1, len(string)))
            return " " * (len(string) - n) + "\u2014" * n
        else:
            return string

    @classmethod
    def from_data(cls, data, **options):
        if NUMPY_LT_2_0:
            from numpy.core.arrayprint import _get_format_function
        else:
            from numpy._core.arrayprint import _get_format_function

        return cls(_get_format_function(data, **options))


def _array2string(a, options, separator=" ", prefix=""):
    # Mostly copied from numpy.core.arrayprint, except:
    # - The format function is wrapped in a mask-aware class;
    # - Arrays scalars are not cast as arrays.
    if NUMPY_LT_2_0:
        from numpy.core.arrayprint import _formatArray, _leading_trailing
    else:
        from numpy._core.arrayprint import _formatArray, _leading_trailing

    data = np.asarray(a)

    if a.size > options["threshold"]:
        summary_insert = "..."
        data = _leading_trailing(data, options["edgeitems"])
    else:
        summary_insert = ""

    # find the right formatting function for the array
    format_function = MaskedFormat.from_data(data, **options)

    # skip over "["
    next_line_prefix = " "
    # skip over array(
    next_line_prefix += " " * len(prefix)

    lst = _formatArray(
        a,
        format_function,
        options["linewidth"],
        next_line_prefix,
        separator,
        options["edgeitems"],
        summary_insert,
        options["legacy"],
    )
    return lst


@dispatched_function
def array2string(
    a,
    max_line_width=None,
    precision=None,
    suppress_small=None,
    separator=" ",
    prefix="",
    style=np._NoValue,
    formatter=None,
    threshold=None,
    edgeitems=None,
    sign=None,
    floatmode=None,
    suffix="",
):
    # Copied from numpy.core.arrayprint, but using _array2string above.
    if NUMPY_LT_2_0:
        from numpy.core.arrayprint import _format_options, _make_options_dict
    else:
        from numpy._core.arrayprint import _format_options, _make_options_dict

    overrides = _make_options_dict(
        precision,
        threshold,
        edgeitems,
        max_line_width,
        suppress_small,
        None,
        None,
        sign,
        formatter,
        floatmode,
    )
    options = _format_options.copy()
    options.update(overrides)

    options["linewidth"] -= len(suffix)

    # treat as a null array if any of shape elements == 0
    if a.size == 0:
        result = "[]"
    else:
        result = _array2string(a, options, separator, prefix)

    return result, None, None


def _array_str_scalar(x):
    # This wraps np.array_str for use as a format function in
    # MaskedFormat. We cannot use it directly as format functions
    # expect numpy scalars, while np.array_str expects an array.
    return np.array_str(np.array(x))


@dispatched_function
def array_str(a, max_line_width=None, precision=None, suppress_small=None):
    # Override to change special treatment of array scalars, since the numpy
    # code turns the masked array scalar into a regular array scalar.
    # By going through MaskedFormat, we can replace the string as needed.
    if a.shape == () and a.dtype.names is None:
        return MaskedFormat(_array_str_scalar)(a), None, None
    else:
        return array2string(a, max_line_width, precision, suppress_small, " ", "")


# For the nanfunctions, we just treat any nan as an additional mask.
_nanfunc_fill_values = {"nansum": 0, "nancumsum": 0, "nanprod": 1, "nancumprod": 1}


def masked_nanfunc(nanfuncname):
    np_func = getattr(np, nanfuncname[3:])
    fill_value = _nanfunc_fill_values.get(nanfuncname, None)

    def nanfunc(a, *args, **kwargs):
        from astropy.utils.masked import Masked

        a, mask = Masked._get_data_and_mask(a)
        if issubclass(a.dtype.type, np.inexact):
            nans = np.isnan(a)
            mask = nans if mask is None else (nans | mask)

        if mask is not None:
            a = Masked(a, mask)
            if fill_value is not None:
                a = a.filled(fill_value)

        return np_func(a, *args, **kwargs), None, None

    doc = f"Like `numpy.{nanfuncname}`, skipping masked values as well.\n\n"
    if fill_value is not None:
        # sum, cumsum, prod, cumprod
        doc += (
            f"Masked/NaN values are replaced with {fill_value}. "
            "The output is not masked."
        )
    elif "arg" in nanfuncname:
        doc += (
            "No exceptions are raised for fully masked/NaN slices.\n"
            "Instead, these give index 0."
        )
    else:
        doc += (
            "No warnings are given for fully masked/NaN slices.\n"
            "Instead, they are masked in the output."
        )

    nanfunc.__doc__ = doc
    nanfunc.__name__ = nanfuncname

    return nanfunc


_nplibnanfunctions = np.lib.nanfunctions if NUMPY_LT_2_0 else np.lib._nanfunctions_impl
for nanfuncname in _nplibnanfunctions.__all__:
    globals()[nanfuncname] = dispatched_function(
        masked_nanfunc(nanfuncname), helps=getattr(np, nanfuncname)
    )


# Add any dispatched or helper function that has a docstring to
# __all__, so they will be typeset by sphinx. The logic is that for
# those presumably the use of the mask is not entirely obvious.
__all__ += sorted(  # noqa: PLE0605
    helper.__name__
    for helper in (
        set(APPLY_TO_BOTH_FUNCTIONS.values()) | set(DISPATCHED_FUNCTIONS.values())
    )
    if helper.__doc__
)
