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

from astropy.units.quantity_helper.function_helpers import (
    FunctionAssigner, IGNORED_FUNCTIONS)


# This module should not really be imported, but we define __all__
# such that sphinx can typeset the functions with docstrings.
# The latter are added to __all__ at the end.
__all__ = ['MASKED_SAFE_FUNCTIONS', 'APPLY_TO_BOTH_FUNCTIONS',
           'DISPATCHED_FUNCTIONS', 'UNSUPPORTED_FUNCTIONS']


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
"""

DISPATCHED_FUNCTIONS = {}
"""Dict of functions that provide the numpy function's functionality.

These are for more complicated versions where the numpy function itself
cannot easily be used.  It should return either the result of the
function, or a tuple consisting of the unmasked result, the mask for the
result and a possible output instance.
"""

UNSUPPORTED_FUNCTIONS = set()
"""Set of numpy functions that are not supported for masked arrays.

For most, masked input simply makes no sense.  If a function is included
for which this does not seem to hold, please raise an issue.
"""

MASKED_SAFE_FUNCTIONS |= {
    # np.core.numeric
    np.isclose, np.allclose,
    # np.lib.fromnumeric (covered by methods)
    np.amax, np.amin, np.sum, np.cumsum, np.any, np.all,
    np.sometrue, np.alltrue, np.prod, np.product, np.cumprod, np.cumproduct,
    # np.lib.function_base
    np.diff, np.extract,
    # np.lib.shape_base
    np.split, np.array_split, np.hsplit, np.vsplit, np.dsplit,
    np.put_along_axis, np.apply_along_axis,
    # np.lib.index_tricks
    np.diag_indices_from, np.triu_indices_from, np.tril_indices_from,
    np.fill_diagonal,
    # np.lib.ufunclike
    np.fix, np.isneginf, np.isposinf,
    # np.lib.function_base
    np.angle, np.i0,
}

# No support for the functions also not supported by Quantity
# (io, polynomial, etc.).
UNSUPPORTED_FUNCTIONS |= IGNORED_FUNCTIONS

# TODO: the following could in principle be supported.
UNSUPPORTED_FUNCTIONS |= {
    np.pad, np.is_busday, np.busday_count, np.busday_offset,
    np.unravel_index, np.ravel_multi_index, np.ix_}

apply_to_both = FunctionAssigner(APPLY_TO_BOTH_FUNCTIONS)
dispatched_function = FunctionAssigner(DISPATCHED_FUNCTIONS)


def _data_masks(*args):
    """Separate out arguments into tuples of data and masks.

    An all-False mask is created if an argument does not have a mask.
    """
    from .core import Masked

    data = []
    masks = []
    for arg in args:
        if isinstance(arg, Masked):
            data.append(arg.unmasked)
            masks.append(arg.mask)
        else:
            data.append(arg)
            masks.append(np.zeros(np.shape(arg), bool))

    return tuple(data), tuple(masks)


@dispatched_function
def datetime_as_string(arr, *args, **kwargs):
    return (np.datetime_as_string(arr.unmasked, *args, **kwargs),
            arr.mask.copy(), None)


@dispatched_function
def sinc(x):
    return np.sinc(x.unmasked), x.mask.copy(), None


@dispatched_function
def iscomplex(x):
    return np.iscomplex(x.unmasked), x.mask.copy(), None


@dispatched_function
def unwrap(p, *args, **kwargs):
    return np.unwrap(p.unmasked, *args, **kwargs), p.mask.copy(), None


@apply_to_both(helps={
    np.copy, np.asfarray, np.real_if_close, np.sort_complex, np.resize,
    np.moveaxis, np.rollaxis, np.expand_dims, np.squeeze, np.roll, np.take,
    np.repeat})
def masked_a_helper(a, *args, **kwargs):
    data, mask = _data_masks(a)
    return data + args, mask + args, kwargs, None


@apply_to_both(helps={np.flip, np.flipud, np.fliplr, np.rot90, np.triu, np.tril})
def masked_m_helper(m, *args, **kwargs):
    data, mask = _data_masks(m)
    return data + args, mask + args, kwargs, None


@apply_to_both(helps={np.diag, np.diagflat})
def masked_v_helper(v, *args, **kwargs):
    data, mask = _data_masks(v)
    return data + args, mask + args, kwargs, None


@apply_to_both(helps={np.tile})
def masked_A_helper(A, *args, **kwargs):
    data, mask = _data_masks(A)
    return data + args, mask + args, kwargs, None


@apply_to_both(helps={np.delete})
def masked_arr_helper(array, *args, **kwargs):
    data, mask = _data_masks(array)
    return data + args, mask + args, kwargs, None


@apply_to_both
def broadcast_to(array,  shape, subok=False):
    """Broadcast array to the given shape.

    Like `numpy.broadcast_to`, and applied to both unmasked data and mask.
    Note that ``subok`` is taken to mean whether or not subclasses of
    the unmasked data and mask are allowed, i.e., for ``subok=False``,
    a `~astropy.utils.masked.MaskedNDArray` will be returned.
    """
    data, mask = _data_masks(array)
    return data, mask, dict(shape=shape, subok=subok), None


@dispatched_function
def empty_like(prototype, dtype=None, order='K', subok=True, shape=None):
    unmasked = np.empty_like(prototype.unmasked, dtype=dtype, order=order,
                             subok=subok, shape=shape)
    if dtype is not None:
        dtype = (np.ma.make_mask_descr(unmasked.dtype)
                 if unmasked.dtype.names else np.dtype('?'))
    mask = np.empty_like(prototype.mask, dtype=dtype, order=order,
                         subok=subok, shape=shape)

    return unmasked, mask, None


@dispatched_function
def zeros_like(a, dtype=None, order='K', subok=True, shape=None):
    unmasked = np.zeros_like(a.unmasked, dtype=dtype, order=order,
                             subok=subok, shape=shape)
    return unmasked, False, None


@dispatched_function
def ones_like(a, dtype=None, order='K', subok=True, shape=None):
    unmasked = np.ones_like(a.unmasked, dtype=dtype, order=order,
                            subok=subok, shape=shape)
    return unmasked, False, None


@dispatched_function
def full_like(a, fill_value, dtype=None, order='K', subok=True, shape=None):
    result = np.zeros_like(a, dtype=dtype, order=order, subok=subok, shape=shape)
    result[...] = fill_value
    return result, None, None


@dispatched_function
def put(a, ind, v, mode='raise'):
    from astropy.utils.masked import Masked
    if isinstance(ind, Masked) or not isinstance(a, Masked):
        return NotImplemented, None, None
    v_data, v_mask = a._data_mask(v)
    if v_data is not None:
        np.put(a.unmasked, ind, v_data, mode=mode)
    np.put(a.mask, ind, v_mask, mode=mode)
    return None, None, None


@dispatched_function
def putmask(a, mask, values):
    from astropy.utils.masked import Masked
    if isinstance(mask, Masked) or not isinstance(a, Masked):
        return NotImplemented, None, None
    values_data, values_mask = a._data_mask(values)
    if values_data is not None:
        np.putmask(a.unmasked, mask, values_data)
    np.putmask(a.mask, mask, values_mask)
    return None, None, None


@dispatched_function
def place(arr, mask, vals):
    from astropy.utils.masked import Masked
    if isinstance(mask, Masked) or not isinstance(arr, Masked):
        return NotImplemented, None, None
    vals_data, vals_mask = arr._data_mask(vals)
    if vals_data is not None:
        np.place(arr.unmasked, mask, vals_data)
    np.place(arr.mask, mask, vals_mask)
    return None, None, None


@dispatched_function
def copyto(dst, src,  casting='same_kind', where=True):
    from astropy.utils.masked import Masked
    if not isinstance(dst, Masked):
        return NotImplemented, None, None
    src_data, src_mask = dst._data_mask(src)

    if src_data is not None:
        np.copyto(dst.unmasked, src_data, casting=casting, where=where)
    if src_mask is not None:
        np.copyto(dst.mask, src_mask, where=where)
    return None, None, None


@dispatched_function
def packbits(a, *args, **kwargs):
    result = np.packbits(a.unmasked, *args, **kwargs)
    mask = np.packbits(a.mask, *args, **kwargs).astype(bool)
    return result, mask, None


@dispatched_function
def unpackbits(a, *args, **kwargs):
    result = np.unpackbits(a.unmasked, *args, **kwargs)
    mask = np.zeros(a.shape, dtype='u1')
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
        weights, w_mask = Masked._data_mask(weights)
        if w_mask is not None:
            mask = np.bincount(x, w_mask.astype(int),
                               minlength=minlength).astype(bool)
    result = np.bincount(x, weights, minlength=0)
    return result, mask, None


@apply_to_both
def concatenate(arrays, axis=0, out=None):
    data, masks = _data_masks(*arrays)
    return (data,), (masks,), dict(axis=axis), out


@apply_to_both
def append(arr, values, axis=None):
    data, masks = _data_masks(arr, values)
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
    (arrays, list_ndim, result_ndim,
     final_size) = np.core.shape_base._block_setup(arrays)
    shape, slices, arrays = np.core.shape_base._block_info_recursion(
        arrays, list_ndim, result_ndim)
    dtype = np.result_type(*[arr.dtype for arr in arrays])
    F_order = all(arr.flags['F_CONTIGUOUS'] for arr in arrays)
    C_order = all(arr.flags['C_CONTIGUOUS'] for arr in arrays)
    order = 'F' if F_order and not C_order else 'C'
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
    data = [(arg.unmasked if is_masked else arg)
            for arg, is_masked in zip(args, are_masked)]
    results = np.broadcast_arrays(*data, subok=subok)

    shape = results[0].shape if isinstance(results, list) else results.shape
    masks = [(np.broadcast_to(arg.mask, shape, subok=subok)
              if is_masked else None)
             for arg, is_masked in zip(args, are_masked)]
    results = [(Masked(result, mask) if mask is not None else result)
               for (result, mask) in zip(results, masks)]
    return (results if len(results) > 1 else results[0]), None, None


@apply_to_both
def insert(arr, obj, values, axis=None):
    data, masks = _data_masks(arr, values)
    return ((data[0], obj, data[1], axis),
            (masks[0], obj, masks[1], axis), {}, None)


@apply_to_both
def where(condition, x=None, y=None):
    """Return elements chosen from ``x`` or ``y`` depending on ``condition``.

    Like `numpy.where`, with any mask on ``condition`` ignored, while
    masks on ``x`` and ``y`` are propagated.

    """
    from astropy.utils.masked import Masked
    if isinstance(condition, Masked):
        condition = condition.unmasked  # simply ignore mask

    if x is None and y is None:
        return condition.nonzero(), None, None

    if x is None or y is None:
        return NotImplemented, None, None

    data, masks = _data_masks(x, y)
    return (condition,) + data, (condition,) + masks, {}, None


@dispatched_function
def choose(a, choices, out=None, **kwargs):
    (a, *data), (mask, *masks) = _data_masks(a, *choices)
    if out is None:
        data_chosen = np.choose(a, data, **kwargs)
        mask_chosen = np.choose(mask, masks, **kwargs) | mask
        return data_chosen, mask_chosen, None
    else:
        from astropy.utils.masked import Masked
        if not isinstance(out, Masked):
            return NotImplemented, None, None
        np.choose(a, data, out=out.unmasked, **kwargs)
        np.choose(mask, masks, out=out.mask, **kwargs)
        return out, None, None


@apply_to_both
def select(condlist, choicelist, default=0):
    from astropy.utils.masked import Masked

    data_list, mask_list = _data_masks(*choicelist)
    default = Masked(default) if default is not np.ma.masked else Masked(0, mask=True)
    return ((condlist, data_list, default.unmasked),
            (condlist, mask_list, default.mask), {}, None)


@dispatched_function
def piecewise(x, condlist, funclist, *args, **kw):
    # Copied implementation from numpy.lib.function_base.piecewise,
    # just to ensure output is Masked.
    n2 = len(funclist)
    # undocumented: single condition is promoted to a list of one condition
    if np.isscalar(condlist) or (
            not isinstance(condlist[0], (list, np.ndarray)) and x.ndim != 0):
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
    xd, xm = Masked._data_mask(x)
    if isinstance(xp, Masked) or isinstance(fp, Masked):
        (xp, fp), (xpm, fpm) = _data_masks(xp, fp)
        if xp.ndim == fp.ndim == 1:
            # Avoid making arrays 1-D; will just raise below.
            m = xpm | fpm
            xp = xp[m]
            fp = fp[m]

    result = np.interp(xd, xp, fp, *args, **kwargs)
    return result, None if xm is None else xm.copy(), None


@dispatched_function
def lexsort(keys, axis=-1):
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
                raise ValueError("function is not returning "
                                 "an array of the correct shape")
    # Returning mask is None to signal nothing should happen to
    # the output.
    return val, None, None
class MaskedFormat:
    def __init__(self, format_function):
        self.format_function = format_function
        # Special case for structured void: we need to make all the
        # format functions for the items masked as well.
        # TODO: maybe is a separate class is more logical?
        ffs = getattr(format_function, 'format_functions', None)
        if ffs:
            self.format_function.format_functions = [MaskedFormat(ff) for ff in ffs]

    def __call__(self, x):
        if x.dtype.names:
            # The replacement of x with a list is needed because the function
            # inside StructuredVoidFormat iterates over x, which works for an
            # np.void but not an array scalar.
            return self.format_function([x[field] for field in x.dtype.names])

        string = self.format_function(x.unmasked[()])
        if x.mask:
            # Strikethrough would be neat, but terminal needs a different
            # formatting than, say, jupyter notebook.
            # return "\x1B[9m"+string+"\x1B[29m"
            # return ''.join(s+'\u0336' for s in string)
            n = min(3, max(1, len(string)))
            return ' ' * (len(string)-n) + '\u2014' * n
        else:
            return string


def _array2string(a, options, separator=' ', prefix=""):
    # Mostly copied from numpy.core.arrayprint, except:
    # - The format function is wrapped in a mask-aware class;
    # - Arrays scalars are not cast as arrays.
    from numpy.core.arrayprint import (_leading_trailing, _get_format_function,
                                       _formatArray)
    data = np.asarray(a)

    if a.size > options['threshold']:
        summary_insert = "..."
        data = _leading_trailing(data, options['edgeitems'])
    else:
        summary_insert = ""

    # find the right formatting function for the array
    format_function = MaskedFormat(_get_format_function(data, **options))

    # skip over "["
    next_line_prefix = " "
    # skip over array(
    next_line_prefix += " "*len(prefix)

    lst = _formatArray(a, format_function, options['linewidth'],
                       next_line_prefix, separator, options['edgeitems'],
                       summary_insert, options['legacy'])
    return lst


@dispatched_function
def array2string(a, max_line_width=None, precision=None,
                 suppress_small=None, separator=' ', prefix="",
                 style=np._NoValue, formatter=None, threshold=None,
                 edgeitems=None, sign=None, floatmode=None, suffix=""):
    # Copied from numpy.core.arrayprint, but using _array2string above.
    from numpy.core.arrayprint import _make_options_dict, _format_options

    overrides = _make_options_dict(precision, threshold, edgeitems,
                                   max_line_width, suppress_small, None, None,
                                   sign, formatter, floatmode)
    options = _format_options.copy()
    options.update(overrides)

    options['linewidth'] -= len(suffix)

    # treat as a null array if any of shape elements == 0
    if a.size == 0:
        return "[]", None, None

    return _array2string(a, options, separator, prefix), None, None


@dispatched_function
def array_str(a, max_line_width=None, precision=None, suppress_small=None):
    # Override to avoid special treatment of array scalars.
    return array2string(a, max_line_width, precision, suppress_small, ' ', "")


# Add any dispatched or helper function that has a docstring to
# __all__, so they will be typeset by sphinx. The logic is that for
# those presumably the use of the mask is not entirely obvious.
__all__ += sorted(helper.__name__ for helper in (
    set(APPLY_TO_BOTH_FUNCTIONS.values())
    | set(DISPATCHED_FUNCTIONS.values())) if helper.__doc__)
