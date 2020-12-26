# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np

from astropy.units.quantity_helper.function_helpers import FunctionAssigner


MASKED_SAFE_FUNCTIONS = set()
APPLY_TO_BOTH_FUNCTIONS = {}
DISPATCHED_FUNCTIONS = {}


MASKED_SAFE_FUNCTIONS |= {np.diff}


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


@apply_to_both
def concatenate(arrays, axis=0, out=None):
    data, masks = _data_masks(*arrays)
    return (data,), (masks,), dict(axis=axis), out


@apply_to_both
def take_along_axis(arr, *args, **kwargs):
    data, mask = _data_masks(arr)
    return data+args, mask+args, kwargs, None


@apply_to_both
def broadcast_to(array, shape, subok=True):
    data, mask = _data_masks(array) if subok else ((array.unmasked,), None)
    return data, mask, dict(shape=shape, subok=subok), None


@dispatched_function
def broadcast_arrays(*args, subok=True):
    from .core import Masked

    are_masked = [isinstance(arg, Masked) for arg in args]
    data = [(arg.unmasked if is_masked else arg)
            for arg, is_masked in zip(args, are_masked)]
    results = np.broadcast_arrays(*data, subok=subok)
    if not subok:
        return results

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
def append(arr, values, *args, **kwargs):
    data, masks = _data_masks(arr, values)
    return data + args, masks + args, kwargs, None


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
                new_key[key.mask] = np.zeros_like(new_key, shape=())
            new_keys.extend([new_key, key.mask])
        else:
            new_keys.append(key)

    return np.lexsort(new_keys, axis=axis), None, None


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
        return "[]"

    return _array2string(a, options, separator, prefix), None, None


@dispatched_function
def array_str(a, max_line_width=None, precision=None, suppress_small=None):
    # Override to avoid special treatment of array scalars.
    return array2string(a, max_line_width, precision, suppress_small, ' ', "")
