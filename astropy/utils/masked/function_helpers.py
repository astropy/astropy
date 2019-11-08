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
