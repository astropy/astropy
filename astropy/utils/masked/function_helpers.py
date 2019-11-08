# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np

from astropy.units.quantity_helper.function_helpers import FunctionAssigner


APPLY_TO_BOTH_FUNCTIONS = {}
DISPATCHED_FUNCTIONS = {}


apply_to_both = FunctionAssigner(APPLY_TO_BOTH_FUNCTIONS)
dispatched_function = FunctionAssigner(DISPATCHED_FUNCTIONS)


def _data_mask(array):
    from .core import Masked

    if isinstance(array, Masked):
        return array.unmasked, array.mask
    else:
        return array, np.zeros(np.shape(array), bool)


def _data_masks(*args):
    data = []
    masks = []
    for arg in args:
        datum, mask = _data_mask(arg)
        data.append(datum)
        masks.append(mask)

    return data, masks


@apply_to_both
def concatenate(arrays, axis=0, out=None):
    data, masks = _data_masks(*arrays)
    return data, masks, (axis,), {}, out


@apply_to_both
def take_along_axis(arr, indices, axis):
    data, mask = _data_mask(arr)
    return data, mask, (indices, axis), {}, None


@apply_to_both
def broadcast_to(array, shape, subok=True):
    data, mask = _data_mask(array) if subok else (array.unmasked, None)
    return data, mask, (shape, subok), {}, None


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
