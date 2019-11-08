# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np

from astropy.units.quantity_helper.function_helpers import FunctionAssigner


APPLY_TO_BOTH_FUNCTIONS = {}


apply_to_both = FunctionAssigner(APPLY_TO_BOTH_FUNCTIONS)


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
