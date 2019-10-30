# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np

from astropy.units.quantity_helper.function_helpers import FunctionAssigner


APPLY_TO_BOTH_FUNCTIONS = {}


apply_to_both = FunctionAssigner(APPLY_TO_BOTH_FUNCTIONS)


def _data_mask(*args):
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

    return data, masks


@apply_to_both
def concatenate(arrays, axis=0, out=None):
    data, masks = _data_mask(*arrays)
    return data, masks, (axis,), {}, out
