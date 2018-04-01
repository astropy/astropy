# coding: utf-8
# Licensed like the corresponding numpy file; see licenses/NUMPY_LICENSE.rst
"""
Utilities that manipulate strides to achieve desirable effects.

An explanation of strides can be found in the "ndarray.rst" file in the
NumPy reference guide.
"""

import warnings

import numpy as np
from numpy.lib.stride_tricks import (
    broadcast_arrays as np_broadcast_arrays,
    broadcast_to as np_broadcast_to)

from ....exceptions import AstropyDeprecationWarning

__all__ = ['broadcast_arrays', 'broadcast_to', 'GE1P10']
__doctest_skip__ = ['*']


def GE1P10(module=np):
    return hasattr(module, 'broadcast_to')


def broadcast_arrays(*args, **kwargs):
    warnings.warn(
        'This function is deprecated, as it is available in all NumPy versions '
        'that this version of Astropy supports. You should use '
        'numpy.broadcast_arrays directly.', AstropyDeprecationWarning)
    return np_broadcast_arrays(*args, **kwargs)


def broadcast_to(*args, **kwargs):
    warnings.warn(
        'This function is deprecated, as it is available in all NumPy versions '
        'that this version of Astropy supports. You should use '
        'numpy.broadcast_to directly.', AstropyDeprecationWarning)
    return np_broadcast_to(*args, **kwargs)
