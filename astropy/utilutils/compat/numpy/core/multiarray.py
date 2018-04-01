# coding: utf-8
# Licensed like numpy; see licenses/NUMPY_LICENSE.rst

import warnings

import numpy as np
from numpy import matmul as np_matmul

from ....exceptions import AstropyDeprecationWarning

__all__ = ['matmul', 'GE1P10']


def GE1P10(module=np):
    return hasattr(module, 'matmul')

def matmul(*args, **kwargs):
    warnings.warn(
        'This function is deprecated, as it is available in all NumPy versions '
        'that this version of Astropy supports. You should use '
        'numpy.matmul directly.', AstropyDeprecationWarning)
    return np_matmul(*args, **kwargs)
