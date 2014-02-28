# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Statistic functions used in `~astropy.modeling.fitting.py`.
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import numpy as np

__all__ = ['leastsquare']

def leastsquare(measured_vals, updated_model, weights, x, y=None):
    """
    Least square statistic with optional weights

    measured_vals : array
        measurd data values 
    updated_model : an instance of `~astropy.modeling.Model
        model with parameters set by the current iteration of the optimizer
    args : list
        other arguments passed through the fitter's objective function
    kwargs : dict
        other keyword arguments passed through the fitter

    """
    if y is None:
        model_vals = updated_model(x)
    else:
        model_vals = updated_model(x, y)
    if weights is None:
        return np.sum((model_vals - measured_vals) ** 2)
    else:
        return np.sum((weights * (model_vals - measured_vals)) ** 2)

