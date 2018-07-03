# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""A few utilities useful for determining whether quantities are close.

Wrappers around numpy equivalents that take units into account.
"""
import numpy as np

from .quantity import Quantity
from ..core import UnitsError, dimensionless_unscaled


__all__ = ['isclose', 'allclose']


def isclose(a, b, rtol=1.e-5, atol=None, **kwargs):
    """
    Notes
    -----
    Returns True if two arrays are element-wise equal within a tolerance.

    This is a :class:`~astropy.units.Quantity`-aware version of
    :func:`numpy.isclose`.
    """
    return np.isclose(*_unquantify_allclose_arguments(a, b, rtol, atol),
                      **kwargs)


def allclose(a, b, rtol=1.e-5, atol=None, **kwargs):
    """
    Notes
    -----
    Returns True if two arrays are element-wise equal within a tolerance.

    This is a :class:`~astropy.units.Quantity`-aware version of
    :func:`numpy.allclose`.
    """
    return np.allclose(*_unquantify_allclose_arguments(a, b, rtol, atol),
                       **kwargs)


def _unquantify_allclose_arguments(actual, desired, rtol, atol):
    actual = Quantity(actual, subok=True, copy=False)

    desired = Quantity(desired, subok=True, copy=False)
    try:
        desired = desired.to(actual.unit)
    except UnitsError:
        raise UnitsError("Units for 'desired' ({0}) and 'actual' ({1}) "
                         "are not convertible"
                         .format(desired.unit, actual.unit))

    if atol is None:
        # by default, we assume an absolute tolerance of 0
        atol = Quantity(0)
    else:
        atol = Quantity(atol, subok=True, copy=False)
        try:
            atol = atol.to(actual.unit)
        except UnitsError:
            raise UnitsError("Units for 'atol' ({0}) and 'actual' ({1}) "
                             "are not convertible"
                             .format(atol.unit, actual.unit))

    rtol = Quantity(rtol, subok=True, copy=False)
    try:
        rtol = rtol.to(dimensionless_unscaled)
    except Exception:
        raise UnitsError("`rtol` should be dimensionless")

    return actual.value, desired.value, rtol.value, atol.value
