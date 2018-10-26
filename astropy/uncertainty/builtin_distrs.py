# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Built-in distribution-creation functions.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from ..units import UnitsError
from .core import Distribution

__all__ = ['normal', 'poisson', 'uniform', 'uniform_center_width']


def normal(center, std=None, var=None, ivar=None, n_samples=1000,
           dcls=Distribution, **kwargs):
    """
    Create a Gaussian/normal distribution.

    Parameters
    ----------
    center : `~astropy.units.Quantity`
        The center of this distribution
    std : `~astropy.units.Quantity` or `None`
        The standard deviation/σ of this distribution. Shape must match and unit
        must be compatible with ``center``, or be `None` (if ``var`` or ``ivar``
        are set).
    var : `~astropy.units.Quantity` or `None`
        The variance of this distribution. Shape must match and unit must be
        compatible with ``center``, or be `None` (if ``std`` or ``ivar`` are set).
    ivar : `~astropy.units.Quantity` or `None`
        The inverse variance of this distribution. Shape must match and unit
        must be compatible with ``center``, or be `None` (if ``std`` or ``var``
        are set).
    n_samples : int
        The number of Monte Carlo samples to use with this distribution
    dcls : class
        The class to use to create this distribution.  Typically a
        `Distribution` subclass.

    Remaining keywords are passed into the constructor of the ``dcls``

    Returns
    -------
    distr : ``dcls``, usually `Distribution`
        The sampled Gaussian distribution.
    """
    center = np.asanyarray(center)
    if var is not None:
        if std is None:
            std = np.asanyarray(var)**0.5
        else:
            raise ValueError('normal cannot take both std and var')
    if ivar is not None:
        if std is None:
            std = np.asanyarray(ivar)**-0.5
        else:
            raise ValueError('normal cannot take both ivar and '
                             'and std or var')
    if std is None:
        raise ValueError('normal requires one of std, var, or ivar')
    else:
        std = np.asanyarray(std)

    randshape = np.broadcast(std, center).shape + (n_samples,)
    samples = center[..., np.newaxis] + np.random.randn(*randshape) * std[..., np.newaxis]
    distr = dcls(samples, **kwargs)
    distr.distr_std = std
    return distr


def poisson(center, n_samples=1000, dcls=Distribution, **kwargs):
    """
    Create a Poisson distribution.

    Parameters
    ----------
    center : `~astropy.units.Quantity`
        The center value of this distribution (i.e., λ).
    n_samples : int
        The number of Monte Carlo samples to use with this distribution
    dcls : class
        The class to use to create this distribution.  Typically a
        `Distribution` subclass.

    Remaining keywords are passed into the constructor of the ``dcls``

    Returns
    -------
    distr : ``dcls``, usually `Distribution`
        The sampled poisson distribution.
    """
    # we convert to arrays because np.random.poisson has trouble with quantities
    has_unit = False
    if hasattr(center, 'unit'):
        has_unit = True
        poissonarr = np.asanyarray(center.value)
    else:
        poissonarr = np.asanyarray(center)
    randshape = poissonarr.shape + (n_samples,)

    samples = np.random.poisson(poissonarr[..., np.newaxis], randshape)
    if has_unit:
        # re-attach the unit
        samples = samples * center.unit

    distr = dcls(samples, **kwargs)
    distr.distr_std = center**0.5
    return distr


def uniform(lower, upper, n_samples=1000, dcls=Distribution, **kwargs):
    """
    Create a Uniform Distribution.

    Parameters
    ----------
    lower : array-like
        The lower edge of this distribution. If a `~astropy.units.Quantity`, the
        distribution will have the same units as ``lower``.
    upper : `~astropy.units.Quantity`
        The upper edge of this distribution. Must match shape and if a
        `~astropy.units.Quantity` must have compatible units with ``lower``.
    n_samples : int
        The number of Monte Carlo samples to use with this distribution
    dcls : class
        The class to use to create this distribution.  Typically a
        `Distribution` subclass.

    Remaining keywords are passed into the constructor of the ``dcls``

    Returns
    -------
    distr : ``dcls``, usually `Distribution`
        The sampled uniform distribution.
    """
    lhasu = hasattr(lower, 'unit')
    uhasu = hasattr(upper, 'unit')
    unit = None
    if lhasu and uhasu:
        if lower.unit != upper.unit:
            upper = upper.to(lower.unit)
        unit = lower.unit

        lowerarr = np.asanyarray(lower.value)
        upperarr = np.asanyarray(upper.value)
    elif not lhasu and not uhasu:
        lowerarr = np.asanyarray(lower)
        upperarr = np.asanyarray(upper)
    else:
        raise UnitsError('lower and upper must have consistent (or no) '
                         'units in uniform.')

    if lowerarr.shape != upperarr.shape:
        raise ValueError('lower and upper must have consistent shapes in '
                         'uniform.')

    newshape = lowerarr.shape + (n_samples,)
    samples = np.random.uniform(lowerarr[..., np.newaxis],
                                upperarr[..., np.newaxis], newshape)
    if unit is not None:
        samples = samples * unit

    return dcls(samples, **kwargs)


def uniform_center_width(center, width, **kwargs):
    """
    Create a uniform distribution from lower/upper bounds (instead of center
    and width as the regular constructor uses).

    Parameters
    ----------
    center : array-like
        The center value of the distribution.
    width : array-like
        The width of the distribution.  Must have the same shape and compatible
        units with ``center`` (if any).

    Remaining keywords are passed into the `uniform` function.

    Returns
    -------
    distr : ``dcls``, usually `Distribution`
        The sampled uniform distribution.
    """
    whalf = width/2
    return uniform(center - whalf, center + whalf, **kwargs)
