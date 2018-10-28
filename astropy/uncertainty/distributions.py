# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Built-in distribution-creation functions.
"""
from warnings import warn

import numpy as np

from .. import units as u
from .core import Distribution

__all__ = ['normal', 'poisson', 'uniform']


def normal(center, *, std=None, var=None, ivar=None, n_samples,
           cls=Distribution, **kwargs):
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
    cls : class
        The class to use to create this distribution.  Typically a
        `Distribution` subclass.

    Remaining keywords are passed into the constructor of the ``cls``

    Returns
    -------
    distr : ``cls``, usually `Distribution`
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
    return cls(samples, **kwargs)


COUNT_UNITS = (u.count, u.electron, u.dimensionless_unscaled, u.chan, u.bin, u.vox, u.bit, u.byte)

def poisson(center, n_samples, cls=Distribution, **kwargs):
    """
    Create a Poisson distribution.

    Parameters
    ----------
    center : `~astropy.units.Quantity`
        The center value of this distribution (i.e., λ).
    n_samples : int
        The number of Monte Carlo samples to use with this distribution
    cls : class
        The class to use to create this distribution.  Typically a
        `Distribution` subclass.

    Remaining keywords are passed into the constructor of the ``cls``

    Returns
    -------
    distr : ``cls``, usually `Distribution`
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
        if center.unit == u.adu:
            warn('ADUs were provided to poisson.  ADUs are not strictly count'
                 'units because they need the gain to be applied. It is '
                 'recommended you apply the gain to convert to e.g. electrons.')
        elif center.unit not in COUNT_UNITS:
            warn('Unit {} was provided to poisson, which is not one of {}, '
                 'and therefore suspect as a "counting" unit.  Ensure you mean '
                 'to use Poisson statistics.'.format(center.unit, COUNT_UNITS))

        # re-attach the unit
        samples = samples * center.unit

    return cls(samples, **kwargs)


def uniform(*, lower=None, upper=None, center=None, width=None, n_samples,
            cls=Distribution, **kwargs):
    """
    Create a Uniform distriution from the lower and upper bounds.

    Note that this function requires keywords to be explicit, and requires
    either ``lower``/``upper`` or ``center``/``width``.

    Parameters
    ----------
    lower : array-like
        The lower edge of this distribution. If a `~astropy.units.Quantity`, the
        distribution will have the same units as ``lower``.
    upper : `~astropy.units.Quantity`
        The upper edge of this distribution. Must match shape and if a
        `~astropy.units.Quantity` must have compatible units with ``lower``.
    center : array-like
        The center value of the distribution. Cannot be provided at the same
        time as ``lower``/``upper``.
    width : array-like
        The width of the distribution.  Must have the same shape and compatible
        units with ``center`` (if any).
    n_samples : int
        The number of Monte Carlo samples to use with this distribution
    cls : class
        The class to use to create this distribution.  Typically a
        `Distribution` subclass.

    Remaining keywords are passed into the constructor of the ``cls``

    Returns
    -------
    distr : ``cls``, usually `Distribution`
        The sampled uniform distribution.
    """
    if center is None and width is None:
        lower = np.asanyarray(lower)
        upper = np.asanyarray(upper)
        if lower.shape != upper.shape:
            raise ValueError('lower and upper must have consistent shapes')
    elif upper is None and lower is None:
        center = np.asanyarray(center)
        width = np.asanyarray(width)
        lower = center - width/2
        upper = center + width/2
    else:
        raise ValueError('either upper/lower or center/width must be given '
                         'to uniform - other combinations are not valid')

    newshape = lower.shape + (n_samples,)
    if lower.shape == tuple() and upper.shape == tuple():
        width = upper - lower  # scalar
    else:
        width = (upper - lower)[:, np.newaxis]
        lower = lower[:, np.newaxis]
    samples = lower + width * np.random.uniform(size=newshape)

    return cls(samples, **kwargs)
