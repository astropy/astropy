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

__all__ = ['NormalDistribution', 'PoissonDistribution',
           'UniformDistribution']


class NormalDistribution(Distribution):
    """
    A Gaussian/normal Distribution

    Parameters
    ----------
    center : `Quantity`
        The center of this `NormalDistribution`
    std : `Quantity` or `None`
        The standard deviation/σ of this distribution. Shape must match and unit
        must be compatible with ``center``, or be `None` (if ``var`` or ``ivar``
        are set).
    var : `Quantity` or `None`
        The variance of this distribution. Shape must match and unit must be
        compatible with ``center``, or be `None` (if ``std`` or ``ivar`` are set).
    ivar : `Quantity` or `None`
        The inverse variance of this distribution. Shape must match and unit
        must be compatible with ``center``, or be `None` (if ``std`` or ``var``
        are set).
    n_samples : int
        The number of Monte Carlo samples to use with this distribution

    Remaining keywords are passed into the `Distribution` constructor
    """
    def __new__(cls, center, std=None, var=None, ivar=None, n_samples=1000, **kwargs):
        center = np.asanyarray(center)
        if var is not None:
            if std is None:
                std = np.asanyarray(var)**0.5
            else:
                raise ValueError('NormalDistribution cannot take both std and var')
        if ivar is not None:
            if std is None:
                std = np.asanyarray(ivar)**-0.5
            else:
                raise ValueError('NormalDistribution cannot take both ivar and '
                                 'and std or var')
        if std is None:
            raise ValueError('NormalDistribution requires one of std, var, or ivar')
        else:
            std = np.asanyarray(std)

        randshape = np.broadcast(std, center).shape + (n_samples,)
        distr = center[..., np.newaxis] + np.random.randn(*randshape) * std[..., np.newaxis]
        self = super(NormalDistribution, cls).__new__(cls, distr, **kwargs)
        self.distr_std = std
        return self


class PoissonDistribution(Distribution):
    """
    A Poisson Distribution

    Parameters
    ----------
    poissonval : `Quantity`
        The center value of this `PoissonDistribution` (i.e., λ).
    n_samples : int
        The number of Monte Carlo samples to use with this distribution

    Remaining keywords are passed into the `Distribution` constructor
    """
    def __new__(cls, poissonval, n_samples=1000, **kwargs):

        # we convert to arrays because np.random.poisson has trouble with quantities
        has_unit = False
        if hasattr(poissonval, 'unit'):
            has_unit = True
            poissonarr = np.asanyarray(poissonval.value)
        else:
            poissonarr = np.asanyarray(poissonval)
        randshape = poissonarr.shape + (n_samples,)

        distr = np.random.poisson(poissonarr[..., np.newaxis], randshape)
        if has_unit:
            # re-attach the unit
            distr = distr * poissonval.unit

        self = super(PoissonDistribution, cls).__new__(cls, distr,
                                                       **kwargs)
        self.distr_std = poissonval**0.5
        return self


class UniformDistribution(Distribution):
    """
    A Uniform Distribution.

    Parameters
    ----------
    lower : array-like
        The lower edge of this distribution. If a `~astropy.units.Quantity`, the
        distribution will have the same units as ``lower``.
    upper : `Quantity`
        The upper edge of this distribution. Must match shape and if a
        `~astropy.units.Quantity` must have compatible units with ``lower``.
    n_samples : int
        The number of Monte Carlo samples to use with this distribution

    Remaining keywords are passed into the `Distribution` constructor.
    """
    def __new__(cls, lower, upper, n_samples=1000, **kwargs):
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
                             'units in UniformDistribution constructor.')

        if lowerarr.shape != upperarr.shape:
            raise ValueError('lower and upper must have consistent shapes in '
                             'UniformDistribution constructor')

        newshape = lowerarr.shape + (n_samples,)
        distr = np.random.uniform(lowerarr[..., np.newaxis],
                                  upperarr[..., np.newaxis], newshape)
        if unit is not None:
            distr = distr * unit

        self = super(UniformDistribution, cls).__new__(cls, distr,
                                                       **kwargs)
        return self

    @classmethod
    def from_center_width(cls, center, width, **kwargs):
        """
        Create a `UniformDistribution` from lower/upper bounds (instead of center
        and width as the regular constructor uses).

        Parameters
        ----------
        center
            The center value of this `UniformDistribution`.
        width
            The width of this `UniformDistribution`.  Must have the same shape and
            compatible units with ``center``.

        Remaining keywords are passed into the `UniformDistribution` constructor.
        """
        whalf = width/2
        return UniformDistribution(center - whalf, center + whalf, **kwargs)
