# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Distribution class and associated machinery.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from . import Quantity

__all__ = ['Distribution', 'NormalDistribution', 'PoissonDistribution',
           'UniformDistribution']


class Distribution(Quantity):
    """
    A unitful value with associated uncertainty distribution.

    Notes
    -----
    This class has an attribute ``stat_view`` which can be set to a an
    `numpy.ndarray` subclass to determine what the summary statistics are
    converted to.  By default this is `~astropy.units.Quantity`, but it is
    available so that subclasses can change it.
    """

    __array_priority__ = Quantity.__array_priority__ + 1

    # this is what the summary statistics get downgraded to. Someone subclassing
    # Distribution might want to use this to get to something other than
    # Quantity.
    stat_view = Quantity

    def __new__(cls, distr, unit=None, *args, **kwargs):
        self = super(Distribution, cls).__new__(cls, distr, unit, *args, **kwargs)

        if len(self.shape) < 1:
            raise TypeError('Attempted to initialize a Distribution with a scalar')

        return self

    @property
    def n_samples(self):
        """
        The number of samples of this distribution.  A single int.
        """
        return self.shape[0]

    @property
    def distr_shape(self):
        """
        The shape of the underlying quantity (i.e., the *non-samples* part of
        this distribution). A tuple (possibly length-0 for scalars).
        """
        return self.shape[1:]

    @property
    def pdf_mean(self):
        """
        The mean of this distribution.
        """
        return self.mean(axis=0).view(self.stat_view)

    @property
    def pdf_std(self):
        """
        The standard deviation of this distribution.
        """
        return self.std(axis=0).view(self.stat_view)

    @property
    def pdf_var(self):
        """
        The variance of this distribution.
        """
        return self.var(axis=0).view(self.stat_view)

    @property
    def pdf_median(self):
        """
        The median of this distribution.
        """
        return np.median(self, axis=0).view(self.stat_view)

    @property
    def pdf_mad(self):
        """
        The median absolute deviation of this distribution.
        """
        return np.median(np.abs(self - self.pdf_median), axis=0).view(self.stat_view)

    # we set this by hand because the symbolic expression (below) requires scipy
    # _smad_scale_factor = 1 / scipy.stats.norm.ppf(0.75)
    _smad_scale_factor = 1.48260221850560203193936104071326553821563720703125

    @property
    def pdf_smad(self):
        """
        The median absolute deviation of this distribution rescaled to match the
        standard deviation for a normal distribution.
        """
        return self.pdf_mad * self._smad_scale_factor

    def percentiles(self, perc, **kwargs):
        """
        Compute percentiles of this Distribution.

        Parameters
        ----------
        perc : float or array of floats
            The desired  precentiles of the distribution (i.e., on [0,100]).
        kwargs
            Additional keywords are passed into `numpy.percentile`.

        Returns
        -------
        percs : Quantity of shape ``distr_shape``
            The ``fracs`` percentiles of this distribution.
        """
        perc = np.percentile(self, perc, axis=0)
        # numpy.percentile strips units for unclear reasons, so we have to make
        # a new object with units
        return self.stat_view(perc, unit=self.unit, copy=False)


class NormalDistribution(Distribution):
    """
    A Gaussian/normal Distribution

    Parameters
    ----------
    center : `Quantity`
        The center of this `NormalDistribution`
    std : `Quantity` or None
        The standard deviation/σ of this distribution. Shape must match and unit
        must be compatible with ``center``, or be None (if ``var`` or ``ivar``
        are set).
    var : `Quantity` or None
        The variance of this distribution. Shape must match and unit must be
        compatible with ``center``, or be None (if ``std`` or ``ivar`` are set).
    var : `Quantity` or None
        The inversev ariance of this distribution. Shape must match and unit
        must be compatible with ``center``, or be None (if ``std`` or ``var``
        are set).
    n_samples : int
        The number of Monte Carlo samples to use with this distribution

    Remaining keywords are passed into the `Quantity` constructor
    """
    def __new__(cls, center, std=None, var=None, ivar=None, n_samples=1000, **kwargs):
        if var is not None:
            if std is None:
                std = var**0.5
            else:
                raise ValueError('NormalDistribution cannot take both std and var')
        if ivar is not None:
            if std is None:
                std = ivar**-0.5
            else:
                raise ValueError('NormalDistribution cannot take both ivar and '
                                 'and std or var')
        if std is None:
            raise ValueError('NormalDistribution requires one of std, var, or ivar')

        randshape = [n_samples] + list(np.broadcast(std, center).shape)
        distr = np.random.randn(*randshape)*std + center
        if distr.unit != center.unit:
            distr = distr.to(center.unit)

        self = super(NormalDistribution, cls).__new__(cls, distr, unit=None,
                                                      **kwargs)
        self.distr_std = std
        self.distr_center = center
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

    Remaining keywords are passed into the `Quantity` constructor
    """
    def __new__(cls, poissonval, n_samples=1000, **kwargs):
        randshape = [n_samples] + list(poissonval.shape)
        distr = np.random.poisson(poissonval.value, randshape)

        self = super(PoissonDistribution, cls).__new__(cls, distr,
                                                       unit=poissonval.unit,
                                                       **kwargs)
        self.distr_std = poissonval**0.5
        return self

class UniformDistribution(Distribution):
    """
    A Uniform Distribution.

    Parameters
    ----------
    lower : `Quantity`
        The lower edge of this distribution.
    upper : `Quantity`
        The upper edge of this distribution. Must match shape and be
        compatible units with ``lower``.
    n_samples : int
        The number of Monte Carlo samples to use with this distribution

    Remaining keywords are passed into the `Quantity` constructor.
    """
    def __new__(cls, lower, upper, n_samples=1000, **kwargs):
        if not hasattr(lower, 'unit'):
            raise TypeError('Inputs must be Quantity')
        if lower.unit != upper.unit:
            upper = upper.to(lower.unit)

        newshape = [n_samples] + list(lower.shape)
        distr = np.random.uniform(lower, upper, newshape)

        self = super(UniformDistribution, cls).__new__(cls, distr,
                                                       unit=lower.unit,
                                                       **kwargs)
        return self

    @classmethod
    def from_center_width(cls, centerq, width, unit=None, **kwargs):
        """
        Create a `UniformDistribution` from lower/upper bounds (instead of center
        and width as the regular constructor uses).

        Parameters
        ----------
        centerq : `Quantity`
            The center value of this `UniformDistribution`.
        width : `Quantity`
            The width of this `UniformDistribution`.  Must have the same shape and
            compatible units with ``centerq``.

        Remaining keywords are passed into the `UniformDistribution` constructor.
        """
        whalf = width/2
        return UniformDistribution(centerq - whalf, centerq + whalf, **kwargs)
