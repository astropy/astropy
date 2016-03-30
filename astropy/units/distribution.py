# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Distribution class and associated machinery.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from . import Quantity

__all__ = ['Distribution']


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
        percs : Quantitiy of shape ``distr_shape``
            The ``fracs`` percentiles of this distribution.
        """
        perc = np.percentile(self, perc, axis=0)
        # numpy.percentile strips units for unclear reasons, so we have to make
        # a new object with units
        return self.stat_view(perc, unit=self.unit, copy=False)
