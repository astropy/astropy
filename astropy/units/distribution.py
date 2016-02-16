# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Distribution class and associated machinery.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ..extern import six

import numpy as np

from . import Quantity, UnitsError

__all__ = ['Distribution']

class Distribution(Quantity):
    """
    A unitful value with associated uncertainty distribution.
    """
    
    __array_priority__ = Quantity.__array_priority__ + 1
    
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
        return self.mean(axis=0)
        
    @property
    def pdf_std(self):
        """
        The standard deviation of this distribution.
        """
        return self.std(axis=0)
        
    @property
    def pdf_var(self):
        """
        The variance of this distribution.
        """
        return self.var(axis=0)
        
    @property
    def pdf_median(self):
        """
        The median of this distribution.
        """
        return np.median(self, axis=0)
        
    @property
    def pdf_mad(self):
        """
        The median absolute deviation of this distribution.
        """
        return np.abs(self - self.pdef_median)
    
    # TODO: decide how to best compute this exactly - it's 1/(phi^-1(0.75)), where phi is the inverse CDF of the normal
    _smad_scale_factor = 1.4826
    
    @property
    def pdf_smad(self):
        """
        The median absolute deviation of this distribution rescaled to match the
        standard deviation for a normal distribution.
        """
        return self.pdf_mad * self._smad_scale_factor
        
    def percentiles(self, fracs):
        """
        Compute percentiles of this Distribution.
        
        Parameters
        ----------
        fracs : the desired  precentiles of the distribution as fractions 
        (i.e., on [0,1]).
        
        Returns
        -------
        percs : Quantitiy of shape ``distr_shape``
            The ``fracs`` percentiles of this distribution.
        """
        raise NotImplementedError
        
    