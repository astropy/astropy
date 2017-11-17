# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Distribution class and associated machinery.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from . import Quantity
from .. import visualization

__all__ = ['Distribution', 'NormalDistribution', 'PoissonDistribution',
           'UniformDistribution']


class Distribution:
    """
    A unitful value with associated uncertainty distribution.  While subclasses
    may have exact analytic forms, this class uses samples to represent the
    distribution.

    Parameters
    ----------
    distr : array-like
        The distribution, with sampling along the *leading* axis.
    distr_center : None or array-like
        The "center" of this distribution.  It must have the same shape as
        ``distr``, aside from the (initial) sample dimension. If None, the
        *median* of ``distr``  along the sampling axis will be used.
    unit : astropy unit
        A unit of the same sort that `~astropy.units.Quantity` accepts.

    Additional argumnets are passed into the `~astropy.units.Quantity`
    constructor.

    Notes
    -----
    This class has an attribute ``stat_view`` which can be set to an
    `numpy.ndarray` subclass to determine what the summary statistics are
    converted to.  By default this is `~astropy.units.Quantity`, but it is
    available so that subclasses can change it.
    """

    def __new__(cls, distr, distr_center=None, unit=None, *args, **kwargs):

        if distr.shape == ():
            raise TypeError('Attempted to initialize a Distribution with a scalar')

        # Not clear whether one wants the distr_center argument still.
        # Would need __array_finalize__ to preserve it through slicing, etc.
        if distr_center is None:
            distr_center = np.median(distr, axis=-1)

        # Not clear why the unit argument is needed.
        if unit is not None:
            distr = Quantity(distr, unit, copy=False, subok=True)

        new_dtype = np.dtype({'names': ['samples'],
                              'formats': [(distr.dtype, (distr.shape[-1],))]})
        distr_cls = type(distr)
        if not issubclass(distr_cls, Distribution):
            # Probably should have a dict on the class with these, so
            # we don't create new classes needlessly.
            new_name = distr_cls.__name__ + cls.__name__
            new_cls = type(new_name, (cls, distr_cls),
                           {'_distr_cls': distr_cls})
        self = distr.view(dtype=new_dtype, type=new_cls)
        # Get rid of trailing dimension of 1.
        self.shape = distr.shape[:-1]
        self.center = distr_center
        return self

    @property
    def distribution(self):
        return self['samples']

    def __getitem__(self, item):
        result = super().__getitem__(item)
        if item == 'samples':
            return super(Distribution, result).view(self._distr_cls)
        else:
            return result

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        converted = []
        outputs = kwargs.pop('out', None)
        if outputs:
            kwargs['out'] = tuple((output.distribution if
                                   isinstance(output, Distribution)
                                   else output) for output in outputs)
        if method in {'reduce', 'accumulate', 'reduceat'}:
            axis = kwargs.get('axis', None)
            if axis is None:
                assert isinstance(inputs[0], Distribution)
                kwargs['axis'] = tuple(range(inputs[0].ndim))

        for input_ in inputs:
            if isinstance(input_, Distribution):
                converted.append(input_.distribution)
            else:
                shape = getattr(input_, 'shape', ())
                if shape:
                    converted.append(input_[..., np.newaxis])
                else:
                    converted.append(input_)

        results = getattr(ufunc, method)(*converted, **kwargs)

        if not isinstance(results, tuple):
            results = (results,)
        if outputs is None:
            outputs = (None,) * len(results)

        finals = []
        for result, output in zip(results, outputs):
            if output is not None:
                finals.append(output)
            else:
                if getattr(result, 'shape', False):
                    finals.append(Distribution(result))
                else:
                    finals.append(result)

        return finals if len(finals) > 1 else finals[0]

    # Override view so that we stay a Distribution version of the new type.
    def view(self, dtype=None, type=None):
        if type is None:
            if issubclass(dtype, np.ndarray):
                type = dtype
                dtype = None
            else:
                raise ValueError('Cannot set just dtype for a Distribution.')

        result = self.distribution.view(dtype, type)
        return Distribution(result)

    def __repr__(self):
        superrepr = super().__repr__()
        toadd = ' with n_samples={}'.format(self.n_samples)
        return superrepr[:-1] + toadd + superrepr[-1]

    def __str__(self):
        superstr = super().__str__()
        toadd = ' with n_samples={}'.format(self.n_samples)
        return superstr + toadd

    def _repr_latex_(self):
        superlatex = super()._repr_latex_()
        toadd = r', \; n_{{\rm samp}}={}'.format(self.n_samples)
        return superlatex[:-1] + toadd + superlatex[-1]

    @property
    def n_samples(self):
        """
        The number of samples of this distribution.  A single int.
        """
        return self.dtype['samples'].shape[0]

    @property
    def pdf_mean(self):
        """
        The mean of this distribution.
        """
        return self.distribution.mean(axis=-1)

    @property
    def pdf_std(self):
        """
        The standard deviation of this distribution.
        """
        return self.distribution.std(axis=-1)

    @property
    def pdf_var(self):
        """
        The variance of this distribution.
        """
        return self.distribution.var(axis=-1)

    @property
    def pdf_median(self):
        """
        The median of this distribution.
        """
        return np.median(self.distribution, axis=-1)

    @property
    def pdf_mad(self):
        """
        The median absolute deviation of this distribution.
        """
        return np.abs(self - self.pdf_median).pdf_median

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
        percs : Quantity
            The ``fracs`` percentiles of this distribution.
        """
        perc = np.percentile(self.distribution, perc, axis=-1)
        # numpy.percentile strips units for unclear reasons, so we have to make
        # a new object with units
        if hasattr(self.distribution, '_new_view'):
            return self.distribution._new_view(perc)
        else:
            return perc

    def hist(self, maxtoshow=10, **kwargs):
        """
        Use `astropy.visualization.hist` (which uses matplotlib's ``hist``
        function) to visualize this distribution.  For N-D distributions, the array
        is flattened following standard numpy rules, and the distributions are shown
        as separate histograms for each element.

        Parameters
        ----------
        maxtoshow : int or None
            The maximum number of distribution elements to show.  If None, there
            will be no limit, but this may overwhelm matplotlib if the distribution
            is large.
        Additional keywords are passed into `astropy.visualization.hist`
        """
        # NOTE: not adjusted to new structure!!
        labelset = 'label' in kwargs
        scalar_distr = len(self.shape) == 1
        if scalar_distr:
            reshaped = [self]
        else:
            reshaped = self.reshape(self.n_samples, self.size//self.n_samples).T

        hists = []
        for i, dat in enumerate(reshaped):
            if i >= maxtoshow:
                break
            if (not scalar_distr) and (not labelset):
                if len(self.shape) == 2:
                    idxstr = str(i)
                else:
                    idxstr = str(np.unravel_index(i, self.shape[1:]))
                kwargs['label'] = 'Distr ' + idxstr
            hists.append(visualization.hist(dat, **kwargs))
        return hists


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

        randshape = np.broadcast(std, center).shape + (n_samples,)
        distr = np.random.randn(*randshape)
        self = super(NormalDistribution, cls).__new__(cls, distr, **kwargs)
        self = center + self * std
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

    Remaining keywords are passed into the `Quantity` constructor
    """
    def __new__(cls, poissonval, n_samples=1000, **kwargs):
        randshape = poissonval.shape + (n_samples,)
        distr = np.random.poisson(poissonval.value[..., np.newaxis], randshape)

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

        newshape = lower.shape + (n_samples,)
        distr = np.random.uniform(lower[..., np.newaxis],
                                  upper[..., np.newaxis], newshape)

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
