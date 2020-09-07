# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Distribution class and associated machinery.
"""

import numpy as np

from astropy import units as u
from astropy import stats

__all__ = ['Distribution']


# we set this by hand because the symbolic expression (below) requires scipy
# SMAD_SCALE_FACTOR = 1 / scipy.stats.norm.ppf(0.75)
SMAD_SCALE_FACTOR = 1.48260221850560203193936104071326553821563720703125


class Distribution:
    """
    A scalar value or array values with associated uncertainty distribution.

    This object will take its exact type from whatever the ``samples`` argument
    is. In general this is expected to be an `~astropy.units.Quantity` or
    `numpy.ndarray`, although anything compatible with `numpy.asanyarray` is
    possible.

    See also: https://docs.astropy.org/en/stable/uncertainty/

    Parameters
    ----------
    samples : array_like
        The distribution, with sampling along the *leading* axis. If 1D, the
        sole dimension is used as the sampling axis (i.e., it is a scalar
        distribution).
    """
    _generated_subclasses = {}

    def __new__(cls, samples):
        if isinstance(samples, Distribution):
            samples = samples.distribution
        else:
            samples = np.asanyarray(samples, order='C')
        if samples.shape == ():
            raise TypeError('Attempted to initialize a Distribution with a scalar')

        new_dtype = np.dtype({'names': ['samples'],
                              'formats': [(samples.dtype, (samples.shape[-1],))]})
        samples_cls = type(samples)
        new_cls = cls._generated_subclasses.get(samples_cls)
        if new_cls is None:
            # Make a new class with the combined name, inserting Distribution
            # itself below the samples class since that way Quantity methods
            # like ".to" just work (as .view() gets intercepted).  However,
            # repr and str are problems, so we put those on top.
            # TODO: try to deal with this at the lower level.  The problem is
            # that array2string does not allow one to override how structured
            # arrays are typeset, leading to all samples to be shown.  It may
            # be possible to hack oneself out by temporarily becoming a void.
            new_name = samples_cls.__name__ + cls.__name__
            new_cls = type(
                new_name,
                (_DistributionRepr, samples_cls, ArrayDistribution),
                {'_samples_cls': samples_cls})
            cls._generated_subclasses[samples_cls] = new_cls

        self = samples.view(dtype=new_dtype, type=new_cls)
        # Get rid of trailing dimension of 1.
        self.shape = samples.shape[:-1]
        return self

    @property
    def distribution(self):
        return self['samples']

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

    @property
    def n_samples(self):
        """
        The number of samples of this distribution.  A single `int`.
        """
        return self.dtype['samples'].shape[0]

    def pdf_mean(self, dtype=None, out=None):
        """
        The mean of this distribution.

        Arguments are as for `numpy.mean`.
        """
        return self.distribution.mean(axis=-1, dtype=dtype, out=out)

    def pdf_std(self, dtype=None, out=None, ddof=0):
        """
        The standard deviation of this distribution.

        Arguments are as for `numpy.std`.
        """
        return self.distribution.std(axis=-1, dtype=dtype, out=out, ddof=ddof)

    def pdf_var(self, dtype=None, out=None, ddof=0):
        """
        The variance of this distribution.

        Arguments are as for `numpy.var`.
        """
        return self.distribution.var(axis=-1, dtype=dtype, out=out, ddof=ddof)

    def pdf_median(self, out=None):
        """
        The median of this distribution.

        Parameters
        ----------
        out : array, optional
            Alternative output array in which to place the result. It must
            have the same shape and buffer length as the expected output,
            but the type (of the output) will be cast if necessary.
        """
        return np.median(self.distribution, axis=-1, out=out)

    def pdf_mad(self, out=None):
        """
        The median absolute deviation of this distribution.

        Parameters
        ----------
        out : array, optional
            Alternative output array in which to place the result. It must
            have the same shape and buffer length as the expected output,
            but the type (of the output) will be cast if necessary.
        """
        median = self.pdf_median(out=out)
        absdiff = np.abs(self - median)
        return np.median(absdiff.distribution, axis=-1, out=median,
                         overwrite_input=True)

    def pdf_smad(self, out=None):
        """
        The median absolute deviation of this distribution rescaled to match the
        standard deviation for a normal distribution.

        Parameters
        ----------
        out : array, optional
            Alternative output array in which to place the result. It must
            have the same shape and buffer length as the expected output,
            but the type (of the output) will be cast if necessary.
        """
        result = self.pdf_mad(out=out)
        result *= SMAD_SCALE_FACTOR
        return result

    def pdf_percentiles(self, percentile, **kwargs):
        """
        Compute percentiles of this Distribution.

        Parameters
        ----------
        percentile : float or array of floats or `~astropy.units.Quantity`
            The desired  precentiles of the distribution (i.e., on [0,100]).
            `~astropy.units.Quantity` will be converted to percent, meaning
            that a ``dimensionless_unscaled`` `~astropy.units.Quantity` will
            be interpreted as a quantile.

        Additional keywords are passed into `numpy.percentile`.

        Returns
        -------
        percentiles : `~astropy.units.Quantity`
            The ``fracs`` percentiles of this distribution.
        """
        percentile = u.Quantity(percentile, u.percent).value
        percs = np.percentile(self.distribution, percentile, axis=-1, **kwargs)
        # numpy.percentile strips units for unclear reasons, so we have to make
        # a new object with units
        if hasattr(self.distribution, '_new_view'):
            return self.distribution._new_view(percs)
        else:
            return percs

    def pdf_histogram(self, **kwargs):
        """
        Compute histogram over the samples in the distribution.

        Parameters
        ----------
        All keyword arguments are passed into `astropy.stats.histogram`. Note
        That some of these options may not be valid for some multidimensional
        distributions.

        Returns
        -------
        hist : array
            The values of the histogram. Trailing dimension is the histogram
            dimension.
        bin_edges : array of dtype float
            Return the bin edges ``(length(hist)+1)``. Trailing dimension is the
            bin histogram dimension.
        """
        distr = self.distribution
        raveled_distr = distr.reshape(distr.size//distr.shape[-1], distr.shape[-1])

        nhists = []
        bin_edges = []
        for d in raveled_distr:
            nhist, bin_edge = stats.histogram(d, **kwargs)
            nhists.append(nhist)
            bin_edges.append(bin_edge)

        nhists = np.array(nhists)
        nh_shape = self.shape + (nhists.size//self.size,)
        bin_edges = np.array(bin_edges)
        be_shape = self.shape + (bin_edges.size//self.size,)
        return nhists.reshape(nh_shape), bin_edges.reshape(be_shape)


class ScalarDistribution(Distribution, np.void):
    """Scalar distribution.

    This class mostly exists to make `~numpy.array2print` possible for
    all subclasses.  It is a scalar element, still with n_samples samples.
    """
    pass


class ArrayDistribution(Distribution, np.ndarray):
    # This includes the important overrride of view and __getitem__
    # which are needed for all ndarray subclass Distributions, but not
    # for the scalar one.
    _samples_cls = np.ndarray

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

    # Override __getitem__ so that 'samples' is returned as the sample class.
    def __getitem__(self, item):
        result = super().__getitem__(item)
        if item == 'samples':
            # Here, we need to avoid our own redefinition of view.
            return super(ArrayDistribution, result).view(self._samples_cls)
        elif isinstance(result, np.void):
            return result.view((ScalarDistribution, result.dtype))
        else:
            return result


class _DistributionRepr:
    def __repr__(self):
        reprarr = repr(self.distribution)
        if reprarr.endswith('>'):
            firstspace = reprarr.find(' ')
            reprarr = reprarr[firstspace+1:-1]  # :-1] removes the ending '>'
            return '<{} {} with n_samples={}>'.format(self.__class__.__name__,
                                                      reprarr, self.n_samples)
        else:  # numpy array-like
            firstparen = reprarr.find('(')
            reprarr = reprarr[firstparen:]
            return '{}{} with n_samples={}'.format(self.__class__.__name__,
                                                   reprarr, self.n_samples)
            return reprarr

    def __str__(self):
        distrstr = str(self.distribution)
        toadd = f' with n_samples={self.n_samples}'
        return distrstr + toadd

    def _repr_latex_(self):
        if hasattr(self.distribution, '_repr_latex_'):
            superlatex = self.distribution._repr_latex_()
            toadd = fr', \; n_{{\rm samp}}={self.n_samples}'
            return superlatex[:-1] + toadd + superlatex[-1]
        else:
            return None


class NdarrayDistribution(_DistributionRepr, ArrayDistribution):
    pass


# Ensure our base NdarrayDistribution is known.
Distribution._generated_subclasses[np.ndarray] = NdarrayDistribution
