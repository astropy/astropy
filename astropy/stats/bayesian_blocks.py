# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Bayesian Blocks for Time Series Analysis
========================================

Dynamic programming algorithm for solving a piecewise-constant model for
various datasets. This is based on the algorithm presented in Scargle
et al 2012 [1]_. This code was ported from the astroML project [2]_.

Applications include:

- finding an optimial histogram with adaptive bin widths
- finding optimal segmentation of time series data
- detecting inflection points in the rate of event data

The primary interface to these routines is the :func:`bayesian_blocks`
function. This module provides fitness functions suitable for three types
of data:

- Irregularly-spaced event data via the :class:`Events` class
- Regularly-spaced event data via the :class:`RegularEvents` class
- Irregularly-spaced point measurements via the :class:`PointMeasures` class

For more fine-tuned control over the fitness functions used, it is possible
to define custom :class:`FitnessFunc` classes directly and use them with
the :func:`bayesian_blocks` routine.

One common application of the Bayesian Blocks algorithm is the determination
of optimal adaptive-width histogram bins. This uses the same fitness function
as for irregularly-spaced time series events. The easiest interface for
creating Bayesian Blocks histograms is the :func:`astropy.stats.histogram`
function.

References
----------
.. [1] http://adsabs.harvard.edu/abs/2012arXiv1207.5578S
.. [2] http://astroML.org/ http://github.com/astroML/astroML/
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
import inspect

# TODO: implement other fitness functions from appendix B of Scargle 2012

__all__ = ['FitnessFunc', 'Events', 'RegularEvents', 'PointMeasures',
           'bayesian_blocks']


def bayesian_blocks(t, x=None, sigma=None,
                    fitness='events', **kwargs):
    """Bayesian Blocks Implementation

    This is a flexible implementation of the Bayesian Blocks algorithm
    described in Scargle 2012 [1]_.

    Parameters
    ----------
    t : array_like
        data times (one dimensional, length N)
    x : array_like (optional)
        data values
    sigma : array_like or float (optional)
        data errors
    fitness : str or object
        the fitness function to use for the model.
        If a string, the following options are supported:

        - 'events' : binned or unbinned event data
            extra arguments are ``p0``, which gives the false alarm probability
            to compute the prior, or ``gamma`` which gives the slope of the
            prior on the number of bins.
        - 'regular_events' : non-overlapping events measured at multiples
            of a fundamental tick rate, ``dt``, which must be specified as an
            additional argument.  The prior can be specified through ``gamma``,
            which gives the slope of the prior on the number of bins.
        - 'measures' : fitness for a measured sequence with Gaussian errors
            The prior can be specified using ``gamma``, which gives the slope
            of the prior on the number of bins.  If ``gamma`` is not specified,
            then a simulation-derived prior will be used.

        Alternatively, the fitness parameter can be an instance of
        :class:`FitnessFunc` or a subclass thereof.

    **kwargs :
        any additional keyword arguments will be passed to the specified
        :class:`FitnessFunc` derived class.

    Returns
    -------
    edges : ndarray
        array containing the (N+1) edges defining the N bins

    Examples
    --------
    Event data:

    >>> t = np.random.normal(size=100)
    >>> edges = bayesian_blocks(t, fitness='events', p0=0.01)

    Event data with repeats:

    >>> t = np.random.normal(size=100)
    >>> t[80:] = t[:20]
    >>> edges = bayesian_blocks(t, fitness='events', p0=0.01)

    Regular event data:

    >>> dt = 0.01
    >>> t = dt * np.arange(1000)
    >>> x = np.zeros(len(t))
    >>> x[np.random.randint(0, len(t), len(t) // 10)] = 1
    >>> edges = bayesian_blocks(t, x, fitness='regular_events', dt=dt)

    Measured point data with errors:

    >>> t = 100 * np.random.random(100)
    >>> x = np.exp(-0.5 * (t - 50) ** 2)
    >>> sigma = 0.1
    >>> x_obs = np.random.normal(x, sigma)
    >>> edges = bayesian_blocks(t, x_obs, sigma, fitness='measures')

    References
    ----------
    .. [1] Scargle, J et al. (2012)
       http://adsabs.harvard.edu/abs/2012arXiv1207.5578S

    See Also
    --------
    astropy.stats.histogram : compute a histogram using bayesian blocks
    """
    FITNESS_DICT = {'events': Events,
                    'regular_events': RegularEvents,
                    'measures': PointMeasures}
    fitness = FITNESS_DICT.get(fitness, fitness)

    if type(fitness) is type and issubclass(fitness, FitnessFunc):
        fitfunc = fitness(**kwargs)
    elif isinstance(fitness, FitnessFunc):
        fitfunc = fitness
    else:
        raise ValueError("fitness parameter not understood")

    return fitfunc.fit(t, x, sigma)


class FitnessFunc(object):
    """Base class for fitness functions

    Derived classes should overload the following method:

    ``fitness(self, **kwargs)``:
      Compute the fitness given a set of named arguments.
      Arguments accepted by fitness must be among ``[T_k, N_k, a_k, b_k, c_k]``
      (See Scargle2012_ for details on the meaning of these parameters).

    Additionally, other methods may be overloaded as well:

    ``__init__(self, **kwargs)``:
      Initialize the fitness function with any parameters beyond the normal
      ``p0`` and ``gamma``.

    ``validate_input(self, t, x, sigma)``:
      Enable specific checks of the input data (``t``, ``x``, ``sigma``)
      to be performed prior to the fit.

    ``prior(self, N, Ntot)``:
      Specify the form of the Bayesian prior on the number of bins ``N`` given
      the total number of points ``Ntot``. By default it will call
      the ``gamma_prior()`` method if ``gamma`` is specified, and otherwise
      call the ``p0_prior()`` method.

    ``p0_prior(self, N, Ntot)``:
      Specify the form of the prior given the false-alarm probability ``p0``
      (See Scargle2012_ for details).

    ``gamma_prior(self, N, Ntot)``:
      Specify the form of the prior given the slope ``gamma`` of the prior
      on the number of bins (See Scargle2012_ for details).

    For examples of implemented fitness functions, see :class:`Events`,
    :class:`RegularEvents`, and :class:`PointMeasures`.

    References
    ----------
    .. [Scargle2012] Scargle, J et al. (2012)
       http://adsabs.harvard.edu/abs/2012arXiv1207.5578S
    """
    def __init__(self, p0=0.05, gamma=None):
        self.p0 = p0
        self.gamma = gamma

    def validate_input(self, t, x=None, sigma=None):
        """Validate inputs to the model.

        Parameters
        ----------
        t : array_like
            times of observations
        x : array_like (optional)
            values observed at each time
        sigma : float or array_like (optional)
            errors in values x

        Returns
        -------
        t, x, sigma : array_like, float or None
            validated and perhaps modified versions of inputs
        """
        # validate array input
        t = np.asarray(t, dtype=float)
        if x is not None:
            x = np.asarray(x)
        if sigma is not None:
            sigma = np.asarray(sigma)

        # find unique values of t
        t = np.array(t)
        if t.ndim != 1:
            raise ValueError("t must be a one-dimensional array")
        unq_t, unq_ind, unq_inv = np.unique(t, return_index=True,
                                            return_inverse=True)

        # if x is not specified, x will be counts at each time
        if x is None:
            if sigma is not None:
                raise ValueError("If sigma is specified, x must be specified")
            else:
                sigma = 1

            if len(unq_t) == len(t):
                x = np.ones_like(t)
            else:
                x = np.bincount(unq_inv)

            t = unq_t

        # if x is specified, then we need to simultaneously sort t and x
        else:
            # TODO: allow broadcasted x?
            x = np.asarray(x)
            if x.shape not in [(), (1,), (t.size,)]:
                raise ValueError("x does not match shape of t")
            x += np.zeros_like(t)

            if len(unq_t) != len(t):
                raise ValueError("Repeated values in t not supported when "
                                 "x is specified")
            t = unq_t
            x = x[unq_ind]

        # verify the given sigma value
        if sigma is None:
            sigma = 1
        else:
            sigma = np.asarray(sigma)
            if sigma.shape not in [(), (1,), (t.size,)]:
                raise ValueError('sigma does not match the shape of x')

        return t, x, sigma

    def fitness(self, **kwargs):
        raise NotImplementedError()

    def prior(self, N, Ntot):
        """
        Compute the prior on the number of blocks ``N``,
        given the maximum possible blocks ``Ntot``.
        """
        if self.gamma is None:
            return self.p0_prior(N, Ntot)
        else:
            return self.gamma_prior(N, Ntot)

    def p0_prior(self, N, Ntot):
        """
        Empirical prior, parametrized by the false alarm probability ``p0``
        See  eq. 21 in Scargle (2012)
        """
        return 4 - np.log(73.53 * self.p0 * (N ** -0.478))

    def gamma_prior(self, N, Ntot):
        """
        Basic prior, parametrized by slope gamma.
        See eq. 3 in Scargle (2012)
        """
        if self.gamma == 1:
            return 0
        else:
            return (np.log(1 - self.gamma)
                    - np.log(1 - self.gamma ** (Ntot + 1))
                    + N * np.log(self.gamma))

    # the fitness_args property will return the list of arguments accepted by
    # the method fitness().  This allows more efficient computation below.
    @property
    def _fitness_args(self):
        return inspect.getargspec(self.fitness)[0][1:]

    def fit(self, t, x=None, sigma=None):
        """Fit the Bayesian Blocks model given the specified fitness function.

        Parameters
        ----------
        t : array_like
            data times (one dimensional, length N)
        x : array_like (optional)
            data values
        sigma : array_like or float (optional)
            data errors

        Returns
        -------
        edges : ndarray
            array containing the (M+1) edges defining the M optimal bins
        """
        t, x, sigma = self.validate_input(t, x, sigma)

        # compute values needed for computation, below
        if 'a_k' in self._fitness_args:
            ak_raw = np.ones_like(x) / sigma ** 2
        if 'b_k' in self._fitness_args:
            bk_raw = x / sigma ** 2
        if 'c_k' in self._fitness_args:
            ck_raw = x * x / sigma ** 2

        # create length-(N + 1) array of cell edges
        edges = np.concatenate([t[:1],
                                0.5 * (t[1:] + t[:-1]),
                                t[-1:]])
        block_length = t[-1] - edges

        # arrays to store the best configuration
        N = len(t)
        best = np.zeros(N, dtype=float)
        last = np.zeros(N, dtype=int)

        #-----------------------------------------------------------------
        # Start with first data cell; add one cell at each iteration
        #-----------------------------------------------------------------
        for R in range(N):
            # Compute fit_vec : fitness of putative last block (end at R)
            kwds = {}

            # T_k: width/duration of each block
            if 'T_k' in self._fitness_args:
                kwds['T_k'] = block_length[:R + 1] - block_length[R + 1]

            # N_k: number of elements in each block
            if 'N_k' in self._fitness_args:
                kwds['N_k'] = np.cumsum(x[:R + 1][::-1])[::-1]

            # a_k: eq. 31
            if 'a_k' in self._fitness_args:
                kwds['a_k'] = 0.5 * np.cumsum(ak_raw[:R + 1][::-1])[::-1]

            # b_k: eq. 32
            if 'b_k' in self._fitness_args:
                kwds['b_k'] = - np.cumsum(bk_raw[:R + 1][::-1])[::-1]

            # c_k: eq. 33
            if 'c_k' in self._fitness_args:
                kwds['c_k'] = 0.5 * np.cumsum(ck_raw[:R + 1][::-1])[::-1]

            # evaluate fitness function
            fit_vec = self.fitness(**kwds)

            A_R = fit_vec - self.prior(R + 1, N)
            A_R[1:] += best[:R]

            i_max = np.argmax(A_R)
            last[R] = i_max
            best[R] = A_R[i_max]

        #-----------------------------------------------------------------
        # Now find changepoints by iteratively peeling off the last block
        #-----------------------------------------------------------------
        change_points = np.zeros(N, dtype=int)
        i_cp = N
        ind = N
        while True:
            i_cp -= 1
            change_points[i_cp] = ind
            if ind == 0:
                break
            ind = last[ind - 1]
        change_points = change_points[i_cp:]

        return edges[change_points]


class Events(FitnessFunc):
    """Fitness for binned or unbinned events

    Parameters
    ----------
    p0 : float (optional)
        False alarm probability, used to compute the prior on N
        (see eq. 21 of Scargle 2012).  Default prior is for p0 = 0.
    gamma : float (optional)
        If specified, then use this gamma to compute the general prior form,
        p ~ gamma^N.  If gamma is specified, p0 is ignored.
    """
    def fitness(self, N_k, T_k):
        # eq. 19 from Scargle 2012
        return N_k * (np.log(N_k) - np.log(T_k))

    def validate_input(self, t, x, sigma):
        t, x, sigma = super(Events, self).validate_input(t, x, sigma)
        if x is not None and np.any(x % 1 > 0):
            raise ValueError("x must be integer counts for fitness='events'")
        return t, x, sigma

    def p0_prior(self, N, Ntot):
        # eq. 21 from Scargle 2012
        return 4 - np.log(73.53 * self.p0 * (N ** -0.478))


class RegularEvents(FitnessFunc):
    """Fitness for regular events

    This is for data which has a fundamental "tick" length, so that all
    measured values are multiples of this tick length.  In each tick, there
    are either zero or one counts.

    Parameters
    ----------
    dt : float
        tick rate for data
    p0 : float (optional)
        False alarm probability, used to compute the prior on N
        (see eq. 21 of Scargle 2012).  Default prior is for p0 = 0.
    gamma : float (optional)
        If specified, then use this gamma to compute the general prior form,
        p ~ gamma^N.  If gamma is specified, p0 is ignored.
    """
    def __init__(self, dt, p0=0.05, gamma=None):
        self.dt = dt
        self.p0 = p0
        self.gamma = gamma

    def validate_input(self, t, x, sigma):
        t, x, sigma = super(RegularEvents, self).validate_input(t, x, sigma)
        if not np.all((x == 0) | (x == 1)):
            raise ValueError("Regular events must have only 0 and 1 in x")
        return t, x, sigma

    def fitness(self, T_k, N_k):
        # Eq. 75 of Scargle 2012
        M_k = T_k / self.dt
        N_over_M = N_k / M_k

        eps = 1E-8
        if np.any(N_over_M > 1 + eps):
            import warnings
            warnings.warn('regular events: N/M > 1.  '
                          'Is the time step correct?')

        one_m_NM = 1 - N_over_M
        N_over_M[N_over_M <= 0] = 1
        one_m_NM[one_m_NM <= 0] = 1

        return N_k * np.log(N_over_M) + (M_k - N_k) * np.log(one_m_NM)


class PointMeasures(FitnessFunc):
    """Fitness for point measures

    Parameters
    ----------
    p0 : float (optional)
        False alarm probability, used to compute the prior on N
        (see eq. 21 of Scargle 2012).  Default prior is for p0 = 0.
    gamma : float (optional)
        If specified, then use this gamma to compute the general prior form,
        p ~ gamma^N.  If gamma is specified, p0 is ignored.

    Notes
    -----
    if neither p0 nor gamma is specified, the prior will be computed using
    the empirical form given in section 3.3 of Scargle (2012).
    """
    def __init__(self, p0=None, gamma=None):
        self.p0 = p0
        self.gamma = gamma

    def fitness(self, a_k, b_k):
        # eq. 41 from Scargle 2012
        return (b_k * b_k) / (4 * a_k)

    def validate_input(self, t, x, sigma):
        if x is None:
            raise ValueError("x must be specified for point measures")
        return super(PointMeasures, self).validate_input(t, x, sigma)

    def prior(self, N, Ntot):
        if self.gamma is not None:
            return self.gamma_prior(N, Ntot)
        elif self.p0 is not None:
            return self.p0_prior(N, Ntot)
        else:
            # eq. at end of sec 3.3 in Scargle 2012
            return 1.32 + 0.577 * np.log10(N)
