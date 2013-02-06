# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""General purpose timer related functions."""

from __future__ import absolute_import, division

# STDLIB
import time
from collections import Iterable
from functools import partial

# THIRD-PARTY
import numpy as np

# LOCAL
from . import OrderedDict
from ..logger import log


__all__ = ['timeit', 'SimpleRunTimePredictor']


def timeit(num_tries=1, verbose=True):
    """Decorator to time a function or method."""
    def real_decorator(function):
        def wrapper(*args, **kwargs):
            ts = time.time()
            for i in xrange(num_tries):
                result = function(*args, **kwargs)
            te = time.time()
            tt = (te - ts) / num_tries
            if verbose:
                log.info('{0} took {1} s on AVERAGE for {2} calls.'.format(function.__name__, tt, num_tries))
            return tt, result
        return wrapper
    return real_decorator


class SimpleRunTimePredictor(object):
    """Class to predict run time.

    .. note: Only predict for single varying numeric input parameter.

    Parameters
    ----------
    func : function
        Function to time.

    args : tuple of numbers
        Fixed positional argument(s) for the function.

    kwargs : dict
        Fixed keyword argument(s) for the function.

    """
    def __init__(self, func, *args, **kwargs):
        self._pfunc = partial(func, *args, **kwargs)
        self._cache_good = OrderedDict()
        self._cache_bad = []
        self._cache_est = OrderedDict()
        self._fit_func = None
        self._power = None

    @timeit(num_tries=1, verbose=False)
    def _timed_pfunc(self, arg):
        """Run partial func once for single arg and time it."""
        return self._pfunc(arg)

    def _cache_time(self, arg):
        """Cache timing results without repetition."""
        if arg not in self._cache_good and arg not in self._cache_bad:
            try:
                result = self._timed_pfunc(arg)
            except Exception as e:
                log.warn(e.message)
                self._cache_bad.append(arg)
            else:
                self._cache_good[arg] = result[0]  # Only run time is kept

    def time_func(self, arglist):
        """Time the partial function for a list of single args
        and store results in a cache. This forms a baseline for
        the prediction.

        Parameters
        ----------
        arglist : list of numbers
            List of input arguments to time.

        """
        if not isinstance(arglist, Iterable):
            arglist = [arglist]
        dummy = map(self._cache_time, arglist)

    def do_fit(self, power=1, deg=1, min_datapoints=3):
        """Fit a function to the lists of arguments and
        their respective run time in the cache.

        Function::

            t = a[0] + a[1] * arg**power + ... + a[deg] * (arg**power)**deg

        Parameters
        ----------
        power : int, optional
            Power of values to fit.

        deg : int, optional
            Degree of polynomial to fit.

        min_datapoints : int, optional
            Minimum number of data points required for fitting.
            They can be built up with `time_func`.

        Returns
        -------
        a : array_like
            Fitted coefficients from `numpy.polyfit`.

        """
        # Reset related attributes
        self._power = power
        self._cache_est = OrderedDict()

        x_arr = np.array(self._cache_good.keys())
        assert x_arr.size >= min_datapoints, \
            'Requires {0} points but has {1}'.format(min_datapoints, x_arr.size)

        a = np.polyfit(x_arr**power, self._cache_good.values(), deg)
        self._fit_func = np.poly1d(a)

        return a

    def predict_time(self, arg):
        """Predict run time for given argument.
        If prediction is already cached, cached value is returned.

        Parameters
        ----------
        arg : number
            Input argument to predict run time for.

        Returns
        -------
        t_est : float
            Estimated run time for given argument.

        """
        if arg in self._cache_est:
            t_est = self._cache_est[arg]
        else:
            assert self._fit_func is not None, 'No fitted data for prediction'
            t_est = self._fit_func(arg**self._power)
            self._cache_est[arg] = t_est
        return t_est

    def plot(self, xlabeltext='args', save_as=''):
        """Plot prediction.

        .. note:: Uses :mod:`matplotlib`.

        Parameters
        ----------
        xlabeltext : str, optional
            Text for X-label.

        save_as : str, optional
            Save plot as given filename.

        """
        import matplotlib.pyplot as plt

        # Actual data
        x_arr = sorted(self._cache_good)
        y_arr = [self._cache_good[x] for x in x_arr]

        assert len(x_arr) > 1, 'Insufficient data for plotting'

        fig, ax = plt.subplots()
        ax.plot(x_arr, y_arr, 'kx-', label='Actual')

        # Fitted data
        if self._fit_func is not None:
            x_est = self._cache_est.keys()
            ax.scatter(x_est, self._cache_est.values(), marker='o', c='r', label='Fitted')

            x_fit = np.array(sorted(x_arr + x_est))
            y_fit = self._fit_func(x_fit**self._power)
            ax.plot(x_fit, y_fit, 'b--', label='Fit')

        ax.set_xlabel(xlabeltext)
        ax.set_ylabel('t (s)')
        ax.legend(loc='best', numpoints=1)

        plt.draw()

        if save_as:  # pragma: no cover
            plt.savefig(save_as)
