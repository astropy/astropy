# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""General purpose timer related functions."""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..extern import six

# STDLIB
import time
import warnings
from collections import Iterable
from functools import partial, wraps

# THIRD-PARTY
import numpy as np

# LOCAL
from . import OrderedDict
from .. import units as u
from .. import log
from .. import modeling
from .exceptions import AstropyUserWarning


__all__ = ['timefunc', 'RunTimePredictor']
__doctest_skip__ = ['timefunc']


def timefunc(num_tries=1, verbose=True):
    """Decorator to time a function or method.

    Parameters
    ----------
    num_tries : int, optional
        Number of calls to make. Timer will take the
        average run time.

    verbose : bool, optional
        Extra log information.

    Returns
    -------
    tt : float
        Average run time in seconds.

    result
        Output(s) from the function.

    Examples
    --------
    To add timer to time `numpy.log` for 100 times with
    verbose output::

        import numpy as np
        from astropy.utils.timer import timefunc

        @timefunc(100)
        def timed_log(x):
            return np.log(x)

    To run the decorated function above:

    >>> t, y = timed_log(100)
    INFO: timed_log took 9.29832458496e-06 s on AVERAGE for 100 call(s). [...]
    >>> t
    9.298324584960938e-06
    >>> y
    4.6051701859880918

    """
    def real_decorator(function):
        @wraps(function)
        def wrapper(*args, **kwargs):
            ts = time.time()
            for i in range(num_tries):
                result = function(*args, **kwargs)
            te = time.time()
            tt = (te - ts) / num_tries
            if verbose:  # pragma: no cover
                log.info('{0} took {1} s on AVERAGE for {2} call(s).'.format(
                    function.__name__, tt, num_tries))
            return tt, result
        return wrapper
    return real_decorator


class RunTimePredictor(object):
    """Class to predict run time.

    .. note:: Only predict for single varying numeric input parameter.

    Parameters
    ----------
    func : function
        Function to time.

    args : tuple
        Fixed positional argument(s) for the function.

    kwargs : dict
        Fixed keyword argument(s) for the function.

    Examples
    --------
    >>> from astropy.utils.timer import RunTimePredictor

    Set up a predictor for :math:`10^{x}`:

    >>> p = RunTimePredictor(pow, 10)

    Give it baseline data to use for prediction and
    get the function output values:

    >>> p.time_func(range(10, 1000, 200))
    >>> for input, result in sorted(p.results.items()):
    ...     print("pow(10, {0})\\n{1}".format(input, result))
    pow(10, 10)
    10000000000
    pow(10, 210)
    10000000000...
    pow(10, 410)
    10000000000...
    pow(10, 610)
    10000000000...
    pow(10, 810)
    10000000000...

    Fit a straight line assuming :math:`\\textnormal{arg}^{1}` relationship
    (coefficients are returned):

    >>> p.do_fit()  # doctest: +SKIP
    array([1.16777420e-05,  1.00135803e-08])

    Predict run time for :math:`10^{5000}`:

    >>> p.predict_time(5000)  # doctest: +SKIP
    6.174564361572262e-05

    Plot the prediction:

    >>> p.plot(xlabeltext='Power of 10')  # doctest: +SKIP

    .. image:: /_static/timer_prediction_pow10.png
        :width: 450px
        :alt: Example plot from `astropy.utils.timer.RunTimePredictor`

    When the changing argument is not the last, e.g.,
    :math:`x^{2}`, something like this might work:

    >>> p = RunTimePredictor(lambda x: pow(x, 2))
    >>> p.time_func([2, 3, 5])
    >>> sorted(p.results.items())
    [(2, 4), (3, 9), (5, 25)]

    """
    def __init__(self, func, *args, **kwargs):
        self._funcname = func.__name__
        self._pfunc = partial(func, *args, **kwargs)
        self._cache_good = OrderedDict()
        self._cache_bad = []
        self._cache_est = OrderedDict()
        self._cache_out = OrderedDict()
        self._fit_func = None
        self._power = None

    @property
    def results(self):
        """Function outputs from `time_func`.

        A dictionary mapping input arguments (fixed arguments
        are not included) to their respective output values.

        """
        return self._cache_out

    @timefunc(num_tries=1, verbose=False)
    def _timed_pfunc(self, arg):
        """Run partial func once for single arg and time it."""
        return self._pfunc(arg)

    def _cache_time(self, arg):
        """Cache timing results without repetition."""
        if arg not in self._cache_good and arg not in self._cache_bad:
            try:
                result = self._timed_pfunc(arg)
            except Exception as e:
                warnings.warn(str(e), AstropyUserWarning)
                self._cache_bad.append(arg)
            else:
                self._cache_good[arg] = result[0]  # Run time
                self._cache_out[arg] = result[1]  # Function output

    def time_func(self, arglist):
        """Time the partial function for a list of single args
        and store run time in a cache. This forms a baseline for
        the prediction.

        This also stores function outputs in `results`.

        Parameters
        ----------
        arglist : list of numbers
            List of input arguments to time.

        """
        if not isinstance(arglist, Iterable):
            arglist = [arglist]

        # Preserve arglist order
        for arg in arglist:
            self._cache_time(arg)

    # FUTURE: Implement N^x * O(log(N)) fancy fitting.
    def do_fit(self, model=None, fitter=None, power=1, min_datapoints=3):
        """Fit a function to the lists of arguments and
        their respective run time in the cache.

        By default, this does a linear least-square fitting
        to a straight line on run time w.r.t. argument values
        raised to the given power, and returns the optimal
        intercept and slope.

        Parameters
        ----------
        model : `astropy.modeling.Model`
            Model for the expected trend of run time (Y-axis)
            w.r.t. :math:`\\textnormal{arg}^{\\textnormal{power}}` (X-axis).
            If `None`, will use `~astropy.modeling.polynomial.Polynomial1D`
            with ``degree=1``.

        fitter : `astropy.modeling.fitting.Fitter`
            Fitter for the given model to extract optimal coefficient values.
            If `None`, will use `~astropy.modeling.fitting.LinearLSQFitter`.

        power : int, optional
            Power of values to fit.

        min_datapoints : int, optional
            Minimum number of data points required for fitting.
            They can be built up with `time_func`.

        Returns
        -------
        a : array_like
            Fitted `~astropy.modeling.FittableModel` parameters.

        Raises
        ------
        AssertionError
            Insufficient data points for fitting.

        ModelsError
            Invalid model or fitter.

        """
        # Reset related attributes
        self._power = power
        self._cache_est = OrderedDict()

        x_arr = np.array(list(six.iterkeys(self._cache_good)))
        assert x_arr.size >= min_datapoints, \
            'Requires {0} points but has {1}'.format(min_datapoints,
                                                     x_arr.size)

        if model is None:
            model = modeling.models.Polynomial1D(1)
        elif not isinstance(model, modeling.core.Model):
            raise modeling.fitting.ModelsError(
                '{0} is not a model.'.format(model))

        if fitter is None:
            fitter = modeling.fitting.LinearLSQFitter()
        elif not isinstance(fitter, modeling.fitting.Fitter):
            raise modeling.fitting.ModelsError(
                '{0} is not a fitter.'.format(fitter))

        self._fit_func = fitter(
            model, x_arr**power, list(six.itervalues(self._cache_good)))

        return self._fit_func.parameters

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

        Raises
        ------
        AssertionError
            No fitted data for prediction.

        """
        if arg in self._cache_est:
            t_est = self._cache_est[arg]
        else:
            assert self._fit_func is not None, 'No fitted data for prediction'
            t_est = self._fit_func(arg**self._power)
            self._cache_est[arg] = t_est
        return t_est

    def plot(self, xscale='linear', yscale='linear', xlabeltext='args',
             save_as=''):  # pragma: no cover
        """Plot prediction.

        .. note:: Uses `matplotlib <http://matplotlib.sourceforge.net/>`_.

        Parameters
        ----------
        xscale, yscale : {'linear', 'log', 'symlog'}
            Scaling for `matplotlib.axes.Axes`.

        xlabeltext : str, optional
            Text for X-label.

        save_as : str, optional
            Save plot as given filename.

        Raises
        ------
        AssertionError
            Insufficient data for plotting.

        """
        import matplotlib.pyplot as plt

        # Actual data
        x_arr = sorted(self._cache_good)
        y_arr = np.array([self._cache_good[x] for x in x_arr])

        assert len(x_arr) > 1, 'Insufficient data for plotting'

        # Auto-ranging
        qmean = y_arr.mean() * u.second
        for cur_u in (u.minute, u.second, u.millisecond, u.microsecond,
                      u.nanosecond):
            val = qmean.to(cur_u).value
            if 1000 > val >= 1:
                break
        y_arr = (y_arr * u.second).to(cur_u).value

        fig, ax = plt.subplots()
        ax.plot(x_arr, y_arr, 'kx-', label='Actual')

        # Fitted data
        if self._fit_func is not None:
            x_est = list(six.iterkeys(self._cache_est))
            y_est = (np.array(list(six.itervalues(self._cache_est))) *
                     u.second).to(cur_u).value
            ax.scatter(x_est, y_est, marker='o', c='r', label='Predicted')

            x_fit = np.array(sorted(x_arr + x_est))
            y_fit = (self._fit_func(x_fit**self._power) *
                     u.second).to(cur_u).value
            ax.plot(x_fit, y_fit, 'b--', label='Fit')

        ax.set_xscale(xscale)
        ax.set_yscale(yscale)

        ax.set_xlabel(xlabeltext)
        ax.set_ylabel('Run time ({})'.format(cur_u.to_string()))
        ax.set_title(self._funcname)
        ax.legend(loc='best', numpoints=1)

        plt.draw()

        if save_as:
            plt.savefig(save_as)
