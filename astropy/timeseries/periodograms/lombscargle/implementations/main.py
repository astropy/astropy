"""
Main Lomb-Scargle Implementation

The ``lombscargle`` function here is essentially a sophisticated switch
statement for the various implementations available in this submodule
"""

__all__ = ['lombscargle', 'available_methods']

import numpy as np

from .slow_impl import lombscargle_slow
from .fast_impl import lombscargle_fast
from .scipy_impl import lombscargle_scipy
from .chi2_impl import lombscargle_chi2
from .fastchi2_impl import lombscargle_fastchi2
from .cython_impl import lombscargle_cython


METHODS = {'slow': lombscargle_slow,
           'fast': lombscargle_fast,
           'chi2': lombscargle_chi2,
           'scipy': lombscargle_scipy,
           'fastchi2': lombscargle_fastchi2,
           'cython': lombscargle_cython}


def available_methods():
    methods = ['auto', 'slow', 'chi2', 'cython', 'fast', 'fastchi2']

    # Scipy required for scipy algorithm (obviously)
    try:
        import scipy
    except ImportError:
        pass
    else:
        methods.append('scipy')
    return methods


def _is_regular(frequency):
    frequency = np.asarray(frequency)

    if frequency.ndim != 1:
        return False
    elif len(frequency) == 1:
        return True
    else:
        diff = np.diff(frequency)
        return np.allclose(diff[0], diff)


def _get_frequency_grid(frequency, assume_regular_frequency=False):
    """Utility to get grid parameters from a frequency array

    Parameters
    ----------
    frequency : array_like or Quantity
        input frequency grid
    assume_regular_frequency : bool (default = False)
        if True, then do not check whether frequency is a regular grid

    Returns
    -------
    f0, df, N : scalars
        Parameters such that all(frequency == f0 + df * np.arange(N))
    """
    frequency = np.asarray(frequency)
    if frequency.ndim != 1:
        raise ValueError("frequency grid must be 1 dimensional")
    elif len(frequency) == 1:
        return frequency[0], frequency[0], 1
    elif not (assume_regular_frequency or _is_regular(frequency)):
        raise ValueError("frequency must be a regular grid")

    return frequency[0], frequency[1] - frequency[0], len(frequency)


def validate_method(method, dy, fit_mean, nterms,
                    frequency, assume_regular_frequency):
    """
    Validate the method argument, and if method='auto'
    choose the appropriate method
    """
    methods = available_methods()
    prefer_fast = (len(frequency) > 200
                   and (assume_regular_frequency or _is_regular(frequency)))
    prefer_scipy = 'scipy' in methods and dy is None and not fit_mean

    # automatically choose the appropriate method
    if method == 'auto':

        if nterms != 1:
            if prefer_fast:
                method = 'fastchi2'
            else:
                method = 'chi2'
        elif prefer_fast:
            method = 'fast'
        elif prefer_scipy:
            method = 'scipy'
        else:
            method = 'cython'

    if method not in METHODS:
        raise ValueError(f"invalid method: {method}")

    return method


def lombscargle(t, y, dy=None,
                frequency=None,
                method='auto',
                assume_regular_frequency=False,
                normalization='standard',
                fit_mean=True, center_data=True,
                method_kwds=None, nterms=1):
    """
    Compute the Lomb-scargle Periodogram with a given method.

    Parameters
    ----------
    t : array_like
        sequence of observation times
    y : array_like
        sequence of observations associated with times t
    dy : float or array_like, optional
        error or sequence of observational errors associated with times t
    frequency : array_like
        frequencies (not angular frequencies) at which to evaluate the
        periodogram. If not specified, optimal frequencies will be chosen using
        a heuristic which will attempt to provide sufficient frequency range
        and sampling so that peaks will not be missed. Note that in order to
        use method='fast', frequencies must be regularly spaced.
    method : str, optional
        specify the lomb scargle implementation to use. Options are:

        - 'auto': choose the best method based on the input
        - 'fast': use the O[N log N] fast method. Note that this requires
          evenly-spaced frequencies: by default this will be checked unless
          ``assume_regular_frequency`` is set to True.
        - `slow`: use the O[N^2] pure-python implementation
        - `chi2`: use the O[N^2] chi2/linear-fitting implementation
        - `fastchi2`: use the O[N log N] chi2 implementation. Note that this
          requires evenly-spaced frequencies: by default this will be checked
          unless `assume_regular_frequency` is set to True.
        - `scipy`: use ``scipy.signal.lombscargle``, which is an O[N^2]
          implementation written in C. Note that this does not support
          heteroskedastic errors.

    assume_regular_frequency : bool, optional
        if True, assume that the input frequency is of the form
        freq = f0 + df * np.arange(N). Only referenced if method is 'auto'
        or 'fast'.
    normalization : str, optional
        Normalization to use for the periodogram.
        Options are 'standard' or 'psd'.
    fit_mean : bool, optional
        if True, include a constant offset as part of the model at each
        frequency. This can lead to more accurate results, especially in the
        case of incomplete phase coverage.
    center_data : bool, optional
        if True, pre-center the data by subtracting the weighted mean
        of the input data. This is especially important if `fit_mean = False`
    method_kwds : dict, optional
        additional keywords to pass to the lomb-scargle method
    nterms : int, optional
        number of Fourier terms to use in the periodogram.
        Not supported with every method.

    Returns
    -------
    PLS : array_like
        Lomb-Scargle power associated with each frequency omega
    """
    # frequencies should be one-dimensional arrays
    output_shape = frequency.shape
    frequency = frequency.ravel()

    # we'll need to adjust args and kwds for each method
    args = (t, y, dy)
    kwds = dict(frequency=frequency,
                center_data=center_data,
                fit_mean=fit_mean,
                normalization=normalization,
                nterms=nterms,
                **(method_kwds or {}))

    method = validate_method(method, dy=dy, fit_mean=fit_mean, nterms=nterms,
                             frequency=frequency,
                             assume_regular_frequency=assume_regular_frequency)

    # scipy doesn't support dy or fit_mean=True
    if method == 'scipy':
        if kwds.pop('fit_mean'):
            raise ValueError("scipy method does not support fit_mean=True")
        if dy is not None:
            dy = np.ravel(np.asarray(dy))
            if not np.allclose(dy[0], dy):
                raise ValueError("scipy method only supports "
                                 "uniform uncertainties dy")
        args = (t, y)

    # fast methods require frequency expressed as a grid
    if method.startswith('fast'):
        f0, df, Nf = _get_frequency_grid(kwds.pop('frequency'),
                                         assume_regular_frequency)
        kwds.update(f0=f0, df=df, Nf=Nf)

    # only chi2 methods support nterms
    if not method.endswith('chi2'):
        if kwds.pop('nterms') != 1:
            raise ValueError("nterms != 1 only supported with 'chi2' "
                             "or 'fastchi2' methods")

    PLS = METHODS[method](*args, **kwds)
    return PLS.reshape(output_shape)
