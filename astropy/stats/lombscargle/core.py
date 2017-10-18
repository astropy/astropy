"""Main Lomb-Scargle Implementation"""
from __future__ import print_function, division

import numpy as np

from .implementations import lombscargle, available_methods
from .implementations.mle import periodic_fit
from ... import units
from ...utils.compat.numpy import broadcast_arrays
from ...extern.six.moves import map


def has_units(obj):
    return hasattr(obj, 'unit')


def get_unit(obj):
    return getattr(obj, 'unit', 1)


def strip_units(*arrs):
    strip = lambda a: None if a is None else np.asarray(a)
    if len(arrs) == 1:
        return strip(arrs[0])
    else:
        return map(strip, arrs)


class LombScargle(object):
    """Compute the Lomb-Scargle Periodogram

    This implementations here are based on code presented in [1]_ and [2]_;
    if you use this functionality in an academic application, citation of
    those works would be appreciated.

    Parameters
    ----------
    t : array_like or Quantity
        sequence of observation times
    y : array_like or Quantity
        sequence of observations associated with times t
    dy : float, array_like or Quantity (optional)
        error or sequence of observational errors associated with times t
    fit_mean : bool (optional, default=True)
        if True, include a constant offset as part of the model at each
        frequency. This can lead to more accurate results, especially in the
        case of incomplete phase coverage.
    center_data : bool (optional, default=True)
        if True, pre-center the data by subtracting the weighted mean
        of the input data. This is especially important if fit_mean = False
    nterms : int (optional, default=1)
        number of terms to use in the Fourier fit

    Examples
    --------
    Generate noisy periodic data:

    >>> rand = np.random.RandomState(42)
    >>> t = 100 * rand.rand(100)
    >>> y = np.sin(2 * np.pi * t) + rand.randn(100)

    Compute the Lomb-Scargle periodogram on an automatically-determined
    frequency grid & find the frequency of max power:

    >>> frequency, power = LombScargle(t, y).autopower()
    >>> frequency[np.argmax(power)]  # doctest: +FLOAT_CMP
    1.0016662310392956

    Compute the Lomb-Scargle periodogram at a user-specified frequency grid:

    >>> freq = np.arange(0.8, 1.3, 0.1)
    >>> LombScargle(t, y).power(freq)  # doctest: +FLOAT_CMP
    array([0.0204304 , 0.01393845, 0.35552682, 0.01358029, 0.03083737])

    If the inputs are astropy Quantities with units, the units will be
    validated and the outputs will also be Quantities with appropriate units:

    >>> from astropy import units as u
    >>> t = t * u.s
    >>> y = y * u.mag
    >>> frequency, power = LombScargle(t, y).autopower()
    >>> frequency.unit
    Unit("1 / s")
    >>> power.unit
    Unit(dimensionless)

    Note here that the Lomb-Scargle power is always a unitless quantity,
    because it is related to the :math:`\\chi^2` of the best-fit periodic
    model at each frequency.

    References
    ----------
    .. [1] Vanderplas, J., Connolly, A. Ivezic, Z. & Gray, A. *Introduction to
        astroML: Machine learning for astrophysics*. Proceedings of the
        Conference on Intelligent Data Understanding (2012)
    .. [2] VanderPlas, J. & Ivezic, Z. *Periodograms for Multiband Astronomical
        Time Series*. ApJ 812.1:18 (2015)
    """
    available_methods = available_methods()

    def __init__(self, t, y, dy=None,
                 fit_mean=True, center_data=True, nterms=1):
        self.t, self.y, self.dy = self._validate_inputs(t, y, dy)
        self.fit_mean = fit_mean
        self.center_data = center_data
        self.nterms = nterms

    def _validate_inputs(self, t, y, dy):
        # Validate shapes of inputs
        if dy is None:
            t, y = broadcast_arrays(t, y, subok=True)
        else:
            t, y, dy = broadcast_arrays(t, y, dy, subok=True)
        if t.ndim != 1:
            raise ValueError("Inputs (t, y, dy) must be 1-dimensional")

        # validate units of inputs if any is a Quantity
        if any(has_units(arr) for arr in (t, y, dy)):
            t, y = map(units.Quantity, (t, y))
            if dy is not None:
                dy = units.Quantity(dy)
                try:
                    dy = units.Quantity(dy, unit=y.unit)
                except units.UnitConversionError:
                    raise ValueError("Units of dy not equivalent "
                                     "to units of y")
        return t, y, dy

    def _validate_frequency(self, frequency):
        frequency = np.asanyarray(frequency)

        if has_units(self.t):
            frequency = units.Quantity(frequency)
            try:
                frequency = units.Quantity(frequency, unit=1./self.t.unit)
            except units.UnitConversionError:
                raise ValueError("Units of frequency not equivalent to "
                                 "units of 1/t")
        else:
            if has_units(frequency):
                raise ValueError("frequency have units while 1/t doesn't.")
        return frequency

    def _validate_t(self, t):
        t = np.asanyarray(t)

        if has_units(self.t):
            t = units.Quantity(t)
            try:
                t = units.Quantity(t, unit=self.t.unit)
            except units.UnitConversionError:
                raise ValueError("Units of t not equivalent to "
                                 "units of input self.t")
        return t

    def _power_unit(self, norm):
        if has_units(self.y):
            if self.dy is None and norm == 'psd':
                return self.y.unit ** 2
            else:
                return units.dimensionless_unscaled
        else:
            return 1

    def autofrequency(self, samples_per_peak=5, nyquist_factor=5,
                      minimum_frequency=None, maximum_frequency=None):
        """Determine a suitable frequency grid for data.

        Note that this assumes the peak width is driven by the observational
        baseline, which is generally a good assumption when the baseline is
        much larger than the oscillation period.
        If you are searching for periods longer than the baseline of your
        observations, this may not perform well.

        Even with a large baseline, be aware that the maximum frequency
        returned is based on the concept of "average Nyquist frequency", which
        may not be useful for irregularly-sampled data. The maximum frequency
        can be adjusted via the nyquist_factor argument, or through the
        maximum_frequency argument.

        Parameters
        ----------
        samples_per_peak : float (optional, default=5)
            The approximate number of desired samples across the typical peak
        nyquist_factor : float (optional, default=5)
            The multiple of the average nyquist frequency used to choose the
            maximum frequency if maximum_frequency is not provided.
        minimum_frequency : float (optional)
            If specified, then use this minimum frequency rather than one
            chosen based on the size of the baseline.
        maximum_frequency : float (optional)
            If specified, then use this maximum frequency rather than one
            chosen based on the average nyquist frequency.

        Returns
        -------
        frequency : ndarray or Quantity
            The heuristically-determined optimal frequency bin
        """
        baseline = self.t.max() - self.t.min()
        n_samples = self.t.size

        df = 1.0 / baseline / samples_per_peak

        if minimum_frequency is None:
            minimum_frequency = 0.5 * df

        if maximum_frequency is None:
            avg_nyquist = 0.5 * n_samples / baseline
            maximum_frequency = nyquist_factor * avg_nyquist

        Nf = 1 + int(np.round((maximum_frequency - minimum_frequency) / df))

        return minimum_frequency + df * np.arange(Nf)

    def autopower(self, method='auto', method_kwds=None,
                  normalization='standard', samples_per_peak=5,
                  nyquist_factor=5, minimum_frequency=None,
                  maximum_frequency=None):
        """Compute Lomb-Scargle power at automatically-determined frequencies

        Parameters
        ----------
        method : string (optional)
            specify the lomb scargle implementation to use. Options are:

            - 'auto': choose the best method based on the input
            - 'fast': use the O[N log N] fast method. Note that this requires
              evenly-spaced frequencies: by default this will be checked unless
              ``assume_regular_frequency`` is set to True.
            - 'slow': use the O[N^2] pure-python implementation
            - 'cython': use the O[N^2] cython implementation. This is slightly
              faster than method='slow', but much more memory efficient.
            - 'chi2': use the O[N^2] chi2/linear-fitting implementation
            - 'fastchi2': use the O[N log N] chi2 implementation. Note that this
              requires evenly-spaced frequencies: by default this will be checked
              unless ``assume_regular_frequency`` is set to True.
            - 'scipy': use ``scipy.signal.lombscargle``, which is an O[N^2]
              implementation written in C. Note that this does not support
              heteroskedastic errors.

        method_kwds : dict (optional)
            additional keywords to pass to the lomb-scargle method
        normalization : string (optional, default='standard')
            Normalization to use for the periodogram.
            Options are 'standard', 'model', or 'psd'.
        samples_per_peak : float (optional, default=5)
            The approximate number of desired samples across the typical peak
        nyquist_factor : float (optional, default=5)
            The multiple of the average nyquist frequency used to choose the
            maximum frequency if maximum_frequency is not provided.
        minimum_frequency : float (optional)
            If specified, then use this minimum frequency rather than one
            chosen based on the size of the baseline.
        maximum_frequency : float (optional)
            If specified, then use this maximum frequency rather than one
            chosen based on the average nyquist frequency.

        Returns
        -------
        frequency, power : ndarrays
            The frequency and Lomb-Scargle power
        """
        frequency = self.autofrequency(samples_per_peak=samples_per_peak,
                                       nyquist_factor=nyquist_factor,
                                       minimum_frequency=minimum_frequency,
                                       maximum_frequency=maximum_frequency)
        power = self.power(frequency,
                           normalization=normalization,
                           method=method, method_kwds=method_kwds,
                           assume_regular_frequency=True)
        return frequency, power

    def power(self, frequency, normalization='standard', method='auto',
              assume_regular_frequency=False, method_kwds=None):
        """Compute the Lomb-Scargle power at the given frequencies

        Parameters
        ----------
        frequency : array_like or Quantity
            frequencies (not angular frequencies) at which to evaluate the
            periodogram. Note that in order to use method='fast', frequencies
            must be regularly-spaced.
        method : string (optional)
            specify the lomb scargle implementation to use. Options are:

            - 'auto': choose the best method based on the input
            - 'fast': use the O[N log N] fast method. Note that this requires
              evenly-spaced frequencies: by default this will be checked unless
              ``assume_regular_frequency`` is set to True.
            - 'slow': use the O[N^2] pure-python implementation
            - 'cython': use the O[N^2] cython implementation. This is slightly
              faster than method='slow', but much more memory efficient.
            - 'chi2': use the O[N^2] chi2/linear-fitting implementation
            - 'fastchi2': use the O[N log N] chi2 implementation. Note that this
              requires evenly-spaced frequencies: by default this will be checked
              unless ``assume_regular_frequency`` is set to True.
            - 'scipy': use ``scipy.signal.lombscargle``, which is an O[N^2]
              implementation written in C. Note that this does not support
              heteroskedastic errors.

        assume_regular_frequency : bool (optional)
            if True, assume that the input frequency is of the form
            freq = f0 + df * np.arange(N). Only referenced if method is 'auto'
            or 'fast'.
        normalization : string (optional, default='standard')
            Normalization to use for the periodogram.
            Options are 'standard', 'model', 'log', or 'psd'.
        fit_mean : bool (optional, default=True)
            if True, include a constant offset as part of the model at each
            frequency. This can lead to more accurate results, especially in
            the case of incomplete phase coverage.
        center_data : bool (optional, default=True)
            if True, pre-center the data by subtracting the weighted mean of
            the input data. This is especially important if fit_mean = False
        method_kwds : dict (optional)
            additional keywords to pass to the lomb-scargle method

        Returns
        -------
        power : ndarray
            The Lomb-Scargle power at the specified frequency
        """
        frequency = self._validate_frequency(frequency)
        power = lombscargle(*strip_units(self.t, self.y, self.dy),
                            frequency=strip_units(frequency),
                            center_data=self.center_data,
                            fit_mean=self.fit_mean,
                            nterms=self.nterms,
                            normalization=normalization,
                            method=method, method_kwds=method_kwds,
                            assume_regular_frequency=assume_regular_frequency)
        return power * self._power_unit(normalization)

    def model(self, t, frequency):
        """Compute the Lomb-Scargle model at the given frequency

        Parameters
        ----------
        t : array_like or Quantity, length n_samples
            times at which to compute the model
        frequency : float
            the frequency for the model

        Returns
        -------
        y : np.ndarray, length n_samples
            The model fit corresponding to the input times
        """
        frequency = self._validate_frequency(frequency)
        t = self._validate_t(t)
        y_fit = periodic_fit(*strip_units(self.t, self.y, self.dy),
                             frequency=strip_units(frequency),
                             t_fit=strip_units(t),
                             center_data=self.center_data,
                             fit_mean=self.fit_mean,
                             nterms=self.nterms)
        return y_fit * get_unit(self.y)
