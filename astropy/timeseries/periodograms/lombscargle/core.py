"""Main Lomb-Scargle Implementation."""

import numpy as np

from astropy import units
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.timeseries.periodograms.base import BasePeriodogram

from . import _statistics
from .implementations import available_methods, lombscargle
from .implementations.mle import design_matrix, periodic_fit


def has_units(obj):
    return hasattr(obj, "unit")


def get_unit(obj):
    return getattr(obj, "unit", 1)


def strip_units(*arrs):
    strip = lambda a: None if a is None else np.asarray(a)
    if len(arrs) == 1:
        return strip(arrs[0])
    else:
        return map(strip, arrs)


class LombScargle(BasePeriodogram):
    """Compute the Lomb-Scargle Periodogram.

    This implementations here are based on code presented in [1]_ and [2]_;
    if you use this functionality in an academic application, citation of
    those works would be appreciated.

    Parameters
    ----------
    t : array-like or `~astropy.units.Quantity` ['time']
        sequence of observation times
    y : array-like or `~astropy.units.Quantity`
        sequence of observations associated with times t
    dy : float, array-like, or `~astropy.units.Quantity`, optional
        error or sequence of observational errors associated with times t
    fit_mean : bool, optional
        if True, include a constant offset as part of the model at each
        frequency. This can lead to more accurate results, especially in the
        case of incomplete phase coverage.
    center_data : bool, optional
        if True, pre-center the data by subtracting the weighted mean
        of the input data. This is especially important if fit_mean = False
    nterms : int, optional
        number of terms to use in the Fourier fit
    normalization : {'standard', 'model', 'log', 'psd'}, optional
        Normalization to use for the periodogram.

    Examples
    --------
    Generate noisy periodic data:

    >>> rand = np.random.default_rng(42)
    >>> t = 100 * rand.random(100)
    >>> y = np.sin(2 * np.pi * t) + rand.standard_normal(100)

    Compute the Lomb-Scargle periodogram on an automatically-determined
    frequency grid & find the frequency of max power:

    >>> frequency, power = LombScargle(t, y).autopower()
    >>> frequency[np.argmax(power)]  # doctest: +FLOAT_CMP
    1.0007641728995051

    Compute the Lomb-Scargle periodogram at a user-specified frequency grid:

    >>> freq = np.arange(0.8, 1.3, 0.1)
    >>> LombScargle(t, y).power(freq)  # doctest: +FLOAT_CMP
    array([0.0792948 , 0.01778874, 0.25328167, 0.01064157, 0.01471387])

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

    def __init__(
        self,
        t,
        y,
        dy=None,
        fit_mean=True,
        center_data=True,
        nterms=1,
        normalization="standard",
    ):
        # If t is a TimeDelta, convert it to a quantity. The units we convert
        # to don't really matter since the user gets a Quantity back at the end
        # so can convert to any units they like.
        if isinstance(t, TimeDelta):
            t = t.to("day")

        # We want to expose self.t as being the times the user passed in, but
        # if the times are absolute, we need to convert them to relative times
        # internally, so we use self._trel and self._tstart for this.

        self.t = t

        if isinstance(self.t, Time):
            self._tstart = self.t[0]
            trel = (self.t - self._tstart).to(u.day)
        else:
            self._tstart = None
            trel = self.t

        self._trel, self.y, self.dy = self._validate_inputs(trel, y, dy)

        self.fit_mean = fit_mean
        self.center_data = center_data
        self.nterms = nterms
        self.normalization = normalization

    def _validate_inputs(self, t, y, dy):
        # Validate shapes of inputs
        if dy is None:
            t, y = np.broadcast_arrays(t, y, subok=True)
        else:
            t, y, dy = np.broadcast_arrays(t, y, dy, subok=True)
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
                    raise ValueError("Units of dy not equivalent to units of y")
        return t, y, dy

    def _validate_frequency(self, frequency):
        frequency = np.asanyarray(frequency)

        if has_units(self._trel):
            frequency = units.Quantity(frequency)
            try:
                frequency = units.Quantity(frequency, unit=1.0 / self._trel.unit)
            except units.UnitConversionError:
                raise ValueError("Units of frequency not equivalent to units of 1/t")
        else:
            if has_units(frequency):
                raise ValueError("frequency have units while 1/t doesn't.")
        return frequency

    def _validate_t(self, t):
        t = np.asanyarray(t)

        if has_units(self._trel):
            t = units.Quantity(t)
            try:
                t = units.Quantity(t, unit=self._trel.unit)
            except units.UnitConversionError:
                raise ValueError("Units of t not equivalent to units of input self.t")
        return t

    def _power_unit(self, norm):
        if has_units(self.y):
            if self.dy is None and norm == "psd":
                return self.y.unit**2
            else:
                return units.dimensionless_unscaled
        else:
            return 1

    def autofrequency(
        self,
        samples_per_peak=5,
        nyquist_factor=5,
        minimum_frequency=None,
        maximum_frequency=None,
        return_freq_limits=False,
    ):
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
        samples_per_peak : float, optional
            The approximate number of desired samples across the typical peak
        nyquist_factor : float, optional
            The multiple of the average nyquist frequency used to choose the
            maximum frequency if maximum_frequency is not provided.
        minimum_frequency : float, optional
            If specified, then use this minimum frequency rather than one
            chosen based on the size of the baseline.
        maximum_frequency : float, optional
            If specified, then use this maximum frequency rather than one
            chosen based on the average nyquist frequency.
        return_freq_limits : bool, optional
            if True, return only the frequency limits rather than the full
            frequency grid.

        Returns
        -------
        frequency : ndarray or `~astropy.units.Quantity` ['frequency']
            The heuristically-determined optimal frequency bin
        """
        baseline = self._trel.max() - self._trel.min()
        n_samples = self._trel.size

        df = 1.0 / baseline / samples_per_peak

        if minimum_frequency is None:
            minimum_frequency = 0.5 * df

        if maximum_frequency is None:
            avg_nyquist = 0.5 * n_samples / baseline
            maximum_frequency = nyquist_factor * avg_nyquist

        Nf = 1 + int(np.round((maximum_frequency - minimum_frequency) / df))

        if return_freq_limits:
            return minimum_frequency, minimum_frequency + df * (Nf - 1)
        else:
            return minimum_frequency + df * np.arange(Nf)

    def autopower(
        self,
        method="auto",
        method_kwds=None,
        normalization=None,
        samples_per_peak=5,
        nyquist_factor=5,
        minimum_frequency=None,
        maximum_frequency=None,
    ):
        """Compute Lomb-Scargle power at automatically-determined frequencies.

        Parameters
        ----------
        method : str, optional
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

        method_kwds : dict, optional
            additional keywords to pass to the lomb-scargle method
        normalization : {'standard', 'model', 'log', 'psd'}, optional
            If specified, override the normalization specified at instantiation.
        samples_per_peak : float, optional
            The approximate number of desired samples across the typical peak
        nyquist_factor : float, optional
            The multiple of the average nyquist frequency used to choose the
            maximum frequency if maximum_frequency is not provided.
        minimum_frequency : float or `~astropy.units.Quantity` ['frequency'], optional
            If specified, then use this minimum frequency rather than one
            chosen based on the size of the baseline. Should be `~astropy.units.Quantity`
            if inputs to LombScargle are `~astropy.units.Quantity`.
        maximum_frequency : float or `~astropy.units.Quantity` ['frequency'], optional
            If specified, then use this maximum frequency rather than one
            chosen based on the average nyquist frequency. Should be `~astropy.units.Quantity`
            if inputs to LombScargle are `~astropy.units.Quantity`.

        Returns
        -------
        frequency, power : ndarray
            The frequency and Lomb-Scargle power
        """
        frequency = self.autofrequency(
            samples_per_peak=samples_per_peak,
            nyquist_factor=nyquist_factor,
            minimum_frequency=minimum_frequency,
            maximum_frequency=maximum_frequency,
        )
        power = self.power(
            frequency,
            normalization=normalization,
            method=method,
            method_kwds=method_kwds,
            assume_regular_frequency=True,
        )
        return frequency, power

    def power(
        self,
        frequency,
        normalization=None,
        method="auto",
        assume_regular_frequency=False,
        method_kwds=None,
    ):
        """Compute the Lomb-Scargle power at the given frequencies.

        Parameters
        ----------
        frequency : array-like or `~astropy.units.Quantity` ['frequency']
            frequencies (not angular frequencies) at which to evaluate the
            periodogram. Note that in order to use method='fast', frequencies
            must be regularly-spaced.
        method : str, optional
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

        assume_regular_frequency : bool, optional
            if True, assume that the input frequency is of the form
            freq = f0 + df * np.arange(N). Only referenced if method is 'auto'
            or 'fast'.
        normalization : {'standard', 'model', 'log', 'psd'}, optional
            If specified, override the normalization specified at instantiation.

        method_kwds : dict, optional
            additional keywords to pass to the lomb-scargle method

        Returns
        -------
        power : ndarray
            The Lomb-Scargle power at the specified frequency
        """
        if normalization is None:
            normalization = self.normalization
        frequency = self._validate_frequency(frequency)
        power = lombscargle(
            *strip_units(self._trel, self.y, self.dy),
            frequency=strip_units(frequency),
            center_data=self.center_data,
            fit_mean=self.fit_mean,
            nterms=self.nterms,
            normalization=normalization,
            method=method,
            method_kwds=method_kwds,
            assume_regular_frequency=assume_regular_frequency,
        )
        return power * self._power_unit(normalization)

    def _as_relative_time(self, name, times):
        """
        Convert the provided times (if absolute) to relative times using the
        current _tstart value. If the times provided are relative, they are
        returned without conversion (though we still do some checks).
        """
        if isinstance(times, TimeDelta):
            times = times.to("day")

        if self._tstart is None:
            if isinstance(times, Time):
                raise TypeError(
                    f"{name} was provided as an absolute time but "
                    "the LombScargle class was initialized "
                    "with relative times."
                )
        else:
            if isinstance(times, Time):
                times = (times - self._tstart).to(u.day)
            else:
                raise TypeError(
                    f"{name} was provided as a relative time but "
                    "the LombScargle class was initialized "
                    "with absolute times."
                )

        return times

    def model(self, t, frequency):
        """Compute the Lomb-Scargle model at the given frequency.

        The model at a particular frequency is a linear model:
        model = offset + dot(design_matrix, model_parameters)

        Parameters
        ----------
        t : array-like or `~astropy.units.Quantity` ['time']
            Times (length ``n_samples``) at which to compute the model.
        frequency : float
            the frequency for the model

        Returns
        -------
        y : np.ndarray
            The model fit corresponding to the input times
            (will have length ``n_samples``).

        See Also
        --------
        design_matrix
        offset
        model_parameters
        """
        frequency = self._validate_frequency(frequency)
        t = self._validate_t(self._as_relative_time("t", t))
        y_fit = periodic_fit(
            *strip_units(self._trel, self.y, self.dy),
            frequency=strip_units(frequency),
            t_fit=strip_units(t),
            center_data=self.center_data,
            fit_mean=self.fit_mean,
            nterms=self.nterms,
        )
        return y_fit * get_unit(self.y)

    def offset(self):
        """Return the offset of the model.

        The offset of the model is the (weighted) mean of the y values.
        Note that if self.center_data is False, the offset is 0 by definition.

        Returns
        -------
        offset : scalar

        See Also
        --------
        design_matrix
        model
        model_parameters
        """
        y, dy = strip_units(self.y, self.dy)
        if dy is None:
            dy = 1
        dy = np.broadcast_to(dy, y.shape)
        if self.center_data:
            w = dy**-2.0
            y_mean = np.dot(y, w) / w.sum()
        else:
            y_mean = 0
        return y_mean * get_unit(self.y)

    def model_parameters(self, frequency, units=True):
        r"""Compute the best-fit model parameters at the given frequency.

        The model described by these parameters is:

        .. math::

            y(t; f, \vec{\theta}) = \theta_0 + \sum_{n=1}^{\tt nterms} [\theta_{2n-1}\sin(2\pi n f t) + \theta_{2n}\cos(2\pi n f t)]

        where :math:`\vec{\theta}` is the array of parameters returned by this function.

        Parameters
        ----------
        frequency : float
            the frequency for the model
        units : bool
            If True (default), return design matrix with data units.

        Returns
        -------
        theta : np.ndarray (n_parameters,)
            The best-fit model parameters at the given frequency.

        See Also
        --------
        design_matrix
        model
        offset
        """
        frequency = self._validate_frequency(frequency)
        t, y, dy = strip_units(self._trel, self.y, self.dy)

        if self.center_data:
            y = y - strip_units(self.offset())

        dy = np.ones_like(y) if dy is None else np.asarray(dy)
        X = self.design_matrix(frequency)
        parameters = np.linalg.solve(np.dot(X.T, X), np.dot(X.T, y / dy))
        if units:
            parameters = get_unit(self.y) * parameters
        return parameters

    def design_matrix(self, frequency, t=None):
        """Compute the design matrix for a given frequency.

        Parameters
        ----------
        frequency : float
            the frequency for the model
        t : array-like, `~astropy.units.Quantity`, or `~astropy.time.Time` (optional)
            Times (length ``n_samples``) at which to compute the model.
            If not specified, then the times and uncertainties of the input
            data are used.

        Returns
        -------
        X : array
            The design matrix for the model at the given frequency.
            This should have a shape of (``len(t)``, ``n_parameters``).

        See Also
        --------
        model
        model_parameters
        offset
        """
        if t is None:
            t, dy = strip_units(self._trel, self.dy)
        else:
            t, dy = strip_units(self._validate_t(self._as_relative_time("t", t)), None)
        return design_matrix(t, frequency, dy, nterms=self.nterms, bias=self.fit_mean)

    def distribution(self, power, cumulative=False):
        """Expected periodogram distribution under the null hypothesis.

        This computes the expected probability distribution or cumulative
        probability distribution of periodogram power, under the null
        hypothesis of a non-varying signal with Gaussian noise. Note that
        this is not the same as the expected distribution of peak values;
        for that see the ``false_alarm_probability()`` method.

        Parameters
        ----------
        power : array-like
            The periodogram power at which to compute the distribution.
        cumulative : bool, optional
            If True, then return the cumulative distribution.

        See Also
        --------
        false_alarm_probability
        false_alarm_level

        Returns
        -------
        dist : np.ndarray
            The probability density or cumulative probability associated with
            the provided powers.
        """
        dH = 1 if self.fit_mean or self.center_data else 0
        dK = dH + 2 * self.nterms
        dist = _statistics.cdf_single if cumulative else _statistics.pdf_single
        return dist(power, len(self._trel), self.normalization, dH=dH, dK=dK)

    def false_alarm_probability(
        self,
        power,
        method="baluev",
        samples_per_peak=5,
        nyquist_factor=5,
        minimum_frequency=None,
        maximum_frequency=None,
        method_kwds=None,
    ):
        """False alarm probability of periodogram maxima under the null hypothesis.

        This gives an estimate of the false alarm probability given the height
        of the largest peak in the periodogram, based on the null hypothesis
        of non-varying data with Gaussian noise.

        Parameters
        ----------
        power : array-like
            The periodogram value.
        method : {'baluev', 'davies', 'naive', 'bootstrap'}, optional
            The approximation method to use.
        maximum_frequency : float
            The maximum frequency of the periodogram.
        method_kwds : dict, optional
            Additional method-specific keywords.

        Returns
        -------
        false_alarm_probability : np.ndarray
            The false alarm probability

        Notes
        -----
        The true probability distribution for the largest peak cannot be
        determined analytically, so each method here provides an approximation
        to the value. The available methods are:

        - "baluev" (default): the upper-limit to the alias-free probability,
          using the approach of Baluev (2008) [1]_.
        - "davies" : the Davies upper bound from Baluev (2008) [1]_.
        - "naive" : the approximate probability based on an estimated
          effective number of independent frequencies.
        - "bootstrap" : the approximate probability based on bootstrap
          resamplings of the input data.

        Note also that for normalization='psd', the distribution can only be
        computed for periodograms constructed with errors specified.

        See Also
        --------
        distribution
        false_alarm_level

        References
        ----------
        .. [1] Baluev, R.V. MNRAS 385, 1279 (2008)
        """
        if self.nterms != 1:
            raise NotImplementedError(
                "false alarm probability is not implemented for multiterm periodograms."
            )
        if not (self.fit_mean or self.center_data):
            raise NotImplementedError(
                "false alarm probability is implemented "
                "only for periodograms of centered data."
            )

        fmin, fmax = self.autofrequency(
            samples_per_peak=samples_per_peak,
            nyquist_factor=nyquist_factor,
            minimum_frequency=minimum_frequency,
            maximum_frequency=maximum_frequency,
            return_freq_limits=True,
        )
        return _statistics.false_alarm_probability(
            power,
            fmax=fmax,
            t=self._trel,
            y=self.y,
            dy=self.dy,
            normalization=self.normalization,
            method=method,
            method_kwds=method_kwds,
        )

    def false_alarm_level(
        self,
        false_alarm_probability,
        method="baluev",
        samples_per_peak=5,
        nyquist_factor=5,
        minimum_frequency=None,
        maximum_frequency=None,
        method_kwds=None,
    ):
        """Level of maximum at a given false alarm probability.

        This gives an estimate of the periodogram level corresponding to a
        specified false alarm probability for the largest peak, assuming a
        null hypothesis of non-varying data with Gaussian noise.

        Parameters
        ----------
        false_alarm_probability : array-like
            The false alarm probability (0 < fap < 1).
        maximum_frequency : float
            The maximum frequency of the periodogram.
        method : {'baluev', 'davies', 'naive', 'bootstrap'}, optional
            The approximation method to use; default='baluev'.
        method_kwds : dict, optional
            Additional method-specific keywords.

        Returns
        -------
        power : np.ndarray
            The periodogram peak height corresponding to the specified
            false alarm probability.

        Notes
        -----
        The true probability distribution for the largest peak cannot be
        determined analytically, so each method here provides an approximation
        to the value. The available methods are:

        - "baluev" (default): the upper-limit to the alias-free probability,
          using the approach of Baluev (2008) [1]_.
        - "davies" : the Davies upper bound from Baluev (2008) [1]_.
        - "naive" : the approximate probability based on an estimated
          effective number of independent frequencies.
        - "bootstrap" : the approximate probability based on bootstrap
          resamplings of the input data. The number of samples can
          be set with the method-specific keyword "n_bootstraps" (default=1000).

        Note also that for normalization='psd', the distribution can only be
        computed for periodograms constructed with errors specified.

        See Also
        --------
        distribution
        false_alarm_probability

        References
        ----------
        .. [1] Baluev, R.V. MNRAS 385, 1279 (2008)
        """
        if self.nterms != 1:
            raise NotImplementedError(
                "false alarm probability is not implemented for multiterm periodograms."
            )
        if not (self.fit_mean or self.center_data):
            raise NotImplementedError(
                "false alarm probability is implemented "
                "only for periodograms of centered data."
            )

        fmin, fmax = self.autofrequency(
            samples_per_peak=samples_per_peak,
            nyquist_factor=nyquist_factor,
            minimum_frequency=minimum_frequency,
            maximum_frequency=maximum_frequency,
            return_freq_limits=True,
        )
        return _statistics.false_alarm_level(
            false_alarm_probability,
            fmax=fmax,
            t=self._trel,
            y=self.y,
            dy=self.dy,
            normalization=self.normalization,
            method=method,
            method_kwds=method_kwds,
        )
