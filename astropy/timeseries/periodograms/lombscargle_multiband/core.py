"""Main Lomb-Scargle Multiband Implementation"""

import numpy as np

from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.timeseries.binned import BinnedTimeSeries
from astropy.timeseries.periodograms.lombscargle import LombScargle
from astropy.timeseries.sampled import TimeSeries

from .implementations import available_methods, lombscargle_multiband
from .implementations.mle import construct_regularization, design_matrix, periodic_fit

__all__ = ["LombScargleMultiband"]


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


class LombScargleMultiband(LombScargle):
    """Compute the Lomb-Scargle Periodogram.

    This implementation is based on code presented in [1]_ and [2]_;
    if you use this functionality in an academic application, citation of
    those works would be appreciated.

    Parameters
    ----------
    t : array-like or `~astropy.units.Quantity` ['time']
        sequence of observation times
    y : array-like or `~astropy.units.Quantity`
        sequence of observations associated with times t
    bands : array-like
        sequence of passband labels associated with times t, each unique label
        defines a single band of data.
    dy : float, array-like, or `~astropy.units.Quantity`, optional
        error or sequence of observational errors associated with times t
    normalization : {'standard', 'model', 'log', 'psd'}, optional
        Normalization to use for the periodogram.
    nterms_base : int, optional
        number of frequency terms to use for the base model common to all bands.
        In the case of the fast algorithm, this parameter is passed along to
        the single band LombScargle method as the ``nterms`` parameter.
    nterms_band : int, optional
        number of frequency terms to use for the residuals between the base
        model and each individual band
    reg_base : float or None (default = None)
        amount of regularization to use on the base model parameters
    reg_band : float or None (default = 1E-6)
        amount of regularization to use on the band model parameters
    regularize_by_trace : bool (default = True)
        if True, then regularization is expressed in units of the trace of
        the normal matrix
    center_data : bool, optional
        if True, pre-center the data by subtracting the weighted mean
        of the input data.
    fit_mean : bool, optional
        if True, include a constant offset as part of the model at each
        frequency. This can lead to more accurate results, especially in the
        case of incomplete phase coverage. Only applicable to the "fast" method

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
        bands,
        dy=None,
        normalization="standard",
        nterms_base=1,
        nterms_band=1,
        reg_base=None,
        reg_band=1e-6,
        regularize_by_trace=True,
        center_data=True,
        fit_mean=True,
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

        self._trel, self.y, self.bands, self.dy = self._validate_inputs(
            trel, y, bands, dy
        )

        self.normalization = normalization
        self.nterms_base = nterms_base
        self.nterms_band = nterms_band
        self.reg_base = reg_base
        self.reg_band = reg_band
        self.regularize_by_trace = regularize_by_trace
        self.center_data = center_data
        self.fit_mean = fit_mean
        self.nterms = self.nterms_base  # determined by the base model params

    @classmethod
    def from_timeseries(
        cls,
        timeseries,
        signal_column=None,
        uncertainty_column=None,
        band_labels=None,
        **kwargs,
    ):
        """
        Initialize a multiband periodogram from a time series object.

        If a binned time series is passed, the time at the center of the bins is
        used. Also note that this method automatically gets rid of NaN/undefined
        values when initializing the periodogram.

        Parameters
        ----------
        signal_column : list
            The names of columns containing the signal values to use.
        uncertainty_column : list, optional
            The names of columns containing the errors on the signal.
        band_labels : list, optional
            The labels for each band, matched by index. If none, uses the
            labels of ``signal_column`` as band names.
        **kwargs
            Additional keyword arguments are passed to the initializer for this
            periodogram class.
        """
        if signal_column is None:
            raise ValueError(
                "signal_column_name should be set to a list of valid column names"
            )

        if band_labels is not None:
            if len(band_labels) != len(signal_column):
                raise ValueError(
                    "band_labels have an equal number of elements to signal_column"
                )
        else:
            band_labels = signal_column  # use the flux labels as band labels

        if isinstance(timeseries, TimeSeries):
            time = timeseries.time
        elif isinstance(timeseries, BinnedTimeSeries):
            time = timeseries.time_bin_center
        else:
            raise TypeError(
                "Input time series should be an instance of "
                "TimeSeries or BinnedTimeSeries"
            )

        # Build lombscargle_multiband inputs
        t = []
        y = []
        dy = []
        band = []
        for i, col in enumerate(signal_column):
            if np.ma.is_masked(timeseries[col]):
                signal_mask = ~timeseries[col].mask
            else:
                signal_mask = ~np.isnan(timeseries[col])

            if uncertainty_column is not None:
                dy_col = timeseries[uncertainty_column[i]]
                if np.ma.is_masked(dy_col):
                    signal_mask &= ~dy_col.mask
                else:
                    signal_mask &= ~np.isnan(dy_col)

            t.append(time[signal_mask].mjd * u.day)
            y.append(timeseries[col][signal_mask])
            band.append([band_labels[i]] * sum(signal_mask))
            dy.append(timeseries[uncertainty_column[i]][signal_mask])

        t = np.hstack(t)
        y = np.hstack(y)
        band = np.hstack(band)
        if uncertainty_column is not None:
            dy = np.hstack(dy)

        if len(dy) == 0:
            dy = None
        return cls(t, y, band, dy=dy, **kwargs)

    def _validate_inputs(self, t, y, bands, dy):
        # Validate shapes of inputs
        if dy is None:
            t, y, bands = np.broadcast_arrays(t, y, bands, subok=True)
        else:
            t, y, bands, dy = np.broadcast_arrays(t, y, bands, dy, subok=True)
        if t.ndim != 1:
            raise ValueError("Inputs (t, y, bands, dy) must be 1-dimensional")

        # validate units of inputs if any is a Quantity
        if any(has_units(arr) for arr in (t, y, bands, dy)):
            t, y = map(u.Quantity, (t, y))
            if dy is not None:
                dy = u.Quantity(dy)
                try:
                    dy = u.Quantity(dy, unit=y.unit)
                except u.UnitConversionError:
                    raise ValueError("Units of dy not equivalent " "to units of y")
        return t, y, bands, dy

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
        if hasattr(self._trel, "unit"):
            unit = self._trel.unit
            trel = self._trel.to(u.day)  # do frequency calculations in days
        else:
            trel = self._trel
        baseline = trel.max() - trel.min()
        n_samples = trel.size

        df = 1.0 / baseline / samples_per_peak

        if minimum_frequency is None:
            minimum_frequency = 0.5 * df

        if maximum_frequency is None:
            avg_nyquist = 0.5 * n_samples / baseline
            maximum_frequency = nyquist_factor * avg_nyquist

        # Convert back to the input units
        if hasattr(self._trel, "unit"):
            df = df.to(1 / unit)
            minimum_frequency = minimum_frequency.to(1 / unit)
            maximum_frequency = maximum_frequency.to(1 / unit)

        n_freq = 1 + int(np.round((maximum_frequency - minimum_frequency) / df))

        if return_freq_limits:
            return minimum_frequency, minimum_frequency + df * (n_freq - 1)
        else:
            return minimum_frequency + df * np.arange(n_freq)

    def autopower(
        self,
        method="flexible",
        sb_method="auto",
        normalization="standard",
        samples_per_peak=5,
        nyquist_factor=5,
        minimum_frequency=None,
        maximum_frequency=None,
    ):
        """Compute Lomb-Scargle power at automatically-determined frequencies.

        Parameters
        ----------
        method : str, optional
            specify the multi-band lomb scargle implementation to use. Options are:

            - 'flexible': Constructs a common model, and an offset model per individual
                band. Applies regularization to the resulting model to constrain
                complexity.
            - 'fast': Passes each band individually through LombScargle (single-band),
                combines periodograms at the end by weight. Speed depends on single-band
                method chosen in 'sb_method'.

        sb_method : str, optional
            specify the single-band lomb scargle implementation to use, only in
            the case of using the 'fast' multiband method. Options are:

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
            frequency, method=method, sb_method=sb_method, normalization=normalization
        )

        return frequency, power

    def power(
        self, frequency, method="flexible", sb_method="auto", normalization="standard"
    ):
        """Compute the Lomb-Scargle power at the given frequencies.

        Parameters
        ----------
        frequency : array-like or `~astropy.units.Quantity` ['frequency']
            frequencies (not angular frequencies) at which to evaluate the
            periodogram. Note that in order to use method='fast', frequencies
            must be regularly-spaced.
        method : str, optional
            specify the multi-band lomb scargle implementation to use. Options are:

            - 'flexible': Constructs a common model, and an offset model per individual
                band. Applies regularization to the resulting model to constrain
                complexity.
            - 'fast': Passes each band individually through LombScargle (single-band),
                combines periodograms at the end by weight. Speed depends on single-band
                method chosen in 'sb_method'.
        sb_method : str, optional
            specify the single-band lomb scargle implementation to use, only in
            the case of using the 'fast' multiband method. Options can be found
            in `~astropy.timeseries.LombScargle`.
        normalization : {'standard', 'model', 'log', 'psd'}, optional
            If specified, override the normalization specified at instantiation.

        Returns
        -------
        power : ndarray
            The Lomb-Scargle power at the specified frequency
        """
        if normalization is None:
            normalization = self.normalization
        frequency = self._validate_frequency(frequency)
        f_shape = np.shape(frequency)
        power = lombscargle_multiband(
            strip_units(self._trel),
            strip_units(self.y),
            strip_units(self.bands),
            dy=strip_units(self.dy),
            frequency=strip_units(np.ravel(frequency)),
            method=method,
            sb_method=sb_method,
            normalization=normalization,
            nterms_base=self.nterms_base,
            nterms_band=self.nterms_band,
            reg_base=self.reg_base,
            reg_band=self.reg_band,
            regularize_by_trace=self.regularize_by_trace,
            center_data=self.center_data,
            fit_mean=self.fit_mean,
        )
        return np.reshape(power * self._power_unit(normalization), f_shape)

    def design_matrix(self, frequency, t_fit=None, bands_fit=None):
        """Compute the design matrix for a given frequency

        Parameters
        ----------
        frequency : float
            the frequency for the model
        t_fit : array-like, `~astropy.units.Quantity`, or `~astropy.time.Time` (optional)
            Times (length ``n_samples``) at which to compute the model.
            If not specified, then the times and uncertainties of the input
            data are used.
        bands_fit : array-like, or str
            Bands to use in fitting, must be subset of bands in input data.

        Returns
        -------
        ndarray
            The design matrix for the model at the given frequency.
            This should have a shape of (``len(t)``, ``n_parameters``).

        See Also
        --------
        model
        model_parameters
        offset
        """
        if t_fit is None:
            t_fit, dy = strip_units(self._trel, self.dy)
        else:
            t_fit, dy = strip_units(
                self._validate_t(self._as_relative_time("t", t_fit)), None
            )

        if bands_fit is None:
            bands_fit = np.unique(self.bands)
        elif type(bands_fit) == str:
            bands_fit = [bands_fit]
        unique_bands = np.unique(bands_fit)
        unique_bands_fit = np.unique(bands_fit)

        bands_fit = bands_fit[:, np.newaxis]

        if not set(unique_bands_fit).issubset(set(unique_bands)):
            raise ValueError(
                "bands_fit does not match training data: "
                f"input: {set(unique_bands_fit)} output: {set(unique_bands)}"
            )

        t_fit, bands_fit = np.broadcast_arrays(t_fit, bands_fit)
        return design_matrix(
            t_fit.ravel(),
            bands_fit.ravel(),
            frequency,
            dy,
            nterms_base=self.nterms_base,
            nterms_band=self.nterms_band,
        )

    def model(self, t, frequency, bands_fit=None):
        """Compute the Lomb-Scargle model at the given frequency.

        The model at a particular frequency is a linear model:
        model = offset + dot(design_matrix, model_parameters)

        Parameters
        ----------
        t : array-like or `~astropy.units.Quantity` ['time']
            Times (length ``n_samples``) at which to compute the model.
        frequency : float
            the frequency for the model
        bands_fit : list or array-like
            the unique bands to fit in the model

        Returns
        -------
        y : np.ndarray
            The model fit corresponding to the input times.
            Will have shape (``n_bands``,``n_samples``).

        See Also
        --------
        design_matrix
        offset
        model_parameters
        """
        if bands_fit is None:
            bands_fit = np.unique(self.bands)

        frequency = self._validate_frequency(frequency)
        t = self._validate_t(self._as_relative_time("t", t))
        y_fit = periodic_fit(
            *strip_units(self._trel, self.y, self.dy),
            bands=self.bands,
            frequency=strip_units(frequency),
            t_fit=strip_units(t),
            bands_fit=bands_fit,
            center_data=self.center_data,
            nterms_base=self.nterms_base,
            nterms_band=self.nterms_band,
            reg_base=self.reg_base,
            reg_band=self.reg_band,
            regularize_by_trace=self.regularize_by_trace,
        )

        return y_fit * get_unit(self.y)

    def offset(self, t_fit=None, bands_fit=None):
        """Return the offset array of the model

        The offset array of the model is the (weighted) mean of the y values in each band.
        Note that if self.center_data is False, the offset is 0 by definition.

        Parameters
        ----------
        t_fit : array-like, `~astropy.units.Quantity`, or `~astropy.time.Time` (optional)
            Times (length ``n_samples``) at which to compute the model.
            If not specified, then the times and uncertainties of the input
            data are used.
        bands_fit : array-like, or str
            Bands to use in fitting, must be subset of bands in input data.

        Returns
        -------
        offset : array

        See Also
        --------
        design_matrix
        model
        model_parameters
        """
        if bands_fit is None:
            bands_fit = np.unique(self.bands)
        if t_fit is None:
            on_fit = False
            t_fit = self.t
        else:
            on_fit = True

        bands_fit = bands_fit[:, np.newaxis]
        unique_bands = np.unique(self.bands)
        unique_bands_fit = np.unique(bands_fit)

        if not set(unique_bands_fit).issubset(set(unique_bands)):
            raise ValueError(
                "bands_fit does not match training data: "
                f"input: {set(unique_bands_fit)} output: {set(unique_bands)}"
            )
        y, dy = strip_units(self.y, self.dy)

        if np.shape(t_fit) != np.shape(
            bands_fit
        ):  # No need to broadcast if bands map to timestamps
            t_fit, bands_fit = np.broadcast_arrays(t_fit, bands_fit)

        # need to make sure all unique filters are represented
        u, i = np.unique(
            np.concatenate([bands_fit.ravel(), unique_bands]), return_inverse=True
        )

        if not self.center_data:
            return 0

        if dy is None:
            dy = 1
        dy = np.broadcast_to(dy, y.shape)

        # Calculate ymeans -- per band
        ymeans = np.zeros(y.shape)  # filter specific means
        for band in unique_bands:
            mask = self.bands == band
            ymeans[mask] = np.average(y[mask], weights=1 / dy[mask] ** 2)
        ymeans_fit = ymeans[i[: -len(unique_bands)]]

        if on_fit:
            return ymeans_fit * get_unit(self.y)
        else:
            return ymeans * get_unit(self.y)

    def model_parameters(self, frequency, units=True):
        r"""Compute the best-fit model parameters at the given frequency.

        The model described by these parameters is:

        .. math::

            y(t; f, \vec{\theta}) = \theta_0 + \sum_{n=1}^{\tt nterms_base} [\theta_{2n-1}\sin(2\pi n f t) + \theta_{2n}\cos(2\pi n f t)]
            + \theta_0^{(k)} + \sum_{n=1}^{\tt nterms_band} [\theta_{2n-1}^{(k)}\sin(2\pi n f t) +

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
        X = design_matrix(
            t,
            self.bands,
            frequency,
            dy=dy,
            nterms_base=self.nterms_base,
            nterms_band=self.nterms_band,
        )

        regularization = construct_regularization(
            self.bands,
            nterms_base=self.nterms_base,
            nterms_band=self.nterms_band,
            reg_base=self.reg_base,
            reg_band=self.reg_band,
        )

        M = np.dot(X.T, X)

        if regularization is not None:
            # M is being affected by operations on diag
            diag = M.ravel(order="K")[:: M.shape[0] + 1]

            if self.regularize_by_trace:
                diag += diag.sum() * np.asarray(regularization)
            else:
                diag += np.asarray(regularization)

        try:
            parameters = np.linalg.solve(M, np.dot(X.T, y / dy))
        except np.linalg.LinAlgError:
            parameters = np.dot(M, np.linalg.lstsq(X.T, y / dy)[0])
        if units:
            parameters = get_unit(self.y) * parameters

        return parameters

    def false_alarm_probability(self):
        """Not Implemented"""
        raise NotImplementedError

    def false_alarm_level(self):
        """Not Implemented"""
        raise NotImplementedError

    def distribution(self):
        """Not Implemented"""
        raise NotImplementedError
