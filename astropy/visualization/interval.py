# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Classes that deal with computing intervals from arrays of values based on
various criteria.
"""

import abc

import numpy as np

from astropy.utils.masked import get_data_and_mask

from .transform import BaseTransform

__all__ = [
    "AsymmetricPercentileInterval",
    "BaseInterval",
    "ManualInterval",
    "MinMaxInterval",
    "PercentileInterval",
    "ZScaleInterval",
]


class BaseInterval(BaseTransform):
    """
    Base class for the interval classes, which, when called with an
    array of values, return an interval computed following different
    algorithms.
    """

    @abc.abstractmethod
    def get_limits(self, values):
        """
        Return the minimum and maximum value in the interval based on
        the values provided.

        Parameters
        ----------
        values : ndarray
            The image values.

        Returns
        -------
        vmin, vmax : float
            The mininium and maximum image value in the interval.
        """
        raise NotImplementedError("Needs to be implemented in a subclass.")

    @staticmethod
    def _process_values(values):
        """
        Process the input values.

        This function filters out masked and/or invalid values (inf,
        nan) and returns a flattened 1D array.

        Parameters
        ----------
        values : array-like
            The input values.

        Returns
        -------
        result : 1D ndarray
            The processed values.
        """
        data, mask = get_data_and_mask(np.asanyarray(values))
        ok = np.isfinite(data)
        if mask is not None:
            ok &= ~mask

        return data[ok]

    def __call__(self, values, clip=True, out=None):
        """
        Transform values using this interval.

        The ``vmin`` and ``vmax`` values are determined by the
        `get_limits` method.

        The following transformation is then applied to the values:

        .. math::

            {\\rm result} = \\frac{{\\rm values} - v_{\\rm min}}
                                   {v_{\\rm max} - v_{\\rm min}}

        If ``clip`` is `True` (default), the result is then clipped to
        the [0:1] range.

        Parameters
        ----------
        values : array-like
            The input values.
        clip : bool, optional
            If `True` (default), values outside the [0:1] range are
            clipped to the [0:1] range.
        out : ndarray, optional
            If specified, the output values will be placed in this array
            (typically used for in-place calculations).

        Returns
        -------
        result : ndarray
            The transformed values.
        """
        vmin, vmax = self.get_limits(values)

        if out is None:
            values = np.subtract(values, float(vmin))
        else:
            if out.dtype.kind != "f":
                raise TypeError(
                    "Can only do in-place scaling for floating-point arrays"
                )
            values = np.subtract(values, float(vmin), out=out)

        if (vmax - vmin) != 0:
            np.true_divide(values, vmax - vmin, out=values)

        if clip:
            np.clip(values, 0.0, 1.0, out=values)

        return values


class ManualInterval(BaseInterval):
    """
    Interval based on user-specified values.

    Parameters
    ----------
    vmin : float, optional
        The minimum value in the scaling.  Defaults to the image
        minimum (ignoring NaNs)
    vmax : float, optional
        The maximum value in the scaling.  Defaults to the image
        maximum (ignoring NaNs)
    """

    def __init__(self, vmin=None, vmax=None):
        self.vmin = vmin
        self.vmax = vmax

    def get_limits(self, values):
        # Avoid overhead of preparing array if both limits have been specified
        # manually, for performance.
        if self.vmin is not None and self.vmax is not None:
            return self.vmin, self.vmax

        values = self._process_values(values)
        vmin = np.min(values) if self.vmin is None else self.vmin
        vmax = np.max(values) if self.vmax is None else self.vmax

        return vmin, vmax


class MinMaxInterval(BaseInterval):
    """
    Interval based on the minimum and maximum values in the data.
    """

    def get_limits(self, values):
        values = self._process_values(values)

        return np.min(values), np.max(values)


class AsymmetricPercentileInterval(BaseInterval):
    """
    Interval based on a keeping a specified fraction of pixels (can be
    asymmetric).

    Parameters
    ----------
    lower_percentile : float or None
        The lower percentile below which to ignore pixels. If None, then
        defaults to 0.
    upper_percentile : float or None
        The upper percentile above which to ignore pixels. If None, then
        defaults to 100.
    n_samples : int, optional
        Maximum number of values to use. If this is specified, and there
        are more values in the dataset as this, then values are randomly
        sampled from the array (with replacement).
    """

    def __init__(self, lower_percentile=None, upper_percentile=None, n_samples=None):
        self.lower_percentile = (
            lower_percentile if lower_percentile is not None else 0.0
        )
        self.upper_percentile = (
            upper_percentile if upper_percentile is not None else 100.0
        )
        self.n_samples = n_samples

    def get_limits(self, values):
        values = self._process_values(values)

        # If needed, limit the number of samples. We sample with replacement
        # since this is much faster.
        if self.n_samples is not None and values.size > self.n_samples:
            values = np.random.choice(values, self.n_samples)

        # Determine values at percentiles
        vmin, vmax = np.percentile(
            values, (self.lower_percentile, self.upper_percentile)
        )

        return vmin, vmax


class PercentileInterval(AsymmetricPercentileInterval):
    """
    Interval based on a keeping a specified fraction of pixels.

    Parameters
    ----------
    percentile : float
        The fraction of pixels to keep. The same fraction of pixels is
        eliminated from both ends.
    n_samples : int, optional
        Maximum number of values to use. If this is specified, and there
        are more values in the dataset as this, then values are randomly
        sampled from the array (with replacement).
    """

    def __init__(self, percentile, n_samples=None):
        lower_percentile = (100 - percentile) * 0.5
        upper_percentile = 100 - lower_percentile
        super().__init__(lower_percentile, upper_percentile, n_samples=n_samples)


class ZScaleInterval(BaseInterval):
    """
    Interval based on IRAF's zscale.

    Original implementation:
    https://github.com/spacetelescope/stsci.numdisplay/blob/master/lib/stsci/numdisplay/zscale.py

    Licensed under a 3-clause BSD style license (see AURA_LICENSE.rst).

    Parameters
    ----------
    n_samples : int, optional
        The number of points in the array to sample for determining
        scaling factors.  Defaults to 1000.

        .. versionchanged:: 7.0
            ``nsamples`` parameter is removed.

    contrast : float, optional
        The scaling factor (between 0 and 1) for determining the minimum
        and maximum value.  Larger values decrease the difference
        between the minimum and maximum values used for display.
        Defaults to 0.25.
    max_reject : float, optional
        If more than ``max_reject * npixels`` pixels are rejected, then
        the returned values are the minimum and maximum of the data.
        Defaults to 0.5.
    min_npixels : int, optional
        If there are less than ``min_npixels`` pixels remaining after
        the pixel rejection, then the returned values are the minimum
        and maximum of the data.  Defaults to 5.
    krej : float, optional
        The number of sigma used for the rejection. Defaults to 2.5.
    max_iterations : int, optional
        The maximum number of iterations for the rejection. Defaults to
        5.
    """

    def __init__(
        self,
        n_samples=1000,
        contrast=0.25,
        max_reject=0.5,
        min_npixels=5,
        krej=2.5,
        max_iterations=5,
    ):
        self.n_samples = n_samples
        self.contrast = contrast
        self.max_reject = max_reject
        self.min_npixels = min_npixels
        self.krej = krej
        self.max_iterations = max_iterations

    def get_limits(self, values):
        values = self._process_values(values)

        # Sample the image
        stride = int(max(1.0, values.size / self.n_samples))
        samples = values[::stride][: self.n_samples]
        samples.sort()

        npix = len(samples)
        vmin = samples[0]
        vmax = samples[-1]

        # Fit a line to the sorted array of samples
        minpix = max(self.min_npixels, int(npix * self.max_reject))
        x = np.arange(npix)
        ngoodpix = npix
        last_ngoodpix = npix + 1

        # Bad pixels mask used in k-sigma clipping
        badpix = np.zeros(npix, dtype=bool)

        # Kernel used to dilate the bad pixels mask
        ngrow = max(1, int(npix * 0.01))
        kernel = np.ones(ngrow, dtype=bool)

        for _ in range(self.max_iterations):
            if ngoodpix >= last_ngoodpix or ngoodpix < minpix:
                break

            fit = np.polyfit(x, samples, deg=1, w=(~badpix).astype(int))
            fitted = np.poly1d(fit)(x)

            # Subtract fitted line from the data array
            flat = samples - fitted

            # Compute the k-sigma rejection threshold
            threshold = self.krej * flat[~badpix].std()

            # Detect and reject pixels further than k*sigma from the
            # fitted line
            badpix[(flat < -threshold) | (flat > threshold)] = True

            # Convolve with a kernel of length ngrow
            badpix = np.convolve(badpix, kernel, mode="same")

            last_ngoodpix = ngoodpix
            ngoodpix = np.sum(~badpix)

        if ngoodpix >= minpix:
            slope, _ = fit

            if self.contrast > 0:
                slope = slope / self.contrast
            center_pixel = (npix - 1) // 2
            median = np.median(samples)
            vmin = max(vmin, median - (center_pixel - 1) * slope)
            vmax = min(vmax, median + (npix - center_pixel) * slope)

        return vmin, vmax
