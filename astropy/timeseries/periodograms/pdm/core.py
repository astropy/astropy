"""
Main Phase Dispersion Minimization (PDM) implementation
"""

from collections.abc import Sequence

import numpy as np

from astropy import units as u
from astropy.timeseries.periodograms.base import BasePeriodogram


class PhaseDispersionMinimization(BasePeriodogram):
    """
    Compute the Phase Dispersion Minimization (PDM) periodogram

    As described by Stellingwerf (1978) in [1]_. Mathematically,
    this is a least-squares fitting technique, but rather than a fit to a given curve (such as Fourier),
    the fit is relative to the mean curve as defined by the means of each bin. Simultaneously, we obtain
    the best least-squares light curve and the best period.

    Parameters
    ----------
    t : array-like
        The time values of the observations.
    y : array-like
        The signal values of the observations.
    dy : array-like, optional
        The uncertainties on the signal values.

    References
    ----------
    .. [1] Stellingwerf, R. F., “Period determination using phase dispersion minimization.”,
    <i>The Astrophysical Journal</i>, vol. 224, IOP, pp. 953–960, 1978. doi:10.1086/156444.

    """

    def __init__(
        self,
        t: Sequence[u.Quantity] | np.ndarray,
        y: Sequence[u.Quantity] | np.ndarray,
        dy=None,
    ):
        self.t = t
        self.y = y
        self.dy = dy

        self._validate_inputs()

        self.global_variance = np.var(y, ddof=1)  # denominator of theta

    def _validate_inputs(self):
        """
        Validate the inputs to the PDM periodogram.

        Raises
        ------
        ValueError
            If the input arrays are not of the same length or if they are empty.
        """
        if len(self.t) != len(self.y):
            raise ValueError("Time and signal arrays must be of the same length.")
        if len(self.t) == 0:
            raise ValueError("Input arrays must not be empty.")
        if len(self.t) < 2:
            raise ValueError(
                "At least two data points are required to compute the periodogram."
            )

    def run(self, periods, bins: int = 10) -> np.ndarray:
        """
        Calculate theta for each period
        """
        if len(periods) == 0:
            raise ValueError("Periods array must not be empty.")

        theta = np.zeros(len(periods))

        for i, period in enumerate(periods):
            # phase folding!
            phases = self.t % period / period

            bin_variances, bin_observations = self._bin_variances(
                phases, self.y, bins=bins
            )

            valid = (
                bin_observations > 1
            )  # valid bins are those with more than 1 observation

            num = np.sum(bin_variances[valid] * (bin_observations[valid] - 1))

            den = np.sum(bin_observations[valid]) - bins

            s2 = np.divide(num, den)  # numerator of theta

            theta[i] = s2 / self.global_variance

        return theta

    @staticmethod
    def _bin_variances(phases, y, bins: int = 10) -> tuple[np.ndarray, np.ndarray]:

        variances = np.zeros(bins)
        observations = np.zeros(bins)
        for i in range(bins):
            # bitwise operation to put as True the values that are in the bin
            bin_mask = (phases >= i / bins) & (phases < (i + 1) / bins)

            bin_y = y[bin_mask]
            observations[i] = len(bin_y)

            if len(bin_y) > 1:
                variances[i] = np.var(bin_y, ddof=1)
            else:
                variances[i] = (
                    0.0  # not necessary but if found a bin with 1 or 0 just put 0.0
                )

        return variances, observations
