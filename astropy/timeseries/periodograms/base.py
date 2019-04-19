import numpy as np
from astropy.timeseries import TimeSeries, BinnedTimeSeries

__all__ = ['BasePeriodogram']


class BasePeriodogram:

    @classmethod
    def from_timeseries(cls, timeseries, *, column=None, error=None, **kwargs):
        """
        Initialize a periodogram from a time series object.

        If a binned time series is passed, the time at the center of the bins is
        used. Also note that this method automatically gets rid of NaN/undefined
        values when initalizing the periodogram.

        Parameters
        ----------
        column : str, optional
            The name of the column containing the y values to use.
        error : str or float or `~astropy.units.Quantity`, optional
            The name of the column containing the y error values, or the value
            to use for the error, if a scalar.
        **kwargs
            Additional keyword arguments are passed to the initializer for this
            periodogram class.
        """

        if column is None:
            raise ValueError('column should be set to a valid column name')

        y = timeseries[column]
        keep = ~np.isnan(y)

        if isinstance(error, str):
            dy = timeseries[error]
            keep &= ~np.isnan(dy)
            dy = dy[keep]
        else:
            dy = error

        if isinstance(timeseries, TimeSeries):
            time = timeseries.time
        elif isinstance(timeseries, BinnedTimeSeries):
            time = timeseries.time_bin_center
        else:
            raise TypeError('Input time series should be an instance of '
                            'TimeSeries or BinnedTimeSeries')

        return cls(time[keep], y[keep], dy=dy, **kwargs)
