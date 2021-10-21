import abc
import numpy as np
from astropy.timeseries import TimeSeries, BinnedTimeSeries

__all__ = ['BasePeriodogram']


class BasePeriodogram:

    @abc.abstractmethod
    def __init__(self, t, y, dy=None):
        pass

    @classmethod
    def from_timeseries(cls, timeseries, signal_column_name=None, uncertainty=None, **kwargs):
        """
        Initialize a periodogram from a time series object.

        If a binned time series is passed, the time at the center of the bins is
        used. Also note that this method automatically gets rid of NaN/undefined
        values when initializing the periodogram.

        Parameters
        ----------
        signal_column_name : str
            The name of the column containing the signal values to use.
        uncertainty : str or float or `~astropy.units.Quantity`, optional
            The name of the column containing the errors on the signal, or the
            value to use for the error, if a scalar.
        **kwargs
            Additional keyword arguments are passed to the initializer for this
            periodogram class.
        """

        if signal_column_name is None:
            raise ValueError('signal_column_name should be set to a valid column name')

        y = timeseries[signal_column_name]
        keep = ~np.isnan(y)

        if isinstance(uncertainty, str):
            dy = timeseries[uncertainty]
            keep &= ~np.isnan(dy)
            dy = dy[keep]
        else:
            dy = uncertainty

        if isinstance(timeseries, TimeSeries):
            time = timeseries.time
        elif isinstance(timeseries, BinnedTimeSeries):
            time = timeseries.time_bin_center
        else:
            raise TypeError('Input time series should be an instance of '
                            'TimeSeries or BinnedTimeSeries')

        return cls(time[keep], y[keep], dy=dy, **kwargs)
