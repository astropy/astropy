# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings

from copy import deepcopy

import numpy as np

from ..table import groups, QTable, Table
from ..time import Time, TimeDelta
from .. import units as u
from ..units import Quantity

from .core import TimeSeries
from .binned import BinnedTimeSeries

__all__ = ['SampledTimeSeries']


def reduceat(array, indices, function):
    """
    Manual reduceat functionality for cases where Numpy functions don't have a reduceat
    """
    result = [function(array[indices[i]:indices[i+1]]) for i in range(len(indices) - 1)]
    result.append(function(array[indices[-1]:]))
    return np.array(result)


class SampledTimeSeries(TimeSeries):

    _require_time_column = False

    def __init__(self, data=None, time=None, time_delta=None, n_samples=None, **kwargs):
        """
        """

        super().__init__(data=data, **kwargs)

        # FIXME: this is because for some operations, an empty time series needs
        # to be created, then columns added one by one. We should check that
        # when columns are added manually, time is added first and is of the
        # right type.
        if data is None and time is None and time_delta is None:
            self._required_columns = ['time']
            return

        # First if time has been given in the table data, we should extract it
        # and treat it as if it had been passed as a keyword argument.

        if data is not None:
            # TODO: raise error if also passed explicily and inconsistent
            n_samples = len(self)

        if 'time' in self.colnames:
            if time is None:
                time = self.columns['time']
                self.remove_column('time')
            else:
                raise TypeError("'time' has been given both in the table and as a keyword argument")

        if time is None:
            raise TypeError("'time' has not been specified")

        if not isinstance(time, Time):
            time = Time(time)

        if time_delta is not None and not isinstance(time_delta, (Quantity, TimeDelta)):
            raise TypeError("'time_delta' should be a Quantity or a TimeDelta")

        if isinstance(time_delta, TimeDelta):
            time_delta = time_delta.sec * u.s

        if time.isscalar:

            # We interpret this as meaning that time is that of the first
            # sample and that the interval is given by time_delta.

            if time_delta is None:
                raise TypeError("'time' is scalar, so 'bin_size' is required")

            if time_delta.isscalar:
                time_delta = np.repeat(time_delta, n_samples)

            time_delta = np.cumsum(time_delta)
            time_delta = np.roll(time_delta, 1)
            time_delta[0] = 0. * u.s

            time = time + time_delta

        else:

            if len(self.colnames) > 0 and len(time) != len(self):
                raise ValueError("Length of 'time' ({0}) should match "
                                 "table length ({1})".format(len(time), n_samples))

            if time_delta is not None:
                raise TypeError("'time_delta' should not be specified since "
                                "'time' is an array")

        self.add_column(time, index=0, name='time')
        self.add_index('time')

    @property
    def time(self):
        return self['time']

    def downsample(self, bin_size, func=None, start_time=None, n_bins=None, ):
        """
        Downsample the time series by binning values into bins with a fixed
        size, and return a :class:`~astropy.timeseries.BinnedTimeSeries`

        Parameters
        ----------
        bin_size : `~astropy.units.Quantity`
            The time interval for the binned time series
        func : callable, optional
            The function to use for combining points in the same bin. Defaults
            to np.nanmean.
        start_time : `~astropy.time.Time`, optional
            The start time for the binned time series. Defaults to the first
            time in the sampled time series.
        n_bins : int, optional
            The number of bins to use. Defaults to the number needed to fit all
            the original points.
        """

        bin_size_sec = bin_size.to_value(u.s)

        # Use the table sorted by time
        sorted = self.iloc[:]

        # Determine start time if needed
        if start_time is None:
            start_time = sorted.time[0]

        # Find the relative time since the start time, in seconds
        relative_time_sec = (sorted.time - start_time).sec

        # Determine the number of bins if needed
        if n_bins is None:
            n_bins = int(np.ceil(relative_time_sec[-1] / bin_size_sec))

        if func is None:
            func = np.nanmedian

        # Determine the bins
        relative_bins_sec = np.cumsum(np.hstack([0, np.repeat(bin_size_sec, n_bins)]))
        bins = start_time + relative_bins_sec * u.s

        # Find the subset of the table that is inside the bins
        keep = ((relative_time_sec >= relative_bins_sec[0]) &
                (relative_time_sec < relative_bins_sec[-1]))
        subset = sorted[keep]

        # Figure out which bin each row falls in - the -1 is because items
        # falling in the first bins will have index 1 but we want that to be 0
        indices = np.searchsorted(relative_bins_sec, relative_time_sec[keep]) - 1

        # Create new binned time series
        binned = BinnedTimeSeries(start_time=bins[:-1], end_time=bins[-1])

        # Determine rows where values are defined
        groups = np.hstack([0, np.nonzero(np.diff(indices))[0] + 1])

        # Find unique indices to determine which rows in the final time series
        # will not be empty.
        unique_indices = np.unique(indices)

        # Add back columns

        for colname in subset.colnames:

            if colname == 'time':
                continue

            values = subset[colname]

            # FIXME: figure out how to avoid the following, if possible
            if not isinstance(values, (np.ndarray, u.Quantity)):
                warnings.warn("Skipping column {0} since it has a mix-in type")
                continue

            data = np.ma.zeros(n_bins, dtype=values.dtype)
            data.mask = 1

            if isinstance(values, u.Quantity):
                data[unique_indices] = u.Quantity(reduceat(values.value, groups, func),
                                                  values.unit, copy=False)
            else:
                data[unique_indices] = reduceat(values, groups, func)

            data.mask[unique_indices] = 0
            binned[colname] = data

        return binned

    def fold(self, period=None, midpoint_epoch=None):
        """
        Return a new SampledTimeSeries folded with a period and midpoint epoch.

        Parameters
        ----------
        period : `~astropy.units.Quantity`
            The period to use for folding
        midpoint_epoch : `~astropy.time.Time`
            The time to use as the midpoint epoch, at which the relative
            time offset will be 0. Defaults to the first time in the time
            series.
        """

        folded = self.copy()
        folded.remove_column('time')

        if midpoint_epoch is None:
            midpoint_epoch = self.time[0]
        else:
            midpoint_epoch = Time(midpoint_epoch)

        period_sec = period.to_value(u.s)
        relative_time_sec = ((self.time - midpoint_epoch).sec + period_sec / 2) % period_sec - period_sec / 2

        folded_time = TimeDelta(relative_time_sec * u.s)

        folded.add_column(folded_time, name='time', index=0)

        return folded

    def __getitem__(self, item):
        if self._is_list_or_tuple_of_str(item):
            if 'time' not in item:
                out = QTable([self[x] for x in item],
                             meta=deepcopy(self.meta),
                             copy_indices=self._copy_indices)
                out._groups = groups.TableGroups(out, indices=self.groups._indices,
                                                 keys=self.groups._keys)
                return out
        return super().__getitem__(item)

    def add_columns(self, *args, **kwargs):
        result = super().add_columns(*args, **kwargs)
        if len(self.indices) == 0 and 'time' in self.colnames:
            self.add_index('time')
        return result

    @classmethod
    def from_pandas(self, df):
        """
        Convert a :class:`~pandas.DataFrame` to a
        :class:`astropy.timeseries.SampledTimeSeries`.
        """

        from pandas import DataFrame, DatetimeIndex

        if not isinstance(df, DataFrame):
            raise TypeError("Input should be a pandas dataframe")

        if not isinstance(df.index, DatetimeIndex):
            raise TypeError("DataFrame does not have a DatetimeIndex")

        # TODO: determine how user can specify time scale
        time = Time(df.index)

        # Create table without the time column
        table = Table.from_pandas(df)

        return SampledTimeSeries(time=time, data=table)

    def to_pandas(self):
        """
        Convert this time series to a :class:`~pandas.DataFrame` with a
        :class:`~pandas.DatetimeIndex` index.
        """

        # Extract table without time column
        table = self[[x for x in self.colnames if x != 'time']]

        # First make a normal pandas dataframe
        df = table.to_pandas()

        # Set index
        df.set_index(self.time.datetime64, inplace=True)

        return df
