# Licensed under a 3-clause BSD style license - see LICENSE.rst

from copy import deepcopy

import numpy as np

from ..table import groups, QTable
from ..time import Time, TimeDelta
from .. import units as u
from ..units import Quantity


__all__ = ['TimeSeries', 'SampledTimeSeries', 'BinnedTimeSeries']


def reduceat(array, indices, function):
    """
    Manual reduceat functionality for cases where Numpy functions don't have a reduceat
    """
    result = [function(array[indices[i]:indices[i+1]]) for i in range(len(indices) - 1)]
    result.append(function(array[indices[-1]:]))
    return np.array(result)


class TimeSeries(QTable):
    pass


class SampledTimeSeries(TimeSeries):

    def __init__(self, data=None, time=None, time_delta=None, n_samples=None, **kwargs):
        """
        """

        super().__init__(data=data, **kwargs)

        # FIXME: this is because for some operations, an empty time series needs
        # to be created, then columns added one by one. We should check that
        # when columns are added manually, time is added first and is of the
        # right type.
        if data is None and time is None and time_delta is None:
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
        size, and return a `BinnedTimeSeries`

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

        # Determine the bins
        relative_bins_sec = np.cumsum(np.hstack([0, np.repeat(bin_size_sec, n_bins)]))
        bins = start_time + relative_bins_sec * u.s

        # Find the subset of the table that is inside the bins
        keep = (relative_time_sec >= relative_bins_sec[0]) & (relative_time_sec < relative_bins_sec[-1])
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
                data[unique_indices] = u.Quantity(reduceat(values.value, groups, func), values.unit, copy=False)
            else:
                data[unique_indices] = reduceat(values, groups, func)

            data.mask[unique_indices] = 0
            binned[colname] = data

        return binned

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


class BinnedTimeSeries(TimeSeries):

    def __init__(self, data=None, start_time=None, end_time=None, bin_size=None, n_bins=None, **kwargs):

        super().__init__(data=data, **kwargs)

        # First if start_time and end_time have been given in the table data, we
        # should extract them and treat them as if they had been passed as
        # keyword arguments.

        if 'start_time' in self.colnames:
            if start_time is None:
                start_time = self.columns['start_time']
                self.remove_column('start_time')
            else:
                raise TypeError("'start_time' has been given both in the table and as a keyword argument")

        if 'end_time' in self.colnames:
            if end_time is None:
                end_time = self.columns['end_time']
                self.remove_column('end_time')
            else:
                raise TypeError("'end_time' has been given both in the table and as a keyword argument")

        if start_time is None:
            raise TypeError("'start_time' has not been specified")

        if end_time is None and bin_size is None:
            raise TypeError("Cannot specify both 'end_time' and 'bin_size'")

        if not isinstance(start_time, Time):
            start_time = Time(start_time)

        if end_time is not None and not isinstance(end_time, Time):
            end_time = Time(end_time)

        if bin_size is not None and not isinstance(bin_size, (Quantity, TimeDelta)):
            raise TypeError("'bin_size' should be a Quantity or a TimeDelta")

        if isinstance(bin_size, TimeDelta):
            bin_size = bin_size.sec * u.s

        if start_time.isscalar:

            # We interpret this as meaning that this is the start of the
            # first bin and that the bins are contiguous. In this case,
            # we require bin_size to be specified.

            if bin_size is None:
                raise TypeError("'start_time' is scalar, so 'bin_size' is required")

            if bin_size.isscalar:

                if data is not None:
                    # TODO: raise error if also passed explicily and inconsistent
                    n_bins = len(self)

                bin_size = np.repeat(bin_size, n_bins)

            time_delta = np.cumsum(bin_size)
            end_time = start_time + time_delta

            # Now shift the array so that the first entry is 0
            time_delta = np.roll(time_delta, 1)
            time_delta[0] = 0. * u.s

            # Make start_time into an array
            start_time = start_time + time_delta

        else:

            if len(self.colnames) > 0 and len(start_time) != len(self):
                raise ValueError("Length of 'start_time' ({0}) should match "
                                 "table length ({1})".format(len(start_time), len(self)))

            if bin_size is not None:
                end_time = start_time + bin_size
            if end_time is not None:
                if end_time.isscalar:
                    times = start_time.copy()
                    times[:-1] = times[1:]
                    times[-1] = end_time
                    end_time = times
            else:
                raise TypeError("Either 'bin_size' or 'end_time' should be specified")

        self.add_column(start_time, index=0, name='start_time')
        self.add_column(end_time, index=1, name='end_time')
        self.add_index('start_time')
        self.add_index('end_time')

    @property
    def start_time(self):
        return self['start_time']

    @property
    def end_time(self):
        return self['end_time']
