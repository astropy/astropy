# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings

import numpy as np
from astropy import units as u
from astropy.utils.exceptions import AstropyUserWarning

from astropy.timeseries.sampled import TimeSeries
from astropy.timeseries.binned import BinnedTimeSeries

__all__ = ['aggregate_downsample']


def reduceat(array, indices, function):
    """
    Manual reduceat functionality for cases where Numpy functions don't have a reduceat.
    It will check if the input function has a reduceat and call that if it does.
    """
    if hasattr(function, 'reduceat'):
        return np.array(function.reduceat(array, indices))
    else:
        result = []
        for i in range(len(indices) - 1):
            if indices[i+1] <= indices[i]+1:
                result.append(function(array[indices[i]]))
            else:
                result.append(function(array[indices[i]:indices[i+1]]))
        result.append(function(array[indices[-1]:]))
        return np.array(result)


def aggregate_downsample(time_series, *, time_bin_size=None, time_bin_start=None,
                      n_bins=None, aggregate_func=None):
    """
    Downsample a time series by binning values into bins with a fixed size,
    using a single function to combine the values in the bin.

    Parameters
    ----------
    time_series : :class:`~astropy.timeseries.TimeSeries`
        The time series to downsample.
    time_bin_size : `~astropy.units.Quantity`
        The time interval for the binned time series.
    time_bin_start : `~astropy.time.Time`, optional
        The start time for the binned time series. Defaults to the first
        time in the sampled time series.
    n_bins : int, optional
        The number of bins to use. Defaults to the number needed to fit all
        the original points.
    aggregate_func : callable, optional
        The function to use for combining points in the same bin. Defaults
        to np.nanmean.

    Returns
    -------
    binned_time_series : :class:`~astropy.timeseries.BinnedTimeSeries`
        The downsampled time series.
    """

    if not isinstance(time_series, TimeSeries):
        raise TypeError("time_series should be a TimeSeries")

    if not isinstance(time_bin_size, u.Quantity):
        raise TypeError("time_bin_size should be a astropy.unit quantity")

    bin_size_sec = time_bin_size.to_value(u.s)

    # Use the table sorted by time
    sorted = time_series.iloc[:]

    # Determine start time if needed
    if time_bin_start is None:
        time_bin_start = sorted.time[0]

    # Find the relative time since the start time, in seconds
    relative_time_sec = (sorted.time - time_bin_start).sec

    # Determine the number of bins if needed
    if n_bins is None:
        n_bins = int(np.ceil(relative_time_sec[-1] / bin_size_sec))

    if aggregate_func is None:
        aggregate_func = np.nanmean

    # Determine the bins
    relative_bins_sec = np.cumsum(np.hstack([0, np.repeat(bin_size_sec, n_bins)]))
    bins = time_bin_start + relative_bins_sec * u.s

    # Find the subset of the table that is inside the bins
    keep = ((relative_time_sec >= relative_bins_sec[0]) &
            (relative_time_sec < relative_bins_sec[-1]))
    subset = sorted[keep]

    # Figure out which bin each row falls in - the -1 is because items
    # falling in the first bins will have index 1 but we want that to be 0
    indices = np.searchsorted(relative_bins_sec, relative_time_sec[keep]) - 1
    # Add back the first time.
    indices[relative_time_sec[keep] == relative_bins_sec[0]] = 0

    # Create new binned time series
    binned = BinnedTimeSeries(time_bin_start=bins[:-1], time_bin_end=bins[-1])

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
            warnings.warn("Skipping column {0} since it has a mix-in type", AstropyUserWarning)
            continue

        if isinstance(values, u.Quantity):
            data = u.Quantity(np.repeat(np.nan,  n_bins), unit=values.unit)
            data[unique_indices] = u.Quantity(reduceat(values.value, groups, aggregate_func),
                                              values.unit, copy=False)
        else:
            data = np.ma.zeros(n_bins, dtype=values.dtype)
            data.mask = 1
            data[unique_indices] = reduceat(values, groups, aggregate_func)
            data.mask[unique_indices] = 0

        binned[colname] = data

    return binned
