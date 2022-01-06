# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings

import numpy as np
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.utils.exceptions import AstropyUserWarning

from astropy.timeseries.sampled import TimeSeries
from astropy.timeseries.binned import BinnedTimeSeries

__all__ = ['aggregate_downsample']


def reduceat(array, indices, function):
    """
    Manual reduceat functionality for cases where Numpy functions don't have a reduceat.
    It will check if the input function has a reduceat and call that if it does.
    """
    if len(indices) == 0:
        return np.array([])
    elif hasattr(function, 'reduceat'):
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
                         time_bin_end=None, n_bins=None, aggregate_func=None):
    """
    Downsample a time series by binning values into bins with a fixed size or
    custom sizes, using a single function to combine the values in the bin.

    Parameters
    ----------
    time_series : :class:`~astropy.timeseries.TimeSeries`
        The time series to downsample.
    time_bin_size : `~astropy.units.Quantity` or `~astropy.time.TimeDelta` ['time'], optional
        The time interval for the binned time series - this is either a scalar
        value (in which case all time bins will be assumed to have the same
        duration) or as an array of values (in which case each time bin can
        have a different duration). If this argument is provided,
        ``time_bin_end`` should not be provided.
    time_bin_start : `~astropy.time.Time` or iterable, optional
        The start time for the binned time series - this can be either given
        directly as a `~astropy.time.Time` array or as any iterable that
        initializes the `~astropy.time.Time` class. This can also be a scalar
        value if ``time_bin_size`` or ``time_bin_end`` is provided.
        Defaults to the first time in the sampled time series.
    time_bin_end : `~astropy.time.Time` or iterable, optional
        The times of the end of each bin - this can be either given directly as
        a `~astropy.time.Time` array or as any iterable that initializes the
        `~astropy.time.Time` class. This can only be given if ``time_bin_start``
        is provided or its default is used. If ``time_bin_end`` is scalar and
        ``time_bin_start`` is an array, time bins are assumed to be contiguous;
        the end of each bin is the start of the next one, and ``time_bin_end`` gives
        the end time for the last bin.  If ``time_bin_end`` is an array and
        ``time_bin_start`` is scalar, bins will be contiguous. If both ``time_bin_end``
        and ``time_bin_start`` are arrays, bins do not need to be contiguous.
        If this argument is provided, ``time_bin_size`` should not be provided.
    n_bins : int, optional
        The number of bins to use. Defaults to the number needed to fit all
        the original points. If both ``time_bin_start`` and ``time_bin_size``
        are provided and are scalar values, this determines the total bins
        within that interval. If ``time_bin_start`` is an iterable, this
        parameter will be ignored.
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

    if time_bin_size is not None and not isinstance(time_bin_size, (u.Quantity, TimeDelta)):
        raise TypeError("'time_bin_size' should be a Quantity or a TimeDelta")

    if time_bin_start is not None and not isinstance(time_bin_start, (Time, TimeDelta)):
        time_bin_start = Time(time_bin_start)

    if time_bin_end is not None and not isinstance(time_bin_end, (Time, TimeDelta)):
        time_bin_end = Time(time_bin_end)

    # Use the table sorted by time
    ts_sorted = time_series.iloc[:]

    # If start time is not provided, it is assumed to be the start of the timeseries
    if time_bin_start is None:
        time_bin_start = ts_sorted.time[0]

    # Total duration of the timeseries is needed for determining either
    # `time_bin_size` or `nbins` in the case of scalar `time_bin_start`
    if time_bin_start.isscalar:
        time_duration = (ts_sorted.time[-1] - time_bin_start).sec

    if time_bin_size is None and time_bin_end is None:
        if time_bin_start.isscalar:
            if n_bins is None:
                raise TypeError("With single 'time_bin_start' either 'n_bins', "
                                "'time_bin_size' or time_bin_end' must be provided")
            else:
                # `nbins` defaults to the number needed to fit all points
                time_bin_size = time_duration / n_bins * u.s
        else:
            time_bin_end = np.maximum(ts_sorted.time[-1], time_bin_start[-1])

    if time_bin_start.isscalar:
        if time_bin_size is not None:
            if time_bin_size.isscalar:
                # Determine the number of bins
                if n_bins is None:
                    bin_size_sec = time_bin_size.to_value(u.s)
                    n_bins = int(np.ceil(time_duration/bin_size_sec))
        elif time_bin_end is not None:
            if not time_bin_end.isscalar:
                # Convert start time to an array and populate using `time_bin_end`
                scalar_start_time = time_bin_start
                time_bin_start = time_bin_end.replicate(copy=True)
                time_bin_start[0] = scalar_start_time
                time_bin_start[1:] = time_bin_end[:-1]

    # Check for overlapping bins, and warn if they are present
    if time_bin_end is not None:
        if (not time_bin_end.isscalar and not time_bin_start.isscalar and
                np.any(time_bin_start[1:] < time_bin_end[:-1])):
            warnings.warn("Overlapping bins should be avoided since they "
                          "can lead to double-counting of data during binning.",
                          AstropyUserWarning)

    binned = BinnedTimeSeries(time_bin_size=time_bin_size,
                              time_bin_start=time_bin_start,
                              time_bin_end=time_bin_end,
                              n_bins=n_bins)

    if aggregate_func is None:
        aggregate_func = np.nanmean

    # Start and end times of the binned timeseries
    bin_start = binned.time_bin_start
    bin_end = binned.time_bin_end

    # Set `n_bins` to match the length of `time_bin_start` if
    # `n_bins` is unspecified or if `time_bin_start` is an iterable
    if n_bins is None or not time_bin_start.isscalar:
        n_bins = len(bin_start)

    # Find the subset of the table that is inside the union of all bins
    keep = ((ts_sorted.time >= bin_start[0]) & (ts_sorted.time <= bin_end[-1]))

    # Find out indices to be removed because of uncontiguous bins
    for ind in range(n_bins-1):
        delete_indices = np.where(np.logical_and(ts_sorted.time > bin_end[ind],
                                                 ts_sorted.time < bin_start[ind+1]))
        keep[delete_indices] = False

    subset = ts_sorted[keep]

    # Figure out which bin each row falls in by sorting with respect
    # to the bin end times
    indices = np.searchsorted(bin_end, ts_sorted.time[keep])

    # For time == bin_start[i+1] == bin_end[i], let bin_start takes precedence
    if len(indices) and np.all(bin_start[1:] >= bin_end[:-1]):
        indices_start = np.searchsorted(subset.time, bin_start[bin_start <= ts_sorted.time[-1]])
        indices[indices_start] = np.arange(len(indices_start))

    # Determine rows where values are defined
    if len(indices):
        groups = np.hstack([0, np.nonzero(np.diff(indices))[0] + 1])
    else:
        groups = np.array([])

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
