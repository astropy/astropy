# Licensed under a 3-clause BSD style license - see LICENSE.rst

from copy import deepcopy

import numpy as np

from astropy.table import groups, QTable, Table
from astropy.time import Time, TimeDelta
from astropy import units as u
from astropy.units import Quantity

from astropy.timeseries.core import BaseTimeSeries, autocheck_required_columns

__all__ = ['TimeSeries']


@autocheck_required_columns
class TimeSeries(BaseTimeSeries):
    """
    A class to represent time series data in tabular form.

    `~astropy.timeseries.TimeSeries` provides a class for representing time
    series as a collection of values of different quantities measured at specific
    points in time (for time series with finite time bins, see the
    `~astropy.timeseries.BinnedTimeSeries` class).
    `~astropy.timeseries.TimeSeries` is a sub-class of `~astropy.table.QTable`
    and thus provides all the standard table maniplation methods available to
    tables, but it also provides additional conveniences for dealing with time
    series, such as a flexible initializer for setting up the times, a method
    for folding time series, and a ``time`` attribute for easy access to the
    time values.

    See also: http://docs.astropy.org/en/stable/timeseries/

    Parameters
    ----------
    data : numpy ndarray, dict, list, `~astropy.table.Table`, or table-like object, optional
        Data to initialize time series. This does not need to contain the times,
        which can be provided separately, but if it does contain the times they
        should be in a column called ``'time'`` to be automatically recognized.
    time : `~astropy.time.Time` or iterable
        The times at which the values are sampled - this can be either given
        directly as a `~astropy.time.Time` array or as any iterable that
        initializes the `~astropy.time.Time` class. If this is given, then
        the remaining time-related arguments should not be used.
    time_start : `~astropy.time.Time` or str
        The time of the first sample in the time series. This is an alternative
        to providing ``time`` and requires that ``time_delta`` is also provided.
    time_delta : `~astropy.time.TimeDelta` or `~astropy.units.Quantity`
        The step size in time for the series. This can either be a scalar if
        the time series is evenly sampled, or an array of values if it is not.
    n_samples : int
        The number of time samples for the series. This is only used if both
        ``time_start`` and ``time_delta`` are provided and are scalar values.
    **kwargs : dict, optional
        Additional keyword arguments are passed to `~astropy.table.QTable`.
    """

    _required_columns = ['time']

    def __init__(self, data=None, *, time=None, time_start=None,
                 time_delta=None, n_samples=None, **kwargs):

        super().__init__(data=data, **kwargs)

        # For some operations, an empty time series needs to be created, then
        # columns added one by one. We should check that when columns are added
        # manually, time is added first and is of the right type.
        if data is None and time is None and time_start is None and time_delta is None:
            self._required_columns_relax = True
            return

        # First if time has been given in the table data, we should extract it
        # and treat it as if it had been passed as a keyword argument.

        if data is not None:
            if n_samples is not None:
                if n_samples != len(self):
                    raise TypeError("'n_samples' has been given both and it is not the "
                                    "same length as the input data.")
            else:
                n_samples = len(self)

        if 'time' in self.colnames:
            if time is None:
                time = self.columns['time']
            else:
                raise TypeError("'time' has been given both in the table and as a keyword argument")

        if time is None and time_start is None:
            raise TypeError("Either 'time' or 'time_start' should be specified")
        elif time is not None and time_start is not None:
            raise TypeError("Cannot specify both 'time' and 'time_start'")

        if time is not None and not isinstance(time, Time):
            time = Time(time)

        if time_start is not None and not isinstance(time_start, Time):
            time_start = Time(time_start)

        if time_delta is not None and not isinstance(time_delta, (Quantity, TimeDelta)):
            raise TypeError("'time_delta' should be a Quantity or a TimeDelta")

        if isinstance(time_delta, TimeDelta):
            time_delta = time_delta.sec * u.s

        if time_start is not None:

            # We interpret this as meaning that time is that of the first
            # sample and that the interval is given by time_delta.

            if time_delta is None:
                raise TypeError("'time' is scalar, so 'time_delta' is required")

            if time_delta.isscalar:
                time_delta = np.repeat(time_delta, n_samples)

            time_delta = np.cumsum(time_delta)
            time_delta = np.roll(time_delta, 1)
            time_delta[0] = 0. * u.s

            time = time_start + time_delta

        elif len(self.colnames) > 0 and len(time) != len(self):
            raise ValueError("Length of 'time' ({0}) should match "
                             "data length ({1})".format(len(time), n_samples))

        elif time_delta is not None:
            raise TypeError("'time_delta' should not be specified since "
                            "'time' is an array")

        with self._delay_required_column_checks():
            if 'time' in self.colnames:
                self.remove_column('time')
            self.add_column(time, index=0, name='time')

    @property
    def time(self):
        """
        The time values.
        """
        return self['time']

    def fold(self, period=None, midpoint_epoch=None):
        """
        Return a new `~astropy.timeseries.TimeSeries` folded with a period and
        midpoint epoch.

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

        if midpoint_epoch is None:
            midpoint_epoch = self.time[0]
        else:
            midpoint_epoch = Time(midpoint_epoch)

        period_sec = period.to_value(u.s)
        relative_time_sec = ((self.time - midpoint_epoch).sec + period_sec / 2) % period_sec - period_sec / 2

        folded_time = TimeDelta(relative_time_sec * u.s)

        with folded._delay_required_column_checks():
            folded.remove_column('time')
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
        """
        See :meth:`~astropy.table.Table.add_columns`.
        """
        # Note that the docstring is inherited from QTable
        result = super().add_columns(*args, **kwargs)
        if len(self.indices) == 0 and 'time' in self.colnames:
            self.add_index('time')
        return result

    @classmethod
    def from_pandas(self, df, time_scale='utc'):
        """
        Convert a :class:`~pandas.DataFrame` to a
        :class:`astropy.timeseries.TimeSeries`.

        Parameters
        ----------
        df : :class:`pandas.DataFrame`
            A pandas :class:`pandas.DataFrame` instance.
        time_scale : str
            The time scale to pass into `astropy.time.Time`.
            Defaults to ``UTC``.

        """
        from pandas import DataFrame, DatetimeIndex

        if not isinstance(df, DataFrame):
            raise TypeError("Input should be a pandas DataFrame")

        if not isinstance(df.index, DatetimeIndex):
            raise TypeError("DataFrame does not have a DatetimeIndex")

        time = Time(df.index, scale=time_scale)
        table = Table.from_pandas(df)

        return TimeSeries(time=time, data=table)

    def to_pandas(self):
        """
        Convert this :class:`~astropy.timeseries.TimeSeries` to a
        :class:`~pandas.DataFrame` with a :class:`~pandas.DatetimeIndex` index.

        Returns
        -------
        dataframe : :class:`pandas.DataFrame`
            A pandas :class:`pandas.DataFrame` instance
        """
        return Table(self).to_pandas(index='time')

    @classmethod
    def read(self, filename, time_column=None, time_format=None, time_scale=None, format=None, *args, **kwargs):
        """
        Read and parse a file and returns a `astropy.timeseries.TimeSeries`.

        This method uses the unified I/O infrastructure in Astropy which makes
        it easy to define readers/writers for various classes
        (http://docs.astropy.org/en/stable/io/unified.html). By default, this
        method will try and use readers defined specifically for the
        `astropy.timeseries.TimeSeries` class - however, it is also
        possible to use the ``format`` keyword to specify formats defined for
        the `astropy.table.Table` class - in this case, you will need to also
        provide the column names for column containing the start times for the
        bins, as well as other column names (see the Parameters section below
        for details)::

            >>> from astropy.timeseries import TimeSeries
            >>> ts = TimeSeries.read('sampled.dat', format='ascii.ecsv',
            ...                      time_column='date')  # doctest: +SKIP

        Parameters
        ----------
        filename : str
            File to parse.
        format : str
            File format specifier.
        time_column : str, optional
            The name of the time column.
        time_format : str, optional
            The time format for the time column.
        time_scale : str, optional
            The time scale for the time column.
        *args : tuple, optional
            Positional arguments passed through to the data reader.
        **kwargs : dict, optional
            Keyword arguments passed through to the data reader.

        Returns
        -------
        out : `astropy.timeseries.sampled.TimeSeries`
            TimeSeries corresponding to file contents.

        Notes
        -----
        """
        try:

            # First we try the readers defined for the BinnedTimeSeries class
            return super().read(filename, format=format, *args, **kwargs)

        except TypeError:

            # Otherwise we fall back to the default Table readers

            if time_column is None:
                raise ValueError("``time_column`` should be provided since the default Table readers are being used.")

            table = Table.read(filename, format=format, *args, **kwargs)

            if time_column in table.colnames:
                time = Time(table.columns[time_column], scale=time_scale, format=time_format)
                table.remove_column(time_column)
            else:
                raise ValueError("Time column '{}' not found in the input data.".format(time_column))

            return TimeSeries(time=time, data=table)
