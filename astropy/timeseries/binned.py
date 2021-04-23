# Licensed under a 3-clause BSD style license - see LICENSE.rst

from copy import deepcopy

import numpy as np

from astropy.table import groups, Table, QTable
from astropy.time import Time, TimeDelta
from astropy import units as u
from astropy.units import Quantity

from astropy.timeseries.core import BaseTimeSeries, autocheck_required_columns

__all__ = ['BinnedTimeSeries']


@autocheck_required_columns
class BinnedTimeSeries(BaseTimeSeries):
    """
    A class to represent binned time series data in tabular form.

    `~astropy.timeseries.BinnedTimeSeries` provides a class for
    representing time series as a collection of values of different
    quantities measured in time bins (for time series with values
    sampled at specific times, see the `~astropy.timeseries.TimeSeries`
    class). `~astropy.timeseries.BinnedTimeSeries` is a sub-class of
    `~astropy.table.QTable` and thus provides all the standard table
    maniplation methods available to tables, but it also provides
    additional conveniences for dealing with time series, such as a
    flexible initializer for setting up the times, and attributes to
    access the start/center/end time of bins.

    See also: https://docs.astropy.org/en/stable/timeseries/

    Parameters
    ----------
    data : numpy ndarray, dict, list, table-like object, optional
        Data to initialize time series. This does not need to contain the
        times, which can be provided separately, but if it does contain the
        times they should be in columns called ``'time_bin_start'`` and
        ``'time_bin_size'`` to be automatically recognized.
    time_bin_start : `~astropy.time.Time` or iterable
        The times of the start of each bin - this can be either given
        directly as a `~astropy.time.Time` array or as any iterable that
        initializes the `~astropy.time.Time` class. If this is given, then
        the remaining time-related arguments should not be used. This can also
        be a scalar value if ``time_bin_size`` is provided.
    time_bin_end : `~astropy.time.Time` or iterable
        The times of the end of each bin - this can be either given directly
        as a `~astropy.time.Time` array or as any value or iterable that
        initializes the `~astropy.time.Time` class. If this is given, then the
        remaining time-related arguments should not be used. This can only be
        given if ``time_bin_start`` is an array of values. If ``time_bin_end``
        is a scalar, time bins are assumed to be contiguous, such that the end
        of each bin is the start of the next one, and ``time_bin_end`` gives
        the end time for the last bin. If ``time_bin_end`` is an array, the
        time bins do not need to be contiguous. If this argument is provided,
        ``time_bin_size`` should not be provided.
    time_bin_size : `~astropy.time.TimeDelta` or `~astropy.units.Quantity`
        The size of the time bins, either as a scalar value (in which case all
        time bins will be assumed to have the same duration) or as an array of
        values (in which case each time bin can have a different duration).
        If this argument is provided, ``time_bin_end`` should not be provided.
    n_bins : int
        The number of time bins for the series. This is only used if both
        ``time_bin_start`` and ``time_bin_size`` are provided and are scalar
        values.
    **kwargs : dict, optional
        Additional keyword arguments are passed to `~astropy.table.QTable`.
    """

    _required_columns = ['time_bin_start', 'time_bin_size']

    def __init__(self, data=None, *, time_bin_start=None, time_bin_end=None,
                 time_bin_size=None, n_bins=None, **kwargs):

        super().__init__(data=data, **kwargs)

        # For some operations, an empty time series needs to be created, then
        # columns added one by one. We should check that when columns are added
        # manually, time is added first and is of the right type.
        if (data is None and time_bin_start is None and time_bin_end is None and
                time_bin_size is None and n_bins is None):
            self._required_columns_relax = True
            return

        # First if time_bin_start and time_bin_end have been given in the table data, we
        # should extract them and treat them as if they had been passed as
        # keyword arguments.

        if 'time_bin_start' in self.colnames:
            if time_bin_start is None:
                time_bin_start = self.columns['time_bin_start']
            else:
                raise TypeError("'time_bin_start' has been given both in the table "
                                "and as a keyword argument")

        if 'time_bin_size' in self.colnames:
            if time_bin_size is None:
                time_bin_size = self.columns['time_bin_size']
            else:
                raise TypeError("'time_bin_size' has been given both in the table "
                                "and as a keyword argument")

        if time_bin_start is None:
            raise TypeError("'time_bin_start' has not been specified")

        if time_bin_end is None and time_bin_size is None:
            raise TypeError("Either 'time_bin_size' or 'time_bin_end' should be specified")

        if not isinstance(time_bin_start, (Time, TimeDelta)):
            time_bin_start = Time(time_bin_start)

        if time_bin_end is not None and not isinstance(time_bin_end, (Time, TimeDelta)):
            time_bin_end = Time(time_bin_end)

        if time_bin_size is not None and not isinstance(time_bin_size, (Quantity, TimeDelta)):
            raise TypeError("'time_bin_size' should be a Quantity or a TimeDelta")

        if isinstance(time_bin_size, TimeDelta):
            time_bin_size = time_bin_size.sec * u.s

        if n_bins is not None and time_bin_size is not None:
            if not (time_bin_start.isscalar and time_bin_size.isscalar):
                raise TypeError("'n_bins' cannot be specified if 'time_bin_start' or "
                                "'time_bin_size' are not scalar'")

        if time_bin_start.isscalar:

            # We interpret this as meaning that this is the start of the
            # first bin and that the bins are contiguous. In this case,
            # we require time_bin_size to be specified.

            if time_bin_size is None:
                raise TypeError("'time_bin_start' is scalar, so 'time_bin_size' is required")

            if time_bin_size.isscalar:
                if data is not None:
                    if n_bins is not None:
                        if n_bins != len(self):
                            raise TypeError("'n_bins' has been given and it is not the "
                                            "same length as the input data.")
                    else:
                        n_bins = len(self)

                time_bin_size = np.repeat(time_bin_size, n_bins)

            time_delta = np.cumsum(time_bin_size)
            time_bin_end = time_bin_start + time_delta

            # Now shift the array so that the first entry is 0
            time_delta = np.roll(time_delta, 1)
            time_delta[0] = 0. * u.s

            # Make time_bin_start into an array
            time_bin_start = time_bin_start + time_delta

        else:

            if len(self.colnames) > 0 and len(time_bin_start) != len(self):
                raise ValueError("Length of 'time_bin_start' ({}) should match "
                                 "table length ({})".format(len(time_bin_start), len(self)))

            if time_bin_end is not None:
                if time_bin_end.isscalar:
                    times = time_bin_start.copy()
                    times[:-1] = times[1:]
                    times[-1] = time_bin_end
                    time_bin_end = times
                time_bin_size = (time_bin_end - time_bin_start).sec * u.s

        if time_bin_size.isscalar:
            time_bin_size = np.repeat(time_bin_size, len(self))

        with self._delay_required_column_checks():

            if 'time_bin_start' in self.colnames:
                self.remove_column('time_bin_start')

            if 'time_bin_size' in self.colnames:
                self.remove_column('time_bin_size')

            self.add_column(time_bin_start, index=0, name='time_bin_start')
            self.add_index('time_bin_start')
            self.add_column(time_bin_size, index=1, name='time_bin_size')

    @property
    def time_bin_start(self):
        """
        The start times of all the time bins.
        """
        return self['time_bin_start']

    @property
    def time_bin_center(self):
        """
        The center times of all the time bins.
        """
        return self['time_bin_start'] + self['time_bin_size'] * 0.5

    @property
    def time_bin_end(self):
        """
        The end times of all the time bins.
        """
        return self['time_bin_start'] + self['time_bin_size']

    @property
    def time_bin_size(self):
        """
        The sizes of all the time bins.
        """
        return self['time_bin_size']

    def __getitem__(self, item):
        if self._is_list_or_tuple_of_str(item):
            if 'time_bin_start' not in item or 'time_bin_size' not in item:
                out = QTable([self[x] for x in item],
                             meta=deepcopy(self.meta),
                             copy_indices=self._copy_indices)
                out._groups = groups.TableGroups(out, indices=self.groups._indices,
                                                 keys=self.groups._keys)
                return out
        return super().__getitem__(item)

    @classmethod
    def read(self, filename, time_bin_start_column=None, time_bin_end_column=None,
             time_bin_size_column=None, time_bin_size_unit=None, time_format=None, time_scale=None,
             format=None, *args, **kwargs):
        """
        Read and parse a file and returns a `astropy.timeseries.BinnedTimeSeries`.

        This method uses the unified I/O infrastructure in Astropy which makes
        it easy to define readers/writers for various classes
        (https://docs.astropy.org/en/stable/io/unified.html). By default, this
        method will try and use readers defined specifically for the
        `astropy.timeseries.BinnedTimeSeries` class - however, it is also
        possible to use the ``format`` keyword to specify formats defined for
        the `astropy.table.Table` class - in this case, you will need to also
        provide the column names for column containing the start times for the
        bins, as well as other column names (see the Parameters section below
        for details)::

            >>> from astropy.timeseries.binned import BinnedTimeSeries
            >>> ts = BinnedTimeSeries.read('binned.dat', format='ascii.ecsv',
            ...                            time_bin_start_column='date_start',
            ...                            time_bin_end_column='date_end')  # doctest: +SKIP

        Parameters
        ----------
        filename : str
            File to parse.
        format : str
            File format specifier.
        time_bin_start_column : str
            The name of the column with the start time for each bin.
        time_bin_end_column : str, optional
            The name of the column with the end time for each bin. Either this
            option or ``time_bin_size_column`` should be specified.
        time_bin_size_column : str, optional
            The name of the column with the size for each bin. Either this
            option or ``time_bin_end_column`` should be specified.
        time_bin_size_unit : `astropy.units.Unit`, optional
            If ``time_bin_size_column`` is specified but does not have a unit
            set in the table, you can specify the unit manually.
        time_format : str, optional
            The time format for the start and end columns.
        time_scale : str, optional
            The time scale for the start and end columns.
        *args : tuple, optional
            Positional arguments passed through to the data reader.
        **kwargs : dict, optional
            Keyword arguments passed through to the data reader.

        Returns
        -------
        out : `astropy.timeseries.binned.BinnedTimeSeries`
            BinnedTimeSeries corresponding to the file.

        """

        try:

            # First we try the readers defined for the BinnedTimeSeries class
            return super().read(filename, format=format, *args, **kwargs)

        except TypeError:

            # Otherwise we fall back to the default Table readers

            if time_bin_start_column is None:
                raise ValueError("``time_bin_start_column`` should be provided since the default Table readers are being used.")
            if time_bin_end_column is None and time_bin_size_column is None:
                raise ValueError("Either `time_bin_end_column` or `time_bin_size_column` should be provided.")
            elif time_bin_end_column is not None and time_bin_size_column is not None:
                raise ValueError("Cannot specify both `time_bin_end_column` and `time_bin_size_column`.")

            table = Table.read(filename, format=format, *args, **kwargs)

            if time_bin_start_column in table.colnames:
                time_bin_start = Time(table.columns[time_bin_start_column],
                                      scale=time_scale, format=time_format)
                table.remove_column(time_bin_start_column)
            else:
                raise ValueError(f"Bin start time column '{time_bin_start_column}' not found in the input data.")

            if time_bin_end_column is not None:

                if time_bin_end_column in table.colnames:
                    time_bin_end = Time(table.columns[time_bin_end_column],
                                        scale=time_scale, format=time_format)
                    table.remove_column(time_bin_end_column)
                else:
                    raise ValueError(f"Bin end time column '{time_bin_end_column}' not found in the input data.")

                time_bin_size = None

            elif time_bin_size_column is not None:

                if time_bin_size_column in table.colnames:
                    time_bin_size = table.columns[time_bin_size_column]
                    table.remove_column(time_bin_size_column)
                else:
                    raise ValueError(f"Bin size column '{time_bin_size_column}' not found in the input data.")

                if time_bin_size.unit is None:
                    if time_bin_size_unit is None or not isinstance(time_bin_size_unit, u.UnitBase):
                        raise ValueError("The bin size unit should be specified as an astropy Unit using ``time_bin_size_unit``.")
                    time_bin_size = time_bin_size * time_bin_size_unit
                else:
                    time_bin_size = u.Quantity(time_bin_size)

                time_bin_end = None

            if time_bin_start.isscalar and time_bin_size.isscalar:
                return BinnedTimeSeries(data=table,
                                    time_bin_start=time_bin_start,
                                    time_bin_end=time_bin_end,
                                    time_bin_size=time_bin_size,
                                    n_bins=len(table))
            else:
                return BinnedTimeSeries(data=table,
                                    time_bin_start=time_bin_start,
                                    time_bin_end=time_bin_end,
                                    time_bin_size=time_bin_size)
