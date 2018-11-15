# Licensed under a 3-clause BSD style license - see LICENSE.rst

from copy import deepcopy

import numpy as np

from astropy.table import groups, QTable, Table
from astropy.time import Time, TimeDelta
from astropy import units as u
from astropy.units import Quantity

from .core import TimeSeries

__all__ = ['SampledTimeSeries']


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
        df = Table(table).to_pandas()

        # Set index
        df.set_index(self.time.datetime64, inplace=True)

        return df
