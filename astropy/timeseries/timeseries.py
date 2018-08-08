# Licensed under a 3-clause BSD style license - see LICENSE.rst

from ..table import QTable, Column, Row, MaskedColumn, TableColumns
from ..time import Time, TimeDelta


__all__ = ['TimeSeries']


class TimeSeriesColumn(Column):
    pass


class TimeSeriesRow(Row):
    pass


class TimeSeriesMaskedColumn(MaskedColumn):
    pass


class TimeSeriesTableColumns(TableColumns):
    pass


class TimeSeries(QTable):

    # Row = TimeSeriesRow
    # Column = TimeSeriesColumn
    # MaskedColumn = TimeSeriesMaskedColumn
    # TableColumns = TimeSeriesTableColumns

    def __init__(self, data=None, time=None, time_delta=None, **kwargs):
        """
        """
        if not (isinstance(time, (Time, TimeDelta)) or not None):
            raise ValueError("'time' should be Time or TimeDelta or provided in 'data'")

        super().__init__(data=data, **kwargs)


        if 'time' in self.colnames and time is not None:
                raise ValueError("'time' is ambiguous, it has been provided both "
                                 "in the data and in the time arguments.")
        if time is None:
            if 'time' in self.colnames:
                # TODO: design decision: is 'time' always the first column?
                time = self.columns['time']
                self.time = time
            else:
                # TODO: figure out what to do with empty TimeSeries, Time() is not an option
                self.time = None
        else:
            self.time = time

            if self.time.info.name is None:
                self.time.info.name = 'time'

            self.columns['time'] = time
