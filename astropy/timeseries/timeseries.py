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

    Row = TimeSeriesRow
    Column = TimeSeriesColumn
    MaskedColumn = TimeSeriesMaskedColumn
    TableColumns = TimeSeriesTableColumns

    def __init__():
        """
        """
        pass
