# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

from ..table import QTable, Column, Row, MaskedColumn, TableColumns
from ..time import Time, TimeDelta
from .. import units as u
from ..units import Quantity


__all__ = ['TimeSeries', 'SampledTimeSeries', 'BinnedTimeSeries']


class TimeSeries(QTable):
    pass


class SampledTimeSeries(TimeSeries):

    def __init__(self, data=None, time=None, start_time=None, end_time=None, time_delta=None, **kwargs):
        """
        """

        super().__init__(data=data, **kwargs)



        if time is None:

            if 'time' not in self.colnames:
                raise ValueError("Expected a 'time' column in the time series")

        else:

            if 'time' in self.colnames:
                raise ValueError("'time' has been given both in the table and "
                                 "as a keyword argument.")

            if time.info.name is None:
                time.info.name = 'time'

            self.columns['time'] = time

        if not isinstance(self['time'], (Time, TimeDelta)):
            raise ValueError("The 'time' column should be a Time or TimeDelta object")

        # TODO: design decision: is 'time' always the first column?
        self.time = self.columns['time']


class BinnedTimeSeries(TimeSeries):

    def __init__(self, data=None, start_time=None, end_time=None, bin_size=None, **kwargs):

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

        if bin_size is not None and not isinstance(bin_size, Quantity):
            raise TypeError("'bin_size' should be a Quantity")

        if start_time.isscalar:

            # We interpret this as meaning that this is the start of the
            # first bin and that the bins are contiguous. In this case,
            # we require bin_size to be specified.

            if bin_size is None:
                raise TypeError("'start_time' is scalar, so 'bin_size' is required")

            if bin_size.isscalar:
                bin_size = np.repeat(bin_size, len(self))

            time_delta = np.cumsum(bin_size)
            end_time = start_time + time_delta

            # Now shift the array so that the first entry is 0
            time_delta = np.roll(time_delta, 1)
            time_delta[0] = 0. * u.s

            # Make start_time into an array
            start_time = start_time + time_delta

        else:

            if len(start_time) != len(self):
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

    @property
    def start_time(self):
        return self['start_time']

    @property
    def end_time(self):
        return self['end_time']
