# Licensed under a 3-clause BSD style license - see LICENSE.rst

from copy import deepcopy

import numpy as np

from ..table import groups, QTable
from ..time import Time, TimeDelta
from .. import units as u
from ..units import Quantity


__all__ = ['TimeSeries', 'SampledTimeSeries', 'BinnedTimeSeries']


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
