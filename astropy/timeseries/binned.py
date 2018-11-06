# Licensed under a 3-clause BSD style license - see LICENSE.rst

from copy import deepcopy

import numpy as np

from astropy.table import groups, QTable
from astropy.time import Time, TimeDelta
from astropy import units as u
from astropy.units import Quantity

from .core import TimeSeries

__all__ = ['BinnedTimeSeries']


class BinnedTimeSeries(TimeSeries):

    _require_time_column = False

    def __init__(self, data=None, start_time=None, end_time=None,
                 bin_size=None, n_bins=None, **kwargs):

        super().__init__(data=data, **kwargs)

        # FIXME: this is because for some operations, an empty time series needs
        # to be created, then columns added one by one. We should check that
        # when columns are added manually, time is added first and is of the
        # right type.
        if (data is None and start_time is None and end_time is None and
                bin_size is None and n_bins is None):
            self._required_columns = ['start_time', 'bin_size']
            return

        # First if start_time and end_time have been given in the table data, we
        # should extract them and treat them as if they had been passed as
        # keyword arguments.

        if 'start_time' in self.colnames:
            if start_time is None:
                start_time = self.columns['start_time']
                self.remove_column('start_time')
            else:
                raise TypeError("'start_time' has been given both in the table "
                                "and as a keyword argument")

        if 'bin_size' in self.colnames:
            if bin_size is None:
                bin_size = self.columns['bin_size']
                self.remove_column('bin_size')
            else:
                raise TypeError("'bin_size' has been given both in the table "
                                "and as a keyword argument")

        if start_time is None:
            raise TypeError("'start_time' has not been specified")

        if end_time is None and bin_size is None:
            raise TypeError("Either 'bin_size' or 'end_time' should be specified")

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

            if end_time is not None:
                if end_time.isscalar:
                    times = start_time.copy()
                    times[:-1] = times[1:]
                    times[-1] = end_time
                    end_time = times
                bin_size = (end_time - start_time).sec * u.s
            elif bin_size is None:
                raise TypeError("Either 'bin_size' or 'end_time' should be specified")

        self.add_column(start_time, index=0, name='start_time')
        self.add_index('start_time')

        if bin_size.isscalar:
            bin_size = np.repeat(bin_size, len(self))

        self.add_column(bin_size, index=1, name='bin_size')

    @property
    def start_time(self):
        return self['start_time']

    @property
    def end_time(self):
        return self['start_time'] + self['bin_size']

    @property
    def centre_time(self):
        return self['start_time'] + self['bin_size'] * 0.5

    def __getitem__(self, item):
        if self._is_list_or_tuple_of_str(item):
            if 'start_time' not in item or 'bin_size' not in item:
                out = QTable([self[x] for x in item],
                             meta=deepcopy(self.meta),
                             copy_indices=self._copy_indices)
                out._groups = groups.TableGroups(out, indices=self.groups._indices,
                                                 keys=self.groups._keys)
                return out
        return super().__getitem__(item)
