# Licensed under a 3-clause BSD style license - see LICENSE.rst

from copy import deepcopy

import numpy as np

from astropy.table import groups, QTable
from astropy.time import Time, TimeDelta
from astropy import units as u
from astropy.units import Quantity

from .core import BaseTimeSeries

__all__ = ['BinnedTimeSeries']


class BinnedTimeSeries(BaseTimeSeries):

    _require_time_column = False

    def __init__(self, data=None, time_bin_start=None, time_bin_end=None,
                 time_bin_size=None, n_bins=None, **kwargs):

        super().__init__(data=data, **kwargs)

        # FIXME: this is because for some operations, an empty time series needs
        # to be created, then columns added one by one. We should check that
        # when columns are added manually, time is added first and is of the
        # right type.
        if (data is None and time_bin_start is None and time_bin_end is None and
                time_bin_size is None and n_bins is None):
            self._required_columns = ['time_bin_start', 'time_bin_size']
            return

        # First if time_bin_start and time_bin_end have been given in the table data, we
        # should extract them and treat them as if they had been passed as
        # keyword arguments.

        if 'time_bin_start' in self.colnames:
            if time_bin_start is None:
                time_bin_start = self.columns['time_bin_start']
                self.remove_column('time_bin_start')
            else:
                raise TypeError("'time_bin_start' has been given both in the table "
                                "and as a keyword argument")

        if 'time_bin_size' in self.colnames:
            if time_bin_size is None:
                time_bin_size = self.columns['time_bin_size']
                self.remove_column('time_bin_size')
            else:
                raise TypeError("'time_bin_size' has been given both in the table "
                                "and as a keyword argument")

        if time_bin_start is None:
            raise TypeError("'time_bin_start' has not been specified")

        if time_bin_end is None and time_bin_size is None:
            raise TypeError("Either 'time_bin_size' or 'time_bin_end' should be specified")

        if not isinstance(time_bin_start, Time):
            time_bin_start = Time(time_bin_start)

        if time_bin_end is not None and not isinstance(time_bin_end, Time):
            time_bin_end = Time(time_bin_end)

        if time_bin_size is not None and not isinstance(time_bin_size, (Quantity, TimeDelta)):
            raise TypeError("'time_bin_size' should be a Quantity or a TimeDelta")

        if isinstance(time_bin_size, TimeDelta):
            time_bin_size = time_bin_size.sec * u.s

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
                raise ValueError("Length of 'time_bin_start' ({0}) should match "
                                 "table length ({1})".format(len(time_bin_start), len(self)))

            if time_bin_end is not None:
                if time_bin_end.isscalar:
                    times = time_bin_start.copy()
                    times[:-1] = times[1:]
                    times[-1] = time_bin_end
                    time_bin_end = times
                time_bin_size = (time_bin_end - time_bin_start).sec * u.s

        self.add_column(time_bin_start, index=0, name='time_bin_start')
        self.add_index('time_bin_start')

        if time_bin_size.isscalar:
            time_bin_size = np.repeat(time_bin_size, len(self))

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
