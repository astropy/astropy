# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

import matplotlib.units as units
import matplotlib.pyplot as plt

from astropy.time import Time
from astropy import units as u
from astropy.visualization.wcsaxes.utils import select_step_hour

from matplotlib.ticker import MaxNLocator, ScalarFormatter

__all__ = ['time_support']

YMDHMS_FORMATS = ('fits', 'iso', 'isot')


def time_support(scale='utc', format='isot'):
    """
    Enable support for plotting `astropy.time.Time` instances in
    matplotlib.

    May be (optionally) used with a ``with`` statement.

      >>> import matplotlib.pyplot as plt
      >>> from astropy import units as u
      >>> from astropy import visualization
      >>> with visualization.time_support():
      ...     plt.figure()
      ...     plt.plot(Time(['2016-03-22T12:30:31', '2016-03-22T12:30:38', '2016-03-22T12:34:40']))
      ...     plt.draw()

    Parameters
    ----------
    format : str, optional
        The time format to use for the times on the axis
    scale : str, optional
        The time scale to use for the times on the axis

    """
    return MplTimeConverter(scale=scale, format=format)


class AstropyTimeLocator(MaxNLocator):

    # Note: we default to AutoLocator since many time formats
    # can just use this.

    def __init__(self, converter, *args, **kwargs):
        kwargs['nbins'] = 4
        super().__init__(*args, **kwargs)
        self._converter = converter

    def tick_values(self, vmin, vmax):

        # Where we put the ticks depends on the format we are using
        if self._converter.format in YMDHMS_FORMATS:

            # If we are here, we need to check what the range of values
            # is and decide how to find tick locations accordingly

            vrange = vmax - vmin

            if vrange > 366:   # greater than a year

                # We need to be careful here since not all years have the
                # same length

                pass

            elif vrange > 31:  # greater than a month

                # We need to be careful here since not all months have the
                # same length

                pass

            elif vrange > 1:  # greater than a day

                self.set_params(steps=[1, 2, 5, 10])
                return super().tick_values(vmin, vmax)

            else:

                # Determine ideal step
                dv = (vmax - vmin) / 3 * 24 << u.hourangle

                # And round to nearest sensible value
                dv = select_step_hour(dv).to_value(u.hourangle) / 24

                # Determine tick locations
                imin = np.ceil(vmin / dv)
                imax = np.floor(vmax / dv)
                return np.arange(imin, imax + 1, dtype=int) * dv

        else:

            return super().tick_values(vmin, vmax)

    def __call__(self):
        vmin, vmax = self.axis.get_view_interval()
        return self.tick_values(vmin, vmax)


class AstropyTimeFormatter(ScalarFormatter):

    def __init__(self, converter, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._converter = converter
        self.set_useOffset(False)
        self.set_scientific(False)

    def format_ticks(self, values):
        if self._converter.format in YMDHMS_FORMATS:
            if len(values) == 0:
                return []
            times = Time(values, format='mjd', scale=self._converter.scale)
            return getattr(times, self._converter.format)
        else:
            return super().format_ticks(values)


class MplTimeConverter(units.ConversionInterface):

    def __init__(self, scale=None, format=None):

        super().__init__()

        self._format = format
        self._scale = scale

        if Time not in units.registry:
            units.registry[Time] = self
            self._remove = True
        else:
            self._remove = False

    def __enter__(self):
        return self

    def __exit__(self, type, value, tb):
        if self._remove:
            del units.registry[Time]

    @property
    def format(self):
        return self._format

    @format.setter
    def format(self, value):
        self._format = value

    @property
    def scale(self):
        return self._scale

    @scale.setter
    def scale(self, value):
        self._scale = value

    def convert(self, value, unit, axis):
        'Convert a Time value to a scalar or array'
        scaled = getattr(value, self.scale)
        if self.format in YMDHMS_FORMATS:
            return scaled.mjd
        else:
            return getattr(scaled, self.format)

    def axisinfo(self, unit, axis):
        'Return major and minor tick locators and formatters'
        if unit != 'date':
            return None
        majloc = AstropyTimeLocator(self)
        majfmt = AstropyTimeFormatter(self)
        return units.AxisInfo(majfmt=majfmt,
                              majloc=majloc,
                              label='blah')

    @staticmethod
    def default_units(x, axis):
        'Return the default unit for x or None'
        return 'date'
