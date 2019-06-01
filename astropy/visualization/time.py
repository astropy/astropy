# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from contextlib import contextmanager

import numpy as np

import matplotlib.units as units
from datetime import datetime

from astropy.time import Time
from astropy import units as u
from astropy.visualization.wcsaxes.utils import select_step_hour, select_step_scalar

from matplotlib.ticker import MaxNLocator, ScalarFormatter

__all__ = ['time_support']

YMDHMS_FORMATS = ('fits', 'iso', 'isot')


@contextmanager
def time_support(scale='utc', format='isot', simplify=True):
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
    simplify : bool, optional
        If possible, simplify labels, e.g. by removing 00:00:00.000 times from
        ISO strings if all labels fall on that time.
    """
    converter = MplTimeConverter(scale=scale, format=format, simplify=simplify)
    units.registry[Time] = converter
    yield converter
    units.registry.pop(Time)


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

            if vrange > 31:  # greater than a month

                # We need to be careful here since not all years and months have
                # the same length

                # Start off by converting the values from the range to
                # datetime objects, so that we can easily extract the year and
                # month.

                tmin = Time(vmin, scale=self._converter.scale, format='mjd').datetime
                tmax = Time(vmax, scale=self._converter.scale, format='mjd').datetime

                # Find the range of years
                ymin = tmin.year
                ymax = tmax.year

                if ymax > ymin + 1:  # greater than a year

                    # Find the step we want to use
                    ystep = int(select_step_scalar(max(1, (ymax - ymin) / 3)))

                    ymin = ystep * (ymin // ystep)

                    # Generate the years for these steps
                    times = []
                    for year in range(ymin, ymax + 1, ystep):
                        times.append(datetime(year=year, month=1, day=1))

                else:  # greater than a month but less than a year

                    mmin = tmin.month
                    mmax = tmax.month + 12 * (ymax - ymin)

                    mstep = int(select_step_scalar(max(1, (mmax - mmin) / 3)))

                    mmin = mstep * max(1, mmin // mstep)

                    # Generate the months for these steps
                    times = []
                    for month in range(mmin, mmax + 1, mstep):
                        times.append(datetime(year=ymin + month // 12,
                                              month=month % 12, day=1))

                # Convert back to MJD
                values = Time(times, scale=self._converter.scale).mjd

                # Get rid of values outside of the input interval
                values = values[(values >= vmin) & (values <= vmax)]

                return values

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
            formatted = getattr(times, self._converter.format)
            if self._converter.simplify:
                if all([x.endswith('00:00:00.000') for x in formatted]):
                    split = 'T' if self._converter.format == 'isot' else ' '
                    formatted = [x.split(split)[0] for x in formatted]
            return formatted
        else:
            return super().format_ticks(values)


class MplTimeConverter(units.ConversionInterface):

    def __init__(self, scale=None, format=None, simplify=None):

        super().__init__()

        self.format = format
        self.scale = scale
        self.simplify = simplify

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
                              majloc=majloc)

    @staticmethod
    def default_units(x, axis):
        'Return the default unit for x or None'
        return 'date'
