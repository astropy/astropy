# Licensed under a 3-clause BSD style license - see LICENSE.rst

from datetime import datetime

import numpy as np

from astropy import units as u
from astropy.time import Time

__all__ = ["time_support"]

__doctest_requires__ = {"time_support": ["matplotlib"]}

UNSUPPORTED_FORMATS = ("datetime", "datetime64")
YMDHMS_FORMATS = ("fits", "iso", "isot", "yday")
STR_FORMATS = YMDHMS_FORMATS + ("byear_str", "jyear_str")


def time_support(*, scale=None, format=None, simplify=True):
    """
    Enable support for plotting `astropy.time.Time` instances in
    matplotlib.

    May be (optionally) used with a ``with`` statement.

      >>> import matplotlib.pyplot as plt
      >>> from astropy import units as u
      >>> from astropy import visualization
      >>> with visualization.time_support():  # doctest: +IGNORE_OUTPUT
      ...     plt.figure()
      ...     plt.plot(Time(['2016-03-22T12:30:31', '2016-03-22T12:30:38', '2016-03-22T12:34:40']))
      ...     plt.draw()

    Parameters
    ----------
    scale : str, optional
        The time scale to use for the times on the axis. If not specified,
        the scale of the first Time object passed to Matplotlib is used.
    format : str, optional
        The time format to use for the times on the axis. If not specified,
        the format of the first Time object passed to Matplotlib is used.
    simplify : bool or str, optional
        If `True`, simplify labels by removing ``00:00:00.000`` time suffixes
        from ISO strings when all labels fall on that time. If ``'concise'``,
        use concise formatting that shows only the part of the label that
        varies across ticks, with the common prefix set as an axis offset label.
    """
    import matplotlib.units as units
    from matplotlib.ticker import MaxNLocator, ScalarFormatter

    from astropy.visualization.wcsaxes.utils import select_step_hour, select_step_scalar

    class AstropyTimeLocator(MaxNLocator):
        # Note: we default to AutoLocator since many time formats
        # can just use this.

        # Parameters controlling tick density.  Subclasses can override these
        # to change tick placement without duplicating tick_values().
        _nbins_default = 4
        # Ranges of at most this many years use monthly ticks instead of yearly.
        _year_threshold = 1
        # Divisor used when estimating the year / month step size.
        _step_divisor = 3
        # Divisor used when estimating the sub-day (hour/minute/second) step.
        _hour_divisor = 3

        def __init__(self, converter, *args, **kwargs):
            kwargs["nbins"] = self._nbins_default
            super().__init__(*args, **kwargs)
            self._converter = converter

        def tick_values(self, vmin, vmax):
            # Where we put the ticks depends on the format we are using
            if self._converter.format in YMDHMS_FORMATS:
                # If we are here, we need to check what the range of values
                # is and decide how to find tick locations accordingly

                vrange = vmax - vmin

                if (
                    self._converter.format != "yday" and vrange > 31
                ) or vrange > 366:  # greater than a month
                    # We need to be careful here since not all years and months have
                    # the same length

                    # Start off by converting the values from the range to
                    # datetime objects, so that we can easily extract the year and
                    # month.

                    tmin = Time(
                        vmin, scale=self._converter.scale, format="mjd"
                    ).datetime
                    tmax = Time(
                        vmax, scale=self._converter.scale, format="mjd"
                    ).datetime

                    # Find the range of years
                    ymin = tmin.year
                    ymax = tmax.year

                    if ymax > ymin + self._year_threshold:  # greater than a year
                        # Find the step we want to use
                        ystep = int(
                            select_step_scalar(
                                max(1, (ymax - ymin) / self._step_divisor)
                            )
                        )

                        ymin = ystep * (ymin // ystep)

                        # Generate the years for these steps
                        times = []
                        for year in range(ymin, ymax + 1, ystep):
                            times.append(datetime(year=year, month=1, day=1))

                    else:  # greater than a month but less than a year
                        mmin = tmin.month
                        mmax = tmax.month + 12 * (ymax - ymin)

                        mstep = int(
                            select_step_scalar(
                                max(1, (mmax - mmin) / self._step_divisor)
                            )
                        )

                        mmin = mstep * max(1, mmin // mstep)

                        # Generate the months for these steps
                        times = []
                        for month in range(mmin, mmax + 1, mstep):
                            times.append(
                                datetime(
                                    year=ymin + (month - 1) // 12,
                                    month=(month - 1) % 12 + 1,
                                    day=1,
                                )
                            )

                    # Convert back to MJD
                    values = Time(times, scale=self._converter.scale).mjd

                elif vrange > 1:  # greater than a day
                    self.set_params(steps=[1, 2, 5, 10])
                    values = super().tick_values(vmin, vmax)

                else:
                    # Determine ideal step
                    dv = (vmax - vmin) / self._hour_divisor * 24 << u.hourangle

                    # And round to nearest sensible value
                    dv = select_step_hour(dv).to_value(u.hourangle) / 24

                    # Determine tick locations
                    imin = np.ceil(vmin / dv)
                    imax = np.floor(vmax / dv)
                    values = np.arange(imin, imax + 1, dtype=np.int64) * dv

            else:
                values = super().tick_values(vmin, vmax)

            # Get rid of values outside of the input interval
            return values[(values >= vmin) & (values <= vmax)]

        def __call__(self):
            vmin, vmax = self.axis.get_view_interval()
            return self.tick_values(vmin, vmax)

    class AstropyConciseTimeLocator(AstropyTimeLocator):
        """Time locator for use with `~astropy.visualization.time_support` and
        ``simplify='concise'``.

        Produces more ticks than `AstropyTimeLocator` (targeting ~8 ticks
        rather than ~4) and subdivides multi-year spans into monthly ticks,
        matching the density of `~matplotlib.dates.AutoDateLocator`.
        """

        _nbins_default = 8
        _year_threshold = 3
        _step_divisor = 6
        _hour_divisor = 6

        def __init__(self, converter, *args, **kwargs):
            # Call MaxNLocator directly so _nbins_default from this subclass
            # is used rather than AstropyTimeLocator's value.
            MaxNLocator.__init__(self, nbins=self._nbins_default)
            self._converter = converter

    class AstropyTimeFormatter(ScalarFormatter):
        _MONTH_ABBR = [
            "Jan",
            "Feb",
            "Mar",
            "Apr",
            "May",
            "Jun",
            "Jul",
            "Aug",
            "Sep",
            "Oct",
            "Nov",
            "Dec",
        ]

        def __init__(self, converter, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self._converter = converter
            self.set_useOffset(False)
            self.set_scientific(False)
            self._offset_string = ""

        def get_offset(self):
            return self._offset_string

        def format_ticks(self, values):
            if len(values) == 0:
                return []
            if self._converter.format in YMDHMS_FORMATS:
                times = Time(values, format="mjd", scale=self._converter.scale)
                formatted = list(getattr(times, self._converter.format))
                if self._converter.simplify == "concise":
                    if self._converter.format == "yday":
                        return self._concise_yday(formatted)
                    return self._concise_calendar(formatted)
                self._offset_string = ""
                if self._converter.simplify:
                    if self._converter.format in ("fits", "iso", "isot"):
                        if all(x.endswith("00:00:00.000") for x in formatted):
                            split = " " if self._converter.format == "iso" else "T"
                            formatted = [x.split(split)[0] for x in formatted]
                    elif self._converter.format == "yday":
                        if all(x.endswith(":001:00:00:00.000") for x in formatted):
                            formatted = [x.split(":", 1)[0] for x in formatted]
                return formatted
            self._offset_string = ""
            if self._converter.format == "byear_str":
                return Time(
                    values, format="byear", scale=self._converter.scale
                ).byear_str
            if self._converter.format == "jyear_str":
                return Time(
                    values, format="jyear", scale=self._converter.scale
                ).jyear_str
            return super().format_ticks(values)

        def _concise_calendar(self, formatted):
            """Produce concise tick labels for iso/isot/fits formats.

            Shows only the component that varies across ticks and sets
            ``_offset_string`` to the common date/time prefix.
            """
            sep = " " if self._converter.format == "iso" else "T"
            date_parts = [f.split(sep)[0] for f in formatted]
            time_parts = [f.split(sep)[1] for f in formatted]

            years = [int(d[:4]) for d in date_parts]
            months = [int(d[5:7]) for d in date_parts]
            days = [int(d[8:10]) for d in date_parts]
            hours = [int(t[:2]) for t in time_parts]
            minutes = [int(t[3:5]) for t in time_parts]
            secs = [t[6:] for t in time_parts]

            at_year_tick = all(
                m == 1 and d == 1 and h == 0 and mi == 0 and float(s) == 0
                for m, d, h, mi, s in zip(months, days, hours, minutes, secs)
            )
            at_month_tick = all(
                d == 1 and h == 0 and mi == 0 and float(s) == 0
                for d, h, mi, s in zip(days, hours, minutes, secs)
            )
            at_day_tick = all(
                h == 0 and mi == 0 and float(s) == 0
                for h, mi, s in zip(hours, minutes, secs)
            )

            if at_year_tick:
                self._offset_string = ""
                return [str(y) for y in years]

            if at_month_tick:
                if len(set(years)) == 1:
                    self._offset_string = str(years[0])
                    return [self._MONTH_ABBR[m - 1] for m in months]
                else:
                    self._offset_string = ""
                    return [
                        f"{y} {self._MONTH_ABBR[m - 1]}" for y, m in zip(years, months)
                    ]

            if at_day_tick:
                if len(set(zip(years, months))) == 1:
                    self._offset_string = (
                        f"{years[0]} {self._MONTH_ABBR[months[0] - 1]}"
                    )
                    return [str(d) for d in days]
                else:
                    self._offset_string = ""
                    return [
                        f"{y}-{m:02d}-{d:02d}" for y, m, d in zip(years, months, days)
                    ]

            # Sub-day ticks
            same_date = len(set(date_parts)) == 1
            self._offset_string = date_parts[0] if same_date else ""
            if not same_date:
                return formatted

            all_zero_secs = all(float(s) == 0 for s in secs)
            all_zero_mins = all_zero_secs and all(mi == 0 for mi in minutes)

            if all_zero_mins:
                return [f"{h:02d}:00" for h in hours]
            if all_zero_secs:
                return [f"{h:02d}:{mi:02d}" for h, mi in zip(hours, minutes)]
            return [
                f"{h:02d}:{mi:02d}:{self._strip_sec(s)}"
                for h, mi, s in zip(hours, minutes, secs)
            ]

        def _concise_yday(self, formatted):
            """Produce concise tick labels for the yday format."""
            years = [int(f[:4]) for f in formatted]
            doys = [int(f[5:8]) for f in formatted]
            hours = [int(f[9:11]) for f in formatted]
            mins = [int(f[12:14]) for f in formatted]
            secs = [f[15:] for f in formatted]

            at_year_tick = all(
                doy == 1 and h == 0 and mi == 0 and float(s) == 0
                for doy, h, mi, s in zip(doys, hours, mins, secs)
            )
            at_day_tick = all(
                h == 0 and mi == 0 and float(s) == 0
                for h, mi, s in zip(hours, mins, secs)
            )

            if at_year_tick:
                self._offset_string = ""
                return [str(y) for y in years]

            if at_day_tick:
                if len(set(years)) == 1:
                    self._offset_string = str(years[0])
                    return [str(doy) for doy in doys]
                else:
                    self._offset_string = ""
                    return [f"{y}:{doy:03d}" for y, doy in zip(years, doys)]

            # Sub-day ticks
            same_date = len(set(zip(years, doys))) == 1
            if not same_date:
                self._offset_string = ""
                return formatted
            self._offset_string = f"{years[0]}:{doys[0]:03d}"

            all_zero_secs = all(float(s) == 0 for s in secs)
            all_zero_mins = all_zero_secs and all(mi == 0 for mi in mins)

            if all_zero_mins:
                return [f"{h:02d}:00" for h in hours]
            if all_zero_secs:
                return [f"{h:02d}:{mi:02d}" for h, mi in zip(hours, mins)]
            return [
                f"{h:02d}:{mi:02d}:{self._strip_sec(s)}"
                for h, mi, s in zip(hours, mins, secs)
            ]

        @staticmethod
        def _strip_sec(s):
            """Strip trailing zeros from a seconds string like '31.000'."""
            return s.rstrip("0").rstrip(".")

    class MplTimeConverter(units.ConversionInterface):
        def __init__(self, scale=None, format=None, simplify=None):
            super().__init__()

            self.format = format
            self.scale = scale
            self.simplify = simplify

            # Keep track of original converter in case the context manager is
            # used in a nested way.
            self._original_converter = units.registry.get(Time)

            units.registry[Time] = self

        @property
        def format(self):
            return self._format

        @format.setter
        def format(self, value):
            if value in UNSUPPORTED_FORMATS:
                raise ValueError(f"time_support does not support format={value}")
            self._format = value

        def __enter__(self):
            return self

        def __exit__(self, type, value, tb):
            if self._original_converter is None:
                del units.registry[Time]
            else:
                units.registry[Time] = self._original_converter

        def default_units(self, x, axis):
            if isinstance(x, tuple):
                x = x[0]
            if self.format is None:
                self.format = x.format
            if self.scale is None:
                self.scale = x.scale
            return "astropy_time"

        def convert(self, value, unit, axis):
            """
            Convert a Time value to a scalar or array.
            """
            scaled = getattr(value, self.scale)
            if self.format in YMDHMS_FORMATS:
                return scaled.mjd
            elif self.format == "byear_str":
                return scaled.byear
            elif self.format == "jyear_str":
                return scaled.jyear
            else:
                return getattr(scaled, self.format)

        def axisinfo(self, unit, axis):
            """
            Return major and minor tick locators and formatters.
            """
            if self.simplify == "concise":
                majloc = AstropyConciseTimeLocator(self)
            else:
                majloc = AstropyTimeLocator(self)
            majfmt = AstropyTimeFormatter(self)
            return units.AxisInfo(
                majfmt=majfmt, majloc=majloc, label=f"Time ({self.scale})"
            )

    return MplTimeConverter(scale=scale, format=format, simplify=simplify)
