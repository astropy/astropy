# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The astropy.utils.iers package provides access to the tables provided by
the International Earth Rotation and Reference Systems Service, in
particular allowing interpolation of published UT1-UTC values for given
times.  These are used in `astropy.time` to provide UT1 values.  The polar
motions are also used for determining earth orientation for
celestial-to-terrestrial coordinate transformations
(in `astropy.coordinates`).
"""

import os
import re
from datetime import datetime, timedelta
from typing import Self
from urllib.parse import urlparse
from warnings import warn

import erfa
import numpy as np
from astropy_iers_data import (
    IERS_A_FILE,
    IERS_A_README,
    IERS_A_URL,
    IERS_A_URL_MIRROR,
    IERS_B_FILE,
    IERS_B_README,
    IERS_B_URL,
)
from astropy_iers_data import IERS_LEAP_SECOND_FILE as _IERS_LEAP_SECOND_FILE
from astropy_iers_data import IERS_LEAP_SECOND_URL as _IERS_LEAP_SECOND_URL
from astropy_iers_data import IERS_LEAP_SECOND_URL_MIRROR as _IETF_LEAP_SECOND_URL

from astropy import config as _config
from astropy import units as u
from astropy import utils
from astropy.table import MaskedColumn, QTable
from astropy.time import Time
from astropy.utils.data import (
    get_readable_fileobj,
)
from astropy.utils.decorators import deprecated
from astropy.utils.exceptions import AstropyDeprecationWarning, AstropyWarning
from astropy.utils.state import ScienceState

__all__ = [
    "FROM_IERS_A",
    "FROM_IERS_A_PREDICTION",
    "FROM_IERS_B",
    "IERS",
    "IERS_A",
    "IERS_A_FILE",
    "IERS_A_README",
    "IERS_A_URL",
    "IERS_A_URL_MIRROR",
    "IERS_B",
    "IERS_B_FILE",
    "IERS_B_README",
    "IERS_B_URL",
    "TIME_BEFORE_IERS_RANGE",
    "TIME_BEYOND_IERS_RANGE",
    "Conf",
    "IERSDegradedAccuracyWarning",
    "IERSRangeError",
    "IERSStaleWarning",
    "IERSWarning",
    "IERS_Auto",
    "LeapSeconds",
    "conf",
    "earth_orientation_table",
]

# Status/source values returned by IERS.ut1_utc
FROM_IERS_B = 0
FROM_IERS_A = 1
FROM_IERS_A_PREDICTION = 2
TIME_BEFORE_IERS_RANGE = -1
TIME_BEYOND_IERS_RANGE = -2

MJD_ZERO = 2400000.5

INTERPOLATE_ERROR = """\
interpolating from IERS_Auto using predictive values that are more
than {0} days old.

Normally you should not see this error because this class
automatically downloads the latest IERS-A table.  Perhaps you are
offline?  If you understand what you are doing then this error can be
suppressed by setting the auto_max_age configuration variable to
``None``:

  from astropy.utils.iers import conf
  conf.auto_max_age = None
"""

MONTH_ABBR = [
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


class IERSWarning(AstropyWarning):
    """
    Generic warning class for IERS.
    """


class IERSDegradedAccuracyWarning(AstropyWarning):
    """
    IERS time conversion has degraded accuracy normally due to setting
    ``conf.auto_download = False`` and ``conf.iers_degraded_accuracy = 'warn'``.
    """


class IERSStaleWarning(IERSWarning):
    """
    Downloaded IERS table may be stale.
    """


def download_file(*args, **kwargs):
    """
    Overload astropy.utils.data.download_file within iers module to use a
    custom (longer) wait time.  This just passes through ``*args`` and
    ``**kwargs`` after temporarily setting the download_file remote timeout to
    the local ``iers.conf.remote_timeout`` value.
    """
    kwargs.setdefault(
        "http_headers",
        {
            "User-Agent": "astropy/iers",
            "Accept": "*/*",
        },
    )

    with utils.data.conf.set_temp("remote_timeout", conf.remote_timeout):
        return utils.data.download_file(*args, **kwargs)


def _none_to_float(value):
    """
    Convert None to a valid floating point value.  Especially
    for auto_max_age = None.
    """
    return value if value is not None else np.finfo(float).max


_DEPRECATED_LEAP_SECOND_CONSTANTS = {
    "IERS_LEAP_SECOND_FILE": _IERS_LEAP_SECOND_FILE,
    "IERS_LEAP_SECOND_URL": _IERS_LEAP_SECOND_URL,
    "IETF_LEAP_SECOND_URL": _IETF_LEAP_SECOND_URL,
}


def __getattr__(name):
    if name in _DEPRECATED_LEAP_SECOND_CONSTANTS:
        warn(
            f"{name} is deprecated since 8.1 and may be removed in a future "
            "version. Leap seconds have been discontinued, so the table built "
            "into ERFA is final; use LeapSeconds.from_erfa() to get it.",
            AstropyDeprecationWarning,
            stacklevel=2,
        )
        return _DEPRECATED_LEAP_SECOND_CONSTANTS[name]

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.utils.iers`.
    """

    auto_download = _config.ConfigItem(
        True,
        "Enable auto-downloading of the latest IERS-A data.  If set to False "
        "then the bundled IERS-A file will be used by default (even if a "
        "newer version of the IERS-A file was previously downloaded and cached). "
        "Default is True.",
    )
    auto_max_age = _config.ConfigItem(
        30.0,
        "Maximum age (days) of predictive data before auto-downloading. "
        'See "Auto refresh behavior" in astropy.utils.iers documentation for details. '
        "Default is 30.",
    )
    iers_auto_url = _config.ConfigItem(
        IERS_A_URL, "URL for auto-downloading IERS file data."
    )
    iers_auto_url_mirror = _config.ConfigItem(
        IERS_A_URL_MIRROR, "Mirror URL for auto-downloading IERS file data."
    )
    remote_timeout = _config.ConfigItem(
        10.0, "Remote timeout downloading IERS file data (seconds)."
    )
    iers_degraded_accuracy = _config.ConfigItem(
        ["error", "warn", "ignore"],
        "IERS behavior if the range of available IERS data does not "
        "cover the times when converting time scales, potentially leading "
        "to degraded accuracy.",
    )
    # The three leap-second items below are deprecated since 8.1 and unused:
    # leap seconds have been discontinued, so the table is never updated.
    system_leap_second_file = _config.ConfigItem(
        "", "Deprecated and unused. System file with leap seconds."
    )
    iers_leap_second_auto_url = _config.ConfigItem(
        _IERS_LEAP_SECOND_URL,
        "Deprecated and unused. URL for auto-downloading leap seconds.",
    )
    ietf_leap_second_auto_url = _config.ConfigItem(
        _IETF_LEAP_SECOND_URL,
        "Deprecated and unused. Alternate URL for auto-downloading leap seconds.",
    )


conf = Conf()


class IERSRangeError(IndexError):
    """
    Any error for when dates are outside of the valid range for IERS.
    """


class IERS(QTable):
    """Generic IERS table class, defining interpolation functions.

    Sub-classed from `astropy.table.QTable`.  The table should hold columns
    'MJD', 'UT1_UTC', 'dX_2000A'/'dY_2000A', and 'PM_x'/'PM_y'.
    """

    iers_table = None
    """Cached table, returned if ``open`` is called without arguments."""

    @classmethod
    def open(cls, file=None, cache=False, **kwargs):
        """Open an IERS table, reading it from a file if not loaded before.

        Parameters
        ----------
        file : str or None
            full local or network path to the ascii file holding IERS data,
            for passing on to the ``read`` class methods (further optional
            arguments that are available for some IERS subclasses can be added).
            If None, use the default location from the ``read`` class method.
        cache : bool
            Whether to use cache. Defaults to False, since IERS files
            are regularly updated.

        Returns
        -------
        IERS
            An IERS table class instance

        Notes
        -----
        On the first call in a session, the table will be memoized (in the
        ``iers_table`` class attribute), and further calls to ``open`` will
        return this stored table if ``file=None`` (the default).

        If a table needs to be re-read from disk, pass on an explicit file
        location or use the (sub-class) close method and re-open.

        If the location is a network location it is first downloaded via
        download_file.

        For the IERS class itself, an IERS_B sub-class instance is opened.

        """
        if file is not None or cls.iers_table is None:
            if file is not None:
                if urlparse(file).netloc:
                    kwargs.update(file=download_file(file, cache=cache))
                else:
                    kwargs.update(file=file)

            # TODO: the below is really ugly and probably a bad idea.  Instead,
            # there should probably be an IERSBase class, which provides
            # useful methods but cannot really be used on its own, and then
            # *perhaps* an IERS class which provides best defaults.  But for
            # backwards compatibility, we use the IERS_B reader for IERS here.
            if cls is IERS:
                cls.iers_table = IERS_B.read(**kwargs)
            else:
                cls.iers_table = cls.read(**kwargs)
        return cls.iers_table

    @classmethod
    def close(cls):
        """Remove the IERS table from the class.

        This allows the table to be re-read from disk during one's session
        (e.g., if one finds it is out of date and has updated the file).
        """
        cls.iers_table = None

    def mjd_utc(self, jd1, jd2=0.0):
        """Turn a time to MJD, returning integer and fractional parts.

        Parameters
        ----------
        jd1 : float, array, or `~astropy.time.Time`
            first part of two-part JD, or Time object
        jd2 : float or array, optional
            second part of two-part JD.
            Default is 0., ignored if jd1 is `~astropy.time.Time`.

        Returns
        -------
        mjd : float or array
            integer part of MJD
        utc : float or array
            fractional part of MJD
        """
        try:  # see if this is a Time object
            jd1, jd2 = jd1.utc.jd1, jd1.utc.jd2
        except Exception:
            pass

        mjd = np.floor(jd1 - MJD_ZERO + jd2)
        utc = jd1 - (MJD_ZERO + mjd) + jd2
        return mjd, utc

    def ut1_utc(self, jd1, jd2=0.0, return_status=False):
        """Interpolate UT1-UTC corrections in IERS Table for given dates.

        Parameters
        ----------
        jd1 : float, array of float, or `~astropy.time.Time` object
            first part of two-part JD, or Time object
        jd2 : float or float array, optional
            second part of two-part JD.
            Default is 0., ignored if jd1 is `~astropy.time.Time`.
        return_status : bool
            Whether to return status values.  If False (default),
            raise ``IERSRangeError`` if any time is out of the range covered
            by the IERS table.

        Returns
        -------
        ut1_utc : float or float array
            UT1-UTC, interpolated in IERS Table
        status : int or int array
            Status values (if ``return_status``=``True``)::
            ``iers.FROM_IERS_B``
            ``iers.FROM_IERS_A``
            ``iers.FROM_IERS_A_PREDICTION``
            ``iers.TIME_BEFORE_IERS_RANGE``
            ``iers.TIME_BEYOND_IERS_RANGE``
        """
        return self._interpolate(
            jd1, jd2, ["UT1_UTC"], self.ut1_utc_source if return_status else None
        )

    def dcip_xy(self, jd1, jd2=0.0, return_status=False):
        """Interpolate CIP corrections in IERS Table for given dates.

        Parameters
        ----------
        jd1 : float, array of float, or `~astropy.time.Time` object
            first part of two-part JD, or Time object
        jd2 : float or float array, optional
            second part of two-part JD (default 0., ignored if jd1 is Time)
        return_status : bool
            Whether to return status values.  If False (default),
            raise ``IERSRangeError`` if any time is out of the range covered
            by the IERS table.

        Returns
        -------
        D_x : `~astropy.units.Quantity` ['angle']
            x component of CIP correction for the requested times.
        D_y : `~astropy.units.Quantity` ['angle']
            y component of CIP correction for the requested times
        status : int or int array
            Status values (if ``return_status``=``True``)::
            ``iers.FROM_IERS_B``
            ``iers.FROM_IERS_A``
            ``iers.FROM_IERS_A_PREDICTION``
            ``iers.TIME_BEFORE_IERS_RANGE``
            ``iers.TIME_BEYOND_IERS_RANGE``
        """
        return self._interpolate(
            jd1,
            jd2,
            ["dX_2000A", "dY_2000A"],
            self.dcip_source if return_status else None,
        )

    def pm_xy(self, jd1, jd2=0.0, return_status=False):
        """Interpolate polar motions from IERS Table for given dates.

        Parameters
        ----------
        jd1 : float, array of float, or `~astropy.time.Time` object
            first part of two-part JD, or Time object
        jd2 : float or float array, optional
            second part of two-part JD.
            Default is 0., ignored if jd1 is `~astropy.time.Time`.
        return_status : bool
            Whether to return status values.  If False (default),
            raise ``IERSRangeError`` if any time is out of the range covered
            by the IERS table.

        Returns
        -------
        PM_x : `~astropy.units.Quantity` ['angle']
            x component of polar motion for the requested times.
        PM_y : `~astropy.units.Quantity` ['angle']
            y component of polar motion for the requested times.
        status : int or int array
            Status values (if ``return_status``=``True``)::
            ``iers.FROM_IERS_B``
            ``iers.FROM_IERS_A``
            ``iers.FROM_IERS_A_PREDICTION``
            ``iers.TIME_BEFORE_IERS_RANGE``
            ``iers.TIME_BEYOND_IERS_RANGE``
        """
        return self._interpolate(
            jd1, jd2, ["PM_x", "PM_y"], self.pm_source if return_status else None
        )

    def _check_interpolate_indices(self, indices_orig, indices_clipped, max_input_mjd):
        """
        Check that the indices from interpolation match those after clipping
        to the valid table range.  This method gets overridden in the IERS_Auto
        class because it has different requirements.
        """
        if np.any(indices_orig != indices_clipped):
            if conf.iers_degraded_accuracy == "error":
                msg = (
                    "(some) times are outside of range covered by IERS table. Cannot"
                    " convert with full accuracy. To allow conversion with degraded"
                    " accuracy set astropy.utils.iers.conf.iers_degraded_accuracy to"
                    ' "warn" or "silent". For more information about setting this'
                    " configuration parameter or controlling its value globally, see"
                    " the Astropy configuration system documentation"
                    " https://docs.astropy.org/en/stable/config/index.html."
                )
                raise IERSRangeError(msg)
            elif conf.iers_degraded_accuracy == "warn":
                # No IERS data covering the time(s) and user requested a warning.
                msg = (
                    "(some) times are outside of range covered by IERS table, "
                    "accuracy is degraded."
                )
                warn(msg, IERSDegradedAccuracyWarning)
            # No IERS data covering the time(s) and user is OK with no warning.

    def _interpolate(self, jd1, jd2, columns, source=None):
        mjd, utc = self.mjd_utc(jd1, jd2)
        # enforce array
        is_scalar = not hasattr(mjd, "__array__") or mjd.ndim == 0
        if is_scalar:
            mjd = np.array([mjd])
            utc = np.array([utc])

        self._refresh_table_as_needed(mjd)

        # For typical format, will always find a match (since MJD are integer)
        # hence, important to define which side we will be; this ensures
        # self['MJD'][i-1]<=mjd<self['MJD'][i]
        i = np.searchsorted(self["MJD"].value, mjd, side="right")

        # Get index to MJD at or just below given mjd, clipping to ensure we
        # stay in range of table (status will be set below for those outside)
        i1 = np.clip(i, 1, len(self) - 1)
        i0 = i1 - 1
        mjd_0, mjd_1 = self["MJD"][i0].value, self["MJD"][i1].value
        results = []
        for column in columns:
            val_0, val_1 = self[column][i0], self[column][i1]
            d_val = val_1 - val_0
            if column == "UT1_UTC":
                # Check & correct for possible leap second (correcting diff.,
                # not 1st point, since jump can only happen right at 2nd point)
                d_val -= d_val.round()
            # Linearly interpolate (which is what TEMPO does for UT1-UTC, but
            # may want to follow IERS gazette #13 for more precise
            # interpolation and correction for tidal effects;
            # https://maia.usno.navy.mil/iers-gaz13)
            val = val_0 + (mjd - mjd_0 + utc) / (mjd_1 - mjd_0) * d_val

            # Do not extrapolate outside range, instead just propagate last values.
            val[i == 0] = self[column][0]
            val[i == len(self)] = self[column][-1]

            if is_scalar:
                val = val[0]

            results.append(val)

        if source:
            # Set status to source, using the routine passed in.
            status = source(i1)
            # Check for out of range
            status[i == 0] = TIME_BEFORE_IERS_RANGE
            status[i == len(self)] = TIME_BEYOND_IERS_RANGE
            if is_scalar:
                status = status[0]
            results.append(status)
            return results
        else:
            # Pass in initial to np.max to allow things to work for empty mjd.
            self._check_interpolate_indices(i1, i, np.max(mjd, initial=50000))
            return results[0] if len(results) == 1 else results

    def _refresh_table_as_needed(self, mjd):
        """
        Potentially update the IERS table in place depending on the requested
        time values in ``mdj`` and the time span of the table.  The base behavior
        is not to update the table.  ``IERS_Auto`` overrides this method.
        """

    def ut1_utc_source(self, i):
        """Source for UT1-UTC.  To be overridden by subclass."""
        return np.zeros_like(i)

    def dcip_source(self, i):
        """Source for CIP correction.  To be overridden by subclass."""
        return np.zeros_like(i)

    def pm_source(self, i):
        """Source for polar motion.  To be overridden by subclass."""
        return np.zeros_like(i)

    @property
    @deprecated(since="8.0", alternative="Time.now()", obj_type="property")
    def time_now(self):
        """
        Property to provide the current time, but also allow for explicitly setting
        the _time_now attribute for testing purposes.
        """
        try:
            return self._time_now
        except Exception:
            return Time.now()

    def _convert_col_for_table(self, col):
        # Fill masked columns with units to avoid dropped-mask warnings
        # when converting to Quantity.
        # TODO: Once we support masked quantities, we can drop this and
        # in the code below replace b_bad with table['UT1_UTC_B'].mask, etc.
        if getattr(col, "unit", None) is not None and isinstance(col, MaskedColumn):
            col = col.filled(np.nan)

        return super()._convert_col_for_table(col)


class IERS_A(IERS):
    """IERS Table class targeted to IERS A, provided by USNO.

    These include rapid turnaround and predicted times.
    See https://datacenter.iers.org/eop.php

    Notes
    -----
    If the package IERS A file (``iers.IERS_A_FILE``) is out of date, a new
    version can be downloaded from ``iers.IERS_A_URL`` or ``iers.IERS_A_URL_MIRROR``.
    See ``iers.__doc__`` for instructions on use in ``Time``, etc.
    """

    iers_table = None

    @classmethod
    def _combine_a_b_columns(cls, iers_a):
        """
        Return a new table with appropriate combination of IERS_A and B columns.
        """
        # IERS A has some rows at the end that hold nothing but dates & MJD
        # presumably to be filled later.  Exclude those a priori -- there
        # should at least be a predicted UT1-UTC and PM!
        table = iers_a[np.isfinite(iers_a["UT1_UTC_A"]) & (iers_a["PolPMFlag_A"] != "")]

        # This does nothing for IERS_A, but allows IERS_Auto to ensure the
        # IERS B values in the table are consistent with the true ones.
        table = cls._substitute_iers_b(table)

        # Combine A and B columns, using B where possible.
        b_bad = np.isnan(table["UT1_UTC_B"])
        table["UT1_UTC"] = np.where(b_bad, table["UT1_UTC_A"], table["UT1_UTC_B"])
        table["UT1Flag"] = np.where(b_bad, table["UT1Flag_A"], "B")
        # Repeat for polar motions.
        b_bad = np.isnan(table["PM_X_B"]) | np.isnan(table["PM_Y_B"])
        table["PM_x"] = np.where(b_bad, table["PM_x_A"], table["PM_X_B"])
        table["PM_y"] = np.where(b_bad, table["PM_y_A"], table["PM_Y_B"])
        table["PolPMFlag"] = np.where(b_bad, table["PolPMFlag_A"], "B")

        b_bad = np.isnan(table["dX_2000A_B"]) | np.isnan(table["dY_2000A_B"])
        table["dX_2000A"] = np.where(b_bad, table["dX_2000A_A"], table["dX_2000A_B"])
        table["dY_2000A"] = np.where(b_bad, table["dY_2000A_A"], table["dY_2000A_B"])
        table["NutFlag"] = np.where(b_bad, table["NutFlag_A"], "B")

        # Get the table index for the first row that has predictive values
        # PolPMFlag_A  IERS (I) or Prediction (P) flag for
        #              Bull. A polar motion values
        # UT1Flag_A    IERS (I) or Prediction (P) flag for
        #              Bull. A UT1-UTC values
        # Since only 'P' and 'I' are possible and 'P' is guaranteed to come
        # after 'I', we can use searchsorted for 100 times speed up over
        # finding the first index where the flag equals 'P'.
        p_index = min(
            np.searchsorted(table["UT1Flag_A"], "P"),
            np.searchsorted(table["PolPMFlag_A"], "P"),
        )
        table.meta["predictive_index"] = p_index
        table.meta["predictive_mjd"] = table["MJD"][p_index].value

        return table

    @classmethod
    def _substitute_iers_b(cls, table):
        # See documentation in IERS_Auto.
        return table

    @classmethod
    def read(
        cls,
        file: str | os.PathLike[str] | None = None,
        readme: str | os.PathLike[str] | None = None,
    ) -> Self:
        """Read IERS-A table from a finals2000a.* file provided by USNO.

        Parameters
        ----------
        file : str or os.PathLike[str]
            full path to ascii file holding IERS-A data.
            Defaults to ``iers.IERS_A_FILE``.
        readme : str or os.PathLike[str]
            full path to ascii file holding CDS-style readme.
            Defaults to package version, ``iers.IERS_A_README``.

        Returns
        -------
        ``IERS_A`` class instance
        """
        if file is None:
            # In prior versions of astropy, read() would only work without a
            # file argument if a finals2000A.all file was present in the
            # current working directory. We can now use the versions from
            # astropy-iers-data but for backward-compatibility we first check
            # if there is a file in the current working directory and use that
            # if so, emitting a deprecation warning
            if os.path.exists("finals2000A.all"):
                file = "finals2000A.all"
                warn(
                    "The file= argument was not specified but "
                    "'finals2000A.all' is present in the current working "
                    "directory, so reading IERS data from that file. To "
                    "continue reading a local file from the current working "
                    "directory, specify file= explicitly otherwise a bundled "
                    "file will be used in future.",
                    AstropyDeprecationWarning,
                )
            else:
                file = IERS_A_FILE
        else:
            file = os.fspath(file)
        if readme is None:
            readme = IERS_A_README
        else:
            readme = os.fspath(readme)

        iers_a = super().read(file, format="cds", readme=readme)

        # Combine the A and B data for UT1-UTC and PM columns
        table = cls._combine_a_b_columns(iers_a)
        table.meta["data_path"] = file
        table.meta["readme_path"] = readme

        return table

    def ut1_utc_source(self, i):
        """Set UT1-UTC source flag for entries in IERS table."""
        ut1flag = self["UT1Flag"][i]
        source = np.ones_like(i) * FROM_IERS_B
        source[ut1flag == "I"] = FROM_IERS_A
        source[ut1flag == "P"] = FROM_IERS_A_PREDICTION
        return source

    def dcip_source(self, i):
        """Set CIP correction source flag for entries in IERS table."""
        nutflag = self["NutFlag"][i]
        source = np.ones_like(i) * FROM_IERS_B
        source[nutflag == "I"] = FROM_IERS_A
        source[nutflag == "P"] = FROM_IERS_A_PREDICTION
        return source

    def pm_source(self, i):
        """Set polar motion source flag for entries in IERS table."""
        pmflag = self["PolPMFlag"][i]
        source = np.ones_like(i) * FROM_IERS_B
        source[pmflag == "I"] = FROM_IERS_A
        source[pmflag == "P"] = FROM_IERS_A_PREDICTION
        return source


class IERS_B(IERS):
    """IERS Table class targeted to IERS B, provided by IERS itself.

    These are final values; see https://www.iers.org/IERS/EN/Home/home_node.html

    Notes
    -----
    If the package IERS B file (``iers.IERS_B_FILE``) is out of date, a new
    version can be downloaded from ``iers.IERS_B_URL``.

    See `~astropy.utils.iers.IERS_B.read` for instructions on how to read
    a pre-2023 style IERS B file (usually named ``eopc04_IAU2000.62-now``).
    """

    iers_table = None

    @classmethod
    def read(
        cls,
        file: str | os.PathLike[str] | None = None,
        readme: str | os.PathLike[str] | None = None,
        data_start: int = 6,
    ) -> Self:
        """Read IERS-B table from a eopc04.* file provided by IERS.

        Parameters
        ----------
        file : str or os.PathLike[str]
            full path to ascii file holding IERS-B data.
            Defaults to package version, ``iers.IERS_B_FILE``.
        readme : str or os.PathLike[str]
            full path to ascii file holding CDS-style readme.
            Defaults to package version, ``iers.IERS_B_README``.
        data_start : int
            Starting row. Default is 6, appropriate for standard IERS files.

        Returns
        -------
        ``IERS_B`` class instance

        Notes
        -----
        To read a pre-2023 style IERS B file (usually named something like
        ``eopc04_IAU2000.62-now``), do something like this example with an
        excerpt that is used for testing::

            >>> from astropy.utils.iers import IERS_B
            >>> from astropy.utils.data import get_pkg_data_filename
            >>> old_style_file = get_pkg_data_filename(
            ...     "tests/data/iers_b_old_style_excerpt",
            ...     package="astropy.utils.iers")
            >>> iers_b = IERS_B.read(
            ...     old_style_file,
            ...     readme=get_pkg_data_filename("data/ReadMe.eopc04_IAU2000",
            ...                                  package="astropy.utils.iers"),
            ...     data_start=14)

        """
        if file is None:
            file = IERS_B_FILE
        else:
            file = os.fspath(file)
        if readme is None:
            readme = IERS_B_README
        else:
            readme = os.fspath(readme)
        table = super().read(file, format="cds", readme=readme, data_start=data_start)

        table.meta["data_path"] = file
        table.meta["readme_path"] = readme
        return table

    def ut1_utc_source(self, i):
        """Set UT1-UTC source flag for entries in IERS table."""
        return np.ones_like(i) * FROM_IERS_B

    def dcip_source(self, i):
        """Set CIP correction source flag for entries in IERS table."""
        return np.ones_like(i) * FROM_IERS_B

    def pm_source(self, i):
        """Set PM source flag for entries in IERS table."""
        return np.ones_like(i) * FROM_IERS_B


class IERS_Auto(IERS_A):
    """
    Provide most-recent IERS data and automatically handle downloading
    of updated values as necessary.

    The returned table combines the IERS-A and IERS-B files, with the data
    in the IERS-B file considered to be official values and thus superseding
    values from the IERS-A file at the same times.
    """

    iers_table = None

    @classmethod
    def open(cls):
        """
        This reads in the bundled IERS_A table (or the version of the file in
        the current working directory, if present - see `IERS_A.read`).

        Subsequently, if values are requested which are more recent than the
        range of validity of the bundled table, an updated version of the table
        will be downloaded if the configuration setting
        ``astropy.utils.iers.conf.auto_download`` is set to True (default) and
        if the available version of this file is older than
        ``astropy.utils.iers.conf.auto_max_age`` days old.

        If the configuration setting ``astropy.utils.iers.conf.auto_download``
        is set to False then the bundled IERS-A table will always be used.

        Returns
        -------
        `~astropy.table.QTable` instance
            With IERS (Earth rotation) data columns

        """
        if cls.iers_table is None:
            cls.iers_table = cls.read()
        return cls.iers_table

    def _check_interpolate_indices(self, indices_orig, indices_clipped, max_input_mjd):
        """Check that the indices from interpolation match those after clipping to the
        valid table range.  The IERS_Auto class is exempted as long as it has
        sufficiently recent available data so the clipped interpolation is
        always within the confidence bounds of current Earth rotation
        knowledge.
        """
        predictive_mjd = self.meta["predictive_mjd"]

        # See explanation in _refresh_table_as_needed for these conditions
        auto_max_age = _none_to_float(conf.auto_max_age)
        if (
            max_input_mjd > predictive_mjd
            and Time.now().mjd - predictive_mjd > auto_max_age
        ):
            # If predictive values are stale, behave like IERS_B does for
            # out-of-range times: honor the iers_degraded_accuracy setting
            # instead of always raising.  Only the default "error" setting
            # raises; "warn" warns and "ignore" is silent.
            if conf.iers_degraded_accuracy == "error":
                raise ValueError(INTERPOLATE_ERROR.format(auto_max_age))
            elif conf.iers_degraded_accuracy == "warn":
                msg = (
                    f"interpolating from IERS_Auto using predictive values that "
                    f"are more than {auto_max_age} days old; accuracy is degraded."
                )
                warn(msg, IERSDegradedAccuracyWarning)

    def _refresh_table_as_needed(self, mjd):
        """Potentially update the IERS table in place depending on the requested
        time values in ``mjd`` and the time span of the table.

        For IERS_Auto the behavior is that the table is refreshed from the IERS
        server if both the following apply:

        - Any of the requested IERS values are predictive.  The IERS-A table
          contains predictive data out for a year after the available
          definitive values.
        - The first predictive values are at least ``conf.auto_max_age days`` old.
          In other words the IERS-A table was created by IERS long enough
          ago that it can be considered stale for predictions.
        """
        # If downloading is disabled, bail out silently.
        # _check_interpolate_indices() will error later if appropriate.
        if not conf.auto_download:
            return

        # Pass in initial to np.max to allow things to work for empty mjd.
        max_input_mjd = np.max(mjd, initial=50000)
        now_mjd = Time.now().mjd

        # IERS-A table contains predictive data out for a year after
        # the available definitive values.
        fpi = self.meta["predictive_index"]
        predictive_mjd = self.meta["predictive_mjd"]

        # Update table in place if necessary
        auto_max_age = _none_to_float(conf.auto_max_age)

        # If auto_max_age is smaller than IERS update time then repeated downloads may
        # occur without getting updated values (giving a IERSStaleWarning).
        if auto_max_age < 10:
            raise ValueError(
                "IERS auto_max_age configuration value must be larger than 10 days"
            )

        if max_input_mjd > predictive_mjd and (now_mjd - predictive_mjd) > auto_max_age:
            all_urls = (conf.iers_auto_url, conf.iers_auto_url_mirror)

            # Get the latest version
            try:
                filename = download_file(all_urls[0], sources=all_urls, cache="update")
            except Exception as err:
                # Issue a warning here, perhaps user is offline.  An exception
                # will be raised downstream when actually trying to interpolate
                # predictive values.
                warn(
                    AstropyWarning(
                        f"failed to download {' and '.join(all_urls)}: {err}.\nA"
                        " coordinate or time-related calculation might be compromised"
                        " or fail because the dates are not covered by the available"
                        ' IERS file.  See the "IERS data access" section of the'
                        " astropy documentation for additional information on working"
                        " offline."
                    )
                )
                return

            try:
                new_table = self.__class__.read(file=filename)
            except Exception as err:
                # The download succeeded but the content could not be parsed as
                # an IERS table (e.g. the server returned an error page).  Keep
                # using the existing table; an exception will be raised
                # downstream when actually trying to interpolate predictive
                # values.
                warn(
                    AstropyWarning(
                        "malformed IERS table downloaded from"
                        f" {' and '.join(all_urls)}: {err}.\nA coordinate or"
                        " time-related calculation might be compromised or fail"
                        " because the dates are not covered by the available IERS"
                        ' file.  See the "IERS data access" section of the astropy'
                        " documentation for additional information on working"
                        " offline."
                    )
                )
                return
            new_table.meta["data_url"] = str(all_urls[0])

            # New table has new values?
            if new_table["MJD"][-1] > self["MJD"][-1]:
                # Replace *replace* current values from the first predictive index through
                # the end of the current table.  This replacement is much faster than just
                # deleting all rows and then using add_row for the whole duration.
                new_fpi = np.searchsorted(
                    new_table["MJD"].value, predictive_mjd, side="right"
                )
                n_replace = len(self) - fpi
                self[fpi:] = new_table[new_fpi : new_fpi + n_replace]

                # Sanity check for continuity
                if new_table["MJD"][new_fpi + n_replace] - self["MJD"][-1] != 1.0 * u.d:
                    raise ValueError("unexpected gap in MJD when refreshing IERS table")

                # Now add new rows in place
                for row in new_table[new_fpi + n_replace :]:
                    self.add_row(row)

                self.meta.update(new_table.meta)
            else:
                warn(
                    IERSStaleWarning(
                        "IERS_Auto predictive values are older than"
                        f" {conf.auto_max_age} days but downloading the latest table"
                        " did not find newer values"
                    )
                )

    @classmethod
    def _substitute_iers_b(cls, table):
        """Substitute IERS B values with those from a real IERS B table.

        IERS-A has IERS-B values included, but for reasons unknown these
        do not match the latest IERS-B values (see comments in #4436).
        Here, we use the bundled astropy IERS-B table to overwrite the values
        in the IERS-A table.
        """
        iers_b = IERS_B.open()
        # Substitute IERS-B values for existing B values in IERS-A table
        mjd_b = table["MJD"][np.isfinite(table["UT1_UTC_B"])]
        i0 = np.searchsorted(iers_b["MJD"], mjd_b[0], side="left")
        i1 = np.searchsorted(iers_b["MJD"], mjd_b[-1], side="right")
        iers_b = iers_b[i0:i1]
        n_iers_b = len(iers_b)
        # If there is overlap then replace IERS-A values from available IERS-B
        if n_iers_b > 0:
            # Sanity check that we are overwriting the correct values
            if not u.allclose(table["MJD"][:n_iers_b], iers_b["MJD"]):
                raise ValueError(
                    "unexpected mismatch when copying IERS-B values into IERS-A table."
                )
            # Finally do the overwrite
            table["UT1_UTC_B"][:n_iers_b] = iers_b["UT1_UTC"]
            table["PM_X_B"][:n_iers_b] = iers_b["PM_x"]
            table["PM_Y_B"][:n_iers_b] = iers_b["PM_y"]
            table["dX_2000A_B"][:n_iers_b] = iers_b["dX_2000A"]
            table["dY_2000A_B"][:n_iers_b] = iers_b["dY_2000A"]

        return table


class earth_orientation_table(ScienceState):
    """Default IERS table for Earth rotation and reference systems service.

    These tables are used to calculate the offsets between ``UT1`` and ``UTC``
    and for conversion to Earth-based coordinate systems.

    The state itself is an IERS table, as an instance of one of the
    `~astropy.utils.iers.IERS` classes.  The default, the auto-updating
    `~astropy.utils.iers.IERS_Auto` class, should suffice for most
    purposes.

    Examples
    --------
    To temporarily use the IERS-B file packaged with astropy::

      >>> from astropy.utils import iers
      >>> from astropy.time import Time
      >>> iers_b = iers.IERS_B.open(iers.IERS_B_FILE)
      >>> with iers.earth_orientation_table.set(iers_b):
      ...     print(Time('2000-01-01').ut1.isot)
      2000-01-01T00:00:00.355

    To use the most recent IERS-A file for the whole session::

      >>> iers_a = iers.IERS_A.open(iers.IERS_A_URL)  # doctest: +SKIP
      >>> iers.earth_orientation_table.set(iers_a)  # doctest: +SKIP
      <ScienceState earth_orientation_table: <IERS_A length=17463>...>

    To go back to the default (of `~astropy.utils.iers.IERS_Auto`)::

      >>> iers.earth_orientation_table.set(None)  # doctest: +SKIP
      <ScienceState earth_orientation_table: <IERS_Auto length=17428>...>
    """

    _value = None

    @classmethod
    def validate(cls, value):
        if value is None:
            value = IERS_Auto.open()
        if not isinstance(value, IERS):
            raise ValueError("earth_orientation_table requires an IERS Table.")
        return value


class LeapSeconds(QTable):
    """Leap seconds class, holding TAI-UTC differences.

    The table holds columns 'year', 'month', 'day', 'mjd' and 'tai_utc', with
    one row for each date on which TAI-UTC changed.  Use
    `~astropy.utils.iers.LeapSeconds.from_erfa` to get the table built into
    ERFA, which is the table `~astropy.time.Time` uses.

    Notes
    -----
    Leap seconds have been discontinued: under Resolution 5 of the 28th CGPM
    (2026), UTC is continuous from 2027 May 20, and IERS Bulletin C 72 had
    already ruled out a leap second before 2027 June 28.  The last leap second
    was inserted on 2017 January 1 (TAI-UTC = 37 s), and the table built into
    ERFA already includes it, so it is final.

    Reading leap-second files (``Leap_Second.dat`` from IERS or
    ``leap-seconds.list`` from IETF) and loading tables into ERFA are
    deprecated since 8.1, since they can no longer give a different table.
    """

    # Note: Time instances in this class should use scale='tai' to avoid
    # needing leap seconds in their creation or interpretation.

    _re_expires = re.compile(r"^#.*File expires on[:\s]+(\d+\s\w+\s\d+)\s*$")
    _expires = None

    @classmethod
    def open(cls, file="erfa", cache=False):
        """Open a leap-second list.

        Parameters
        ----------
        file : str
            'erfa' (default), to get the table used by the ERFA library.
            Passing anything else, i.e., a local or network path to a file
            holding leap-second data, or `None`, is deprecated since 8.1.
        cache : bool
            Whether to use the download cache for a network path.  Defaults
            to False.

        Returns
        -------
        leap_seconds : `~astropy.utils.iers.LeapSeconds`
            Table with 'year', 'month', and 'tai_utc' columns, plus possibly
            others.
        """
        if isinstance(file, str) and file.lower() == "erfa":
            return cls.from_erfa()

        warn(
            "Calling LeapSeconds.open() with anything other than 'erfa' is "
            "deprecated since 8.1 and may be removed in a future version. Leap "
            "seconds have been discontinued, so the table built into ERFA is "
            "final; use LeapSeconds.from_erfa() instead.",
            AstropyDeprecationWarning,
            stacklevel=2,
        )
        return cls._open(file, cache=cache)

    @classmethod
    def _open(cls, file=None, cache=False):
        # Implementation of the deprecated part of open, to be removed with it.
        if file is None:
            return cls._auto_open()

        if isinstance(file, str) and file.lower() == "erfa":
            return cls.from_erfa()

        file = os.fspath(file)
        if urlparse(file).netloc:
            file = download_file(file, cache=cache)

        # Just try both reading methods.
        try:
            return cls._from_iers_leap_seconds(file)
        except Exception:
            return cls._from_leap_seconds_list(file)

    @classmethod
    @deprecated(
        since="8.1",
        message=(
            "LeapSeconds.auto_open is deprecated and may be removed in a future "
            "version. Leap seconds have been discontinued, so the table built "
            "into ERFA is final; use LeapSeconds.from_erfa() instead."
        ),
    )
    def auto_open(cls, files=None):
        """Get the leap-second list with the latest expiration date.

        Parameters
        ----------
        files : list of path-like, optional
            List of files/URLs to attempt to open.  By default, uses the
            leap-second file bundled with ``astropy-iers-data``.

        Returns
        -------
        leap_seconds : `~astropy.utils.iers.LeapSeconds`
            The table from the readable file with the latest expiration date.
        """
        return cls._auto_open(files)

    @classmethod
    def _auto_open(cls, files=None):
        # Implementation of the deprecated auto_open, to be removed with it.
        if files is None:
            files = [_IERS_LEAP_SECOND_FILE]

        self = None
        err_list = []
        for f in (os.fspath(f) for f in files if f):
            try:
                trial = cls._open(f, cache=True)
            except Exception as exc:
                err_list.append(exc)
                continue

            if self is None or trial.expires > self.expires:
                self = trial
                self.meta["data_url"] = f

        if self is None:
            raise ValueError(
                "none of the files could be read. The "
                f"following errors were raised:\n {err_list}"
            )

        return self

    @property
    def expires(self):
        """The limit of validity of the table."""
        return self._expires

    @classmethod
    def _read_leap_seconds(cls, file, **kwargs):
        """Read a file, identifying expiration by matching 'File expires'."""
        expires = None
        # Find expiration date.
        with get_readable_fileobj(file) as fh:
            lines = fh.readlines()
            for line in lines:
                match = cls._re_expires.match(line)
                if match:
                    day, month, year = match.groups()[0].split()
                    month_nb = MONTH_ABBR.index(month[:3]) + 1
                    expires = Time(
                        f"{year}-{month_nb:02d}-{day}", scale="tai", out_subfmt="date"
                    )
                    break
            else:
                raise ValueError(f"did not find expiration date in {file}")

        self = cls.read(lines, format="ascii.no_header", **kwargs)
        self._expires = expires
        return self

    @classmethod
    @deprecated(
        since="8.1",
        message=(
            "LeapSeconds.from_iers_leap_seconds is deprecated and may be removed "
            "in a future version. Leap seconds have been discontinued, so the "
            "table built into ERFA is final; use LeapSeconds.from_erfa() instead."
        ),
    )
    def from_iers_leap_seconds(cls, file=_IERS_LEAP_SECOND_FILE):
        """Create a table from a file like the IERS ``Leap_Second.dat``.

        Parameters
        ----------
        file : path-like, optional
            Full local or network path to the file holding leap-second data
            in a format consistent with that used by IERS.  By default, uses
            the file bundled with ``astropy-iers-data``.

        Notes
        -----
        The file *must* contain the expiration date in a comment line, like
        '#  File expires on 28 June 2020'
        """
        return cls._from_iers_leap_seconds(file)

    @classmethod
    def _from_iers_leap_seconds(cls, file=_IERS_LEAP_SECOND_FILE):
        # Implementation of the deprecated from_iers_leap_seconds.
        return cls._read_leap_seconds(
            file, names=["mjd", "day", "month", "year", "tai_utc"]
        )

    @classmethod
    @deprecated(
        since="8.1",
        message=(
            "LeapSeconds.from_leap_seconds_list is deprecated and may be removed "
            "in a future version. Leap seconds have been discontinued, so the "
            "table built into ERFA is final; use LeapSeconds.from_erfa() instead."
        ),
    )
    def from_leap_seconds_list(cls, file):
        """Create a table from a file like the IETF ``leap-seconds.list``.

        Parameters
        ----------
        file : path-like, optional
            Full local or network path to the file holding leap-second data
            in a format consistent with that used by IETF.

        Notes
        -----
        The file *must* contain the expiration date in a comment line, like
        '# File expires on:  28 June 2020'
        """
        return cls._from_leap_seconds_list(file)

    @classmethod
    def _from_leap_seconds_list(cls, file):
        # Implementation of the deprecated from_leap_seconds_list.
        from astropy.io.ascii import convert_numpy  # Here to avoid circular import

        names = ["ntp_seconds", "tai_utc", "comment", "day", "month", "year"]
        # Note: ntp_seconds does not fit in 32 bit, so causes problems on
        # 32-bit systems without the np.int64 converter.
        self = cls._read_leap_seconds(
            file,
            names=names,
            include_names=names[:2],
            converters={"ntp_seconds": [convert_numpy(np.int64)]},
        )
        self["mjd"] = (self["ntp_seconds"] / 86400 + 15020).round()
        # Note: cannot use Time.ymdhms, since that might require leap seconds.
        isot = Time(self["mjd"], format="mjd", scale="tai").isot
        ymd = np.array(
            [[int(part) for part in t.partition("T")[0].split("-")] for t in isot]
        )
        self["year"], self["month"], self["day"] = ymd.T
        return self

    @classmethod
    def from_erfa(cls, built_in=False):
        """Create table from the leap-second list in ERFA.

        Parameters
        ----------
        built_in : bool
            If `False` (default), retrieve the list currently used by ERFA,
            which may have been updated.  If `True`, retrieve the list shipped
            with erfa.

        Notes
        -----
        The 'day' column is always 1 and 'mjd' is the corresponding Modified
        Julian Date, as in tables read from leap-second files.  Unlike those,
        the ERFA table also has entries from before 1972, when UTC was not yet
        an integer number of seconds from TAI.  For those entries, 'tai_utc'
        is only the offset at the given date, since ERFA adds a drift with
        time that is not part of this table.
        """
        current = cls(erfa.leap_seconds.get())
        current["day"] = np.ones_like(current["month"])
        current["mjd"] = erfa.cal2jd(current["year"], current["month"], 1)[1].astype(
            int
        )
        # Leap seconds were discontinued (CGPM 2026, Resolution 5), so the leap-second
        # table built into ERFA is final and cannot expire.  pyerfa versions that
        # include liberfa/pyerfa#359 report this themselves; for older ones, which
        # consider the table expired 180 days after the last leap second, set the same
        # expiration as pyerfa#359 (180 days ahead) so erfa.leap_seconds.expired is
        # False either way.
        if erfa.leap_seconds.expired:
            # Setting this private module attribute is equivalent to the previous
            # auto-update machinery. This code can be removed when #359 is required
            # by astropy.
            erfa.leap_seconds._expires = datetime.now().replace(
                hour=0, minute=0, second=0, microsecond=0
            ) + timedelta(180)

        expires = erfa.leap_seconds.expires
        current._expires = Time(
            f"{expires.year:04d}-{expires.month:02d}-{expires.day:02d}",
            scale="tai",
        )

        if not built_in:
            return current

        try:
            erfa.leap_seconds.set(None)  # reset to defaults
            return cls.from_erfa(built_in=False)
        finally:
            erfa.leap_seconds.set(current)

    @deprecated(
        since="8.1",
        message=(
            "LeapSeconds.update_erfa_leap_seconds is deprecated and may be removed "
            "in a future version. Leap seconds have been discontinued, so the "
            "table built into ERFA is final; if you need to change it anyway, "
            "use erfa.leap_seconds.set() or erfa.leap_seconds.update()."
        ),
    )
    def update_erfa_leap_seconds(self, initialize_erfa=False):
        """Add any leap seconds not already present to the ERFA table.

        This method matches leap seconds with those present in the ERFA table,
        and extends the latter as necessary.

        Parameters
        ----------
        initialize_erfa : bool, or 'only', or 'empty'
            Initialize the ERFA leap second table to its built-in value before
            trying to expand it.  This is generally not needed but can help
            in case it somehow got corrupted.  If equal to 'only', the ERFA
            table is reinitialized and no attempt it made to update it.
            If 'empty', the leap second table is emptied before updating, i.e.,
            it is overwritten altogether (note that this may break things in
            surprising ways, as most leap second tables do not include pre-1970
            pseudo leap-seconds; you were warned).

        Returns
        -------
        n_update : int
            Number of items updated.

        Raises
        ------
        ValueError
            If the leap seconds in the table are not on 1st of January or July,
            or if the matches are inconsistent.  This would normally suggest
            a corrupted leap second table, but might also indicate that the
            ERFA table was corrupted.  If needed, the ERFA table can be reset
            by calling this method with an appropriate value for
            ``initialize_erfa``.
        """
        if initialize_erfa == "empty":
            # Initialize to empty and update is the same as overwrite.
            erfa.leap_seconds.set(self)
            return len(self)

        if initialize_erfa:
            erfa.leap_seconds.set()
            if initialize_erfa == "only":
                return 0

        return erfa.leap_seconds.update(self)
