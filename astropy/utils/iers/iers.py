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

from __future__ import annotations

import os
import re
from datetime import UTC, datetime
from typing import TYPE_CHECKING
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
    IERS_LEAP_SECOND_FILE,
    IERS_LEAP_SECOND_URL,
)
from astropy_iers_data import IERS_LEAP_SECOND_URL_MIRROR as IETF_LEAP_SECOND_URL

from astropy import config as _config
from astropy import units as u
from astropy import utils
from astropy.table import MaskedColumn, QTable
from astropy.time import Time, TimeDelta
from astropy.utils.data import (
    clear_download_cache,
    get_readable_fileobj,
    is_url_in_cache,
)
from astropy.utils.exceptions import AstropyDeprecationWarning, AstropyWarning
from astropy.utils.state import ScienceState

if TYPE_CHECKING:
    from typing import Self

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
    "IERS_LEAP_SECOND_FILE",
    "IERS_LEAP_SECOND_URL",
    "IETF_LEAP_SECOND_URL",
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


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.utils.iers`.
    """

    auto_download = _config.ConfigItem(
        True,
        "Enable auto-downloading of the latest IERS-A data.  If set to False "
        "then the bundled IERS-A file will be used by default (even if a "
        "newer version of the IERS-A file was previously downloaded and cached). "
        "This parameter also controls whether internet resources will be "
        "queried to update the leap second table if the installed version is "
        "out of date. Default is True.",
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
        "to degraded accuracy.  Applies only when using IERS-B data on its own.",
    )
    system_leap_second_file = _config.ConfigItem("", "System file with leap seconds.")
    iers_leap_second_auto_url = _config.ConfigItem(
        IERS_LEAP_SECOND_URL, "URL for auto-downloading leap seconds."
    )
    ietf_leap_second_auto_url = _config.ConfigItem(
        IETF_LEAP_SECOND_URL, "Alternate URL for auto-downloading leap seconds."
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
        """If the configuration setting ``astropy.utils.iers.conf.auto_download``
        is set to True (default), then open a recent version of the IERS-A
        table with predictions for UT1-UTC and polar motion out to
        approximately one year from now.  If the available version of this file
        is older than ``astropy.utils.iers.conf.auto_max_age`` days old
        (or non-existent) then it will be downloaded over the network and cached.

        If the configuration setting ``astropy.utils.iers.conf.auto_download``
        is set to False then the bundled IERS-A table will be used rather than
        any downloaded version of the IERS-A table.

        On the first call in a session, the table will be memoized (in the
        ``iers_table`` class attribute), and further calls to ``open`` will
        return this stored table.

        Returns
        -------
        `~astropy.table.QTable` instance
            With IERS (Earth rotation) data columns

        """
        if not conf.auto_download:
            # If auto_download is changed to False mid-session, iers_table may have already been
            # made from non-bundled files, so it should be remade from bundled files
            if not hasattr(cls, "_iers_table_bundled"):
                cls._iers_table_bundled = cls.read()
            cls.iers_table = cls._iers_table_bundled
            return cls.iers_table

        all_urls = (conf.iers_auto_url, conf.iers_auto_url_mirror)

        if cls.iers_table is not None:
            # If the URL has changed, we need to redownload the file, so we
            # should ignore the internally cached version.

            if cls.iers_table.meta.get("data_url") in all_urls:
                return cls.iers_table

        for url in all_urls:
            try:
                filename = download_file(url, cache=True)
            except Exception as err:
                warn(f"failed to download {url}: {err}", IERSWarning)
                continue

            try:
                cls.iers_table = cls.read(file=filename)
            except Exception as err:
                warn(f"malformed IERS table from {url}: {err}", IERSWarning)
                continue
            cls.iers_table.meta["data_url"] = url
            break

        else:
            # Issue a warning here, perhaps user is offline.  An exception
            # will be raised downstream if actually trying to interpolate
            # predictive values.
            warn(
                "unable to download valid IERS file, using bundled IERS-A", IERSWarning
            )
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
            and self.time_now.mjd - predictive_mjd > auto_max_age
        ):
            raise ValueError(INTERPOLATE_ERROR.format(auto_max_age))

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
        now_mjd = self.time_now.mjd

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

            new_table = self.__class__.read(file=filename)
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

    The table should hold columns 'year', 'month', 'tai_utc'.

    Methods are provided to initialize the table from IERS ``Leap_Second.dat``,
    IETF/ntp ``leap-seconds.list``, or built-in ERFA/SOFA, and to update the
    list used by ERFA.

    Notes
    -----
    Astropy has a built-in ``iers.IERS_LEAP_SECONDS_FILE``. Up to date versions
    can be downloaded from ``iers.IERS_LEAP_SECONDS_URL`` or
    ``iers.LEAP_SECONDS_LIST_URL``.  Many systems also store a version
    of ``leap-seconds.list`` for use with ``ntp`` (e.g., on Debian/Ubuntu
    systems, ``/usr/share/zoneinfo/leap-seconds.list``).

    To prevent querying internet resources if the available local leap second
    file(s) are out of date, set ``iers.conf.auto_download = False``. This
    must be done prior to performing any ``Time`` scale transformations related
    to UTC (e.g. converting from UTC to TAI).
    """

    # Note: Time instances in this class should use scale='tai' to avoid
    # needing leap seconds in their creation or interpretation.

    _re_expires = re.compile(r"^#.*File expires on[:\s]+(\d+\s\w+\s\d+)\s*$")
    _expires = None
    _auto_open_files = [
        "erfa",
        IERS_LEAP_SECOND_FILE,
        "system_leap_second_file",
        "iers_leap_second_auto_url",
        "ietf_leap_second_auto_url",
    ]
    """Files or conf attributes to try in auto_open."""

    @classmethod
    def open(cls, file=None, cache=False):
        """Open a leap-second list.

        Parameters
        ----------
        file : path-like or None
            Full local or network path to the file holding leap-second data,
            for passing on to the various ``from_`` class methods.
            If 'erfa', return the data used by the ERFA library.
            If `None`, use default locations from file and configuration to
            find a table that is not expired.
        cache : bool
            Whether to use cache. Defaults to False, since leap-second files
            are regularly updated.

        Returns
        -------
        leap_seconds : `~astropy.utils.iers.LeapSeconds`
            Table with 'year', 'month', and 'tai_utc' columns, plus possibly
            others.

        Notes
        -----
        Bulletin C is released about 10 days after a possible leap second is
        introduced, i.e., mid-January or mid-July.  Expiration days are thus
        generally at least 150 days after the present.  For the auto-loading,
        a list comprised of the table shipped with astropy, and files and
        URLs in `~astropy.utils.iers.Conf` are tried, returning the first
        that is sufficiently new, or the newest among them all.
        """
        if file is None:
            return cls.auto_open()

        if file.lower() == "erfa":
            return cls.from_erfa()

        if urlparse(file).netloc:
            file = download_file(file, cache=cache)

        # Just try both reading methods.
        try:
            return cls.from_iers_leap_seconds(file)
        except Exception:
            return cls.from_leap_seconds_list(file)

    @staticmethod
    def _today():
        # Get current day in scale='tai' without going through a scale change
        # (so we do not need leap seconds).
        s = "{0.year:04d}-{0.month:02d}-{0.day:02d}".format(datetime.now(tz=UTC))
        return Time(s, scale="tai", format="iso", out_subfmt="date")

    @classmethod
    def auto_open(cls, files=None):
        """Attempt to get an up-to-date leap-second list.

        The routine will try the files in sequence until it finds one
        whose expiration date is "good enough" (see below).  If none
        are good enough, it returns the one with the most recent expiration
        date, warning if that file is expired.

        For remote files that are cached already, the cached file is tried
        first before attempting to retrieve it again.

        Parameters
        ----------
        files : list of path-like, optional
            List of files/URLs to attempt to open.  By default, uses
            ``cls._auto_open_files``.

        Returns
        -------
        leap_seconds : `~astropy.utils.iers.LeapSeconds`
            Up to date leap-second table

        Notes
        -----
        Bulletin C is released about 10 days after a possible leap second is
        introduced, i.e., mid-January or mid-July.  Expiration days are thus
        generally at least 150 days after the present.  We look for a file
        that expires more than 180 - `~astropy.utils.iers.Conf.auto_max_age`
        after the present.
        """
        offset = 180 - (30 if conf.auto_max_age is None else conf.auto_max_age)
        good_enough = cls._today() + TimeDelta(offset, format="jd")

        if files is None:
            # Basic files to go over (entries in _auto_open_files can be
            # configuration items, which we want to be sure are up to date).
            files = [getattr(conf, f, f) for f in cls._auto_open_files]

        # Remove empty entries and normalize Path objects to string
        files = [os.fspath(f) for f in files if f]

        # Our trials start with normal files and remote ones that are
        # already in cache.  The bools here indicate that the cache
        # should be used.
        trials = [
            (f, True) for f in files if not urlparse(f).netloc or is_url_in_cache(f)
        ]
        # If we are allowed to download, we try downloading new versions
        # if none of the above worked.
        if conf.auto_download:
            trials += [(f, False) for f in files if urlparse(f).netloc]

        self = None
        err_list = []
        # Go through all entries, and return the first one that
        # is not expired, or the most up to date one.
        for f, allow_cache in trials:
            if not allow_cache:
                clear_download_cache(f)

            try:
                trial = cls.open(f, cache=True)
            except Exception as exc:
                err_list.append(exc)
                continue

            if self is None or trial.expires > self.expires:
                self = trial
                self.meta["data_url"] = str(f)
                if self.expires > good_enough:
                    break

        if self is None:
            raise ValueError(
                "none of the files could be read. The "
                f"following errors were raised:\n {err_list}"
            )

        if self.expires < self._today() and conf.auto_max_age is not None:
            warn("leap-second file is expired.", IERSStaleWarning)

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
    def from_iers_leap_seconds(cls, file=IERS_LEAP_SECOND_FILE):
        """Create a table from a file like the IERS ``Leap_Second.dat``.

        Parameters
        ----------
        file : path-like, optional
            Full local or network path to the file holding leap-second data
            in a format consistent with that used by IERS.  By default, uses
            ``iers.IERS_LEAP_SECOND_FILE``.

        Notes
        -----
        The file *must* contain the expiration date in a comment line, like
        '#  File expires on 28 June 2020'
        """
        return cls._read_leap_seconds(
            file, names=["mjd", "day", "month", "year", "tai_utc"]
        )

    @classmethod
    def from_leap_seconds_list(cls, file):
        """Create a table from a file like the IETF ``leap-seconds.list``.

        Parameters
        ----------
        file : path-like, optional
            Full local or network path to the file holding leap-second data
            in a format consistent with that used by IETF.  Up to date versions
            can be retrieved from ``iers.IETF_LEAP_SECOND_URL``.

        Notes
        -----
        The file *must* contain the expiration date in a comment line, like
        '# File expires on:  28 June 2020'
        """
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
        """
        current = cls(erfa.leap_seconds.get())
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
