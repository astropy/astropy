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
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from warnings import warn

try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse

import numpy as np

from ... import config as _config
from ... import units as u
from ...table import Table, QTable
from ...utils.data import get_pkg_data_filename, clear_download_cache
from ... import utils
from ...utils.exceptions import AstropyWarning

__all__ = ['Conf', 'conf',
           'IERS', 'IERS_B', 'IERS_A', 'IERS_Auto',
           'FROM_IERS_B', 'FROM_IERS_A', 'FROM_IERS_A_PREDICTION',
           'TIME_BEFORE_IERS_RANGE', 'TIME_BEYOND_IERS_RANGE',
           'IERS_A_FILE', 'IERS_A_URL', 'IERS_A_README',
           'IERS_B_FILE', 'IERS_B_URL', 'IERS_B_README',
           'IERSRangeError', 'IERSStaleWarning']

# IERS-A default file name, URL, and ReadMe with content description
IERS_A_FILE = 'finals2000A.all'
IERS_A_URL = 'http://maia.usno.navy.mil/ser7/finals2000A.all'
IERS_A_README = get_pkg_data_filename('data/ReadMe.finals2000A')

# IERS-B default file name, URL, and ReadMe with content description
IERS_B_FILE = get_pkg_data_filename('data/eopc04_IAU2000.62-now')
IERS_B_URL = 'http://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now'
IERS_B_README = get_pkg_data_filename('data/ReadMe.eopc04_IAU2000')

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


def download_file(*args, **kwargs):
    """
    Overload astropy.utils.data.download_file within iers module to use a
    custom (longer) wait time.  This just passes through ``*args`` and
    ``**kwargs`` after temporarily setting the download_file remote timeout to
    the local ``iers.conf.remote_timeout`` value.
    """
    with utils.data.conf.set_temp('remote_timeout', conf.remote_timeout):
        return utils.data.download_file(*args, **kwargs)


class IERSStaleWarning(AstropyWarning):
    pass


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.utils.iers`.
    """
    auto_download = _config.ConfigItem(
        True,
        'Enable auto-downloading of the latest IERS data.  If set to False '
        'then the local IERS-B file will be used by default. Default is True.')
    auto_max_age = _config.ConfigItem(
        30.0,
        'Maximum age (days) of predictive data before auto-downloading. Default is 30.')
    iers_auto_url = _config.ConfigItem(
        IERS_A_URL,
        'URL for auto-downloading IERS file data.')
    remote_timeout = _config.ConfigItem(
        10.0,
        'Remote timeout downloading IERS file data (seconds).')


conf = Conf()


class IERSRangeError(IndexError):
    """
    Any error for when dates are outside of the valid range for IERS
    """


class IERS(QTable):
    """Generic IERS table class, defining interpolation functions.

    Sub-classed from `astropy.table.QTable`.  The table should hold columns
    'MJD', 'UT1_UTC', and 'PM_x'/'PM_y'.
    """

    iers_table = None

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
            cls.iers_table = cls.read(**kwargs)
        return cls.iers_table

    @classmethod
    def close(cls):
        """Remove the IERS table from the class.

        This allows the table to be re-read from disk during one's session
        (e.g., if one finds it is out of date and has updated the file).
        """
        cls.iers_table = None

    def mjd_utc(self, jd1, jd2=0.):
        """Turn a time to MJD, returning integer and fractional parts.

        Parameters
        ----------
        jd1 : float, array, or Time
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
        utc = jd1 - (MJD_ZERO+mjd) + jd2
        return mjd, utc

    def ut1_utc(self, jd1, jd2=0., return_status=False):
        """Interpolate UT1-UTC corrections in IERS Table for given dates.

        Parameters
        ----------
        jd1 : float, float array, or Time object
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
        return self._interpolate(jd1, jd2, ['UT1_UTC'],
                                 self.ut1_utc_source if return_status else None)

    def pm_xy(self, jd1, jd2=0., return_status=False):
        """Interpolate polar motions from IERS Table for given dates.

        Parameters
        ----------
        jd1 : float, float array, or Time object
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
        PM_x : Quantity with angle units
            x component of polar motion for the requested times
        PM_y : Quantity with angle units
            y component of polar motion for the requested times
        status : int or int array
            Status values (if ``return_status``=``True``)::
            ``iers.FROM_IERS_B``
            ``iers.FROM_IERS_A``
            ``iers.FROM_IERS_A_PREDICTION``
            ``iers.TIME_BEFORE_IERS_RANGE``
            ``iers.TIME_BEYOND_IERS_RANGE``
        """
        return self._interpolate(jd1, jd2, ['PM_x', 'PM_y'],
                                 self.pm_source if return_status else None)

    def _check_interpolate_indices(self, indices_orig, indices_clipped, max_input_mjd):
        """
        Check that the indices from interpolation match those after clipping
        to the valid table range.  This method gets overridden in the IERS_Auto
        class because it has different requirements.
        """
        if np.any(indices_orig != indices_clipped):
            raise IERSRangeError('(some) times are outside of range covered '
                                 'by IERS table.')

    def _interpolate(self, jd1, jd2, columns, source=None):
        mjd, utc = self.mjd_utc(jd1, jd2)
        # enforce array
        is_scalar = not hasattr(mjd, '__array__') or mjd.ndim == 0
        if is_scalar:
            mjd = np.array([mjd])
            utc = np.array([utc])

        self._refresh_table_as_needed(mjd)

        # For typical format, will always find a match (since MJD are integer)
        # hence, important to define which side we will be; this ensures
        # self['MJD'][i-1]<=mjd<self['MJD'][i]
        i = np.searchsorted(self['MJD'].value, mjd, side='right')

        # Get index to MJD at or just below given mjd, clipping to ensure we
        # stay in range of table (status will be set below for those outside)
        i1 = np.clip(i, 1, len(self) - 1)
        i0 = i1 - 1
        mjd_0, mjd_1 = self['MJD'][i0].value, self['MJD'][i1].value
        results = []
        for column in columns:
            val_0, val_1 = self[column][i0], self[column][i1]
            d_val = val_1 - val_0
            if column == 'UT1_UTC':
                # Check & correct for possible leap second (correcting diff.,
                # not 1st point, since jump can only happen right at 2nd point)
                d_val -= d_val.round()
            # Linearly interpolate (which is what TEMPO does for UT1-UTC, but
            # may want to follow IERS gazette #13 for more precise
            # interpolation and correction for tidal effects;
            # http://maia.usno.navy.mil/iers-gaz13)
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
            self._check_interpolate_indices(i1, i, np.max(mjd))
            return results[0] if len(results) == 1 else results

    def _refresh_table_as_needed(self, mjd):
        """
        Potentially update the IERS table in place depending on the requested
        time values in ``mdj`` and the time span of the table.  The base behavior
        is not to update the table.  ``IERS_Auto`` overrides this method.
        """
        pass

    def ut1_utc_source(self, i):
        """Source for UT1-UTC.  To be overridden by subclass."""
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
        from astropy.time import Time
        try:
            return self._time_now
        except Exception:
            return Time.now()


class IERS_A(IERS):
    """IERS Table class targeted to IERS A, provided by USNO.

    These include rapid turnaround and predicted times.
    See http://maia.usno.navy.mil/

    Notes
    -----
    The IERS A file is not part of astropy.  It can be downloaded from
    ``iers.IERS_A_URL``.  See ``iers.__doc__`` for instructions on how to use
    it in ``Time``, etc.
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
        table = iers_a[~iers_a['UT1_UTC_A'].mask &
                       ~iers_a['PolPMFlag_A'].mask]

        # This does nothing for IERS_A, but allows IERS_Auto to ensure the
        # IERS B values in the table are consistent with the true ones.
        table = cls._substitute_iers_b(table)

        # Run np.where on the data from the table columns, since in numpy 1.9
        # it otherwise returns an only partially initialized column.
        table['UT1_UTC'] = np.where(table['UT1_UTC_B'].mask,
                                    table['UT1_UTC_A'].data,
                                    table['UT1_UTC_B'].data)
        # Ensure the unit is correct, for later column conversion to Quantity.
        table['UT1_UTC'].unit = table['UT1_UTC_A'].unit
        table['UT1Flag'] = np.where(table['UT1_UTC_B'].mask,
                                    table['UT1Flag_A'].data,
                                    'B')
        # Repeat for polar motions.
        table['PM_x'] = np.where(table['PM_X_B'].mask,
                                 table['PM_x_A'].data,
                                 table['PM_X_B'].data)
        table['PM_x'].unit = table['PM_x_A'].unit
        table['PM_y'] = np.where(table['PM_Y_B'].mask,
                                 table['PM_y_A'].data,
                                 table['PM_Y_B'].data)
        table['PM_y'].unit = table['PM_y_A'].unit
        table['PolPMFlag'] = np.where(table['PM_X_B'].mask,
                                      table['PolPMFlag_A'].data,
                                      'B')

        # Get the table index for the first row that has predictive values
        # PolPMFlag_A  IERS (I) or Prediction (P) flag for
        #              Bull. A polar motion values
        # UT1Flag_A    IERS (I) or Prediction (P) flag for
        #              Bull. A UT1-UTC values
        is_predictive = (table['UT1Flag_A'] == 'P') | (table['PolPMFlag_A'] == 'P')
        table.meta['predictive_index'] = np.min(np.flatnonzero(is_predictive))
        table.meta['predictive_mjd'] = table['MJD'][table.meta['predictive_index']]

        return table

    @classmethod
    def _substitute_iers_b(cls, table):
        # See documentation in IERS_Auto.
        return table

    @classmethod
    def read(cls, file=None, readme=None):
        """Read IERS-A table from a finals2000a.* file provided by USNO.

        Parameters
        ----------
        file : str
            full path to ascii file holding IERS-A data.
            Defaults to ``iers.IERS_A_FILE``.
        readme : str
            full path to ascii file holding CDS-style readme.
            Defaults to package version, ``iers.IERS_A_README``.

        Returns
        -------
        ``IERS_A`` class instance
        """
        if file is None:
            file = IERS_A_FILE
        if readme is None:
            readme = IERS_A_README

        # Read in as a regular Table, including possible masked columns.
        # Columns will be filled and converted to Quantity in cls.__init__.
        iers_a = Table.read(file, format='cds', readme=readme)

        # Combine the A and B data for UT1-UTC and PM columns
        table = cls._combine_a_b_columns(iers_a)
        table.meta['data_path'] = file
        table.meta['readme_path'] = readme

        # Fill any masked values, and convert to a QTable.
        return cls(table.filled())

    def ut1_utc_source(self, i):
        """Set UT1-UTC source flag for entries in IERS table"""
        ut1flag = self['UT1Flag'][i]
        source = np.ones_like(i) * FROM_IERS_B
        source[ut1flag == 'I'] = FROM_IERS_A
        source[ut1flag == 'P'] = FROM_IERS_A_PREDICTION
        return source

    def pm_source(self, i):
        """Set polar motion source flag for entries in IERS table"""
        pmflag = self['PolPMFlag'][i]
        source = np.ones_like(i) * FROM_IERS_B
        source[pmflag == 'I'] = FROM_IERS_A
        source[pmflag == 'P'] = FROM_IERS_A_PREDICTION
        return source


class IERS_B(IERS):
    """IERS Table class targeted to IERS B, provided by IERS itself.

    These are final values; see http://www.iers.org/

    Notes
    -----
    If the package IERS B file (```iers.IERS_B_FILE``) is out of date, a new
    version can be downloaded from ``iers.IERS_B_URL``.
    """

    iers_table = None

    @classmethod
    def read(cls, file=None, readme=None, data_start=14):
        """Read IERS-B table from a eopc04_iau2000.* file provided by IERS.

        Parameters
        ----------
        file : str
            full path to ascii file holding IERS-B data.
            Defaults to package version, ``iers.IERS_B_FILE``.
        readme : str
            full path to ascii file holding CDS-style readme.
            Defaults to package version, ``iers.IERS_B_README``.
        data_start : int
            starting row. Default is 14, appropriate for standard IERS files.

        Returns
        -------
        ``IERS_B`` class instance
        """
        if file is None:
            file = IERS_B_FILE
        if readme is None:
            readme = IERS_B_README

        # Read in as a regular Table, including possible masked columns.
        # Columns will be filled and converted to Quantity in cls.__init__.
        iers_b = Table.read(file, format='cds', readme=readme,
                            data_start=data_start)
        return cls(iers_b.filled())

    def ut1_utc_source(self, i):
        """Set UT1-UTC source flag for entries in IERS table"""
        return np.ones_like(i) * FROM_IERS_B

    def pm_source(self, i):
        """Set PM source flag for entries in IERS table"""
        return np.ones_like(i) * FROM_IERS_B


class IERS_Auto(IERS_A):
    """
    Provide most-recent IERS data and automatically handle downloading
    of updated values as necessary.
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
        is set to False then ``astropy.utils.iers.IERS()`` is returned.  This
        is normally the IERS-B table that is supplied with astropy.

        On the first call in a session, the table will be memoized (in the
        ``iers_table`` class attribute), and further calls to ``open`` will
        return this stored table.

        Returns
        -------
        `~astropy.table.QTable` instance with IERS (Earth rotation) data columns

        """
        if not conf.auto_download:
            cls.iers_table = IERS.open()
            return cls.iers_table

        if cls.iers_table is not None:

            # If the URL has changed, we need to redownload the file, so we
            # should ignore the internally cached version.

            if cls.iers_table.meta.get('data_url') == conf.iers_auto_url:
                return cls.iers_table

        try:
            filename = download_file(conf.iers_auto_url, cache=True)
        except Exception as err:
            # Issue a warning here, perhaps user is offline.  An exception
            # will be raised downstream when actually trying to interpolate
            # predictive values.
            warn(AstropyWarning('failed to download {}, using local IERS-B: {}'
                                .format(conf.iers_auto_url, str(err))))
            cls.iers_table = IERS.open()
            return cls.iers_table

        cls.iers_table = cls.read(file=filename)
        cls.iers_table.meta['data_url'] = str(conf.iers_auto_url)

        return cls.iers_table

    def _check_interpolate_indices(self, indices_orig, indices_clipped, max_input_mjd):
        """Check that the indices from interpolation match those after clipping to the
        valid table range.  The IERS_Auto class is exempted as long as it has
        sufficiently recent available data so the clipped interpolation is
        always within the confidence bounds of current Earth rotation
        knowledge.
        """
        predictive_mjd = self.meta['predictive_mjd']

        # See explanation in _refresh_table_as_needed for these conditions
        auto_max_age = (conf.auto_max_age if conf.auto_max_age is not None
                        else np.finfo(np.float).max)
        if (max_input_mjd > predictive_mjd and
                self.time_now.mjd - predictive_mjd > auto_max_age):
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
        max_input_mjd = np.max(mjd)
        now_mjd = self.time_now.mjd

        # IERS-A table contains predictive data out for a year after
        # the available definitive values.
        fpi = self.meta['predictive_index']
        predictive_mjd = self.meta['predictive_mjd']

        # Update table in place if necessary
        auto_max_age = (conf.auto_max_age if conf.auto_max_age is not None
                        else np.finfo(np.float).max)

        # If auto_max_age is smaller than IERS update time then repeated downloads may
        # occur without getting updated values (giving a IERSStaleWarning).
        if auto_max_age < 10:
            raise ValueError('IERS auto_max_age configuration value must be larger than 10 days')

        if (max_input_mjd > predictive_mjd and
               now_mjd - predictive_mjd > auto_max_age):

            # Get the latest version
            try:
                clear_download_cache(conf.iers_auto_url)
                filename = download_file(conf.iers_auto_url, cache=True)
            except Exception as err:
                # Issue a warning here, perhaps user is offline.  An exception
                # will be raised downstream when actually trying to interpolate
                # predictive values.
                warn(AstropyWarning('failed to download {}: {}.\nA coordinate or time-related '
                                    'calculation might be compromised or fail because the dates are '
                                    'not covered by the available IERS file.  See the '
                                    '"IERS data access" section of the astropy documentation '
                                    'for additional information on working offline.'
                                    .format(conf.iers_auto_url, str(err))))
                return

            new_table = self.__class__.read(file=filename)

            # New table has new values?
            if new_table['MJD'][-1] > self['MJD'][-1]:
                # Replace *replace* current values from the first predictive index through
                # the end of the current table.  This replacement is much faster than just
                # deleting all rows and then using add_row for the whole duration.
                new_fpi = np.searchsorted(new_table['MJD'].value, predictive_mjd, side='right')
                n_replace = len(self) - fpi
                self[fpi:] = new_table[new_fpi:new_fpi + n_replace]

                # Sanity check for continuity
                if new_table['MJD'][new_fpi + n_replace] - self['MJD'][-1] != 1.0 * u.d:
                    raise ValueError('unexpected gap in MJD when refreshing IERS table')

                # Now add new rows in place
                for row in new_table[new_fpi + n_replace:]:
                    self.add_row(row)

                self.meta.update(new_table.meta)
            else:
                warn(IERSStaleWarning(
                    'IERS_Auto predictive values are older than {} days but downloading '
                    'the latest table did not find newer values'.format(conf.auto_max_age)))

    @classmethod
    def _substitute_iers_b(cls, table):
        """Substitute IERS B values with those from a real IERS B table.

        IERS-A has IERS-B values included, but for reasons unknown these
        do not match the latest IERS-B values (see comments in #4436).
        Here, we use the bundled astropy IERS-B table to overwrite the values
        in the downloaded IERS-A table.
        """
        iers_b = IERS_B.open()
        # Substitute IERS-B values for existing B values in IERS-A table
        mjd_b = table['MJD'][~table['UT1_UTC_B'].mask]
        i0 = np.searchsorted(iers_b['MJD'].value, mjd_b[0], side='left')
        i1 = np.searchsorted(iers_b['MJD'].value, mjd_b[-1], side='right')
        iers_b = iers_b[i0:i1]
        n_iers_b = len(iers_b)
        # If there is overlap then replace IERS-A values from available IERS-B
        if n_iers_b > 0:
            # Sanity check that we are overwriting the correct values
            if not np.allclose(table['MJD'][:n_iers_b], iers_b['MJD'].value):
                raise ValueError('unexpected mismatch when copying '
                                 'IERS-B values into IERS-A table.')
            # Finally do the overwrite
            table['UT1_UTC_B'][:n_iers_b] = iers_b['UT1_UTC'].value
            table['PM_X_B'][:n_iers_b] = iers_b['PM_x'].value
            table['PM_Y_B'][:n_iers_b] = iers_b['PM_y'].value

        return table


# by default for IERS class, read IERS-B table
IERS.read = IERS_B.read
