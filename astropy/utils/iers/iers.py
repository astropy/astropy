# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The astropy.utils.iers package provides access to the tables provided by
the International Earth Rotation and Reference Systems Service, in
particular allowing interpolation of published UT1-UTC values for given
times.  These are used in `astropy.time` to provide UT1 values.  The polar
motions are also used for determining earth orientation for
celestional-to-terrestrial coordinate transformations
(in `astropy.coordinates`).

By default, IERS B values provided as part of astropy are used, but
user-downloaded files can be substituted.

Generally, there is no need to invoke the iers classes oneself.  E.g.,
the default IERS B table is loaded as necessary in `Time`::

    >>> from astropy.time import Time
    >>> t = Time(['2012-06-30 12:00:00', '2012-06-30 23:59:59',
    ...           '2012-06-30 23:59:60', '2012-07-01 00:00:00',
    ...           '2012-07-01 12:00:00'], scale='utc')
    >>> t.ut1
    <Time object: scale='ut1' format='iso' value=['2012-06-30 11:59:59.413'
     '2012-06-30 23:59:58.413' '2012-06-30 23:59:59.413'
     '2012-07-01 00:00:00.413' '2012-07-01 12:00:00.413']>

But if one is dealing with very recent observations, this does not work::

    >>> t2 = Time.now()
    >>> t2.ut1  # doctest: +SKIP
    Traceback (most recent call last):
    ...
    IERSRangeError: (some) times are outside of range covered by IERS table.

In this case, one needs to update the IERS B table or use IERS A instead
(which also has predictions).  In future versions, this may become
configurable or automatic, but currently it requires handiwork.  For
`Time`, easiest is to set the `delta_ut1_utc` property directly::

    >>> from astropy.utils.iers import IERS_A
    >>> iers_a = IERS_A.open('finals2000A.all')    # doctest: +SKIP
    >>> iers_a.ut1_utc(t2)                         # doctest: +SKIP
    <Quantity -0.5631187323911956 s>
    >>> t2.delta_ut1_utc = iers_a.ut1_utc(t2)      # doctest: +SKIP
    >>> t2.ut1.iso                                 # doctest: +SKIP
    '2015-03-24 03:33:57.288'

Instead of local copies of IERS files, one can also download them, using
`iers.IERS_A_URL` and `iers.IERS_B_URL`::

    >>> from astropy.utils.iers import IERS_A, IERS_A_URL
    >>> iers_a = IERS_A.open(IERS_A_URL)                    # doctest: +SKIP
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse

import numpy as np

from ...table import Table, QTable
from ...utils.data import get_pkg_data_filename, download_file

__all__ = ['IERS', 'IERS_B', 'IERS_A',
           'FROM_IERS_B', 'FROM_IERS_A', 'FROM_IERS_A_PREDICTION',
           'TIME_BEFORE_IERS_RANGE', 'TIME_BEYOND_IERS_RANGE',
           'IERS_A_FILE', 'IERS_A_URL', 'IERS_A_README',
           'IERS_B_FILE', 'IERS_B_URL', 'IERS_B_README',
           'IERSRangeError']

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
            for passing on to the `read` class methods (further optional
            arguments that are available for some IERS subclasses can be added).
            If None, use the default location from the `read` class method.
        cache : bool
            Whether to use cache. Defaults to False, since IERS files
            are regularly updated.

        Returns
        -------
        An IERS table class instance

        Notes
        -----
        On the first call in a session, the table will be memoized (in
        `cls.iers_table`), and further calls to `open` will return this stored
        table if `file=None` (the default).

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
            second part of two-part JD (default: 0., ignored if jd1 is Time)

        Returns
        -------
        mjd : float or array
            integer part of MJD
        utc : float or array
            fractional part of MJD
        """
        try:  # see if this is a Time object
            jd1, jd2 = jd1.utc.jd1, jd1.utc.jd2
        except:
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
            second part of two-part JD (default 0., ignored if jd1 is Time)
        return_status : bool
            Whether to return status values.  If False (default),
            raise `IERSRangeError` if any time is out of the range covered
            by the IERS table.

        Returns
        -------
        ut1_utc : float or float array
            UT1-UTC, interpolated in IERS Table
        status : int or int array
            Status values (if `return_status`=`True`)::
            `iers.FROM_IERS_B`
            `iers.FROM_IERS_A`
            `iers.FROM_IERS_A_PREDICTION`
            `iers.TIME_BEFORE_IERS_RANGE`
            `iers.TIME_BEYOND_IERS_RANGE`
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
            second part of two-part JD (default 0., ignored if jd1 is Time)
        return_status : bool
            Whether to return status values.  If False (default),
            raise `IERSRangeError` if any time is out of the range covered
            by the IERS table.

        Returns
        -------
        PM_x : Quantity with angle units
            x component of polar motion for the requested times
        PM_y : Quantity with angle units
            y component of polar motion for the requested times
        status : int or int array
            Status values (if `return_status`=`True`)::
            `iers.FROM_IERS_B`
            `iers.FROM_IERS_A`
            `iers.FROM_IERS_A_PREDICTION`
            `iers.TIME_BEFORE_IERS_RANGE`
            `iers.TIME_BEYOND_IERS_RANGE`
        """
        return self._interpolate(jd1, jd2, ['PM_x', 'PM_y'],
                                 self.pm_source if return_status else None)

    def _interpolate(self, jd1, jd2, columns, source=None):
        mjd, utc = self.mjd_utc(jd1, jd2)
        # enforce array
        is_scalar = not hasattr(mjd, '__array__') or mjd.ndim == 0
        if is_scalar:
            mjd = np.array([mjd])
            utc = np.array([utc])
        # For typical format, will always find a match (since MJD are integer)
        # hence, important to define which side we will be; this ensures
        # self['MJD'][i-1]<=mjd<self['MJD'][i]
        i = np.searchsorted(self['MJD'].value, mjd, side='right')
        # Get index to MJD at or just below given mjd, clipping to ensure we
        # stay in range of table (status will be set below for those outside)
        i1 = np.clip(i, 1, len(self)-1)
        i0 = i1-1
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
            # Not returning status, so raise an exception for out-of-range
            if np.any(i1 != i):
                raise IERSRangeError('(some) times are outside of range covered '
                                 'by IERS table.')
            return results[0] if len(results) == 1 else results

    def ut1_utc_source(self, i):
        """Source for UT1-UTC.  To be overridden by subclass."""
        return np.zeros_like(i)

    def pm_source(self, i):
        """Source for polar motion.  To be overridden by subclass."""
        return np.zeros_like(i)


class IERS_A(IERS):
    """IERS Table class targeted to IERS A, provided by USNO.

    These include rapid turnaround and predicted times.
    See http://maia.usno.navy.mil/

    Notes
    -----
    The IERS A file is not part of astropy.  It can be downloaded from
    `iers.IERS_A_URL`.  See `iers.__doc__` for instructions on how to use
    it in `Time`, etc.
    """

    iers_table = None

    def __init__(self, table):
        """Initialize an IERS-A table that is already read in.
        Use read or open class methods to read it from disk.

        Combines UT1-UTC values, taking UT1_UTC_B if available, else UT1_UTC_A
        """
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
        # Fill any masked values, and convert to a QTable.
        super(IERS_A, self).__init__(table.filled())

    @classmethod
    def read(cls, file=IERS_A_FILE, readme=None):
        """Read IERS-A table from a finals2000a.* file provided by USNO.

        Parameters
        ----------
        file : str
            full path to ascii file holding IERS-A data
            (default: `iers.IERS_A_FILE`)
        readme : str
            full path to ascii file holding CDS-style readme
            (default: package version, `iers.IERS_A_README`)

        Returns
        -------
        `IERS_A` class instance
        """
        if readme is None:
            readme = IERS_A_README
        # Read in as a regular Table, including possible masked columns.
        # Columns will be filled and converted to Quantity in cls.__init__.
        iers_a = Table.read(file, format='cds', readme=readme)
        # IERS A has some rows at the end that hold nothing but dates & MJD
        # presumably to be filled later.  Exclude those a priori -- there
        # should at least be a predicted UT1-UTC and PM!
        return cls(iers_a[~iers_a['UT1_UTC_A'].mask &
                          ~iers_a['PolPMFlag_A'].mask])

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
    If the package IERS B file (`iers.IERS_B_FILE`) is out of date, a new
    version can be downloaded from `iers.IERS_B_URL`.  See `iers.__doc__`
    for instructions on how to use it in `Time`, etc.
    """

    iers_table = None

    @classmethod
    def read(cls, file=None, readme=None, data_start=14):
        """Read IERS-B table from a eopc04_iau2000.* file provided by IERS.

        Parameters
        ----------
        file : str
            full path to ascii file holding IERS-B data
            (default: package version, `iers.IERS_B_FILE`)
        readme : str
            full path to ascii file holding CDS-style readme
            (default: package version, `iers.IERS_B_README`)
        data_start : int
            starting row (default: 14, appropriate for standard IERS files)

        Returns
        -------
        `IERS_B` class instance
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

# by default for IERS class, read IERS-B table
IERS.read = IERS_B.read
