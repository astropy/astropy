# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The astropy.utils.iers package provides access to the tables provided by the
International Earth Rotation Service, in particular allowing interpolation of
published UT1-UTC values for given times.  These are used in astropy.time to
provide UT1 values.  By default, IERS B values provided as part of astropy are
used, but user-downloaded files can be substituted.

Generally, there is no need to invoke the iers classes oneself.  E.g., the
default IERS B table is loaded as necessary in `Time`::
    >>> t = Time(['2012-06-30 12:00:00', '2012-06-30 23:59:59',
                  '2012-06-30 23:59:60', '2012-07-01 00:00:00',
                  '2012-07-01 12:00:00'], scale='utc')
    >>> t.ut1
    <Time object: scale='ut1' format='iso' vals=['2012-06-30 11:59:59.413'
     '2012-06-30 23:59:58.413' '2012-06-30 23:59:59.413'
     '2012-07-01 00:00:00.413' '2012-07-01 12:00:00.413']>

But if one is dealing with very recent observations, this does now work::
    >>> t2 = Time.now()
    >>> t2.ut1
    ERROR: ValueError: (some) times are beyond range covered by IERS table.
    [astropy.time.core]

In this case, one needs to update the IERS B table or use IERS A instead
(which also has predictions).  In future versions, this may become configurable
or automatic, but currently it requires some handiwork.  For `Time`, the
easiest option is to set the `delta_ut1_utc` property directly::
    >>> from astropy.utils.iers import IERS_A
    >>> a = IERS_A.open('finals2000A.all')
    >>> a.ut1_utc(t2)
    (0.069727551794218745, 2)
    >>> t2.delta_ut1_utc = a.ut1_utc(t2)[0]
    >>> t2.ut1
    <Time object: scale='ut1' format='datetime' vals=2013-06-14 02:31:40.441858>

Note that the status returned by ut1_utc should be checked: negative values
indicate a time out of range.  The 2 above is `iers.FROM_IERS_A_PREDICTION`.

(The IERS-A file `finals2000A.all` can be downloaded from `iers.IERS_A_URL`)
"""

from __future__ import division
import numpy as np

from ...table import Table
from ...utils.data import get_pkg_data_filename

__all__ = ['IERS', 'IERS_B', 'IERS_A',
           'FROM_IERS_B', 'FROM_IERS_A', 'FROM_IERS_A_PREDICTION',
           'TIME_BEFORE_IERS_RANGE', 'TIME_BEYOND_IERS_RANGE',
           'IERS_A_FILE', 'IERS_A_URL', 'IERS_A_README',
           'IERS_B_FILE', 'IERS_B_URL', 'IERS_B_README']

# IERS-A default file name, URL, and ReadMe with content description
IERS_A_FILE = 'finals2000A.all'
IERS_A_URL = 'http://maia.usno.navy.mil/ser7/finals2000A.all'
IERS_A_README = get_pkg_data_filename('data/ReadMe.finals2000A')
# IERS-B default file name, URL, and ReadMe with content description
IERS_B_FILE = get_pkg_data_filename('data/eopc04_IAU2000.62-now.gz')
IERS_B_URL = 'http://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now'
IERS_B_README = get_pkg_data_filename('data/ReadMe.eopc04_IAU2000')
# Status/source values returned by IERS.ut1_utc
FROM_IERS_B = 0
FROM_IERS_A = 1
FROM_IERS_A_PREDICTION = 2
TIME_BEFORE_IERS_RANGE = -1
TIME_BEYOND_IERS_RANGE = -2


MJD_ZERO = 2400000.5


class IERS(Table):
    """Generic IERS table class, defining interpolation functions.

    Sub-classed from `astropy.table.Table`.  The table should hold columns
    'MJD' and 'UT1_UTC'.
    """

    iers_table = None

    @classmethod
    def open(cls, *args, **kwargs):
        """Open an IERS table, reading it from a file if not loaded before.

        Returns
        -------
        An IERS table class instance

        Notes
        -----
        If it exists, a previously stored table is returned (`cls.iers_table`)
        If not, the `read` method is called, passing on all parameters.

        If a table needs to be re-read from disk, use the (sub-class) close
        method and re-open.

        For the IERS class itself, an IERS_B sub-class instance is opened.
        """
        if cls.iers_table is None:
            cls.iers_table = cls.read(*args, **kwargs)
        return cls.iers_table

    @classmethod
    def close(cls):
        """Remove the IERS table from the class.

        This allows the table to be re-read from disk during one's session
        (e.g., if one finds it is out of date and has updated the file.
        """
        cls.iers_table = None

    def mjd_utc(self, jd1, jd2=0.):
        """Turn a time to MJD, returning integer and fractional parts.

        Parameters
        ----------
        jd1: float, array, or Time
            first part of two-part JD, or Time object
        jd2: float or array, optional
            second part of two-part JD (default: 0., ignored if jd1 is Time)

        Returns
        -------
        mjd: float or array
            integer part of MJD
        utc: float or array
            fractional part of MJD
        """
        try:  # see if this is a Time object
            jd1, jd2 = jd1.utc.jd1, jd1.utc.jd2
        except:
            pass
        mjd = np.floor(jd1 - MJD_ZERO + jd2)
        utc = jd1 - (MJD_ZERO+mjd) + jd2
        return mjd, utc

    def ut1_utc(self, jd1, jd2=0.):
        """Interpolate UT1-UTC corrections in IERS Table for given dates.

        Parameters
        ----------
        jd1: float, float array, or Time object
            first part of two-part JD, or Time object
        jd2: float or float array, optional
            second part of two-part JD (default: 0., ignored if jd1 is Time)

        Returns
        -------
        ut1_utc: float or float array
            UT1-UTC, interpolated in IERS Table
        status: int or int array
            Status values, as follows::
            `iers.FROM_IERS_B`
            `iers.FROM_IERS_A`
            `iers.FROM_IERS_A_PREDICTION`
            `iers.TIME_BEFORE_IERS_RANGE`
            `iers.TIME_BEYOND_IERS_RANGE`

        The status values are defined as 0, 1, 2, -1, -2, respectively, but
        but this may change. Always, zero means a definitive result, positive
        values a preliminary result or a prediction, and negative values
        indicate a problem.
        """

        mjd, utc = self.mjd_utc(jd1, jd2)
        # enforce array
        is_scalar = not hasattr(mjd, '__array__') or mjd.ndim == 0
        if is_scalar:
            mjd = np.array([mjd])
            utc = np.array([utc])
        # For typical format, will always find a match (since MJD are integer)
        # hence, important to define which side we will be; this ensures
        # self['MJD'][i-1]<=mjd<self['MJD'][i]
        i = np.searchsorted(self['MJD'], mjd, side='right')
        # Get index to MJD at or just below given mjd, clipping to ensure we
        # stay in range of table (status will be set below for those outside)
        i1 = np.clip(i, 1, len(self)-1)
        i0 = i1-1
        mjd_0, mjd_1 = self['MJD'][i0], self['MJD'][i1]
        ut1_utc_0, ut1_utc_1 = self['UT1_UTC'][i0], self['UT1_UTC'][i1]
        # Check and correct for possible leap second (correcting difference,
        # not first point, since jump can only happen right at second point)
        d_ut1_utc = ut1_utc_1 - ut1_utc_0
        d_ut1_utc -= np.round(d_ut1_utc)
        # Linearly interpolate to get UT1-UTC
        # (which is what TEMPO does, but may want to follow IERS gazette #13
        # for more precise interpolation and correction for tidal effects;
        # http://maia.usno.navy.mil/iers-gaz13)
        ut1_utc = ut1_utc_0 + (mjd-mjd_0+utc)/(mjd_1-mjd_0)*d_ut1_utc

        # Set status to source, possibly using routine provided by subclass
        status = self.ut1_utc_source(i1)
        # Check for out of range - more important than above, so OK to override
        # also reset any extrapolated values - better safe than sorry
        if np.any(i1 != i):
            status[i == 0] = TIME_BEFORE_IERS_RANGE
            ut1_utc[i == 0] = 0.
            status[i == len(self)] = TIME_BEYOND_IERS_RANGE
            ut1_utc[i == len(self)] = 0.

        if is_scalar:
            ut1_utc = ut1_utc[0]
            status = status[0]

        return ut1_utc, status

    def ut1_utc_source(self, i):
        """Source for UT1-UTC.  To be overridden by subclass."""
        return np.zeros_like(i)


class IERS_A(IERS):
    """IERS Table class targeted to IERS A, provided by USNO.

    These include rapid turnaround and predicted times.
    See http://maia.usno.navy.mil/

    Notes
    -----
    The IERS A file is not part of astropy.  It can be downloaded from
    `iers.IERS_A_URL`.  See `iers.__doc__` for instructions on how to enable
    its use in `Time`, etc.
    """

    iers_table = None

    def __init__(self, table):
        """Initialize an IERS-A table that is already read in.
        Use read or open class methods to read it from disk.

        Combines UT1-UTC values, taking UT1_UTC_B if available, else UT1_UTC_A
        """
        table['UT1_UTC'] = np.where(table['UT1_UTC_B'].mask,
                                    table['UT1_UTC_A'],
                                    table['UT1_UTC_B'])
        table['UT1Flag'] = np.where(table['UT1_UTC_B'].mask,
                                    table['UT1Flag_A'],
                                    'B')
        super(IERS_A, self).__init__(table.filled())

    @classmethod
    def read(cls, file=IERS_A_FILE, readme=IERS_A_README):
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
        iers_a = Table.read(file, format='cds', readme=readme)
        # IERS A has some rows at the end that hold nothing but dates & MJD
        # presumably to be filled later.  Exclude those a priori -- there
        # should at least be a predicted UT1-UTC!
        return cls(iers_a[~iers_a['UT1_UTC_A'].mask])

    def ut1_utc_source(self, i):
        """Set UT1-UTC source flag for entries in IERS table"""
        ut1flag = self['UT1Flag'][i]
        source = np.ones_like(i) * FROM_IERS_B
        source[ut1flag == 'I'] = FROM_IERS_A
        source[ut1flag == 'P'] = FROM_IERS_A_PREDICTION
        return source


class IERS_B(IERS):
    """IERS Table class targeted to IERS B, provided by IERS itself.

    These are final values; see http://www.iers.org/

    Notes
    -----
    If the package IERS B file (`iers.IERS_B_FILE`) is out of date, a new
    version can be downloaded from `iers.IERS_B_URL`.  See `iers.__doc__`
    for instructions on how to enable its use in `Time`, etc.
    """

    iers_table = None

    @classmethod
    def read(cls, file=IERS_B_FILE, readme=IERS_B_README, data_start=14):
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
        # can this be done more elegantly, initialising directly, without
        # passing a Table to the Table initialiser?
        iers_b = Table.read(file, format='cds', readme=readme,
                            data_start=data_start)
        return cls(iers_b)

    def ut1_utc_source(self, i):
        """Set UT1-UTC source flag for entries in IERS table"""
        return np.ones_like(i) * FROM_IERS_B

# by default for IERS class, read IERS-B table
IERS.read = IERS_B.read
