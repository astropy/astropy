from __future__ import division, print_function
import numpy as np

from ...table import Table
from ...time import Time
from ...utils.data import get_pkg_data_filename

FINALS2000A = 'finals2000A.all'
FINALS2000A_URL = 'http://maia.usno.navy.mil/ser7/finals2000A.all'
README_FINALS2000A = get_pkg_data_filename('ReadMe.finals2000A')

EOPC04_IAU2000 = 'eopc04_IAU2000.62-now'
EOPC04_IAU2000_URL = \
    'http://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now'
README_EOPC04_IAU2000 = get_pkg_data_filename('ReadMe.eopc04_IAU2000')

# Status/source values
FROM_IERS_B = 0
FROM_IERS_A = 1
FROM_IERS_A_PREDICTION = 2
TIME_BEFORE_IERS_RANGE = -1
TIME_BEYOND_IERS_RANGE = -2

MJD_ZERO = 2400000.5


class IERS(Table):
    """Generic IERS table class, defining interpolation functions.

    Should hold columns 'MJD' and 'UT1_UTC'
    """

    iers_table = None

    @classmethod
    def open(cls, *args, **kwargs):
        if cls.iers_table is None:
            cls.iers_table = cls.read(*args, **kwargs)
        return cls.iers_table

    @classmethod
    def close(cls):
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
        if isinstance(jd1, Time):
            jd1, jd2 = jd1.utc.jd1, jd1.utc.jd2
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

             0 : Interpolated in IERS B values
             1 : Interpolated in IERS A values
             2 : Interpolated in IERS A predictions
            -1 : Time falls before IERS table range
            -2 : Time falls beyond IERS table range
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

    # this should be overridden by subclasses
    def ut1_utc_source(self, i):
        return np.zeros_like(i)


class IERS_A(IERS):
    """IERS Table class targeted to IERS A, provided by USNO.

    These include rapid turnaround and predicted times.
    See http://maia.usno.navy.mil/
    """
    def __init__(self, table):
        # combine UT1_UTC, taking UT1_UTC_B if available, else UT1_UTC_A
        table['UT1_UTC'] = np.where(table['UT1_UTC_B'].mask,
                                    table['UT1_UTC_A'],
                                    table['UT1_UTC_B'])
        table['UT1Flag'] = np.where(table['UT1_UTC_B'].mask,
                                    table['UT1Flag_A'],
                                    'B')
        super(IERS_A, self).__init__(table.filled())

    @classmethod
    def read(cls, file=FINALS2000A, readme=README_FINALS2000A):
        iers_a = Table.read(file, format='cds', readme=readme)
        # IERS A has some rows at the end that hold nothing but dates & MJD
        # presumably to be filled later.  Exclude those a priori -- there
        # should at least be a predicted UT1-UTC!
        return cls(iers_a[~iers_a['UT1_UTC_A'].mask])

    def ut1_utc_source(self, i):
        ut1flag = self['UT1Flag'][i]
        source = np.ones_like(i) * FROM_IERS_B
        source[ut1flag == 'I'] = FROM_IERS_A
        source[ut1flag == 'P'] = FROM_IERS_A_PREDICTION
        return source


class IERS_B(IERS):
    """IERS Table class targeted to IERS B, provided by IERS itself.

    These are final values; see http://www.iers.org/
    """
    @classmethod
    def read(cls, file=EOPC04_IAU2000, readme=README_EOPC04_IAU2000):
        # can this be done more elegantly, initialising directly, without
        # passing a Table to the Table initialiser?
        iers_b = Table.read(file, format='cds', readme=readme, data_start=14)
        return cls(iers_b)

    def ut1_utc_source(self, i):
        return np.ones_like(i) * FROM_IERS_B
