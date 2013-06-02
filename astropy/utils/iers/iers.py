from __future__ import division, print_function
import warnings
import numpy as np

from ...table import Table
from ...time import Time
from ...config import ConfigurationItem
from ...utils.data import get_pkg_data_filename

FINALS2000A = ConfigurationItem('IERS_A_DATA', 'finals2000A.all',
                                'Path to IERS A data file (IAU2000 model)')
FINALS2000A_URL = ConfigurationItem(
    'IERS_A_URL', 'http://maia.usno.navy.mil/ser7/finals2000A.all',
    'URL to IERS A final2000A.all')
README_FINALS2000A = get_pkg_data_filename('ReadMe.finals2000A')

EOPC04_IAU2000 = ConfigurationItem('IERS_B_DATA', 'eopc04_IAU2000.62-now',
                                   'Path to IERS B data file (C04 series)')
EOPC04_IAU2000_URL = ConfigurationItem(
    'IERS_B_URL',
    'http://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now',
    'URL to IERS B C04 series eopc04_IAU2000.62-now')
README_EOPC04_IAU2000 = get_pkg_data_filename('ReadMe.eopc04_IAU2000')

# Status values
FROM_IERS_B = 0
FROM_IERS_A = 1
FROM_IERS_A_PREDICTION = 2
TIME_BEFORE_IERS_RANGE = -1
TIME_AFTER_IERS_RANGE = -2

MJD_ZERO = 2400000.5


class IERS(Table):
    """Generic IERS table, which defines the functions to return interpolated
    UT1-UTC values, etc.  Should hold columns 'MJD' and 'UT1_UTC'
    """
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

    def ut1_utc(self, jd1, jd2=0., return_status=False):
        """Interpolate UT1-UTC corrections in IERS Table for given dates.

        Parameters
        ----------
        jd1: float, array, or Time
            first part of two-part JD, or Time object
        jd2: float or array, optional
            second part of two-part JD (default: 0., ignored if jd1 is Time)
        return_status: bool, optional
            whether to return status accompanying each result (default: False)

        Returns
        -------
        ut1_utc: float (array)
            UT1-UTC, interpolated in IERS Table
        status: int (array), if return_status is True

        Note
        ----
        Status values are as follows.  If return_status is False, warnings
        are given for values not equal to 0 or 1.
            0: Interpolated in IERS B values
            1: Interpolated in IERS A values
            2: Interpolated in IERS A predictions
           -1: Time falls before IERS table range
           -2: Time falls after IERS table range

        For status -2 and 2, the IERS table may need to be updated.
        """

        mjd, utc = self.mjd_utc(jd1, jd2)
        # for typical format, will always find a match (since MJD are integer)
        # hence, important to define which side we will be; this ensures
        # mjd<self['MJD'][i]
        i = np.searchsorted(self['MJD'], mjd, side='right')
        # get index to MJD at or just below given mjd, clipping to ensure
        # we stay in range of table (status will be set below)
        i0 = np.clip(i-1, 0, len(self)-2)
        mjd_0, mjd_1 = self['MJD'][i0], self['MJD'][i0+1]
        ut1_utc_0, ut1_utc_1 = self['UT1_UTC'][i0], self['UT1_UTC'][i0+1]
        # check and correct for possible leap second (correcting difference,
        # not first point, since jump happens right at second point)
        d_ut1_utc = ut1_utc_1 - ut1_utc_0
        d_ut1_utc += d_ut1_utc < 0.5
        d_ut1_utc -= d_ut1_utc > 0.5
        # and linearly interpolate to get UT1-UTC
        ut1_utc = ut1_utc_0 + (mjd-mjd_0+utc)/(mjd_1-mjd_0)*d_ut1_utc

        status = (TIME_BEFORE_IERS_RANGE * (i == 0) +
                  TIME_AFTER_IERS_RANGE * (i == len(self)))
        # better safe than sorry - no extrapolation
        ut1_utc *= (status == 0)

        if 'UT1Flag' in self.colnames:
            ut1flag = self['UT1Flag'][i0]
            status += (status == 0) * \
                      (FROM_IERS_B * (ut1flag == 'B') +
                       FROM_IERS_A * (ut1flag == 'I') +
                       FROM_IERS_A_PREDICTION * (ut1flag == 'P'))

        if return_status:
            return ut1_utc, status
        else:
            if np.any(status == TIME_BEFORE_IERS_RANGE):
                warnings.warn('some times fall before the IERS table range.')
            if np.any(status == TIME_AFTER_IERS_RANGE):
                warnings.warn('some times fall after the IERS table range. ' +
                              'Is your table out of date?')
            if np.any(status == FROM_IERS_A_PREDICTION):
                warnings.warn('for some times, only predicted values are ' +
                              'available.  Is your IERS table out of date?')
            return ut1_utc


class IERS_A(IERS):
    """IERS table targeted to files provided by USNO, which include rapid
    turnaround and predicted times (IERS A)
    """
    def __init__(self, table, **kwargs):
        # combine UT1_UTC, taking UT1_UTC_B if available, else UT1_UTC_A
        Table.__init__(self, table)
        self['UT1_UTC'] = np.where(self['UT1_UTC_B'].mask,
                                   self['UT1_UTC_A'],
                                   self['UT1_UTC_B'])
        self['UT1Flag'] = np.where(self['UT1_UTC_B'].mask,
                                   self['UT1Flag_A'],
                                   np.array(['B']*len(self)))

    @classmethod
    def read(cls, file=FINALS2000A(), readme=README_FINALS2000A, **kwargs):
        iers_a = Table.read(file, format='cds', readme=readme)
        return cls(iers_a[~iers_a['UT1_UTC_A'].mask], **kwargs)


class IERS_B(IERS):
    @classmethod
    def read(cls, file=EOPC04_IAU2000(), readme=README_EOPC04_IAU2000):
        # must be possible to do this more elegantly!
        iers_b = Table.read(file, format='cds', readme=readme, data_start=14)
        return cls(iers_b)
