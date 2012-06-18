# Licensed under a 3-clause BSD style license - see LICENSE.rst

import re
import datetime
import types
import numpy as np

from ..utils import vectorize
from ..io.fits.util import lazyproperty, isiterable

# Defining zeropoints where TCB, TCG and TT are linked to ephemeris time
# wikipedia http://en.wikipedia.org/wiki/International_Atomic_Time
SECONDS_PER_DAY = 86400.0
CALENDAR_ZEROPOINT = datetime.datetime(1977, 1, 1, 0, 0, 32, 184000)
JD_ZEROPOINT = 2443144.5003725
JD_ZEROPOINT_SECONDS = JD_ZEROPOINT * SECONDS_PER_DAY

TAI2TT = 32.184  # TT - TAI
MJD2JD = 2400000.5 # JD - MJD

# JD to UTC leap seconds
# http://maia.usno.navy.mil/ser7/tai-utc.dat
# TODO: Consider automatic leapsecs lookup for futures times that we forget to
# add to this table
LEAPSECS = {
    2441317.5: 10.0,
    2441499.5: 11.0,
    2441683.5: 12.0,
    2442048.5: 13.0,
    2442413.5: 14.0,
    2442778.5: 15.0,
    2443144.5: 16.0,
    2443509.5: 17.0,
    2443874.5: 18.0,
    2444239.5: 19.0,
    2444786.5: 20.0,
    2445151.5: 21.0,
    2445516.5: 22.0,
    2446247.5: 23.0,
    2447161.5: 24.0,
    2447892.5: 25.0,
    2448257.5: 26.0,
    2448804.5: 27.0,
    2449169.5: 28.0,
    2449534.5: 29.0,
    2450083.5: 30.0,
    2450630.5: 31.0,
    2451179.5: 32.0,
    2453736.5: 33.0,
    2454832.5: 34.0,
    2456109.5: 35.0
}


ISO8601_RE = re.compile("(?P<year>\d{4})"
           "(-?(?P<month>\d{1,2})"
           "(-?(?P<day>\d{1,2})"
           "((?P<separator>.)"
           "(?P<hour>\d{2})"
           ":?(?P<minute>\d{2})"
           "(:?(?P<second>\d{2})"
           "(\.(?P<fraction>\d+))?)?"
           "(?P<timezone>Z|(([-+])(\d{2}):(\d{2})))?)?)?)?")


@vectorize(otypes=(np.int64,))
def get_total_seconds(time1, time2):
    """Calculate the total seconds between two `datetime.datetime` objects."""

    td = time1 - time2
    return td.seconds + td.days * 24 * 3600


@vectorize(otypes=(np.float64,))
def get_subseconds(time1, time2):
    """
    Return only the fraction of seconds difference between two
    `datetime.datetime` objects.
    """

    td = time1 - time2
    return td.microseconds / 1e6


@vectorize
def convert_seconds_to_timedelta(seconds, s):
    return datetime.timedelta(seconds=int(seconds),
                              microseconds=s * 1e6)


class Time(object):
    """
    Class representing a time scalar.

    The internal format uses JD 2443144.5003725 (1 January 1977 00:00:32.184)
    as the zeropoint (the instant where TCB, TCG and TT were the same) and
    stores the date as days to this zeropoint in `decimal.Decimal`

    Initialize an AstroTime-object with seconds from JD 2443144.5003725
    Parameters
    ----------
    seconds : int
        The number of seconds since 1 January 1977 00:00:32.184

    subseconds : float
        Fraction of a second for higher precision.
    """

    class __metaclass__(type):
        def __new__(mcs, name, bases, d):
            for attr, item in d.items():
                if not (isinstance(item, types.FunctionType) and
                        attr.startswith('_tai2')):
                    continue

                system = attr[5:]
                if not '_%s2tai' % system in d:
                    continue

                def prop(self, system=system):
                    return self._convert_system(self, system)

                prop.__name__ = system
                prop.__doc__ = ('Convert to a new %s object using the %s '
                                'system.' % (name, system.upper()))

                d[system] = lazyproperty(prop)

            return type.__new__(mcs, name, bases, d)

    def __init__(self, seconds, subseconds=None, system='tai'):
        self.seconds = np.int64(seconds)
        if subseconds is not None:
            self.subseconds = np.float64(subseconds)
        else:
            self.subseconds = None

        system = system.lower()
        if not hasattr(self, '_tai2' + system):
            raise ValueError('Unknown/unsupported time system: %r' % system)
        self.system = system

    def __repr__(self):
        # TODO: Return a useful repr() of the time including its system
        return super(Time, self).__repr__()

    # TODO: These from_ methods could probably also be generated given the
    # correct conversion functions
    @classmethod
    def from_tai(cls, tai, precise=True):
        """
        Initialize from a TAI date and time (using `datetime.datetime`).

        Parameters
        ----------
        calendar_date : `datetime.datetime`, str, or list
            String values should be given in an ISO8601 conforming format.

        Examples
        --------

        >>> from astropy import time
        >>> import datetime
        >> dt = datetime.datetime(1546, 12, 14, 12, 0, 0)
        >>> mytime = time.Time.from_tai(dt)
        >>> mytime.jd
        2286072.0

        References
        ----------
        http://asa.usno.navy.mil/SecM/Glossary.html
        http://en.wikipedia.org/wiki/Julian_day#Converting_Gregorian_calendar_date_to_Julian_Day_Number
        """

        if isinstance(tai, basestring):
            #iso8601 parsing
            isarray = False
            raise NotImplementedError('ISO8601 parsing not available yet')

        elif isinstance(tai, datetime.datetime):
            time_to_zeropoint = (tai - CALENDAR_ZEROPOINT)
            seconds = get_total_seconds(tai, CALENDAR_ZEROPOINT)

            if precise:
                subseconds = time_to_zeropoint.microseconds / 1e6
            else:
                subseconds = None

            return cls(seconds=seconds, subseconds=subseconds)

        elif isiterable(tai):
            tais = np.array(tai)
            if isinstance(tais[0], basestring):
                NotImplementedError('ISO8601 parsing not available yet')

            elif isinstance(tais[0], datetime.datetime):
                # TODO: I'm pretty sure Numpy caches calls to vectorize, but it
                # might make more sense to do this at the module level anyways
                seconds = get_total_seconds(tais, CALENDAR_ZEROPOINT)
                if precise:
                    subseconds = get_subseconds(tais, CALENDAR_ZEROPOINT)
                else:
                    subseconds = None

                return cls(seconds=seconds, subseconds=subseconds)

    # TODO: Allow from_jd and from_mjd to be linked to an arbitrary time system
    @classmethod
    def from_jd(cls, jd_time, precise=True):
        """
        Instantiate a Time-object with Julian Date (linked to TAI)

        Parameters
        ----------
        jd_time : float or list of floats
        A Julian date
        :param jd_time:
            A float object
        """

        # TODO: Support conversion from a string as well

        if isiterable(jd_time):
            jd_time = np.array(jd_time)

        seconds = jd_time * SECONDS_PER_DAY - JD_ZEROPOINT_SECONDS
        int_seconds = np.int64(seconds)
        if precise:
            subseconds = np.float64(seconds - int_seconds)
        else:
            subseconds = None

        return cls(seconds=int_seconds, subseconds=subseconds)

    @classmethod
    def from_mjd(cls, mjd_time, low_precision=False):
        """
        Instantiate an AstroTime-object with Modified Julian Date

        Parameters
        ----------
        mjd_time : float
            A Modified Julian date

        """
        return cls.from_jd(mjd_time + MJD2JD)

    @lazyproperty
    def jd(self):
        """Return the date as JD (system) in a float64"""

        if self.subseconds is None:
            secs = JD_ZEROPOINT_SECONDS + self.seconds
        else:
            secs = JD_ZEROPOINT_SECONDS + self.seconds + self.subseconds

        return secs / SECONDS_PER_DAY

    @lazyproperty
    def mjd(self):
        return self.jd - MJD2JD

    # TODO: Make this not linked exclusively to TAI
    @lazyproperty
    def datetime(self):
        """
        Returns the TAI date in a `datetime.datetime` object

        References
        ----------
        http://www.usno.navy.mil/USNO/astronomical-applications/astronomical-information-center/julian-date-form
        """

        if self.subseconds is None:
            subseconds = 0.0
        else:
            subseconds = self.subseconds

        if isinstance(self.seconds, np.ndarray):
            td = convert_seconds_to_timedelta(self.tai.seconds,
                                              self.tai.subseconds)
        else:
            td = datetime.timedelta(seconds=int(self.tai.seconds),
                                    microseconds=self.tai.subseconds * 1e6)

        return CALENDAR_ZEROPOINT + td

    @classmethod
    def _convert_system(cls, timeobj, new_system):
        new_system = new_system.lower()

        if timeobj.system == new_system:
            # TODO: Maybe consider making a copy instead of returning the same
            # object?  If Time objects are meant to be immutable, is there any
            # point to copying?
            return timeobj

        try:
            if new_system == 'tai':
                conv_meth = getattr(timeobj, '_%s2tai' % timeobj.system)
            else:
                conv_meth = getattr(timeobj.tai, '_tai2' + new_system)
        except AttributeError:
            raise ValueError('%r cannot be converted to %r; the required '
                             'conversion methods do not exist' %
                             (timeobj, new_system))

        return cls(*conv_meth(), system=new_system)

    def _tai2tai(self):
        return self.seconds, self.subseconds

    def _tai2tt(self):
        ttsecs = np.int64(TAI2TT)
        if self.subseconds is None:
            return (self.seconds + ttsecs, None)
        else:
            ttsubsecs = np.float64(TAI2TT - ttsecs)
            return (self.seconds + ttsecs, self.subseconds + ttsubsecs)

    def _tt2tai(self):
        ttsecs = np.int64(TAI2TT)
        if self.subseconds is None:
            return (self.seconds - ttsecs, None)
        else:
            ttsubsecs = np.float64(TAI2TT - ttsecs)
            return (self.seconds - ttsecs, self.subseconds - ttsubsecs)
