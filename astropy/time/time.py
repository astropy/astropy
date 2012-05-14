# Licensed under a 3-clause BSD style license - see LICENSE.rst

import re
import datetime
import decimal
import numpy as np

#Defining zeropoints where TCB, TCG and TT are linked to ephemeris time
#wikipedia http://en.wikipedia.org/wiki/International_Atomic_Time

calendar_zeropoint = datetime.datetime(1977, 1, 1, 0, 0, 32, 184000)
jd_zeropoint = 2443144.5003725
jd_zeropoint_seconds = jd_zeropoint * 86400


iso8601_re = re.compile("(?P<year>\d{4})"
           "(-?(?P<month>\d{1,2})"
           "(-?(?P<day>\d{1,2})"
           "((?P<separator>.)"
           "(?P<hour>\d{2})"
           ":?(?P<minute>\d{2})"
           "(:?(?P<second>\d{2})"
           "(\.(?P<fraction>\d+))?)?"
           "(?P<timezone>Z|(([-+])(\d{2}):(\d{2})))?)?)?)?")




def get_total_seconds(time1, time2):
    #calculate the total seconds between two `datetime.datetime` objects
    td = time1 - time2
    return td.seconds + td.days * 24 * 3600
    
def get_fraction_of_seconds(time1, time2):
    #return only the fraction of seconds difference between two `datetime.datetime` objects
    td = time1 - time2
    return td.microseconds / 1e6

def convert_seconds_to_timedelta(seconds, fraction_of_seconds):
    return datetime.timedelta(seconds=seconds, microseconds=fraction_of_seconds * 1e6)

class Time(object):
    """Class to store a time variable.
    The internal format uses JD 2443144.5003725 (1 January 1977 00:00:32.184) as the zeropoint
    (the instant where TCB, TCG and TT were the same)
    and stores the date as days to this zeropoint in `decimal.Decimal`
    
    Initialize an AstroTime-object with seconds from JD 2443144.5003725
    Parameters
    ----------
    seconds : `decimal.Decimal`
        The number of seconds since 1 January 1977 00:00:32.184
    """
    
    @classmethod
    def from_tai(cls, tai, low_precision=False):
        """
        Initialize from a TAI date and time (using `datetime.datetime`).
        
        Parameters
        ----------
        calendar_date : `datetime.datetime` object or a ISO8601 conforming string or list of either
        
        Examples
        --------
        
        >>> from astropy import time
        >>> import datetime
        >>> mytime = time.Time.from_tai(datetime.datetime(1546, 12, 14, 12, 0, 0))
        >>> mytime.to_jd()
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
            time_to_zeropoint = (tai - calendar_zeropoint)
            return cls(seconds=time_to_zeropoint.seconds,
                       fraction_of_seconds=time_to_zeropoint.microseconds/1e6,
                       isarray=False,
                       low_precision=low_precision)
            
        elif np.iterable(tai):
            tais = np.array(tai)
            if isinstance(tais[0], basestring):
                NotImplementedError('ISO8601 parsing not available yet')
            
            elif isinstance(tais[0], datetime.datetime):
                vector_total_seconds = np.vectorize(get_total_seconds, otypes=np.int64)
                seconds = vector_total_seconds(tais)
                if low_precision:
                    cls(seconds, None, isarray=True, low_precision=True)
                else:
                    vector_fraction_of_seconds = np.vectorize(get_fraction_of_seconds, otypes=np.float64)
                    fraction_of_seconds = vector_fraction_of_seconds(tais)
                    return cls(seconds, fraction_of_seconds, isarray=True, low_precision=False)
                
            
    
    @classmethod
    def from_jd(cls, jd_time):
        """
        Instantiate a Time-object with Julian Date (linked to TAI)

        Parameters
        ----------
        jd_time : float or list of floats
        A Julian date 
        :param jd_time:
            A float object 
        """
        
        if np.iterable(jd_time):
            jd_times = np.array(jd_time)
            seconds = jd * 86400 - jd_seconds
            int_seconds = np.int64(seconds)
            if low_precision:
                return cls(int_seconds, fraction_of_seconds=None, isarray=True, low_precision=True)
            else:
                fraction_of_seconds = np.float64(seconds - int_seconds)
                return cls(int_seconds, fraction_of_seconds, isarray=True, low_precision=False)
        else:
            seconds = jd * 86400 - jd_seconds
            int_seconds = np.int64(seconds)
            if low_precision:
                return cls(int_seconds, fraction_of_seconds=None, isarray=False, low_precision=True)
            else:
                fraction_of_seconds = np.float64(seconds - int_seconds)
                return cls(int_seconds, fraction_of_seconds, isarray=False, low_precision=False)
        
    
    @classmethod    
    def from_mjd(cls, mjd_time):
        """
        Instantiate an AstroTime-object with Modified Julian Date

        Parameters
        ----------
        mjd_time : float
            A Modified Julian date 
        
        """
        return cls.from_jd(mjd_time + 2400000.5)    
    
    def __init__(self, seconds, fraction_of_seconds, isarray=False, low_precision=False):
        self.seconds = np.int64(seconds)
        if low_precision:
            self.fraction_of_seconds = None
        else:
            self.fraction_of_seconds = np.float64(fraction_of_seconds)
        self.isarray = isarray
        
    
    def get_jd(self):
        """return the date as JD in a float64"""
        if self.fraction_of_seconds is None:
            return (jd_zeropoint_seconds + self.seconds) / 86400.
        else:
            return (jd_zeropoint_seconds + self.seconds + self.fraction_of_seconds) / 86400. 
        
    jd = property(get_jd)
    
    def get_tai(self):
        """
        returns the TAI date in a `datetime.datetime` object
        
        References
        ----------
        http://www.usno.navy.mil/USNO/astronomical-applications/astronomical-information-center/julian-date-form
        """
        if self.isarray:
            vector_convert_seconds_to_timedelta = np.vectorize(convert_seconds_to_timedelta)
            
            if self.fraction_of_seconds is None:
                td = vector_convert_seconds_to_timedelta(self.seconds, 0)
            else:
                td = vector_convert_seconds_to_timedelta(self.seconds, self.fraction_of_seconds)
        
        else:
            if self.fraction_of_seconds is None:
                td = datetime.timedelta(seconds=self.seconds)
            else:
                td = datetime.timedelta(seconds=self.seconds, microseconds = self.fraction_of_seconds*1e6)
        
        return calendar_zeropoint + td
 
    tai = property(get_tai)