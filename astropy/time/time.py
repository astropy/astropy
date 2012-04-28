# Licensed under a 3-clause BSD style license - see LICENSE.rst


import datetime
import decimal
import numpy as np

#Defining zeropoints where TCB, TCG and TT are linked to ephemeris time
#wikipedia http://en.wikipedia.org/wiki/International_Atomic_Time

calendar_zeropoint = datetime.datetime(1977, 1, 1, 0, 0, 32, 184000)
jd_zeropoint = decimal.Decimal('2443144.5003725')

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
    def from_jd(cls, jd_time):
        """
        Instantiate a Time-object with Julian Date

        Parameters
        ----------
        jd_time : float
        A Julian date 
        :param jd_time:
            A float object
        """
        
        return cls((decimal.Decimal(jd_time) - decimal.Decimal(jd_zeropoint)) * decimal.Decimal(86400))
    
    @classmethod    
    def from_mjd(cls, mjd_time):
        """
        Instantiate an AstroTime-object with Modified Julian Date

        Parameters
        ----------
        mjd_time : float
            A Modified Julian date 
        
        """
        return cls.from_jd(decimal.Decimal(mjd_time) + decimal.Decimal(2400000.5))
        
    @classmethod
    def from_utc(cls, calendar_date):
        """
        Initialize from a UTC date and time (using `datetime.datetime`).
        
        Parameters
        ----------
        calendar_date : `datetime.datetime` object
        
        Examples
        --------
        
        >>> from astropy import astrotime
        >>> import datetime
        >>> mytime = time.Time.from_utc(datetime.datetime(1546, 12, 14, 12, 0, 0))
        >>> mytime.to_jd()
        2286072.0
        
        References
        ----------
        http://asa.usno.navy.mil/SecM/Glossary.html
        http://en.wikipedia.org/wiki/Julian_day#Converting_Gregorian_calendar_date_to_Julian_Day_Number
        """
        
        return cls((calendar_date - calendar_zeropoint).total_seconds())
        
    
    def __init__(self, seconds):
        self.seconds = decimal.Decimal(seconds)    
    
    def to_jd(self):
        """return the date as JD in a float64"""
        
        return decimal.Decimal(jd_zeropoint) + (self.seconds / decimal.Decimal(86400))
        
    jd = property(to_jd)
    
    def to_utc(self):
        """
        returns the UTC date in a `datetime.datetime` object
        
        References
        ----------
        http://www.usno.navy.mil/USNO/astronomical-applications/astronomical-information-center/julian-date-form
        """
        
        days = (self.seconds / decimal.Decimal(86400))
        return calendar_zeropoint + datetime.timedelta(np.float64(days))
 
    utc = property(to_utc)