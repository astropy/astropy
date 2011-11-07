from datetime import datetime
import numpy as np


    
    
    


class AstroTime(object):
    """Class to store a time variable.
    The internal format uses the January 1, 4713 BC Greenwich noon as the zeropoint (Julian Day)
    and stores the date as days to this zeropoint in ~numpy.int64 and the fraction of a day in ~numpy.float64."""
    
    @classmethod
    def from_jd(cls, jd_time):
        """Instantiate an AstroTime-object with Julian Date"""
        days = np.int64(jd_time)
        fraction_of_day = np.fmod(jd_time, 1.0)
        return cls(days, fraction_of_day)
        
    @classmethod
    def from_gregorian_date(cls, gregorian_datetime):
        #def calendar_to_jd(caltime,tz=None,gregorian=True,mjd=False):
        """
        Convert a calendar date and time to julian date.
        
        :param gregorian_datetime: 
            A :class:`datetime.datetime` object                    
        
        
        
        
        **Examples**
        
        >>> from astropy import astrotime
        >>> astrotime.AstroTime.from_jd(2455872.3907755)
        """
        
        # Reference: http://asa.usno.navy.mil/SecM/Glossary.html
        # Reference: http://en.wikipedia.org/wiki/Julian_day#Converting_Gregorian_calendar_date_to_Julian_Day_Number
        
        a = (14 - gregorian_datetime.month) / 12
        y = gregorian_datetime.year + 4800 - a
        m = gregorian_datetime.month + 12*a - 3
        jdn = gregorian_datetime.day + ((153*m + 2) / 5) + 365*y + y/4 - y/100 + y/400 - 32045
        jd = jdn + (gregorian_datetime.hour - 12) / 24. + \
            gregorian_datetime.minute / 1440. + \
            gregorian_datetime.second / 86400. + \
            gregorian_datetime.microsecond / (86400.*1e6) 
        
        return cls.from_jd(jd)
        
    
    def __init__(self, days, fraction_of_day):
        """Initialize an AstroTime-object with the days and fraction of days from the January 1, 4713 BC Greenwich noon
        
        :param days:
        :param fraction_days:
        """
        self.days = np.int64(days)
        self.fraction_of_day = np.float64(fraction_of_day)
    
    
    def to_jd(self):
        """return the date as JD in a float64"""
        return np.float64(self.days+self.fraction_of_day)
    
    def to_gregorian_date(self):
        """
        
        **Examples**
        
        >>> jd_to_calendar(2451545)
        datetime.datetime(2000, 1, 1, 12, 0, tzinfo=tzutc())
        >>> jd_to_calendar(2305812.5)
        datetime.datetime(1600, 12, 31, 0, 0, tzinfo=tzutc())
        >>> jd_to_calendar([2415020.5,2305447.5],output='array')
        array([[1900,    1,    1,    0,    0,    0,    0],
               [1600,    1,    1,    0,    0,    0,    0]])
        >>> jd_to_calendar(0.0,output='fracarray')
        array([[ -4.71200000e+03,   1.00000000e+00,   1.50000000e+00]])
            
        """
        # Reference: http://www.usno.navy.mil/USNO/astronomical-applications/astronomical-information-center/julian-date-form
        
        jd_fraction_of_days = self.fraction_of_day + 0.5    
        jd_days = np.int64(self.days + jd_fraction_of_days // 1.0)
        jd_fraction_of_days = np.fmod(jd_fraction_of_days, 1.0)
            
        L = jd_days + 68569
        N = 4 * L / 146097
        L = L - (146097 * N + 3) / 4
        I = 4000 * (L+1) / 1461001
        L = L - 1461 * I / 4 + 31
        J = 80 * L / 2447
        days = L - 2447 * J / 80
        L = J / 11
        months = J + 2-12*L
        years = 100 * (N - 49) + I + L

        
        
        fraction = 24*jd_fraction_of_days
        hours = np.int64(fraction)
        fraction = np.fmod(fraction, 1.0)*60
        minutes = np.int64(fraction)
        fraction = np.fmod(fraction, 1.0)*60
        seconds = np.int64(fraction)
        fraction = np.fmod(fraction, 1.0)*1e6
        microseconds = np.int64(fraction)
        
        return datetime(years, months, days, hours, minutes, seconds, microseconds)