# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division

__all__ = ['AstroTime']

from datetime import datetime

_mjdoffset = 2400000.5 #JD at MJD=0

class AstroTime(object):
    """
    The astrotime class is the base representation of time for astropy.
    
    Parameters
    ----------
    time :
        A representation of the time to use to generate this `AstroTime` object.
        Can be in any of the following forms:
        
            * An `AstroTime` object
            * A `datetime.datetime` object
            * A string of the form 'J2000' or 'B1950'
                Indicates an epoch.
            * An integer or float
                The number of seconds from the `jd0` parameter.
    secprecision : float
        The precision for this time in fractions of a second.
    jd0 : float
        The zero-point for this `AstroTime` object in Julian Dates.
            
    
    """
    
    def __init__(self,time,secprecision=1e-9,jd0=_jd2000):
        
        if isinstance(time,AstroTime):
            self._jd0 = time._jd0
            self._secprec = time._secprec
            self._val = time._val
            
        else:
            self._jd0 = jd0
            self._secprec = secprecision
            
            if isinstance(time,datetime):
                self._set_from_jd(self._datetime_to_jd(time))
            elif isinstance(time,basestring):
                self._set_from_jd(self._epoch_to_jd(float(time[1:]), time[0]))
            else:
                self._val = long(time)
    
    @classmethod
    def from_jd(cls,jd,secprecision=1e-9,jd0=_jd2000):
        """
        Creates a new `AstroTime` object from a given Julian Date.
        
        Parameters
        ----------
        jd : float
            Julian date to use for the new `AstroTime` object.
        secprecision : float
            See `AstroTime` docstrings.
        jd0 : float
            See `AstroTime` docstrings.
            
        Returns
        -------
        astrotime : `AstroTime`
        
        """
        res = cls(0,secprecision,jd0)
        res._val = res._jd_to_val(jd)
        return res
        
    @classmethod
    def from_mjd(cls,mjd,secprecision=1e-9,jd0=_jd2000):
        """
        Creates a new `AstroTime` object from a given Modified Julian Date.
        
        Parameters
        ----------
        mjd : float
            Modified Julian Date to use for the new `AstroTime` object.
        secprecision : float
            See `AstroTime` docstrings.
        jd0 : float
            See `AstroTime` docstrings.
            
        Returns
        -------
        astrotime : `AstroTime`
        
        """
        res = cls(0,secprecision,jd0)
        res._val = res._jd_to_val(mjd + _mjdoffset)
        return res
        
    def _jd_to_val(self,jd):
        return long((jd - self._jd0) * 86400 / self._secprec)
        
    def _epoch_to_jd(self,epoch,epochtype='J'):
        if epochtype=='J':
            return (epoch - 2000) * 365.25 + 2451545.0
        elif epochtype=='B':
            return (epoch - 1900) * 365.242198781 + 2415020.31352
        else:
            raise ValueError('Invalid epoch string - must start with J or B')
            
    def __add__(self,other):
        #bypass if they precision and jd0 match
        if self._secprec == other._secprec and self._jd0 == other._jd0:
            return self.__class__(self._val + other._val,
                                  self._secprec, self._jd0)
        
        # use self's jd0 for the new object
        jdoffset = self._jd0 - other._jd0
        oval = other._val + jdoffset*86400/other._secprec
        
        #use best precision
        secprec = min(self._secprec,other._secprec)
        
        if secprec == self._secprec:
            newval = self._val + oval * other._secprec / secprec
        else:
            newval = self._val * self._secprec / secprec + oval
        
        return self.__class__(newval, secprec, self._jd0)
        
    def __sub__(self,other):
        #bypass if they precision and jd0 match
        if self._secprec == other._secprec and self._jd0 == other._jd0:
            return self.__class__(self._val - other._val,
                                  self._secprec, self._jd0)
        
        # use self's jd0 for the new object
        jdoffset = self._jd0 - other._jd0
        oval = other._val + jdoffset*86400/other._secprec
        
        #use best precision
        secprec = min(self._secprec,other._secprec)
        
        if secprec == self._secprec:
            newval = self._val - oval * other._secprec / secprec
        else:
            newval = self._val * self._secprec / secprec - oval
        
        return self.__class__(newval, secprec, self._jd0)
            
    def __eq__(self,other):
        if self._secprec == other._secprec and self._jd0 == other._jd0:
            return self._val == other._val
        
        return (self - other)._val == 0
            
    @property
    def bit_length(self):
        """
        The number of bits this `AstroTime` object uses to store the time.
        """
        return self._val.bit_length()
        
    @property
    def jd(self):
        """
        Julian Date of this `AstroTime` object.
        """
        return self._val * self._secprec / 86400 + self._jd0
    
    @property
    def mjd(self):
        """
        Modified Julian Date of this `AstroTime` object.
        """
        return self.jd - _mjdoffset
    
    @property
    def jdn(self):
        """
        Julian Day Number of this `AstroTime` object.
        """
        from math import floor
        
        return floor(self.jd)
        
    @property
    def jepoch(self):
        """
        Julian epoch (as a float)
        """
        return 2000.0 + (self.jd - 2451545.0)/365.25
    
    @property
    def jepochstr(self):
        """
        Julian epoch (as a string)
        """
        from math import floor
        
        je = self.jepoch
        
        if floor(je)==je:
            return 'J{0.0f}'.format(je)
        else:
            
            return 'J{0qf}'.format(je)
            
    @property
    def bepoch(self):
        """
        Besselian epoch (as a float)
        """
        return 1900 + (self.jd - 2415020.31352)/365.242198781
    
    @property
    def bepochstr(self):
        """
        Besselian epoch (as a string)
        """
        from math import floor
        
        be = self.bepoch
        
        if floor(be)==be:
            return 'B{0.0f}'.format(be)
        else:
            
            return 'B{0qf}'.format(be)