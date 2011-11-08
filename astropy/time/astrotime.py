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
        
            * None
                Takes the current time when the object is created.
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
    
    def __init__(self,time,secprecision=1e-9,jd0=2451545.):
        
        if isinstance(time,AstroTime):
            self._jd0 = time._jd0
            self._secprec = time._secprec
            self._val = time._val
        else:
            self._jd0 = jd0
            self._secprec = secprecision
            
            if time is None:
                self._set_from_jd(self._datetime_to_jd(None))
            elif isinstance(time,datetime):
                self._set_from_jd(self._datetime_to_jd(time))
            elif isinstance(time,basestring):
                self._set_from_jd(self._epoch_to_jd(float(time[1:]), time[0]))
            else:
                self._val = long(time)
    
    @classmethod
    def from_jd(cls,jd,secprecision=1e-9,jd0=2451545.):
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
    def from_mjd(cls,mjd,secprecision=1e-9,jd0=2451545.):
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
            
    def _datetime_to_jd(self, caltime,tz=None,gregorian=True):
        """Converts a datetime object into a julian date 
        
        `caltime` is a datetime.datetime or None to do the current instant.
        `tz` can be used to override the timezone
        `gregorian` determines which type of date scheme to use. If True, the
        will be interpreted as in the Gregorian calendar. Otherwise, it will be
        Julian. If None, it will be assumed to switch over on October 4/15,
        1582.
        """
        from datetime import datetime,date,tzinfo
        
        if caltime is None:
            from dateutil.tz import tzlocal
            datetimes = [datetime.now(tzlocal())]
        elif isinstance(caltime,datetime) or isinstance(caltime,date):
            datetimes = [caltime]
        elif all([isinstance(ct,datetime) or isinstance(ct,date) for ct in caltime]):
            datetimes = caltime
        else:
            datetimes = None
            caltime = list(caltime)
            if not (3 <= len(caltime) < 8):
                raise ValueError('caltime input sequence is invalid size')
            while len(caltime) < 7:
                if len(caltime) == 3:
                    #make hours 12
                    caltime.append(12*np.ones_like(caltime[-1]))
                else:
                    caltime.append(np.zeros_like(caltime[-1]))
            yr,month,day,hr,min,sec,msec = caltime
            
        #if input objects are datetime objects, generate arrays
        if datetimes is not None:
            yr,month,day,hr,min,sec,msec = [],[],[],[],[],[],[]
            for dt in datetimes:
                if not hasattr(dt,'hour'):
                    dt = datetime(dt.year,dt.month,dt.day,12)
                
                if tz is None:
                    off = dt.utcoffset()
                    if off is not None:
                        dt = dt - off
                    
                yr.append(dt.year)
                month.append(dt.month)
                day.append(dt.day)
                hr.append(dt.hour)
                min.append(dt.minute)
                sec.append(dt.second)
                msec.append(dt.microsecond)
                    
                    
        
        yr = np.array(yr,dtype='int64',copy=False).ravel()
        month = np.array(month,dtype='int64',copy=False).ravel()
        day = np.array(day,dtype='int64',copy=False).ravel()
        hr = np.array(hr,dtype=float,copy=False).ravel()
        min = np.array(min,dtype=float,copy=False).ravel()
        sec = np.array(sec,dtype=float,copy=False).ravel()
        msec = np.array(msec,dtype=float,copy=False).ravel()
        
        #do tz conversion if tz is provided  
        if isinstance(tz,basestring) or isinstance(tz,tzinfo):
            if isinstance(tz,basestring):
                from dateutil import tz
                tzi = tz.gettz(tz)
            else:
                tzi = tz
            
            utcoffset = []
            for t in zip(yr,month,day,hr,min,sec,msec):
                #microsecond from float component of seconds
                
                dt = datetime(*[int(ti) for ti in t],**dict(tzinfo=tzi))
                utcdt = dt.utcoffset()
                if utcdt is None:
                    utcoffset.append(0)
                else:
                    utcoffset.append(utcdt.days*24 + (utcdt.seconds + utcdt.microseconds*1e-6)/3600)
        else:
            utcoffset = tz
                
    #    ly = ((month-14)/12).astype(int) #In leap years, -1 for Jan, Feb, else 0
    #    jdn = day - 32075l + 1461l*(yr+4800l+ly)//4
        
    #    jdn += 367l*(month - 2-ly*12)//12 - 3*((yr+4900l+ly)//100)//4
        
    #    res = jdn + (hr/24.0) + min/1440.0 + sec/86400.0 - 0.5

        #this algorithm from meeus 2ed
        m3 = month < 3
        yr[m3] -= 1
        month[m3] += 12
            
        cen = yr//100
        
        if gregorian is None:
            gregorian = (1582,10,4)
        if gregorian is True:
            gregoffset = 2 - cen + cen//4
        elif gregorian is False:
            gregoffset = 0
        else:
            gregoffset = 2 - cen + cen//4
            gmask = (yr>gregorian[0])&(month>gregorian[1])&(day>gregorian[2])
            gregoffset[~gmask] = 0
        
            
        jdn = (365.25*(yr+4716)).astype(int) + \
              (30.6001*(month + 1)).astype(int) + \
                   day + gregoffset - 1524.5
        res = jdn + hr/24.0 + min/1440.0 + sec/86400.0
        
        if np.any(utcoffset):
            res -= np.array(utcoffset)/24.0
        
        return res
            
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
    
    @property
    def datetime(self):
        """
        A `datetime.datetime` object representing this time, rounded to the 
        nearest second. Timezone is UTC.
        """        
        jd = self.jd
        
        rounding = 1000000
        #If non-0, Performs a fix for floating-point errors. It specifies the
        #number of milliseconds by which to round the result to the nearest
        #second. If 1000000 (one second), no milliseconds are recorded. If
        #larger, a ValueError is raised.
        gregorian = None
        #If True, the input will be interpreted as in the Gregorian calendar.
        #Otherwise, it will be Julian. If None, it will be assumed to switch over
        #on October 4/15, 1582.
    
        
        if rounding > 1000000:
            raise ValueError('rounding cannot exceed a second')
        elif rounding <= 0:
            jd += .5 
        else:
            rounding = int(rounding)
            roundingfrac = rounding/86400000000
            jd += .5 + roundingfrac 
            
        z = np.floor(jd).astype(int) 
        dec = jd - z #fractional piece
        
        #fix slight floating-point errors if they hapepn TOOD:check
        dgtr1 = dec>=1.0
        dec[dgtr1] -= 1.0
        z[dgtr1] += 1
        
        
        if gregorian is None:
            gregorian = 2299161
            
        if gregorian is True:
            alpha = ((z-1867216.25)/36524.25).astype(int)
            z += 1 + alpha - alpha//4
        elif gregorian is False:
            pass
        else:
            gmask = z >= gregorian
            alpha = ((z[gmask]-1867216.25)/36524.25).astype(int)
            z[gmask] += 1 + alpha - alpha//4
        
        b = z + 1524
        c = ((b-122.1)/365.25).astype(int)
        d = (365.25*c).astype(int)
        e = ((b-d)/30.6001).astype(int)
        
        day = b - d - (30.6001*e).astype(int)
        
        mmask = e<14
        month = e
        month[mmask] -= 1
        month[~mmask] -= 13
        year = c
        year[month>2] -= 4716
        year[month<=2] -= 4715
        
        
        if rounding == 1000000:
            secdec = dec*86400
            sec = secdec.astype(int)
            min = sec//60
            sec -= 60*min
            hr = min//60
            min -= 60*hr
            #sec[sec==secdec] -= 1
            msec = None
        else:
            msec = (dec*86400000000.).astype('int64') 
            if rounding > 0:
                div = (msec//1000000)*1000000
                toround = (msec - div)<(2*rounding)
                msec[toround] = div + rounding
                msec  -= rounding

            sec = msec//1000000
            msec -= 1000000*sec
            
            min = sec//60
            sec -= 60*min
            hr = min//60
            min -= 60*hr
            
        
        if msec is None:
            return datetime.datetime(year,month,day,hr%24,min%60,sec%60)
        else:
            return datetime.datetime(year,month,day,hr%24,min%60,sec%60,msec%1000000)
        
        