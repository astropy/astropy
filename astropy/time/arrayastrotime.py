# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division
from .astrotime import AstroTime
import numpy as np

__all__ = ['ArrayAstroTime']

class ArrayAstroTime(AstroTime):
    """ An object that extends `~astropy.time.astrotime.AstroTime`  to support 
    ndarrays of times.
    
    This class is not as flexible as `~astropy.time.astrotime.AstroTime` in that
    it does not allow arbitrary precision - it is fixed to 
    
    Parameters
    ----------
    days : int64 array
        Julian Days; days from January 1, 4713 BC noon UTC.
    fraction_of_day : float64 array
        Fraction of day - must match the shape of `days`
    """
    def __init__(self, days, fraction_of_day):
        days = np.array(days,copy=False,dtype=np.int64)
        frac = np.array(fraction_of_day,copy=False,dtype=np.float64)
        
        if days.shape != frac.shape:
            if days.shape == tuple(): 
                #scalar
                days = days *np.ones_like(frac).astype(np.int64)
            elif frac.shape == tuple():
                #scalar
                frac = frac * np.ones_like(dats).astype(np.float64)
            else:
                msg = 'days and fraction_of_day do not have matching shapes:{0} and {1}'
                raise ValueError(msg.format(days.shape,frac.shape))
        
        self._day = days
        self._frac = fraction_of_day
        
        #days+frac=0 -> JD of 0
        self._jd0 = 0 
        
        #effective precision set by float64 -  if the fraction is .9xxx, there 
        #can only be as many digits as a float64 accepts.
        #self._secprec = 1.1102230246251578e-16 #10**(-53*math.log10(2))
        # but need to set it to a day to get the superclass behavior correct
        self._secprec = 86400
        
    def _set_longval(self,val):
        self._day = day = np.floor(val).astype(np.int64)
        self._frac = val - day
        
    def _get_longval(self):
        return self._day  + self._frac
    longval = property(_get_longval,_set_longval)