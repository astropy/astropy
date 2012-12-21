"""
This module provides utility functions for the models package
"""
from __future__ import division, print_function
import numpy as np

class InputParametersException(Exception):
    pass

__all__ = ['pmapdomain', 'comb', 'InputParametersException']
          

def pmapdomain(oldx, domain, window):
    """
    Map domain into window by shifting and scaling
    
    Parameters
    ----------
    oldx: array
          original coordinates
    domain: list or tuple of length 2
          function domain
    window: list or tuple of length 2
          range into which to map the domain 
    """
    domain = np.array(domain, dtype=np.float64)
    window = np.array(window, dtype=np.float64)
    scl = (window[1]-window[0])/(domain[1]-domain[0])
    off = (window[0]*domain[1] - window[1]*domain[0])/(domain[1]-domain[0])
    return off + scl*oldx

def comb(N, k):
    """
    The number of combinations of N things taken k at a time.
    
    Parameters
    ----------
    N : int, array
        Number of things.
    k : int, array
        Number of elements taken.
        
    """
    if (k > N) or (N < 0) or (k < 0):
        return 0L
    val = 1L
    for j in xrange(min(k, N-k)):
        val = (val*(N-j))/(j+1)
    return val
