# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains initial mass function distributions with the standard
default parameters.
"""
try:
    from scipy.stats import rv_continuous
except ImportError:
    from warnings import warn
    warn('scipy must be present to use astropy.stats.imfs')
    rv_continuous = object

import numpy as np

class SalpeterGen(rv_continuous):
    """
    Power law distribution:
    p(x) = C * x**(-p)
    where
    C = (p-1) / x_0**(1-p)

    The default parameters are the Salpeter minimum mass of a=0.03 Msun and
    b=inf.  To change the limits of the mass function, do:
    imf_limited = SalpeterGen(a=0.5,b=120)
    """
    def _pdf(self, x, p=2.35):
        return x**(-p) * (p-1) / self.a**(1-p)

    def _cdf(self, x, p=2.35):
        return 1-(x/self.a)**(1-p) 

    def __init__(self, a=0.03, name='Salpeter', **kwargs):
        super(SalpeterGen,self).__init__(a=a, name=name, **kwargs)

Salpeter = SalpeterGen()

class KroupaGen(rv_continuous):
    """
    Power law distribution:
    p(x) = b*x**(-0.3) | x < 0.08
           c*x**(-1.3) | 0.08 < x < 0.5
           d*x**(-2.3) | 0.5 < x
    """
    def _pdf(self, m, p1=0.3, p2=1.3, p3=2.3, break1=0.08, break2=0.5):
        """
        """

        m = np.array(m)

        binv = ((break1**(-(p1-1)) - self.a**(-(p1-1)))/(1-p1) +
                (break2**(-(p2-1)) - break1**(-(p2-1))) * (break1**(p2-p1))/(1-p2) +
                (- break2**(-(p3-1))) * (break1**(p2-p1)) * (break2**(p3-p2))/(1-p3))
        b = 1./binv
        c = b * break1**(p2-p1)
        d = c * break2**(p3-p2)

        zeta = (b*(m**(-(p1))) * (m<break1) +
                c*(m**(-(p2))) * (m>=break1) * (m<break2) +
                d*(m**(-(p3))) * (m>=break2))
        return zeta

    def __init__(self, a=0.03, name='Kroupa', **kwargs):
        super(KroupaGen,self).__init__(a=a, name=name, **kwargs)

Kroupa = KroupaGen()
