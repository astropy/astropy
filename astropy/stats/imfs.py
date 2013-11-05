# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains initial mass function distributions with the standard
default parameters.
"""
import scipy.stats

class Salpeter_gen(scipy.stats.rv_continuous):
    """
    Power law distribution:
    p(x) = C * x**(-p-1)
    where
    C = x_0**(-p)
    """
    def _pdf(self, x, p=1.35):
        return x**(-p-1) / self.a**p

    def __init__(self, a=0.03, **kwargs):
        super(Salpeter_gen,self).__init__(a=a,**kwargs)
