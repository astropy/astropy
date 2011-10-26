"""
Contains astronomical and physical constants for use in Astropy or other places.

The package contains a `~astropy.constants.cgs` and `~astropy.constants.si` 
module that define constants in CGS and SI units, respectively.  A typical use
case might be::

    from astropy.constants.cgs import c
    
    ... define the mass of something you want the rest energy of as m ...
    E = m*c**2

"""

from . import cgs
from . import si