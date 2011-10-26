"""
Contains astronomical and physical constants for use in Astropy or other
places.

The package contains a `~astropy.constants.cgs` and `~astropy.constants.si`
module that define constants in CGS and SI units, respectively.  A typical use
case might be::

    from astropy.constants.cgs import c

    ... define the mass of something you want the rest energy of as m ...
    E = m*c**2

"""

from . import cgs
from . import si
from .constant import Constant


# Update the docstring to include a list of units from the si module
__doc__ += """
The following constants are defined in `~astropy.constants.cgs` and
`~astropy.constants.si` .

"""

for nm, val in si.__dict__.iteritems():
    if isinstance(val, Constant):
        __doc__ += '    * ' + nm + '\n        ' + val.name + '\n'
del nm, val
__doc__ += '\n'
