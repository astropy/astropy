# Licensed under a 3-clause BSD style license - see LICENSE.rst
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

import itertools

# Hack to make circular imports with units work
try:
    from .. import units
    del units
except ImportError:
    pass

from .constant import Constant, EMConstant
from . import si
from . import cgs

for const in itertools.chain(vars(si).values(), vars(cgs).values()):
    if isinstance(const, Constant) and const.abbrev not in locals():
        locals()[const.abbrev] = const.__class__(const.abbrev, const.name,
                                                 const.value, const._unit,
                                                 const.uncertainty,
                                                 const.reference)

# update the constants module docstring
_lines = [
    'The following constants are available:\n',
    '========== ============== ================ ========================='
    '   Name        Value            Unit       Description'
    '========== ============== ================ ========================='
]


#for nm, val in sorted(locals().items()):
#    if not isinstance(val, Constant):
#        continue
#    _lines.append('{0:^10} {1:^14.9g} {2:^16} {3}\n'.format(
#                  nm + '.' + val.system, val.value, val.unit, val.name))

#_lines.append(_lines[1])

#__doc__ += '\n'.join(_lines)

del _lines, const
