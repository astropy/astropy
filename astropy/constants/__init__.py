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

# for updating the constants module docstring
_lines = [
    'The following constants are available:\n',
    '========== ============== ================ =========================',
    '   Name        Value            Unit       Description',
    '========== ============== ================ =========================',
]

for _nm, _c in itertools.chain(sorted(vars(si).items()),
                               sorted(vars(cgs).items())):
    if isinstance(_c, Constant) and _c.abbrev not in locals():
        locals()[_c.abbrev] = _c.__class__(_c.abbrev, _c.name, _c.value,
                                           _c._unit, _c.uncertainty,
                                           _c.reference)

        _lines.append('{0:^10} {1:^14.9g} {2:^16} {3}'.format(
            _c.abbrev, _c.value, _c._unit, _c.name))

_lines.append(_lines[1])

__doc__ += '\n'.join(_lines)

del _lines, _nm, _c
