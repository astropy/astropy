# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains astronomical and physical constants for use in Astropy or other
places.

A typical use case might be::

    >>> from astropy.constants import c, m_e
    >>> # ... define the mass of something you want the rest energy of as m ...
    >>> m = m_e
    >>> E = m * c**2
    >>> E.to('MeV')  # doctest: +FLOAT_CMP
    <Quantity 0.510998927603161 MeV>

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
from . import codata2014, iau2015

# for updating the constants module docstring
_lines = [
    'The following constants are available:\n',
    '========== ============== ================ =========================',
    '   Name        Value            Unit       Description',
    '========== ============== ================ =========================',
]

for _nm, _c in itertools.chain(sorted(vars(codata2014).items()),
                               sorted(vars(iau2015).items())):
    if isinstance(_c, Constant) and _c.abbrev not in locals():
        locals()[_c.abbrev] = _c.__class__(_c.abbrev, _c.name, _c.value,
                                           _c._unit_string, _c.uncertainty,
                                           _c.reference)

        _lines.append('{0:^10} {1:^14.9g} {2:^16} {3}'.format(
            _c.abbrev, _c.value, _c._unit_string, _c.name))

_lines.append(_lines[1])

if __doc__ is not None:
    __doc__ += '\n'.join(_lines)

del _lines, _nm, _c
