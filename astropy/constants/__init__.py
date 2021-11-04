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
import warnings

from astropy.utils import find_current_module

# Hack to make circular imports with units work
# isort: split
from astropy import units

del units

from . import cgs  # noqa
from . import si  # noqa
from . import utils as _utils  # noqa
from .config import codata, iaudata  # noqa
from .constant import Constant, EMConstant  # noqa

# for updating the constants module docstring
_lines = [
    'The following constants are available:\n',
    '========== ============== ================ =========================',
    '   Name        Value            Unit       Description',
    '========== ============== ================ =========================',
]

# Catch warnings about "already has a definition in the None system"
with warnings.catch_warnings():
    warnings.filterwarnings('ignore', 'Constant .*already has a definition')
    _utils._set_c(codata, iaudata, find_current_module(),
                  not_in_module_only=True, doclines=_lines, set_class=True)

_lines.append(_lines[1])

if __doc__ is not None:
    __doc__ += '\n'.join(_lines)


# Clean up namespace
del find_current_module
del warnings
del _utils
del _lines
