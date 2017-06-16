# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants for Astropy v1.3 and earlier.
See :mod:`astropy.constants` for a complete listing of constants
defined in Astropy.
"""


from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import itertools

from .constant import Constant
from . import codata2010, iau2012

for _nm, _c in itertools.chain(sorted(vars(codata2010).items()),
                               sorted(vars(iau2012).items())):
    if (isinstance(_c, Constant) and _c.abbrev not in locals()):
        locals()[_c.abbrev] = _c
