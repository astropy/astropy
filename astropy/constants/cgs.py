# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in cgs units.  See :mod:`astropy.constants`
for a complete listing of constants defined in Astropy.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import itertools

from .constant import Constant
from . import codata2014, iau2015

for _nm, _c in itertools.chain(sorted(vars(codata2014).items()),
                               sorted(vars(iau2015).items())):
    if (isinstance(_c, Constant) and _c.abbrev not in locals()
         and _c.system in ['esu', 'gauss', 'emu']):
        locals()[_c.abbrev] = _c
