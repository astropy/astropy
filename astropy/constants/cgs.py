# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants in cgs units.  See :mod:`astropy.constants`
for a complete listing of constants defined in Astropy.
"""

from .config import codata, iaudata
from .constant import Constant

for _nm, _c in (*sorted(vars(codata).items()), *sorted(vars(iaudata).items())):
    if (
        isinstance(_c, Constant)
        and _c.abbrev not in locals()
        and _c.system in ["esu", "gauss", "emu"]
    ):
        locals()[_c.abbrev] = _c
