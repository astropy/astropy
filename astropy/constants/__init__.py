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

from .._constants.definition import ConstantDefinition

from .constant import Constant, EMConstant

# Define actual Quantity-based Constants. Need to use _c in the loop instead
# of c to avoid overwriting the speed of light constant.

# Import SI constants
from .._constants import si as _si
for nm, val in sorted(vars(_si).iteritems()):
    if isinstance(val, ConstantDefinition):
        if val.em:
            _c = EMConstant(val.value, val.units, val.uncertainty, val.name,
                            val.reference)
        else:
            _c = Constant(val.value, val.units, val.uncertainty, val.name,
                          val.reference)
        locals()[nm] = _c

del _si, nm, val, _c

# Import cgs constants that don't exist in S.I.
from .._constants import cgs as _cgs
for nm, val in sorted(vars(_cgs).iteritems()):
    if isinstance(val, ConstantDefinition):
        if val.em:
            # EM constants have _gauss, _esu, etc. appended
            nm, system = nm.rsplit('_', 1)
            _c = EMConstant(val.value, val.units, val.uncertainty, val.name,
                            val.reference)
            if nm in locals():
                _c_si = locals()[nm]
                setattr(_c_si, system, _c)
            else:
                locals()[nm] = _c
        else:
            _c = Constant(val.value, val.units, val.uncertainty, val.name,
                          val.reference)
            if nm in locals():
                raise ValueError(
                    "Constant {0} has already been imported".format(nm))
            locals()[nm] = _c

del _cgs, nm, val, _c

# update the si cand cgs module doctrings.
__doc__ += """
The following constants are available:

========== ============== ================ =========================
   Name        Value            Unit       Description
========== ============== ================ =========================
"""
for nm, val in sorted(locals().items()):
    if isinstance(val, Constant):
        __doc__ += '{0:^10} {1:^14.9g} {2:^16} {3}\n'.format(
            nm, val.value, val.unit, val.name)

__doc__ += """\
========== ============== ================ =========================
"""

del nm, val
