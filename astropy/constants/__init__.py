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

from . import cgs
from . import si
from .constant import Constant


# Update the docstring to include a list of units from the si
# module. The rows with lots of '=' signs are to tell Sphinx to
# display a table in the documentation.

__doc__ += """
The following constants are defined in `~astropy.constants.cgs` and
`~astropy.constants.si`. The `si` and `cgs` docstrings list the units
and values in each system.

========== ==============================
"""

for nm, val in sorted(si.__dict__.items()):
    if isinstance(val, Constant):
        __doc__ += '{0:^10} {1}\n'.format(nm, val.name)

__doc__ += """\
========== ==============================
"""

# update the si cand cgs module doctrings.
for module in si, cgs:
    module.__doc__ += """
========== ============== ================ =========================
   Name        Value            Unit       Description
========== ============== ================ =========================
"""
    for nm, val in sorted(module.__dict__.items()):
        if isinstance(val, Constant):
            module.__doc__ += '{0:^10} {1:^14.9g} {2:^16} {3}\n'.format(
                nm, val.real, val._units, val.name)

    module.__doc__ += """\
========== ============== ================ =========================
"""

del nm, val

from .._constants.definition import ConstantDefinition

# Define actual Quantity-based Constants
for nm, val in sorted(si.__dict__.items()):
    if isinstance(val, ConstantDefinition):
        c = Constant(val.value, val.units, val.error, val.name, val.origin)
        locals()[nm] = c