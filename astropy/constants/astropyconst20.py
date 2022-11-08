# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astronomical and physics constants for Astropy v2.0.
See :mod:`astropy.constants` for a complete listing of constants defined
in Astropy.
"""
import warnings

from astropy.utils import find_current_module

from . import codata2014, iau2015
from . import utils as _utils

codata = codata2014
iaudata = iau2015

_utils._set_c(codata, iaudata, find_current_module())

# Overwrite the following for consistency.
# https://github.com/astropy/astropy/issues/8920
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", "Constant .*already has a definition")

    # Solar mass (derived from mass parameter and gravitational constant)
    M_sun = iau2015.IAU2015(
        "M_sun",
        "Solar mass",
        iau2015.GM_sun.value / codata2014.G.value,
        "kg",
        (
            (codata2014.G.uncertainty / codata2014.G.value)
            * (iau2015.GM_sun.value / codata2014.G.value)
        ),
        f"IAU 2015 Resolution B 3 + {codata2014.G.reference}",
        system="si",
    )

    # Jupiter mass (derived from mass parameter and gravitational constant)
    M_jup = iau2015.IAU2015(
        "M_jup",
        "Jupiter mass",
        iau2015.GM_jup.value / codata2014.G.value,
        "kg",
        (
            (codata2014.G.uncertainty / codata2014.G.value)
            * (iau2015.GM_jup.value / codata2014.G.value)
        ),
        f"IAU 2015 Resolution B 3 + {codata2014.G.reference}",
        system="si",
    )

    # Earth mass (derived from mass parameter and gravitational constant)
    M_earth = iau2015.IAU2015(
        "M_earth",
        "Earth mass",
        iau2015.GM_earth.value / codata2014.G.value,
        "kg",
        (
            (codata2014.G.uncertainty / codata2014.G.value)
            * (iau2015.GM_earth.value / codata2014.G.value)
        ),
        f"IAU 2015 Resolution B 3 + {codata2014.G.reference}",
        system="si",
    )

# Clean up namespace
del warnings
del find_current_module
del _utils
