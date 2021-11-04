# Licensed under a 3-clause BSD style license - see LICENSE.rst
import warnings

from erfa import core, helpers, ufunc  # noqa
from erfa.core import *  # noqa
from erfa.helpers import leap_seconds  # noqa
from erfa.ufunc import (dt_dmsf, dt_eraASTROM, dt_eraLDBODY,  # noqa
                        dt_eraLEAPSECOND, dt_hmsf, dt_pv, dt_sign, dt_type, dt_ymdf)

from astropy.utils.exceptions import AstropyDeprecationWarning

warnings.warn('The private astropy._erfa module has been made into its '
              'own package, pyerfa, which is a dependency of '
              'astropy and can be imported directly using "import erfa"',
              AstropyDeprecationWarning)
