# Licensed under a 3-clause BSD style license - see LICENSE.rst
import warnings

from erfa import core, ufunc, helpers  # noqa
from erfa.core import *  # noqa
from erfa.ufunc import (dt_eraASTROM, dt_eraLDBODY, dt_eraLEAPSECOND, dt_pv,  # noqa
                        dt_sign, dt_type, dt_ymdf, dt_hmsf, dt_dmsf)
from erfa.helpers import leap_seconds  # noqa

from astropy.utils.exceptions import AstropyDeprecationWarning


warnings.warn('The private astropy._erfa module has been made into its '
              'own package, pyerfa, which is a dependency of '
              'astropy and can be imported directly using "import erfa"',
              AstropyDeprecationWarning)
