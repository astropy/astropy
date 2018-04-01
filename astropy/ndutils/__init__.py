# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
The `astropy.nddata` subpackage provides the `~astropy.nddata.NDData`
class and related tools to manage n-dimensional array-based data (e.g.
CCD images, IFU Data, grid-based simulation data, ...). This is more than
just `numpy.ndarray` objects, because it provides metadata that cannot
be easily provided by a single array.
"""

from .nddata import *
from .nddata_base import *
from .nddata_withmixins import *
from .nduncertainty import *
from .flag_collection import *

from .decorators import *

from .mixins.ndarithmetic import *
from .mixins.ndslicing import *
from .mixins.ndio import *

from .compat import *
from .utils import *
from .ccddata import *

from .. import config as _config


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.nddata`.
    """

    warn_unsupported_correlated = _config.ConfigItem(
        True,
        'Whether to issue a warning if `~astropy.nddata.NDData` arithmetic '
        'is performed with uncertainties and the uncertainties do not '
        'support the propagation of correlated uncertainties.'
    )

    warn_setting_unit_directly = _config.ConfigItem(
        True,
        'Whether to issue a warning when the `~astropy.nddata.NDData` unit '
        'attribute is changed from a non-``None`` value to another value '
        'that data values/uncertainties are not scaled with the unit change.'
    )


conf = Conf()
