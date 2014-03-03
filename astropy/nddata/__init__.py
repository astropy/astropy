# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
The `nddata` subpackage provides the `~astropy.nddata.NDData`
class and related tools to manage n-dimensional array-based data (e.g.
CCD images, IFU Data, grid-based simulation data, ...). This is more than
just `numpy.ndarray` objects, because it provides metadata that cannot
be easily provided by a single array.
"""

from .nddata import *
from .nduncertainty import *
from .flag_collection import *

from .. import config as _config


class _Conf(_config.ConfigNamespace):
    warn_unsupported_correlated = _config.ConfigItem(
        True,
        'Whether to issue a warning if NDData arithmetic is performed with '
        'uncertainties and the uncertainties do not support the propagation '
        'of correlated uncertainties.'
    )
conf = _Conf()
