# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
.. warning::
    ``astropy.samp`` was deprecated in version 9.0 and will be removed in a future version;
    please use ``pyvo.astropy_samp`` instead.

This subpackage provides classes to communicate with other applications via the
`Simple Application Messaging Protocol (SAMP)
<http://www.ivoa.net/documents/SAMP/>`_.

Before integration into Astropy it was known as
`SAMPy <https://pypi.org/project/sampy/>`_, and was developed by Luigi Paioro
(INAF - Istituto Nazionale di Astrofisica).
"""

import warnings

from astropy import config as _config
from astropy.utils.exceptions import AstropyDeprecationWarning

from .client import *
from .constants import *
from .errors import *
from .hub import *
from .hub_proxy import *
from .integrated_client import *
from .utils import *

warnings.warn(
    "astropy.samp was deprecated in version 9.0 "
    "and will be removed in a future version; "
    "please use pyvo.astropy_samp instead.",
    AstropyDeprecationWarning,
)

# Clean up namespace
del warnings
del AstropyDeprecationWarning


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.samp`.
    """

    use_internet = _config.ConfigItem(
        True,
        "Whether to allow `astropy.samp` to use the internet, if available.",
        aliases=["astropy.samp.utils.use_internet"],
    )

    n_retries = _config.ConfigItem(
        10, "How many times to retry communications when they fail"
    )


conf = Conf()
