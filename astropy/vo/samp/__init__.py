# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This subpackage provides classes to communicate with other applications via the
`Simple Application Messaging Protocal (SAMP)
<http://www.ivoa.net/documents/SAMP/>`_.

Before integration into Astropy it was known as
`SAMPy <https://pypi.python.org/pypi/sampy/>`_, and was developed by Luigi Paioro
(INAF - Istituto Nazionale di Astrofisica).
"""

from .constants import *
from .errors import *
from .utils import *
from .hub import *
from .client import *
from .integrated_client import *
from .hub_proxy import *


from ... import config as _config


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.vo.samp`.
    """

    use_internet = _config.ConfigItem(
        True,
        "Whether to allow `astropy.vo.samp` to use "
        "the internet, if available.",
        aliases=['astropy.vo.samp.utils.use_internet'])

    n_retries = _config.ConfigItem(10,
        "How many times to retry communications when they fail")

conf = Conf()
