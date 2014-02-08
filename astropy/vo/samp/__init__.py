# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This subpackage provides classes to communicate with other applications via the
`Simple Application Messaging Protocal (SAMP) <www.ivoa.net/samp>`_.

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
