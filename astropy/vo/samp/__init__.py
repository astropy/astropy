# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The `vo.samp` subpackage provides classes to communicate with other
applications via the `Simple Application Messaging Protocal (SAMP) <www.ivoa.net/samp>`_.

Before integration into `astropy` it was known as `SAMPy <https://pypi.python.org/pypi/sampy/>`_.

Autor: Luigi Paioro (INAF - Istituto Nazionale di Astrofisica)
"""
from .constants import *
from .errors import *
from .utils import *
from .hub import *
from .client import *
from .integrated_client import *
from .hub_proxy import *
