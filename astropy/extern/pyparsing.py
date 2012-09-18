# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module just pulls in the appropriate `pyparsing` module, depending on the
currently installed version of python.
"""
from sys import version_info

if version_info[0] > 2:
    from .pyparsing_py3.pyparsing_py3 import *
    from .pyparsing_py3.pyparsing_py3 import __doc__
else:
    from .pyparsing_py2.pyparsing_py2 import *
    from .pyparsing_py2.pyparsing_py2 import __doc__

del version_info #cleans up the namespace
