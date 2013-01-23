# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module just pulls in the appropriate `configobj` package, depending on the
currently installed version of python.

Also, this should actually never actually show up as a docstring, because
it should get overwritten by the appropriate configobj docstring.
"""
from sys import version_info

if version_info[0] > 2:
    from .configobj_py3 import configobj, validate, __doc__
else:
    from .configobj_py2 import configobj, validate, __doc__

del version_info #cleans up the namespace
