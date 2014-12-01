# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module just pulls in the appropriate `configobj` package, depending on the
currently installed version of python.

Also, this should actually never actually show up as a docstring, because
it should get overwritten by the appropriate configobj docstring.
"""

from .configobj import configobj, validate, __doc__
