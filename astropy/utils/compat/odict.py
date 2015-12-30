# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import
from pprint import pformat

try:
    from collections import OrderedDict as NativeOrderedDict
except ImportError:
    from ._odict_py2 import OrderedDict as NativeOrderedDict

class OrderedDict(NativeOrderedDict):
    r"""Dictionary that remembers insertion order

    The native OrderedDict class is extended here in order to
    provide a more pleasing printed or string representation.

    """
    def __init__(self, *args, **kwds):
        r"""
        Initialize an ordered dictionary. The signature is the same as
        regular dictionaries, but keyword arguments are not recommended
        because their insertion order is arbitrary.
        """
        NativeOrderedDict.__init__(self, *args, **kwds)
    def __str__(self):
        r"""
        Pretty-printed string representation of OrderedDict.
        x.__str__() <==> str(x)
        """
        return('OrderedDict('+ pformat(self.items(),width=80)+')')

