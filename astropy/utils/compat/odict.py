# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import
from pprint import pformat

try:
    from collections import OrderedDict as NativeOrderedDict
except ImportError:
    from ._odict_py2 import OrderedDict as NativeOrderedDict

class OrderedDict(NativeOrderedDict):
    def __init__(self, *args, **kwds):
        NativeOrderedDict.__init__(self, *args, **kwds)
    def __str__(self):
        return('OrderedDict('+ pformat(dict(self),width=1)[1:-1]+')')
