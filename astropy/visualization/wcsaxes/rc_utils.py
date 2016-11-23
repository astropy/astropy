# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import print_function, division, absolute_import

"""
This is a backport of the rc_context class from matplotlib 1.2.
"""

from matplotlib import rcParams

try:
    from matplotlib import rc_context
except ImportError:
    class rc_context(object):

        def __init__(self, rc=None):
            self.rcdict = rc
            self._rcparams = rcParams.copy()
            if self.rcdict:
                rcParams.update(self.rcdict)

        def __enter__(self):
            return self

        def __exit__(self, type, value, tb):
            rcParams.update(self._rcparams)
