# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Handles backports of the standard library's `fractions.py`.

The fractions module in 2.6 does not handle being instantiated using a
float and then calculating an approximate fraction based on that.
This functionality is required by the FITS unit format generator,
since the FITS unit format handles only rational, not decimal point,
powers.
"""

from __future__ import absolute_import

import sys
if sys.version_info[:2] == (2, 6):
    from ._fractions_py2 import *
else:
    from fractions import *
