# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Handles backports of the standard library's `gzip.py`.

Python 2.6 has a `gzip.py` without support for the `with` statement.
Here, the version that ships with Python 2.7 is used instead.
"""

from __future__ import absolute_import

import sys
if sys.version_info[:2] == (2, 6):
    from ._gzip_py2 import *
else:
    from gzip import *
