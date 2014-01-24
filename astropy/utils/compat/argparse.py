from __future__ import absolute_import

import sys

if sys.version_info[:2] <= (2, 6):
    from ._argparse_py2 import *
else:
    from argparse import *
