from __future__ import absolute_import

import sys

if sys.version_info[:2] <= (2, 6):
    from ._argparse import *
elif sys.version_info[0] == 3 and sys.version_info[:2] <= (3, 1):
    from ._argparse import *
else:
    from argparse import *
