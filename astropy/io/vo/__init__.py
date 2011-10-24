from __future__ import division, absolute_import

import os
import sys

if sys.hexversion < 0x02060000:
    raise RuntimeError("vo requires at least Python 2.6")

from .voexceptions import *

__version__ = "0.7.2"
