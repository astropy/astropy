from __future__ import absolute_import

try:
    from argparse import *
except ImportError:
    from .argparse_py2 import *
