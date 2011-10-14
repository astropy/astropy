#!/usr/bin/env python
from __future__ import division
try:
    from astropy.version import version as __version__
except ImportError:
    # TODO: Issue a warning...
    __version__ = ''
# The version number can be found in the "version" variable of version.py

"""_Astropy is a package intended to contain core functionality and some
common tools needed for performing astronomy and astrophysics research with
Python. It also provides an index for other astronomy packages and tools for
managing them."""

from tests.helper import runtests as test
