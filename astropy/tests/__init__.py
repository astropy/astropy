# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package contains utilities to run the astropy test suite, tools
for writing tests, and general tests that are not associated with a
particular package.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from . import helper
from .disable_internet import turn_off_internet, turn_on_internet
