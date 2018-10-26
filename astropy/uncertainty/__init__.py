# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This sub-package contains classes and functions for creating distributions that
work similar to `~astropy.units.Quantity` or array objects, but can propogate
uncertinties.
"""


from .core import *
from .builtin_distrs import *
