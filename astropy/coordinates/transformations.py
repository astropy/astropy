# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the framework for transforming points from
one coordinate system to another (e.g. equatorial to galactic). The
implementation is actually in individual coordinates in the
`builtin_systems` module, while this module provides the framework and
related utilities.
"""

import math

import numpy as np

pi = math.pi

__all__ = []
