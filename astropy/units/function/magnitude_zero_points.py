# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module has been deprecated and moved to astropy.units.photometric.  The
names remain here for backwards compatibility.
"""
from warnings import warn

from ..photometric import AB, ST

_ns = globals()


def enable():
    """
    Enable magnitude zero point units so they appear in results of
    `~astropy.units.UnitBase.find_equivalent_units` and
    `~astropy.units.UnitBase.compose`.

    This may be used with the ``with`` statement to enable these
    units only temporarily.
    """
    warn('The magnitude_zero_points module has been deprecated, and moved to '
         'astropy.units.photometric and are enabled by default. '
         'magnitude_zero_points is retained as aliases to the new units.')
