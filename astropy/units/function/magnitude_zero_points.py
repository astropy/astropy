# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module has been deprecated and moved to astropy.units.photometric.  The
names remain here for backwards compatibility.
"""
from warnings import warn

from ..photometric import AB, ST
from ...utils import deprecated

_ns = globals()


@deprecated(since='3.1', alternative='astropy.units.photometric',
            message='The magnitude_zero_points module has been deprecated, and'
                    ' moved to astropy.units.photometric and are enabled by '
                    'default. magnitude_zero_points is retained as aliases to '
                    'the new units.')
def enable():
    """
    Enable magnitude zero point units so they appear in results of
    `~astropy.units.UnitBase.find_equivalent_units` and
    `~astropy.units.UnitBase.compose`.

    This may be used with the ``with`` statement to enable these
    units only temporarily.
    """
    # While it may seem like the below can be removed, in fact it needs to
    # remain as long as this function is around so that enable acts as a context
    # manager

    # Local import to avoid cyclical import
    from ..core import add_enabled_units
    # Local import to avoid polluting namespace
    import inspect
    return add_enabled_units(inspect.getmodule(enable))
