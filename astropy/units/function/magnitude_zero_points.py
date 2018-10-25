# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module has been deprecated and moved to astropy.units.photometric.  The
names remain here for backwards compatibility.
"""
from warnings import warn

from ..core import def_unit
from ..photometric import ABflux, STflux

_ns = globals()

def_unit(['AB'], ABflux, namespace=_ns)
def_unit(['ST'], STflux, namespace=_ns)


def enable():
    """
    Enable magnitude zero point units so they appear in results of
    `~astropy.units.UnitBase.find_equivalent_units` and
    `~astropy.units.UnitBase.compose`.

    This may be used with the ``with`` statement to enable these
    units only temporarily.
    """
    warn('The magnitude_zero_points module has been deprecated, and moved to '
         'astropy.units.photometric (with more appropriate default names). '
         'the magnitude_zero_points are retained as aliases to the new units.')

    # Local import to avoid cyclical import
    from ..core import add_enabled_units
    # Local import to avoid polluting namespace
    import inspect
    return add_enabled_units(inspect.getmodule(enable))
