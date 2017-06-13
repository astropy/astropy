# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package defines SI prefix ed units that are required by the VOUnit standard
but that are rarely used in practice.  E.g. ``Unit('msolMass')`` will yield
milli-solar mass, but are in a separate module so they must be accessed as
``astropy.units.required_by_vounit.msolMass`` instead of ``astropy.units.msolMass``.

These units are in a separate modeule from `astropy.units.deprecated` because
they need to be enabled by default for `astropy.units` to be able to compliantly
parse VOUnit strings.

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

_ns = globals()


def _initialize_module():
    # Local imports to avoid polluting top-level namespace
    from . import cgs
    from . import astrophys
    from .core import def_unit, _add_prefixes

    _add_prefixes(astrophys.solMass, namespace=_ns, prefixes=True)
    _add_prefixes(astrophys.solRad, namespace=_ns, prefixes=True)
    _add_prefixes(astrophys.solLum, namespace=_ns, prefixes=True)


_initialize_module()


###########################################################################
# DOCSTRING

# This generates a docstring for this module that describes all of the
# standard units defined here.
from .utils import generate_unit_summary as _generate_unit_summary  # noqa
from .utils import generate_prefixonly_unit_summary as _generate_prefixonly_unit_summary  # noqa
if __doc__ is not None:
    __doc__ += _generate_unit_summary(globals())
    __doc__ += _generate_prefixonly_unit_summary(globals())


def _enable():
    """
    Enable deprecated units so they appear in results of
    `~astropy.units.UnitBase.find_equivalent_units` and
    `~astropy.units.UnitBase.compose`.

    This may be used with the ``with`` statement to enable deprecated
    units only temporarily.
    """
    # Local import to avoid cyclical import
    from .core import add_enabled_units
    # Local import to avoid polluting namespace
    import inspect
    return add_enabled_units(inspect.getmodule(_enable))


# Because these are VOUnit mandated units, they start enabled (which is why the
# function is hidden).
_enable()
