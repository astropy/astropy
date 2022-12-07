# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package defines SI prefixed units that are required by the VOUnit standard
but that are rarely used in practice and liable to lead to confusion (such as
``msolMass`` for milli-solar mass). They are in a separate module from
`astropy.units.deprecated` because they need to be enabled by default for
`astropy.units` to parse compliant VOUnit strings. As a result, e.g.,
``Unit('msolMass')`` will just work, but to access the unit directly, use
``astropy.units.required_by_vounit.msolMass`` instead of the more typical idiom
possible for the non-prefixed unit, ``astropy.units.solMass``.
"""

_ns = globals()


def _initialize_module():
    # Local imports to avoid polluting top-level namespace
    from . import astrophys
    from .core import _add_prefixes

    _add_prefixes(astrophys.solMass, namespace=_ns, prefixes=True)
    _add_prefixes(astrophys.solRad, namespace=_ns, prefixes=True)
    _add_prefixes(astrophys.solLum, namespace=_ns, prefixes=True)


_initialize_module()


###########################################################################
# DOCSTRING

if __doc__ is not None:
    # This generates a docstring for this module that describes all of the
    # standard units defined here.
    from .utils import (
        generate_prefixonly_unit_summary as _generate_prefixonly_unit_summary,
    )
    from .utils import generate_unit_summary as _generate_unit_summary

    __doc__ += _generate_unit_summary(globals())
    __doc__ += _generate_prefixonly_unit_summary(globals())


def _enable():
    """
    Enable the VOUnit-required extra units so they appear in results of
    `~astropy.units.UnitBase.find_equivalent_units` and
    `~astropy.units.UnitBase.compose`, and are recognized in the ``Unit('...')``
    idiom.
    """
    # Local import to avoid cyclical import
    # Local import to avoid polluting namespace
    import inspect

    from .core import add_enabled_units

    return add_enabled_units(inspect.getmodule(_enable))


# Because these are VOUnit mandated units, they start enabled (which is why the
# function is hidden).
_enable()
