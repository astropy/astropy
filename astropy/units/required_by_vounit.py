# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package defines SI prefixed units that are required by the VOUnit standard
but that are rarely used in practice and liable to lead to confusion (such as
``msolMass`` for milli-solar mass). The units here are enabled so, e.g.,
``Unit('msolMass')`` will just work, but to access the unit directly, use
``astropy.units.required_by_vounit.msolMass`` instead of the more typical idiom
possible for the non-prefixed unit, ``astropy.units.solMass``.
"""

from . import astrophys
from .core import UnitBase, _add_prefixes
from .utils import generate_prefixonly_unit_summary, generate_unit_summary

_add_prefixes(astrophys.solMass, namespace=globals(), prefixes=True)
_add_prefixes(astrophys.solRad, namespace=globals(), prefixes=True)
_add_prefixes(astrophys.solLum, namespace=globals(), prefixes=True)

__all__ = [name for (name, member) in globals().items() if isinstance(member, UnitBase)]

if __doc__ is not None:
    # This generates a docstring for this module that describes all of the
    # standard units defined here.
    __doc__ += generate_unit_summary(globals())
    __doc__ += generate_prefixonly_unit_summary(globals())
