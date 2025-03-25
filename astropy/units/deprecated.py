# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package defines deprecated units.

These units are not available in the top-level `astropy.units`
namespace. To use these units, you must import the `astropy.units.deprecated`
module::

    >>> from astropy.units import deprecated
    >>> q = 10. * deprecated.emu  # doctest: +SKIP

To include them in `~astropy.units.UnitBase.compose` and the results of
`~astropy.units.UnitBase.find_equivalent_units`, do::

    >>> from astropy.units import deprecated
    >>> deprecated.enable()  # doctest: +SKIP

"""

import warnings

from astropy.utils.decorators import deprecated
from astropy.utils.exceptions import AstropyDeprecationWarning

from . import astrophys, cgs
from .core import UnitBase, _add_prefixes, add_enabled_units, def_unit
from .utils import generate_prefixonly_unit_summary, generate_unit_summary

local_units = {}

def_unit(["emu"], cgs.Bi, namespace=local_units, doc="Biot: CGS (EMU) unit of current")
# Add only some *prefixes* as deprecated units.
_add_prefixes(astrophys.jupiterMass, namespace=local_units, prefixes=True)
_add_prefixes(astrophys.earthMass, namespace=local_units, prefixes=True)
_add_prefixes(astrophys.jupiterRad, namespace=local_units, prefixes=True)
_add_prefixes(astrophys.earthRad, namespace=local_units, prefixes=True)

__all__ = [
    name for (name, member) in local_units.items() if isinstance(member, UnitBase)
]
__all__ += ["enable"]

if __doc__ is not None:
    # This generates a docstring for this module that describes all of the
    # standard units defined here.
    __doc__ += generate_unit_summary(local_units)
    __doc__ += generate_prefixonly_unit_summary(local_units)


def __getattr__(name):
    if unit := local_units.get(name):
        warnings.warn(
            f"{name!r} is deprecated since version 7.1", AstropyDeprecationWarning
        )
        return unit
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


@deprecated(since="7.1")
def enable():
    """
    Enable deprecated units so they appear in results of
    `~astropy.units.UnitBase.find_equivalent_units` and
    `~astropy.units.UnitBase.compose`.

    This may be used with the ``with`` statement to enable deprecated
    units only temporarily.
    """
    return add_enabled_units(local_units)
