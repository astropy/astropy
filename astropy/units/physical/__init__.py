# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Defines the physical types that correspond to different units.

The classes and functions defined here are also available in
(and should be used through) the `astropy.units` namespace.
"""

from . import _core
from ._core import *
from ._definitions import *

__all__ = list(_core.__all__)

_attrname_physical_mapping = _core._attrname_physical_mapping
# sphinx-astropy needs _name_physical_mapping for building docs
_name_physical_mapping = _core._name_physical_mapping
_physical_unit_mapping = _core._physical_unit_mapping  # Needed for tests
_unit_physical_mapping = _core._unit_physical_mapping  # Needed for tests


def __getattr__(name: str) -> _core.PhysicalType:
    try:
        return _attrname_physical_mapping[name]
    except KeyError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None


def __dir__() -> list[str]:
    return list(set(__all__) | set(_attrname_physical_mapping))


# This generates a docstring addition for this module that describes all of the
# standard physical types defined here.
if __doc__ is not None:
    doclines = [
        ".. list-table:: Defined Physical Types",
        "    :header-rows: 1",
        "    :widths: 30 10 50",
        "",
        "    * - Physical type",
        "      - Unit",
        "      - Other physical type(s) with same unit",
    ]

    for name in sorted(_name_physical_mapping):
        ptype = _name_physical_mapping[name]
        doclines += [
            f"    * - _`{name}`",
            f"      - :math:`{ptype._unit.to_string('latex')[1:-1]}`",
            f"      - {', '.join([n for n in ptype if n != name])}",
        ]

    __doc__ += "\n\n" + "\n".join(doclines)
