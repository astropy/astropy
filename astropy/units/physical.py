# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Defines the physical types that correspond to different units.

The classes and functions defined here are also available in
(and should be used through) the `astropy.units` namespace.
"""

from . import _physical_core
from ._physical_core import *  # noqa: F403
from ._physical_types import *  # noqa: F403

__all__ = list(_physical_core.__all__)

_attrname_physical_mapping = _physical_core._attrname_physical_mapping
# sphinx-astropy needs _name_physical_mapping for building docs
_name_physical_mapping = _physical_core._name_physical_mapping
_physical_unit_mapping = _physical_core._physical_unit_mapping  # Needed for tests
_unit_physical_mapping = _physical_core._unit_physical_mapping  # Needed for tests


# For getting the physical types.
def __getattr__(name):
    """Checks for physical types using lazy import.

    This also allows user-defined physical types to be accessible from the
    :mod:`astropy.units.physical` module.
    See `PEP 562 <https://www.python.org/dev/peps/pep-0562/>`_

    Parameters
    ----------
    name : str
        The name of the attribute in this module. If it is already defined,
        then this function is not called.

    Returns
    -------
    ptype : `~astropy.units.physical.PhysicalType`

    Raises
    ------
    AttributeError
        If the ``name`` does not correspond to a physical type
    """
    if name in _attrname_physical_mapping:
        return _attrname_physical_mapping[name]

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    """Return contents directory (__all__ + all physical type names)."""
    return list(set(__all__) | set(_attrname_physical_mapping.keys()))


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

    for name in sorted(_name_physical_mapping.keys()):
        ptype = _name_physical_mapping[name]
        doclines += [
            f"    * - _`{name}`",
            f"      - :math:`{ptype._unit.to_string('latex')[1:-1]}`",
            f"      - {', '.join([n for n in ptype if n != name])}",
        ]

    __doc__ += "\n\n" + "\n".join(doclines)
