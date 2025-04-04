# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Defines the physical types that correspond to different units.

The classes and functions defined here are also available in
(and should be used through) the `astropy.units` namespace.
"""

from . import core
from .core import *
from .core import _name_physical_mapping  # sphinx-astropy needs this when building docs
from .definitions import *

__all__ = list(core.__all__)


def __getattr__(name):
    try:
        return getattr(core, name)
    except AttributeError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None


def __dir__() -> list[str]:
    return dir(core)


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
