# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Utilities for generating documentation for unit definition modules.

None of the functions in the module are meant for use outside of the
package.
"""

from __future__ import annotations

import io
import re
from typing import TYPE_CHECKING

from .core import NamedUnit, PrefixUnit, Unit, UnitBase

if TYPE_CHECKING:
    from collections.abc import Generator, Mapping
    from typing import Literal


def _get_first_sentence(s: str) -> str:
    """
    Get the first sentence from a string and remove any carriage
    returns.
    """
    x = re.match(r".*?\S\.\s", s)
    if x is not None:
        s = x.group(0)
    return s.replace("\n", " ")


def _iter_unit_summary(
    namespace: Mapping[str, object],
) -> Generator[tuple[NamedUnit, str, str, str, Literal["Yes", "No"]], None, None]:
    """
    Generates the ``(unit, doc, represents, aliases, prefixes)``
    tuple used to format the unit summary docs in `generate_unit_summary`.
    """
    # Get all of the units, and keep track of which ones have SI
    # prefixes
    units = []
    has_prefixes = set()
    for key, val in namespace.items():
        # Skip non-unit items
        if not isinstance(val, UnitBase):
            continue

        if not isinstance(val, NamedUnit):
            raise TypeError(f"{key!r} must be defined with 'def_unit()'")

        # Skip aliases
        if key != val.name:
            continue

        if isinstance(val, PrefixUnit):
            # This will return the root unit that is scaled by the prefix
            # attached to it
            has_prefixes.add(val._represents.bases[0].name)
        else:
            units.append(val)

    # Sort alphabetically, case insensitive
    units.sort(key=lambda x: x.name.lower())

    for unit in units:
        doc = _get_first_sentence(unit.__doc__).strip()
        represents = ""
        if isinstance(unit, Unit):
            represents = f":math:`{unit._represents.to_string('latex')[1:-1]}`"
        aliases = ", ".join(f"``{x}``" for x in unit.aliases)

        yield (
            unit,
            doc,
            represents,
            aliases,
            "Yes" if unit.name in has_prefixes else "No",
        )


def generate_unit_summary(namespace: Mapping[str, object]) -> str:
    """
    Generates a summary of units from a given namespace.  This is used
    to generate the docstring for the modules that define the actual
    units.

    Parameters
    ----------
    namespace : dict
        A namespace containing units.

    Returns
    -------
    docstring : str
        A docstring containing a summary table of the units.
    """
    docstring = io.StringIO()

    docstring.write(
        """
.. list-table:: Available Units
   :header-rows: 1
   :widths: 10 20 20 20 1

   * - Unit
     - Description
     - Represents
     - Aliases
     - SI Prefixes
"""
    )
    template = """
   * - ``{}``
     - {}
     - {}
     - {}
     - {}
"""
    for unit_summary in _iter_unit_summary(namespace):
        docstring.write(template.format(*unit_summary))

    return docstring.getvalue()


def generate_prefixonly_unit_summary(namespace: Mapping[str, object]) -> str:
    """
    Generates table entries for units in a namespace that are just prefixes
    without the base unit.  Note that this is intended to be used *after*
    `generate_unit_summary` and therefore does not include the table header.

    Parameters
    ----------
    namespace : dict
        A namespace containing units that are prefixes but do *not* have the
        base unit in their namespace.

    Returns
    -------
    docstring : str
        A docstring containing a summary table of the units.
    """
    faux_namespace = {}
    for unit in namespace.values():
        if isinstance(unit, PrefixUnit):
            base_unit = unit.represents.bases[0]
            faux_namespace[base_unit.name] = base_unit

    docstring = io.StringIO()
    template = """
   * - Prefixes for ``{}``
     - {} prefixes
     - {}
     - {}
     - Only
"""
    for unit_summary in _iter_unit_summary(faux_namespace):
        docstring.write(template.format(*unit_summary))

    return docstring.getvalue()
