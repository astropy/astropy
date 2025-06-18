# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Utilities for generating documentation for unit definition modules.

None of the functions in the module are meant for use outside of the
package.
"""

import re
from collections.abc import Iterable, Mapping
from io import StringIO

from .core import NamedUnit, PrefixUnit, Unit, UnitBase


def _get_first_sentence(s: str) -> str:
    """
    Get the first sentence from a string and remove any carriage
    returns.
    """
    x = re.match(r".*?\S\.\s", s)
    if x is not None:
        s = x.group(0)
    return s.replace("\n", " ")


def _summarize_units(
    units: Iterable[NamedUnit],
    have_prefixes: set[str],
    docstring: StringIO,
    template: str,
) -> str:
    for unit in sorted(units, key=lambda x: x.name.lower()):
        represents = ""
        if isinstance(unit, Unit):
            represents = f":math:`{unit._represents.to_string('latex')[1:-1]}`"
        docstring.write(
            template.format(
                unit,
                _get_first_sentence(unit.__doc__).strip(),
                represents,
                ", ".join(f"``{x}``" for x in unit.aliases),
                "Yes" if unit.name in have_prefixes else "No",
            )
        )
    return docstring.getvalue()


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
    units = []
    have_prefixes = set()
    for key, val in namespace.items():
        if not isinstance(val, UnitBase):
            continue
        if not isinstance(val, NamedUnit):
            raise TypeError(f"{key!r} must be defined with 'def_unit()'")
        if key == val.name:  # Skip aliases
            if isinstance(val, PrefixUnit):
                # This will return the root unit that is scaled by the prefix
                # attached to it
                have_prefixes.add(val._represents.bases[0].name)
            else:
                units.append(val)
    docstring = StringIO()
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
    return _summarize_units(units, have_prefixes, docstring, template)


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
    non_prefixed_units = {
        unit.represents.bases[0]
        for unit in namespace.values()
        if isinstance(unit, PrefixUnit)
    }
    template = """
   * - Prefixes for ``{}``
     - {} prefixes
     - {}
     - {}
     - Only
"""
    return _summarize_units(non_prefixed_units, set(), StringIO(), template)
