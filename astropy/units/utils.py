# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Miscellaneous utilities for `astropy.units`.

None of the functions in the module are meant for use outside of the
package.
"""
from __future__ import absolute_import, division, print_function, unicode_literals

from . import core


def _get_first_sentence(s):
    """
    Get the first sentence from a string and remove any carriage
    returns.
    """
    import re
    x = re.match(".*?\S\.\s", s)
    if x is not None:
        s = x.group(0)
    return s.replace('\n', ' ')


def generate_unit_summary(namespace):
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
    from . import core
    import io

    # Get all of the units, and keep track of which ones have SI
    # prefixes
    units = []
    has_prefixes = set()
    for key, val in namespace.items():
        # Skip non-unit items
        if not isinstance(val, core.UnitBase):
            continue

        # Skip aliases
        if key != val.name:
            continue

        if isinstance(val, core.PrefixUnit):
            has_prefixes.add(val.decompose().bases[0].name)
        else:
            units.append(val)

    # Sort alphabetically, case insensitive
    units.sort(key=lambda x: x.name.lower())

    docstring = io.StringIO()

    docstring.write("""
.. list-table:: Available Units
   :header-rows: 1
   :widths: 10 20 20 20 1

   * - Unit
     - Description
     - Represents
     - Aliases
     - SI Prefixes
""")

    for unit in units:
        if unit.name in has_prefixes:
            unit_has_prefixes = 'Y'
        else:
            unit_has_prefixes = 'N'
        doc = _get_first_sentence(unit.__doc__).strip()
        represents = ''
        if isinstance(unit, core.Unit):
            represents = ":math:`{0}`".format(
                unit._represents.to_string('latex')[1:-1])
        aliases = ', '.join('``{0}``'.format(x) for x in unit.aliases)
        docstring.write("""
   * - ``{0}``
     - {1}
     - {2}
     - {3}
     - {4}
""".format(unit, doc, represents, aliases, unit_has_prefixes))

    return docstring.getvalue()
