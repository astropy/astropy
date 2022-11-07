# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Define a constant to know if the entry points are installed, since this impacts
# whether we can run the tests.

from importlib.metadata import entry_points

import pytest

# TODO: Exclusively use select when Python minversion is 3.10
eps = entry_points()
if hasattr(eps, "select"):
    ep = [entry.name for entry in eps.select(group="asdf_extensions")]
else:
    ep = [entry.name for entry in eps.get("asdf_extensions", [])]
ASDF_ENTRY_INSTALLED = "astropy" in ep and "astropy-asdf" in ep

del entry_points, eps, ep

if not ASDF_ENTRY_INSTALLED:
    pytest.skip(
        "The astropy asdf entry points are not installed", allow_module_level=True
    )
