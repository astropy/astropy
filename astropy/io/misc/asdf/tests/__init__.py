# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

# Define a constant to know if the entry points are installed, since this impacts
# whether we can run the tests.

from importlib.metadata import entry_points
import pytest

ep = [entry.name for entry in entry_points().get('asdf_extensions', [])]
ASDF_ENTRY_INSTALLED = 'astropy' in ep and 'astropy-asdf' in ep

del entry_points, ep

if not ASDF_ENTRY_INSTALLED:
    pytest.skip('The astropy asdf entry points are not installed',
                allow_module_level=True)
