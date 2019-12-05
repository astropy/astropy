# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

# Define a constant to know if the entry points are installed, since this impacts
# whether we can run the tests.

import pytest
from pkg_resources import iter_entry_points

entry_points = [entry.name for entry in iter_entry_points('asdf_extensions')]
ASDF_ENTRY_INSTALLED = 'astropy' in entry_points and 'astropy-asdf' in entry_points

del entry_points, iter_entry_points

if not ASDF_ENTRY_INSTALLED:
    pytest.skip('The astropy asdf entry points are not installed',
                allow_module_level=True)
