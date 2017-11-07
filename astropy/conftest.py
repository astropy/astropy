# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This file contains pytest configuration settings that are astropy-specific
(i.e.  those that would not necessarily be shared by affiliated packages
making use of astropy's test runner).
"""
from astropy.tests.plugins.display import PYTEST_HEADER_MODULES
from astropy.tests.helper import enable_deprecations_as_exceptions

enable_deprecations_as_exceptions(include_astropy_deprecations=False)

try:
    import matplotlib
except ImportError:
    pass
else:
    matplotlib.use('Agg')

# This is astropy-specific so should not be included in a generic plugin that
# could potentially be used by other projects
PYTEST_HEADER_MODULES['Cython'] = 'cython'
