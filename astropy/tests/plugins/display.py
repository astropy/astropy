# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This plugin now lives in pytest-astropy, but we keep the code below during
# a deprecation phase.

import warnings
from astropy.utils.exceptions import AstropyDeprecationWarning

try:
    from pytest_astropy_header.display import (PYTEST_HEADER_MODULES,
                                               TESTED_VERSIONS)
except ImportError:
    PYTEST_HEADER_MODULES = {}
    TESTED_VERSIONS = {}

warnings.warn('The astropy.tests.plugins.display plugin has been deprecated. '
              'See the pytest-astropy-header documentation for information on '
              'migrating to using pytest-astropy-header to customize the '
              'pytest header.', AstropyDeprecationWarning)
