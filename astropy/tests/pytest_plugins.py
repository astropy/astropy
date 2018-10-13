# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module is included only for backwards compatibility. Packages that
want to use these variables should now import them directly from
`astropy.tests.plugins.display` (although eventually it may be possible to
configure them within setup.cfg).

TODO: This entire module should eventually be removed once backwards
compatibility is no longer supported.
"""
import builtins
import warnings
from ..utils.exceptions import AstropyDeprecationWarning

from .helper import enable_deprecations_as_exceptions
from .plugins.display import PYTEST_HEADER_MODULES, TESTED_VERSIONS

# This makes sure that this module is not collected when running the test
# suite. This is necessary in order to get the test suite to run without errors
# using pytest>=3.7
if getattr(builtins, '_pytest_running', False):
    import pytest
    pytest.skip()

_warning_message = "The module `astropy.tests.pytest_plugins has been " \
    "deprecated. The variables `PYTEST_HEADER_MODULES` and `TESTED_VERSIONS`" \
    "should now be imported from `astropy.tests.plugins.display`. The function " \
    "`enable_deprecations_as_exceptions` should be imported from " \
    "`astropy.tests.helper`"
# Unfortunately, pytest does not display warning messages that occur within
# conftest files, which is where these variables are imported by most packages.
warnings.warn(_warning_message, AstropyDeprecationWarning)
