# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module is included only for backwards compatibility. Packages that
want to use these variables should now import them directly from
`astropy.tests.plugins.display` (although eventually it may be possible to
configure them within setup.cfg).

TODO: This entire module should eventually be removed once backwards
compatibility is no longer supported.
"""
import warnings

from .helper import enable_deprecations_as_exceptions
from .plugins.display import TESTED_VERSIONS, PYTEST_HEADER_MODULES
from ..utils.exceptions import AstropyDeprecationWarning

_warning_message = "The module `astropy.tests.pytest_plugins has been " \
    "deprecated. The variables `PYTEST_HEADER_MODULES` and `TESTED_VERSIONS`" \
    "should now be imported from `astropy.tests.plugins.display`. The function " \
    "`enable_deprecations_as_exceptions` should be imported from " \
    "`astropy.tests.helper`"
# Unfortunately, pytest does not display warning messages that occur within
# conftest files, which is where these variables are imported by most packages.
warnings.warn(_warning_message, AstropyDeprecationWarning)
