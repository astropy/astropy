# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This file contains pytest configuration settings that are astropy-specific
(i.e.  those that would not necessarily be shared by affiliated packages
making use of astropy's test runner).
"""
from importlib.util import find_spec

from astropy.tests.plugins.display import PYTEST_HEADER_MODULES
from astropy.tests.helper import enable_deprecations_as_exceptions

if find_spec('asdf') is not None:
    pytest_plugins = ['asdf.tests.schema_tester']


enable_deprecations_as_exceptions(
    include_astropy_deprecations=False,
    # This is a workaround for the OpenSSL deprecation warning that comes from
    # the `requests` module. It only appears when both asdf and sphinx are
    # installed. This can be removed once pyopenssl 1.7.20+ is released.
    modules_to_ignore_on_import=['requests'])

try:
    import matplotlib
except ImportError:
    pass
else:
    matplotlib.use('Agg')

PYTEST_HEADER_MODULES['Cython'] = 'cython'
