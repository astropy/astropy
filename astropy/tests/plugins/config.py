# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This plugin provides customization of configuration used by pytest.
"""
from astropy.tests.helper import treat_deprecations_as_exceptions


def pytest_configure(config):
    treat_deprecations_as_exceptions()
