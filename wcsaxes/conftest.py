# this contains imports plugins that configure py.test for astropy tests.
# by importing them here in conftest.py they are discoverable by py.test
# no matter how it is invoked within the source tree.

from astropy.tests.pytest_plugins import *
from astropy.tests.pytest_plugins import pytest_addoption as astropy_pytest_addoption

# Uncomment the following line to treat all DeprecationWarnings as
# exceptions
enable_deprecations_as_exceptions()

import os
from astropy.tests.helper import pytest


def pytest_addoption(parser):
    parser.addoption('--generate-images-path', help="directory to generate reference images in", action='store')
    return astropy_pytest_addoption(parser)

@pytest.fixture
def generate(request):
    return request.config.getoption("--generate-images-path")
