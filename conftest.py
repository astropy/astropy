# this contains imports plugins that configure py.test for astropy tests.
# by importing them here in conftest.py they are discoverable by py.test
# no matter how it is invoked within the astropy tree.

from astropy.tests.helper import pytest_addoption, pytest_runtest_setup
from astropy.tests.helper import pytest_report_header
