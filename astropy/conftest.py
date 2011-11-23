# this contains plugins that configure py.test for astropy tests

from .tests.helper import pytest

# these pytest hooks allow us to mark tests and run the marked tests with
# specific command line options.
def pytest_addoption(parser):
    parser.addoption("--remotedata", action="store_true",
        help="run tests with online data")


def pytest_runtest_setup(item):
    if ('remote_data' in item.keywords and
        not item.config.getvalue("remotedata")):
        pytest.skip("need --remotedata option to run")
