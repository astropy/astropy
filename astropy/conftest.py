# this contains imports plugins that configure py.test for astropy tests.
# by importing them here in conftest.py they are discoverable by py.test
# no matter how it is invoked within the astropy tree.

from .tests.pytest_plugins import *

from .utils import turn_off_internet,turn_on_internet

# pytest magic:
# http://pytest.org/latest/plugins.html#_pytest.hookspec.pytest_configure
# use pytest.set_trace() to interactively inspect config's features
def pytest_configure(config):
    if config.getoption('remote_data'):
        pass
        #turn_on_internet(verbose=config.option.verbose)
    else:
        turn_off_internet(verbose=config.option.verbose)

    try:
        from astropy.tests.pytest_plugins import pytest_configure

        pytest_configure(config)
    except ImportError:
        # assume astropy v<0.3
        pass
