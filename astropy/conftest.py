# this contains imports plugins that configure py.test for astropy tests.
# by importing them here in conftest.py they are discoverable by py.test
# no matter how it is invoked within the astropy tree.

from .tests.pytest_plugins import *

from .utils import turn_off_internet,turn_on_internet

# need to rename pytest_configure as we are about to redefine/overload it
_pytest_configure = pytest_configure

# pytest magic:
# http://pytest.org/latest/plugins.html#_pytest.hookspec.pytest_configure
# use pytest.set_trace() to interactively inspect config's features
def pytest_configure(config):
    if config.getoption('remote_data'):
        pass
        # This isn't needed; just don't turn off the internet if remote_data is on
        #turn_on_internet(verbose=config.option.verbose)
    else:
        turn_off_internet(verbose=config.option.verbose)

    # imported from .test.pytest_plugins; see above
    _pytest_configure(config)
