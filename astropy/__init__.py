# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astropy is a package intended to contain core functionality and some
common tools needed for performing astronomy and astrophysics research with
Python. It also provides an index for other astronomy packages and tools for
managing them.
"""


try:
    from .version import version as __version__
except ImportError:
    # TODO: Issue a warning...
    __version__ = ''
# The version number can be found in the "version" variable of version.py

# set up the test command
from . import __path__ as base_path
from .tests import helper
test_runner = helper.TestRunner(base_path[0])
test = test_runner.run_tests

# make the setup.py test command
class setup_test_command(helper.astropy_test):
    package_name = 'astropy'
