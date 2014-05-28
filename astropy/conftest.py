# this contains imports plugins that configure py.test for astropy tests.
# by importing them here in conftest.py they are discoverable by py.test
# no matter how it is invoked within the astropy tree.

from .tests.pytest_plugins import *

enable_deprecations_as_exceptions()
