# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package contains utilities to run the astropy test suite, tools
for writing tests, and general tests that are not associated with a
particular package.
"""

# NOTE: This is retained only for backwards compatibility. Affiliated packages
# should no longer import `disable_internet` from `astropy.tests`. It is now
# available from `pytest_remotedata`. However, this is not the recommended
# mechanism for controlling access to remote data in tests. Instead, packages
# should make use of decorators provided by the pytest_remotedata plugin:
# - `@pytest.mark.remote_data` for tests that require remote data access
# - `@pytest.mark.internet_off` for tests that should only run when remote data
#       access is disabled.
# Remote data access for the test suite is controlled by the `--remote-data`
# command line flag. This is either passed to `pytest` directly or to the
# `setup.py test` command.
#
# TODO: This import should eventually be removed once backwards compatibility
# is no longer supported.

from pkgutil import find_loader

if find_loader('pytest_remotedata') is not None:
    from pytest_remotedata import disable_internet
else:
    from ..extern.plugins.pytest_remotedata import disable_internet
