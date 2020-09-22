# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This is retained only for backwards compatibility. Affiliated packages
should no longer import ``disable_internet`` from ``astropy.tests``. It is
now available from ``pytest_remotedata``. However, this is not the
recommended mechanism for controlling access to remote data in tests.
Instead, packages should make use of decorators provided by the
pytest_remotedata plugin: - ``@pytest.mark.remote_data`` for tests that
require remote data access - ``@pytest.mark.internet_off`` for tests that
should only run when remote data access is disabled.  Remote data access for
the test suite is controlled by the ``--remote-data`` command line flag. This
is passed to ``pytest`` directly.

TODO: This module should eventually be removed once backwards compatibility
is no longer supported.
"""
from warnings import warn
from astropy.utils.exceptions import AstropyDeprecationWarning


warn("The ``disable_internet`` module is no longer provided by astropy. It "
     "is now available as ``pytest_remotedata.disable_internet``. However, "
     "developers are encouraged to avoid using this module directly. See "
     "<https://docs.astropy.org/en/latest/whatsnew/3.0.html#pytest-plugins> "
     "for more information.", AstropyDeprecationWarning)


try:
    # This should only be necessary during testing, in which case the test
    # package must be installed anyway.
    from pytest_remotedata.disable_internet import *
except ImportError:
    pass
