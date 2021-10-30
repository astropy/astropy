# Licensed under a 3-clause BSD style license - see LICENSE.rst

# THIRD PARTY
import pytest

# LOCAL
from mypackage.io import ASTROPY_GE_5


@pytest.fixture(scope="session", autouse=True)
def teardown():
    """Clean up module after tests."""

    yield  # to let all tests within the scope run

    if ASTROPY_GE_5:
        from astropy.cosmology import Cosmology
        from astropy.io import registry as io_registry

        io_registry.unregister_reader("myformat", Cosmology)
        io_registry.unregister_writer("myformat", Cosmology)
        io_registry.unregister_identifier("myformat", Cosmology)

        io_registry.unregister_reader("mypackage", Cosmology)
        io_registry.unregister_writer("mypackage", Cosmology)
        io_registry.unregister_identifier("mypackage", Cosmology)
