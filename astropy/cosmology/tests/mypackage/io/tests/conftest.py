# Licensed under a 3-clause BSD style license - see LICENSE.rst

# THIRD PARTY
import pytest

try:
    import astropy
    from astropy.utils.introspection import minversion
except ImportError:
    ASTROPY_GE_5 = False
else:
    ASTROPY_GE_5 = minversion(astropy, "5.0")


@pytest.fixture(scope="session", autouse=True)
def teardown_mypackage():
    """Clean up module after tests."""

    yield  # to let all tests within the scope run

    if ASTROPY_GE_5:
        from astropy.cosmology import Cosmology
        from astropy.cosmology.connect import convert_registry, readwrite_registry

        readwrite_registry.unregister_reader("myformat", Cosmology)
        readwrite_registry.unregister_writer("myformat", Cosmology)
        readwrite_registry.unregister_identifier("myformat", Cosmology)

        convert_registry.unregister_reader("mypackage", Cosmology)
        convert_registry.unregister_writer("mypackage", Cosmology)
        convert_registry.unregister_identifier("mypackage", Cosmology)
