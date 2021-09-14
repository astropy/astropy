# -*- coding: utf-8 -*-

"""
Test that the conversion interface with Astropy works as expected.
"""

# STDLIB
import os

# THIRD PARTY
import pytest

# LOCAL
from mypackage.cosmology import MyCosmology, myplanck

# skip all tests in module if import is missing
astropy = pytest.importorskip("astropy", minversion="4.3")  # isort: skip
# can now import freely from astropy
from astropy import cosmology
from astropy.cosmology import Cosmology
from astropy.io import registry as io_registry
from astropy.utils.compat.optional_deps import HAS_SCIPY

# cosmology instances to test reading and writing
astropy_cosmos = [
    getattr(cosmology.realizations, n) for n in cosmology.parameters.available
]

mypackage_cosmos = [myplanck]  # list of ``mypackage`` cosmology realizations


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy.")
class TestAstropyCosmologyConvert:
    """Test conversion to/from Astropy."""

    @pytest.mark.parametrize("expected", astropy_cosmos)
    def test_roundtrip_from_astropy(self, expected):
        # convert to ``mypackage``
        mycosmo = expected.to_format("mypackage")

        # Read back
        got = Cosmology.from_format(mycosmo, format="mypackage")

        # test round-tripped as expected
        assert got == expected  # tests immutable parameters, e.g. H0

        # NOTE: if your package's cosmology supports metadata
        # assert got.meta == expected.meta  # (metadata not tested above)

    @pytest.mark.parametrize("expected", mypackage_cosmos)
    def test_roundtrip_from_mypackage(self, expected):
        # convert to Astropy
        acosmo = Cosmology.from_format(expected, format="mypackage")

        # convert back to ``mypackage```
        got = acosmo.to_format("mypackage")

        # test round-tripped as expected
        assert isinstance(got, MyCosmology)
        assert got == expected  # assuming ``MyCosmology`` has an __eq__ method

        ...  # more equality tests

    @pytest.mark.parametrize("acosmo", astropy_cosmos)
    def test_package_equality(self, acosmo):
        """
        The most important test: are the ``mypackage`` and ``astropy``
        cosmologies equivalent!? They should be if read from the same source
        file.
        """
        # convert to ``mypacakge``
        mycosmo = acosmo.to_format("mypackage")

        # ------------
        # test that the mypackage and astropy cosmologies are equivalent

        assert mycosmo.name == acosmo.name
        assert mycosmo.hubble_parameter == acosmo.H0.value

        ...  # continue with the tests

        # e.g. test some methods
        # assert mycosmo.age_in_Gyr(1100) == acosmo.age(1100).to_value(u.Gyr)
