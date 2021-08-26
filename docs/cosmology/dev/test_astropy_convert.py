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
astropy = pytest.importorskip("astropy", minversion="5.0")  # isort: skip
# can now import freely from astropy
from astropy import cosmology

# cosmology instances to test reading and writing
astropy_cosmos = [
    getattr(cosmology.realizations, n) for n in cosmology.parameters.available
]

mypackage_cosmos = [myplanck]  # list of ``mypackage`` cosmology realizations


class TestAstropyCosmologyConvert:
    @pytest.mark.parametrize("expected", astropy_cosmos)
    def test_roundtrip_from_astropy(self, expected):
        # convert to ``mypackage``
        mycosmo = expected.to_format("mypackage")

        # Read back, letting autodetection do it's work
        got = Cosmology.from_format(mycosmo)

        # test round-tripped as expected
        assert got == expected  # tests immutable parameters, e.g. H0
        assert got.meta == expected.meta  # (metadata not tested above)

    @pytest.mark.parametrize("expected", mypackage_cosmos)
    def test_roundtrip_from_mypackage(self, expected):
        # convert to Astropy
        acosmo = Cosmology.from_format(expected)

        # convert back to ``mypackage```
        got = acosmo.to_format("mypackage")

        # test round-tripped as expected
        assert isinstance(got, MyCosmology)
        assert got == expected  # assuming ``MyCosmology`` has an __eq__ method

        ...  # more equality tests

    @pytest.mark.parametrize("expected", astropy_cosmos)
    def test_package_equality(self, expected):
        """
        The most important test: are the ``mypackage`` and ``astropy``
        cosmologies equivalent!? They should be if read from the same source
        file.
        """
        # convert to ``mypacakge``
        got = expected.to_format("mypackage")

        # ------------
        # test that the mypackage and astropy cosmologies are equivalent
        assert got.hubble_parameter == expected.H0
        assert got.age(1100) == expected.age(1100)

        ...  # continue with the tests
