# -*- coding: utf-8 -*-

"""
Test that the interface with Astropy works as expected.
"""

# STDLIB
import os

# THIRD PARTY
import pytest

# LOCAL
from mypackage.cosmology import myplanck
from mypackage.io import file_reader

# skip all tests in module if import is missing
astropy = pytest.importorskip("astropy", minversion="4.3")  # isort: skip
# can now import freely from astropy
import astropy.units as u
from astropy import cosmology
from astropy.cosmology import Cosmology
from astropy.io import registry as io_registry
from astropy.utils.compat.optional_deps import HAS_SCIPY

# the formats being registered with Astropy
readwrite_formats = ["myformat"]
# cosmology instances to test reading and writing
astropy_cosmos = cosmology.parameters.available

# list of ``mypackage`` cosmology realizations
mypackage_cosmos = [myplanck]


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy.")
class TestAstropyCosmologyIO:
    """Test read/write interoperability with Astropy."""

    @pytest.mark.parametrize("format", readwrite_formats)
    @pytest.mark.parametrize("instance", astropy_cosmos)
    def test_roundtrip_from_astropy(self, tmp_path, instance, format):
        cosmo = getattr(cosmology.realizations, instance)
        fname = tmp_path / f"{instance}.{format}"

        # write to file
        cosmo.write(str(fname), format=format)

        # also test kwarg "overwrite"
        assert os.path.exists(str(fname))  # file exists
        with pytest.raises(IOError):
            cosmo.write(str(fname), format=format, overwrite=False)

        assert os.path.exists(str(fname))  # overwrite file existing file
        cosmo.write(str(fname), format=format, overwrite=True)

        # Read back
        got = Cosmology.read(fname, format=format)

        # test round-tripped as expected
        assert got == cosmo  # tests immutable parameters, e.g. H0

        # NOTE: if your package's cosmology supports metadata
        # assert got.meta == expected.meta  # (metadata not tested above)

    @pytest.mark.parametrize("format", readwrite_formats)
    @pytest.mark.parametrize("instance", astropy_cosmos)
    def test_package_equality(self, tmp_path, instance, format):
        """
        The most important test: are the ``mypackage`` and ``astropy``
        cosmologies equivalent!? They should be if read from the same source
        file.
        """
        original = getattr(cosmology.realizations, instance)
        fname = tmp_path / f"{instance}.{format}"

        # write with Astropy
        original.write(str(fname), format=format)

        # Read back with ``myformat``
        cosmo = file_reader(str(fname))  # read to instance from mypackage

        # and a read comparison from Astropy
        cosmo2 = cosmology.Cosmology.read(str(fname), format=format)

        # ------------
        # test that the mypackage and astropy cosmologies are equivalent
        assert original.H0.value == cosmo.hubble_parameter
        assert cosmo2.H0.value == cosmo.hubble_parameter

        ...  # continue with the tests

        # e.g. test some methods
        # assert original.age(1100).to_value(u.Gyr) == cosmo.age_in_Gyr(1100)
        # assert cosmo2.age(1100) == cosmo.age(1100)
