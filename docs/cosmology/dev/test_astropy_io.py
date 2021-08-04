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
from mypackage.io import file_reader, file_writer

# skip all tests in module if import is missing
astropy = pytest.importorskip("astropy", minversion="5.0")  # isort: skip
# can now import freely from astropy
from astropy import cosmology

# the formats being registered with Astropy
readwrite_formats = ["myformat"]
# cosmology instances to test reading and writing
astropy_cosmos = cosmology.parameters.available

# list of ``mypackage`` cosmology realizations
mypackage_cosmos = [myplanck]


class TestAstropyCosmologyIO:
    @pytest.mark.parametrize("format", readwrite_formats)
    @pytest.mark.parametrize("instance", astropy_cosmos)
    def test_roundtrip_from_astropy(self, tmpdir, instance, format):
        cosmo = getattr(cosmology.realizations, instance)
        fname = tmpdir / f"{instance}.{format}"

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
        assert got.meta == cosmo.meta  # (metadata not tested in got == cosmo)

    @pytest.mark.parametrize("format", readwrite_formats)
    @pytest.mark.parametrize("instance", astropy_cosmos)
    def test_package_equality(self, tmpdir, instance, format):
        """
        The most important test: are the ``mypackage`` and ``astropy``
        cosmologies equivalent!? They should be if read from the same source
        file.
        """
        original = getattr(cosmology.realizations, instance)
        fname = tmpdir / f"{instance}.{format}"

        # write with Astropy
        original.write(str(fname), format=format)

        # Read back with ``mypackage``
        cosmo = file_reader(str(fname))  # read to instance from mypackage

        # and a read comparison from Astropy
        cosmo2 = cosmology.Cosmology.read(str(fname), format=format)

        # ------------
        # test that the mypackage and astropy cosmologies are equivalent
        assert original.H0 == cosmo.hubble_parameter
        assert cosmo2.H0 == cosmo.hubble_parameter

        assert original.age(1100) == cosmo.age(1100)
        assert cosmo2.age(1100) == cosmo.age(1100)

        # ...  continue with the tests
