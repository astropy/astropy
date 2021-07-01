# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

import pytest

from astropy import cosmology
from astropy.cosmology import Cosmology
from astropy import table
from astropy.io.registry import get_formats
from astropy.utils.compat import optional_deps

cosmo_instances = cosmology.parameters.available
table_formats = ["json", "ascii.ecsv"]

# TODO! remove in astropy v5.0
if not getattr(optional_deps, "HAS_YAML"):
    table_formats.remove("ascii.ecsv")


# make a common directory for reading / writing cosmologies
@pytest.fixture(scope="session")
def cosmo_dir(tmpdir_factory):
    drct = tmpdir_factory.mktemp("cosmo")
    return drct


# -----------------------------------------------------------------------------


class TestWriteCosmology:
    @pytest.mark.parametrize("format", table_formats)
    @pytest.mark.parametrize("instance", cosmo_instances)
    def test_write_table_file(self, cosmo_dir, instance, format):
        """Read tests happen later."""
        cosmo = getattr(cosmology.realizations, instance)
        fname = cosmo_dir / f"{instance}.{format}"

        cosmo.write(str(fname), format=format)

        # Also test kwarg "overwrite"
        assert os.path.exists(str(fname))  # file exists
        with pytest.raises(IOError):
            cosmo.write(str(fname), format=format, overwrite=False)

        assert os.path.exists(str(fname))  # file exists
        cosmo.write(str(fname), format=format, overwrite=True)


class TestReadCosmology:
    @pytest.mark.parametrize("instance", cosmo_instances)
    def test_read_mapping_instance(self, instance):
        expected = getattr(cosmology.realizations, instance)
        params = getattr(cosmology.parameters, instance)

        # read from most general class
        got = Cosmology.read.from_mapping(params)
        assert got.name == expected.name

        # read from expected class
        got = expected.__class__.read.from_mapping(params)
        assert got.name == expected.name
        # assert got == expected  # FIXME! no __eq__ on cosmo

    @pytest.mark.parametrize("instance", cosmo_instances)
    def test_read_mapping_instance(self, instance):
        expected = getattr(cosmology.realizations, instance)
        params = getattr(cosmology.parameters, instance)

        # read from most general class
        got = Cosmology.read.from_mapping(params)
        assert got.name == expected.name

        # read from expected class
        got = expected.__class__.read.from_mapping(params)
        assert got.name == expected.name
        # assert got == expected  # FIXME! no __eq__ on cosmo

    @pytest.mark.parametrize("format", table_formats)
    @pytest.mark.parametrize("instance", cosmo_instances)
    def test_read_table_file(self, cosmo_dir, instance, format):
        expected = getattr(cosmology.realizations, instance)
        got = Cosmology.read(cosmo_dir / f"{instance}.{format}", format=format)

        assert got.name == expected.name  # FIXME!
        # assert got == expected  # FIXME! no __eq__ on cosmo

    @pytest.mark.parametrize("format", table_formats[1:])  # skip json
    @pytest.mark.parametrize("instance", cosmo_instances)
    def test_read_table_instance(self, cosmo_dir, instance, format):
        expected = getattr(cosmology.realizations, instance)

        tbl = table.QTable.read(cosmo_dir / f"{instance}.{format}", format=format)
        got = Cosmology.read.from_table(tbl)

        assert got.name == expected.name  # FIXME!
        # assert got == expected  # FIXME! no __eq__ on cosmo
