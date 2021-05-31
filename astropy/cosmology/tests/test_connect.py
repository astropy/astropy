# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

import pytest

from astropy import cosmology
from astropy.cosmology import Cosmology
from astropy import table
from astropy.io.registry import get_formats
from astropy.utils.compat import optional_deps
from astropy.cosmology.connect import CosmologyRead, CosmologyWrite
from astropy.utils.exceptions import AstropyUserWarning

cosmo_instances = cosmology.parameters.available
save_formats = ["json", "ascii.ecsv"]

# TODO! remove in astropy v5.0
if not getattr(optional_deps, "HAS_YAML"):
    save_formats.remove("ascii.ecsv")


# make a common directory for reading / writing cosmologies
@pytest.fixture(scope="session")
def cosmo_dir(tmpdir_factory):
    drct = tmpdir_factory.mktemp("cosmo")
    return drct


# -----------------------------------------------------------------------------


class TestReadWriteCosmology:

    def test_instantiate_read(self):
        # no error on base class
        assert isinstance(Cosmology.read, CosmologyRead)

        # Error on calling from subclasses. Warns when initiated.
        with pytest.warns(AstropyUserWarning):
            with pytest.raises(TypeError, match="``Cosmology`` base class"):
                cosmology.realizations.Planck18.read()

    @pytest.mark.parametrize("format", save_formats)
    @pytest.mark.parametrize("instance", cosmo_instances)
    def test_write_then_read_file(self, cosmo_dir, instance, format):
        """Read tests happen later."""
        cosmo = getattr(cosmology.realizations, instance)
        fname = cosmo_dir / f"{instance}.{format}"

        cosmo.write(str(fname), format=format)

        # Also test kwarg "overwrite"
        assert os.path.exists(str(fname))  # file exists
        with pytest.raises(IOError):
            cosmo.write(str(fname), format=format, overwrite=False)

        assert os.path.exists(str(fname))  # overwrite file existing file
        cosmo.write(str(fname), format=format, overwrite=True)

        # Read back
        got = Cosmology.read(cosmo_dir / f"{instance}.{format}", format=format)

        assert got.name == cosmo.name
        # assert got == expected  # FIXME! no __eq__ on cosmo

    @pytest.mark.parametrize("instance", cosmo_instances)
    def test_to_mapping_instance(self, instance):
        instance = getattr(cosmology.realizations, instance)
        m = instance.write.to_mapping()

        assert isinstance(m, dict)
        assert "cosmology" in m

    @pytest.mark.parametrize("instance", cosmo_instances)
    def test_to_table_instance(self, instance):
        instance = getattr(cosmology.realizations, instance)
        t = instance.write.to_table()

        assert isinstance(t, table.QTable)
        assert "cosmology" in t.meta

    @pytest.mark.parametrize("instance", cosmo_instances)
    def test_from_mapping_instance(self, instance):
        expected = getattr(cosmology.realizations, instance)
        params = getattr(cosmology.parameters, instance)

        # read with mismatching parameters errors
        with pytest.raises(TypeError, match="There are unused parameters"):
            Cosmology.read.from_mapping(params)

        # unless mismatched are moved to meta
        got = Cosmology.read.from_mapping(params, move_to_meta=True)
        assert got.name == expected.name
        # assert got == expected  # FIXME! no __eq__ on cosmo

    @pytest.mark.parametrize("format", save_formats[1:])  # skip json
    @pytest.mark.parametrize("instance", cosmo_instances)
    def test_from_table_instance(self, cosmo_dir, instance, format):
        expected = getattr(cosmology.realizations, instance)

        tbl = table.QTable.read(cosmo_dir / f"{instance}.{format}", format=format)
        got = Cosmology.read.from_table(tbl)

        assert got.name == expected.name
        # assert got == expected  # FIXME! no __eq__ on cosmo


# -----------------------------------------------------------------------------


@pytest.mark.parametrize("instance", cosmo_instances)
def test_round_trip_of_mapping_instance(instance):
    expected = getattr(cosmology.realizations, instance)

    got = Cosmology.read.from_mapping(expected.write.to_mapping())
    assert got.name == expected.name
    # assert got == expected  # FIXME! no __eq__ on cosmo


@pytest.mark.parametrize("instance", cosmo_instances)
class Test_round_trip_of_table_instance:

    def test_1D(self, instance):
        expected = getattr(cosmology.realizations, instance)

        got = Cosmology.read.from_table(expected.write.to_table())
        assert got.name == expected.name
        # assert got == expected  # FIXME! no __eq__ on cosmo

    def test_ND(self, instance):
        expected = getattr(cosmology.realizations, instance)
        t = table.vstack([expected.write.to_table(), expected.write.to_table()])
        t["name"][1] = "Other"

        got = Cosmology.read.from_table(t, index=0)
        assert got.name == expected.name
        # assert got == expected  # FIXME! no __eq__ on cosmo

        got = Cosmology.read.from_table(t, index=1)
        assert got.name == "Other"
        # assert got == expected  # FIXME! no __eq__ on cosmo

        # Now test index is string
        got = Cosmology.read.from_table(t, index="Other")
        assert got.name == "Other"
        # assert got == expected  # FIXME! no __eq__ on cosmo
