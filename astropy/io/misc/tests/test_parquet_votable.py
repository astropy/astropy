import copy

import numpy as np
import pytest

import astropy.units as u
from astropy.io.misc.parquet import read_parquet_votable, write_parquet_votable
from astropy.table import Table
from astropy.utils.data import get_pkg_data_filename
from astropy.utils.exceptions import AstropyUserWarning
from astropy.utils.misc import _NOT_OVERWRITING_MSG_MATCH

# Skip all tests in this file if we cannot import pyarrow
pyarrow = pytest.importorskip("pyarrow")


number_of_objects = 10
ids = [f"COSMOS_{ii}" for ii in range(number_of_objects)]
redshift = np.random.uniform(low=0, high=3, size=number_of_objects)
mass = np.random.uniform(low=1e8, high=1e10, size=number_of_objects)
sfr = np.random.uniform(low=1, high=100, size=number_of_objects)

input_table = Table(
    [ids, redshift, mass, sfr],
    names=["id", "z", "mass", "sfr"],
)

column_metadata = {
    "id": {
        "unit": "",
        "ucd": "meta.id",
        "utype": "none",
        "description": "The ID of the galaxy.",
    },
    "z": {
        "unit": "",
        "ucd": "src.redshift",
        "utype": "none",
        "description": "The redshift of the galaxy.",
    },
    "mass": {
        "unit": "solMass",
        "ucd": "phys.mass",
        "utype": "none",
        "description": "The stellar mass of the galaxy.",
    },
    "sfr": {
        "unit": "solMass / yr",
        "ucd": "phys.SFR",
        "utype": "none",
        "description": "The star formation rate of the galaxy.",
    },
}


def test_parquet_votable(tmp_path):
    """Test writing and reading a votable into a parquet file"""
    filename = tmp_path / "test_votable.parq"

    write_parquet_votable(input_table, filename, metadata=column_metadata)

    with pytest.warns(AstropyUserWarning, match="No table::len"):
        loaded_table = read_parquet_votable(filename)

    # Compare data content, but not metadata (input had none)
    assert np.all(input_table == loaded_table)


@pytest.mark.parametrize("overwrite_metadata", [True, False])
def test_parquet_votable_input_column_unit(overwrite_metadata, tmp_path):
    """Test round trip of a table that has column units"""
    filename = tmp_path / "test_votable.parq"

    modified_input = copy.deepcopy(input_table)
    modified_input["mass"].unit = u.kg

    write_parquet_votable(
        modified_input,
        filename,
        metadata=column_metadata,
        overwrite_metadata=overwrite_metadata,
    )

    with pytest.warns(AstropyUserWarning, match="No table::len"):
        loaded_table = read_parquet_votable(filename)

    # Compare data content
    assert np.all(modified_input == loaded_table)

    assert (
        modified_input["mass"].unit == loaded_table["mass"].unit
    ) is not overwrite_metadata
    assert loaded_table["sfr"].unit == u.solMass / u.yr


@pytest.mark.xfail(reason="TODO fix: existing column metadata is ignored")
def test_parquet_votable_input_column_metadata(tmp_path):
    """Test preservation of column metadata"""
    filename = tmp_path / "test_votable.parq"

    modified_input = copy.deepcopy(input_table)
    modified_input["mass"].meta["foo"] = "bar"
    write_parquet_votable(modified_input, filename, metadata=column_metadata)

    with pytest.warns(AstropyUserWarning, match="No table::len"):
        loaded_table = read_parquet_votable(filename)

    assert "foo" in loaded_table["mass"].meta
    assert loaded_table["mass"].meta["foo"] == "bar"


@pytest.mark.xfail(reason="TODO fix: existing metadata is ignored")
def test_parquet_votable_input_metadata(tmp_path):
    """Test preservation of table metadata"""
    filename = tmp_path / "test_votable.parq"

    modified_input = copy.deepcopy(input_table)
    modified_input.meta["table_foo"] = "table_bar"
    write_parquet_votable(modified_input, filename, metadata=column_metadata)

    with pytest.warns(AstropyUserWarning, match="No table::len"):
        loaded_table = read_parquet_votable(filename)

    assert "table_foo" in loaded_table.meta
    assert loaded_table.meta["table_foo"] == "table_bar"


def test_compare_parquet_votable(tmp_path):
    """'parquet.votable' preserves column units and metadata unlike 'parquet'"""

    filename = tmp_path / "test_votable.parq"
    write_parquet_votable(input_table, filename, metadata=column_metadata)

    with pytest.warns(AstropyUserWarning, match="No table::len"):
        parquet_table = Table.read(filename, format="parquet")
        parquet_votable = Table.read(filename, format="parquet.votable")

    assert len(parquet_table["sfr"].meta) == 0
    assert parquet_table["sfr"].unit is None

    assert parquet_votable["sfr"].meta["ucd"] == "phys.SFR"
    assert parquet_votable["sfr"].unit == u.solMass / u.yr


def test_write_from_votable_existing_metadata(tmp_path):
    """Read a votable into a Table and write it out as 'parquet.votable' preserving metadata"""
    output_filename = tmp_path / "test_votable.parq"
    input_data = get_pkg_data_filename("data/gaia_source_dr3_select_1_result.vot")

    input_t = Table.read(input_data)

    with pytest.raises(NotImplementedError):
        # Write out the VOTable into a parquet, preserving the metadata
        # without providing an extra dictionary
        input_t.write(output_filename, format="parquet.votable")

        loaded_table = Table.read(output_filename, format="parquet.votable")

        assert np.all(input_t == loaded_table)
        for attr in ["meta", "unit", "dtype", "description"]:
            assert getattr(input_t["ra"], attr) == getattr(loaded_table["ra"], attr)


# Not sure this use case should be supported at all
@pytest.mark.xfail(reason="TODO fix: metadata could be overwritten")
def test_write_from_votable_overriding_metadata(tmp_path):
    """Read a votable into a Table and write it out as 'parquet.votable' preserving metadata"""
    output_filename = tmp_path / "test_votable.parq"
    input_data = get_pkg_data_filename("data/gaia_source_dr3_select_1_result.vot")

    input_t = Table.read(input_data)

    new_metadata = {
        "SOURCE_ID": {
            "description": "Updated description",
        },
    }
    # Write out the VOTable into a parquet, preserving the metadata
    # without providing an extra dictionary
    input_t.write(
        output_filename,
        format="parquet.votable",
        metadata=new_metadata,
        overwrite_metadata=True,
    )

    loaded_table = Table.read(output_filename, format="parquet.votable")

    assert loaded_table["SOURCE_ID"].description == "Updated description"


def test_read_write_existing(tmp_path):
    """Test writing an existing file without overwriting."""
    filename = tmp_path / "test_votable.parquet"
    with open(filename, "w"):
        # create empty file
        pass

    with pytest.raises(OSError, match=_NOT_OVERWRITING_MSG_MATCH):
        write_parquet_votable(input_table, filename, metadata=column_metadata)
