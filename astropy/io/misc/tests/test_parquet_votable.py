import numpy as np
import pytest

from astropy.io.misc.parquet import read_parquet_votable, write_parquet_votable
from astropy.table import Table
import astropy.units as u
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

    assert np.all(input_table == loaded_table)


def test_compare_parquet_votable(tmp_path):
    """ Parquet votable preserves column units and metadata"""

    filename = tmp_path / "test_votable.parq"
    write_parquet_votable(input_table, filename, metadata=column_metadata)

    with pytest.warns(AstropyUserWarning, match="No table::len"):
        parquet = Table.read(filename)
        parquet_votable = Table.read(filename, format='parquet.votable')

    assert len(parquet['sfr'].meta) == 0
    assert parquet['sfr'].unit is None

    assert parquet_votable['sfr'].meta['ucd'] == 'phys.SFR'
    assert parquet_votable['sfr'].unit == u.solMass / u.yr


def test_read_write_existing(tmp_path):
    """Test writing an existing file without overwriting."""
    filename = tmp_path / "test_votable.parquet"
    with open(filename, "w"):
        # create empty file
        pass

    with pytest.raises(OSError, match=_NOT_OVERWRITING_MSG_MATCH):
        write_parquet_votable(input_table, filename, metadata=column_metadata)
