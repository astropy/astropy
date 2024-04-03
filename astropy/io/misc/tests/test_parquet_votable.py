import numpy as np
import pytest

from astropy.io.misc.parquet import read_parquet_vot, write_parquet_vot
from astropy.table import Table
from astropy.utils.exceptions import AstropyUserWarning

# Skip all tests in this file if we cannot import pyarrow
pyarrow = pytest.importorskip("pyarrow")


def test_parquet_votable(tmp_path):
    """Test writing and reading a votable into a parquet file"""
    filename = tmp_path / "test_votable.parq"

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

    write_parquet_vot(input_table, column_metadata, filename)

    with pytest.warns(AstropyUserWarning, match="No table::len"):
        loaded_table = read_parquet_vot(filename)

    assert np.all(input_table == loaded_table)
