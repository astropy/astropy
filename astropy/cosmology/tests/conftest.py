# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Configure the tests for :mod:`astropy.cosmology`."""

##############################################################################
# IMPORTS

# STDLIB
import json
import os

import pytest

import astropy.cosmology.units as cu
import astropy.units as u
from astropy.cosmology import core
from astropy.cosmology.core import Cosmology
from astropy.utils.misc import isiterable

###############################################################################
# FUNCTIONS


def read_json(filename, **kwargs):
    """Read JSON.

    Parameters
    ----------
    filename : str
    **kwargs
        Keyword arguments into :meth:`~astropy.cosmology.Cosmology.from_format`

    Returns
    -------
    `~astropy.cosmology.Cosmology` instance

    """
    # read
    if isinstance(filename, (str, bytes, os.PathLike)):
        with open(filename, "r") as file:
            data = file.read()
    else:  # file-like : this also handles errors in dumping
        data = filename.read()

    mapping = json.loads(data)  # parse json mappable to dict

    # deserialize Quantity
    with u.add_enabled_units(cu.redshift):
        for k, v in mapping.items():
            if isinstance(v, dict) and "value" in v and "unit" in v:
                mapping[k] = u.Quantity(v["value"], v["unit"])
        for k, v in mapping.get("meta", {}).items():  # also the metadata
            if isinstance(v, dict) and "value" in v and "unit" in v:
                mapping["meta"][k] = u.Quantity(v["value"], v["unit"])

    return Cosmology.from_format(mapping, **kwargs)


def write_json(cosmology, file, *, overwrite=False):
    """Write Cosmology to JSON.

    Parameters
    ----------
    cosmology : `astropy.cosmology.Cosmology` subclass instance
    file : path-like or file-like
    overwrite : bool (optional, keyword-only)
    """
    data = cosmology.to_format("mapping")  # start by turning into dict
    data["cosmology"] = data["cosmology"].__qualname__

    # serialize Quantity
    for k, v in data.items():
        if isinstance(v, u.Quantity):
            data[k] = {"value": v.value.tolist(), "unit": str(v.unit)}
    for k, v in data.get("meta", {}).items():  # also serialize the metadata
        if isinstance(v, u.Quantity):
            data["meta"][k] = {"value": v.value.tolist(), "unit": str(v.unit)}

    # check that file exists and whether to overwrite.
    if os.path.exists(file) and not overwrite:
        raise IOError(f"{file} exists. Set 'overwrite' to write over.")
    with open(file, "w") as write_file:
        json.dump(data, write_file)


def json_identify(origin, filepath, fileobj, *args, **kwargs):
    return filepath is not None and filepath.endswith(".json")


###############################################################################
# FIXTURES

@pytest.fixture
def clean_registry():
    # TODO! with monkeypatch instead for thread safety.
    ORIGINAL_COSMOLOGY_CLASSES = core._COSMOLOGY_CLASSES
    core._COSMOLOGY_CLASSES = {}  # set as empty dict

    yield core._COSMOLOGY_CLASSES

    core._COSMOLOGY_CLASSES = ORIGINAL_COSMOLOGY_CLASSES
