# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy
import os
import json

import pytest

import numpy as np

import astropy.units as u
from astropy import cosmology
from astropy.cosmology import Cosmology
from astropy.cosmology.connect import CosmologyRead
from astropy.cosmology.core import _COSMOLOGY_CLASSES, Cosmology
from astropy.io import registry as io_registry
from astropy.utils.exceptions import AstropyUserWarning

cosmo_instances = cosmology.parameters.available
save_formats = ["json"]


###############################################################################
# Setup

def from_mapping(map, **kwargs):
    params = copy.deepcopy(map)  # so can pop

    # get cosmology. if string, parse to class
    cosmology = params.pop("cosmology")
    if isinstance(cosmology, str):
        cosmology = _COSMOLOGY_CLASSES[cosmology]

    # select arguments from mapping that are in the cosmo's signature.
    ba = cosmology._init_signature.bind_partial()  # blank set of args
    ba.apply_defaults()  # fill in the defaults
    for k in cosmology._init_signature.parameters.keys():  # iter thru sig
        if k in params:  # transfer argument, if in params
            ba.arguments[k] = params.pop(k)

    return cosmology(*ba.args, **ba.kwargs)


def to_mapping(cosmology, *args, **kwargs):
    m = {}
    m["cosmology"] = cosmology.__class__
    m.update({k: v for k, v in cosmology._init_arguments.items()
              if k not in ("meta",)})
    m["meta"] = copy.deepcopy(cosmology.meta)

    return m


def read_json(filename, key=None, **kwargs):
    with open(filename, "r") as file:
        data = file.read()
    mapping = json.loads(data)  # parse json mappable to dict
    # deserialize Quantity
    for k, v in mapping.items():
        if isinstance(v, dict) and "value" in v and "unit" in v:
            mapping[k] = u.Quantity(v["value"], v["unit"])
    return from_mapping(mapping)


def write_json(cosmology, file, *, overwrite=False, **kwargs):
    data = to_mapping(cosmology)  # start by turning into dict
    data["cosmology"] = data["cosmology"].__name__  # change class field to str
    # serialize Quantity
    for k, v in data.items():
        if isinstance(v, u.Quantity):
            data[k] = {"value": v.value.tolist(),
                       "unit": str(v.unit)}

    # check that file exists and whether to overwrite.
    if os.path.exists(file) and not overwrite:
        raise IOError(f"{file} exists. Set 'overwrite' to write over.")
    with open(file, "w") as write_file:
        json.dump(data, write_file)


def json_identify(origin, filepath, fileobj, *args, **kwargs):
    return filepath is not None and filepath.endswith(".json")
    

def setup_module(module):
    """Setup module for tests."""
    io_registry.register_reader("json", Cosmology, read_json)
    io_registry.register_writer("json", Cosmology, write_json)
    io_registry.register_identifier("json", Cosmology, json_identify)


def teardown_module(module):
    """clean up module after tests."""
    with pytest.warns(FutureWarning):  # idk
        io_registry.unregister_reader("json", Cosmology)
        io_registry.unregister_writer("json", Cosmology)
        io_registry.unregister_identifier("json", Cosmology)


###############################################################################
# Tests

class TestReadWriteCosmology:

    def test_instantiate_read(self):
        # no error on base class
        assert isinstance(Cosmology.read, CosmologyRead)

        # Warns when initialized. Cannot be used.
        with pytest.warns(AstropyUserWarning):
            assert cosmology.realizations.Planck18.read is NotImplemented

    @pytest.mark.parametrize("format", save_formats)
    @pytest.mark.parametrize("instance", cosmo_instances)
    def test_write_then_read_file(self, tmpdir, instance, format):
        cosmo = getattr(cosmology.realizations, instance)
        fname = tmpdir / f"{instance}.{format}"

        cosmo.write(str(fname), format=format)

        # Also test kwarg "overwrite"
        assert os.path.exists(str(fname))  # file exists
        with pytest.raises(IOError):
            cosmo.write(str(fname), format=format, overwrite=False)

        assert os.path.exists(str(fname))  # overwrite file existing file
        cosmo.write(str(fname), format=format, overwrite=True)

        # Read back
        got = Cosmology.read(tmpdir / f"{instance}.{format}", format=format)

        assert got == cosmo
        assert got.meta == cosmo.meta
