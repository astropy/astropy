# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy
import json
import os

import pytest

import numpy as np

import astropy.units as u
from astropy import cosmology
from astropy.cosmology import Cosmology, w0wzCDM
from astropy.cosmology.connect import CosmologyRead
from astropy.cosmology.core import _COSMOLOGY_CLASSES, Cosmology
from astropy.io import registry as io_registry

cosmo_instances = cosmology.parameters.available
save_formats = ["json"]


###############################################################################
# Setup

def from_mapping(map, **kwargs):
    params = copy.deepcopy(map)  # so can pop

    # get cosmology. if string, parse to class
    # 1st from 'kwargs'. Allows for override of the cosmology, if on file.
    # 2nd from params. This MUST have the cosmology if 'kwargs' did not.
    if "cosmology" in kwargs:
        cosmology = kwargs.pop("cosmology")
    else:
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
    return from_mapping(mapping, **kwargs)


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

    @pytest.mark.parametrize("format", save_formats)
    @pytest.mark.parametrize("instance", cosmo_instances)
    def test_complete_info(self, tmpdir, instance, format):
        """
        Test writing from an instance and reading from the base class.
        This requires full information.
        """
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
        got = Cosmology.read(fname, format=format)

        assert got == cosmo
        assert got.meta == cosmo.meta

    @pytest.mark.parametrize("format", save_formats)
    @pytest.mark.parametrize("instance", cosmo_instances)
    def test_from_subclass_complete_info(self, tmpdir, instance, format):
        """
        Test writing from an instance and reading from that class, when there's
        full information saved.
        """
        cosmo = getattr(cosmology.realizations, instance)
        fname = tmpdir / f"{instance}.{format}"
        cosmo.write(str(fname), format=format)

        # read with the same class that wrote.
        got = cosmo.__class__.read(fname, format=format)
        assert got == cosmo
        assert got.meta == cosmo.meta

        # this should be equivalent to
        got = Cosmology.read(fname, format=format, cosmology=cosmo.__class__)
        assert got == cosmo
        assert got.meta == cosmo.meta

        # and also
        got = Cosmology.read(fname, format=format, cosmology=cosmo.__class__.__qualname__)
        assert got == cosmo
        assert got.meta == cosmo.meta

    @pytest.mark.parametrize("instance", cosmo_instances)
    def test_from_subclass_partial_info(self, tmpdir, instance):
        """
        Test writing from an instance and reading from that class.
        This requires partial information.

        .. todo::

            generalize over all save formats for this test.
        """
        format = "json"
        cosmo = getattr(cosmology.realizations, instance)
        fname = tmpdir / f"{instance}.{format}"

        cosmo.write(str(fname), format=format)

        # partial information
        with open(fname, "r") as file:
            L = file.readlines()
        L[0] = L[0][:L[0].index('"cosmology":')]+L[0][L[0].index(', ')+2:]
        i = L[0].index('"Tcmb0":')  # delete Tcmb0
        L[0] = L[0][:i] + L[0][L[0].index(', ', i)+2:]

        tempfname = tmpdir / f"{instance}_temp.{format}"
        with open(tempfname, "w") as file:
            file.writelines(L)

        # read with the same class that wrote fills in the missing info with
        # the default value
        got = cosmo.__class__.read(tempfname, format=format)
        got2 = Cosmology.read(tempfname, format=format, cosmology=cosmo.__class__)
        got3 = Cosmology.read(tempfname, format=format, cosmology=cosmo.__class__.__qualname__)

        assert (got == got2) and (got2 == got3)  # internal consistency

        # not equal, because Tcmb0 is changed
        assert got != cosmo
        assert got.Tcmb0 == cosmo.__class__._init_signature.parameters["Tcmb0"].default
        assert got.clone(name=cosmo.name, Tcmb0=cosmo.Tcmb0.value) == cosmo
        # but the metadata is the same
        assert got.meta == cosmo.meta

    @pytest.mark.parametrize("format", save_formats)
    @pytest.mark.parametrize("instance", cosmo_instances)
    def test_reader_class_mismatch(self, tmpdir, instance, format):
        """Test when the reader class doesn't match the file."""
        cosmo = getattr(cosmology.realizations, instance)
        fname = tmpdir / f"{instance}.{format}"
        cosmo.write(str(fname), format=format)

        # class mismatch
        # when reading directly
        with pytest.raises(TypeError, match="missing 1 required"):
            w0wzCDM.read(fname, format=format)

        with pytest.raises(TypeError, match="missing 1 required"):
            Cosmology.read(fname, format=format, cosmology=w0wzCDM)

        # when specifying the class
        with pytest.raises(ValueError, match="`cosmology` must be either"):
            w0wzCDM.read(fname, format=format, cosmology="FlatLambdaCDM")
