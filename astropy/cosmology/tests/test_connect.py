# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy
import inspect
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
readwrite_formats = ["json"]
tofrom_formats = [("mapping", dict)]  # (format, data type)


###############################################################################
# Setup

def read_json(filename, **kwargs):
    with open(filename, "r") as file:
        data = file.read()
    mapping = json.loads(data)  # parse json mappable to dict
    # deserialize Quantity
    for k, v in mapping.items():
        if isinstance(v, dict) and "value" in v and "unit" in v:
            mapping[k] = u.Quantity(v["value"], v["unit"])
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
    io_registry.unregister_reader("json", Cosmology)
    io_registry.unregister_writer("json", Cosmology)
    io_registry.unregister_identifier("json", Cosmology)


###############################################################################
# Tests

class TestReadWriteCosmology:

    @pytest.mark.parametrize("format", readwrite_formats)
    def test_write_methods_have_explicit_kwarg_overwrite(self, format):
        writer = io_registry.get_writer(format, Cosmology)
        # test in signature
        sig = inspect.signature(writer)
        assert "overwrite" in sig.parameters

        # also in docstring
        assert "overwrite : bool" in writer.__doc__

    @pytest.mark.parametrize("format", readwrite_formats)
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

    @pytest.mark.parametrize("format", readwrite_formats)
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

    @pytest.mark.parametrize("format", readwrite_formats)
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


class TestCosmologyToFromFormat:
    """Test methods ``astropy.cosmology.Cosmology.to/from_format``."""

    @pytest.mark.parametrize("format_type", tofrom_formats)
    @pytest.mark.parametrize("instance", cosmo_instances)
    def test_format_complete_info(self, instance, format_type):
        """Read tests happen later."""
        format, objtype = format_type
        cosmo = getattr(cosmology.realizations, instance)

        # test to_format
        obj = cosmo.to_format(format)
        assert isinstance(obj, objtype)

        # test from_format
        got = Cosmology.from_format(obj, format=format)
        # and autodetect
        got2 = Cosmology.from_format(obj)

        assert got2 == got  # internal consistency
        assert got == cosmo  # external consistency
        assert got.meta == cosmo.meta

    @pytest.mark.parametrize("format_type", tofrom_formats)
    @pytest.mark.parametrize("instance", cosmo_instances)
    def test_from_subclass_complete_info(self, instance, format_type):
        """
        Test transforming an instance and parsing from that class, when there's
        full information available.
        """
        format, objtype = format_type
        cosmo = getattr(cosmology.realizations, instance)

        # test to_format
        obj = cosmo.to_format(format)
        assert isinstance(obj, objtype)

        # read with the same class that wrote.
        got = cosmo.__class__.from_format(obj, format=format)
        got2 = Cosmology.from_format(obj)  # and autodetect

        assert got2 == got  # internal consistency
        assert got == cosmo  # external consistency
        assert got.meta == cosmo.meta

        # this should be equivalent to
        got = Cosmology.from_format(obj, format=format, cosmology=cosmo.__class__)
        assert got == cosmo
        assert got.meta == cosmo.meta

        # and also
        got = Cosmology.from_format(obj, format=format, cosmology=cosmo.__class__.__qualname__)
        assert got == cosmo
        assert got.meta == cosmo.meta

    @pytest.mark.parametrize("instance", cosmo_instances)
    def test_from_subclass_partial_info(self, instance):
        """
        Test writing from an instance and reading from that class.
        This requires partial information.

        .. todo::

            generalize over all formats for this test.
        """
        format, objtype = ("mapping", dict)
        cosmo = getattr(cosmology.realizations, instance)

        # test to_format
        obj = cosmo.to_format(format)
        assert isinstance(obj, objtype)

        # partial information
        tempobj = copy.deepcopy(obj)
        del tempobj["cosmology"]
        del tempobj["Tcmb0"]

        # read with the same class that wrote fills in the missing info with
        # the default value
        got = cosmo.__class__.from_format(tempobj, format=format)
        got2 = Cosmology.from_format(tempobj, format=format, cosmology=cosmo.__class__)
        got3 = Cosmology.from_format(tempobj, format=format, cosmology=cosmo.__class__.__qualname__)

        assert (got == got2) and (got2 == got3)  # internal consistency

        # not equal, because Tcmb0 is changed
        assert got != cosmo
        assert got.Tcmb0 == cosmo.__class__._init_signature.parameters["Tcmb0"].default
        assert got.clone(name=cosmo.name, Tcmb0=cosmo.Tcmb0.value) == cosmo
        # but the metadata is the same
        assert got.meta == cosmo.meta

    @pytest.mark.parametrize("format_type", tofrom_formats)
    @pytest.mark.parametrize("instance", cosmo_instances)
    def test_reader_class_mismatch(self, instance, format_type):
        """Test when the reader class doesn't match the object."""
        format, objtype = format_type
        cosmo = getattr(cosmology.realizations, instance)

        # test to_format
        obj = cosmo.to_format(format)
        assert isinstance(obj, objtype)

        # class mismatch
        with pytest.raises(TypeError, match="missing 1 required"):
            w0wzCDM.from_format(obj, format=format)

        with pytest.raises(TypeError, match="missing 1 required"):
            Cosmology.from_format(obj, format=format, cosmology=w0wzCDM)

        # when specifying the class
        with pytest.raises(ValueError, match="`cosmology` must be either"):
            w0wzCDM.from_format(obj, format=format, cosmology="FlatLambdaCDM")
