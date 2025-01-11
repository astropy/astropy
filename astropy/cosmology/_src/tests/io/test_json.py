# Licensed under a 3-clause BSD style license - see LICENSE.rst

import json
import os
from pathlib import Path

import pytest

import astropy.units as u
from astropy.cosmology import Cosmology
from astropy.cosmology import units as cu
from astropy.cosmology.io import readwrite_registry

from .base import ReadWriteDirectTestBase, ReadWriteTestMixinBase

###############################################################################


def read_json(filename, **kwargs):
    """Read JSON.

    Parameters
    ----------
    filename : str | bytes | os.PathLike
    **kwargs
        Keyword arguments into :meth:`~astropy.cosmology.Cosmology.from_format`

    Returns
    -------
    `~astropy.cosmology.Cosmology` instance
    """
    # read
    if isinstance(filename, (str, bytes, os.PathLike)):
        data = Path(filename).read_text()
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

    return Cosmology.from_format(mapping, format="mapping", **kwargs)


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
    file = Path(file)
    if file.exists() and not overwrite:
        raise OSError(f"{file} exists. Set 'overwrite' to write over.")
    with file.open("w") as write_file:
        json.dump(data, write_file)


def json_identify(origin, filepath, fileobj, *args, **kwargs):
    return filepath is not None and filepath.endswith(".json")


###############################################################################


class ReadWriteJSONTestMixin(ReadWriteTestMixinBase):
    """
    Tests for a Cosmology[Read/Write] with ``format="json"``.
    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmology`` for an example.
    """

    @pytest.fixture(scope="class", autouse=True)
    def register_and_unregister_json(self):
        """Setup & teardown for JSON read/write tests."""
        # Register
        readwrite_registry.register_reader("json", Cosmology, read_json, force=True)
        readwrite_registry.register_writer("json", Cosmology, write_json, force=True)
        readwrite_registry.register_identifier(
            "json", Cosmology, json_identify, force=True
        )

        yield  # Run all tests in class

        # Unregister
        readwrite_registry.unregister_reader("json", Cosmology)
        readwrite_registry.unregister_writer("json", Cosmology)
        readwrite_registry.unregister_identifier("json", Cosmology)

    # ========================================================================

    def test_readwrite_json_subclass_partial_info(
        self, cosmo_cls, cosmo, read, write, tmp_path, add_cu
    ):
        """
        Test writing from an instance and reading from that class.
        This works with missing information.
        """
        fp = tmp_path / "test_readwrite_json_subclass_partial_info.json"

        # test write
        cosmo.write(fp, format="json")

        # partial information
        with fp.open() as file:
            L = file.readline()
        L = (
            L[: L.index('"cosmology":')] + L[L.index(", ") + 2 :]
        )  # remove cosmology  : #203
        i = L.index('"Tcmb0":')  # delete Tcmb0
        L = (
            L[:i] + L[L.index(", ", L.index(", ", i) + 1) + 2 :]
        )  # second occurrence  : #203

        tempfname = tmp_path / f"{cosmo.name}_temp.json"
        tempfname.write_text("".join(L))
        # read with the same class that wrote fills in the missing info with
        # the default value
        got = cosmo_cls.read(tempfname, format="json")
        got2 = read(tempfname, format="json", cosmology=cosmo_cls)
        got3 = read(tempfname, format="json", cosmology=cosmo_cls.__qualname__)

        assert (got == got2) and (got2 == got3)  # internal consistency

        # not equal, because Tcmb0 is changed, which also changes m_nu
        assert got != cosmo
        assert got.Tcmb0 == cosmo_cls.parameters["Tcmb0"].default
        assert got.clone(name=cosmo.name, Tcmb0=cosmo.Tcmb0, m_nu=cosmo.m_nu) == cosmo
        # but the metadata is the same
        assert got.meta == cosmo.meta


class TestReadWriteJSON(ReadWriteDirectTestBase, ReadWriteJSONTestMixin):
    """
    Directly test ``read/write_json``.
    These are not public API and are discouraged from use, in favor of
    ``Cosmology.read/write(..., format="json")``, but should be
    tested regardless b/c they are used internally.
    """

    def setup_class(self):
        self.functions = {"read": read_json, "write": write_json}
