# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from astropy.cosmology._src.core import _COSMOLOGY_CLASSES
from astropy.cosmology._src.io.builtin.ecsv import read_ecsv, write_ecsv
from astropy.table import QTable, Table, vstack

from .base import ReadWriteDirectTestBase, ReadWriteTestMixinBase

###############################################################################


class ReadWriteECSVTestMixin(ReadWriteTestMixinBase):
    """
    Tests for a Cosmology[Read/Write] with ``format="ascii.ecsv"``.
    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmology`` for an example.
    """

    def test_to_ecsv_bad_index(self, read, write, tmp_path):
        """Test if argument ``index`` is incorrect"""
        fp = tmp_path / "test_to_ecsv_bad_index.ecsv"

        write(fp, format="ascii.ecsv")

        # single-row table and has a non-0/None index
        with pytest.raises(IndexError, match="index 2 out of range"):
            read(fp, index=2, format="ascii.ecsv")

        # string index where doesn't match
        with pytest.raises(KeyError, match="No matches found for key"):
            read(fp, index="row 0", format="ascii.ecsv")

    # -----------------------

    def test_to_ecsv_failed_cls(self, write, tmp_path):
        """Test failed table type."""
        fp = tmp_path / "test_to_ecsv_failed_cls.ecsv"

        with pytest.raises(TypeError, match="'cls' must be"):
            write(fp, format="ascii.ecsv", cls=list)

    @pytest.mark.parametrize("tbl_cls", [QTable, Table])
    def test_to_ecsv_cls(self, write, tbl_cls, tmp_path):
        fp = tmp_path / "test_to_ecsv_cls.ecsv"
        write(fp, format="ascii.ecsv", cls=tbl_cls)

    # -----------------------

    @pytest.mark.parametrize("in_meta", [True, False])
    def test_to_ecsv_in_meta(self, cosmo_cls, write, in_meta, tmp_path, add_cu):
        """Test where the cosmology class is placed."""
        fp = tmp_path / "test_to_ecsv_in_meta.ecsv"
        write(fp, format="ascii.ecsv", cosmology_in_meta=in_meta)

        # if it's in metadata, it's not a column. And vice versa.
        tbl = QTable.read(fp)
        if in_meta:
            assert tbl.meta["cosmology"] == cosmo_cls.__qualname__
            assert "cosmology" not in tbl.colnames  # not also a column
        else:
            assert tbl["cosmology"][0] == cosmo_cls.__qualname__
            assert "cosmology" not in tbl.meta

    # -----------------------

    def test_readwrite_ecsv_instance(
        self, cosmo_cls, cosmo, read, write, tmp_path, add_cu
    ):
        """Test cosmology -> ascii.ecsv -> cosmology."""
        fp = tmp_path / "test_readwrite_ecsv_instance.ecsv"

        # ------------
        # To Table

        write(fp, format="ascii.ecsv")

        # some checks on the saved file
        tbl = QTable.read(fp)
        assert tbl.meta["cosmology"] == cosmo_cls.__qualname__
        assert tbl["name"] == cosmo.name

        # ------------
        # From Table

        tbl["mismatching"] = "will error"
        tbl.write(fp, format="ascii.ecsv", overwrite=True)

        # tests are different if the last argument is a **kwarg
        if cosmo._init_has_kwargs:
            got = read(fp, format="ascii.ecsv")

            assert got.__class__ is cosmo_cls
            assert got.name == cosmo.name
            assert "mismatching" not in got.meta

            return  # don't continue testing

        # read with mismatching parameters errors
        with pytest.raises(TypeError, match="there are unused parameters"):
            read(fp, format="ascii.ecsv")

        # unless mismatched are moved to meta
        got = read(fp, format="ascii.ecsv", move_to_meta=True)
        assert got == cosmo
        assert got.meta["mismatching"] == "will error"

        # it won't error if everything matches up
        tbl.remove_column("mismatching")
        tbl.write(fp, format="ascii.ecsv", overwrite=True)
        got = read(fp, format="ascii.ecsv")
        assert got == cosmo

        # and it will also work if the cosmology is a class
        # Note this is not the default output of ``write``.
        tbl.meta["cosmology"] = _COSMOLOGY_CLASSES[tbl.meta["cosmology"]]
        got = read(fp, format="ascii.ecsv")
        assert got == cosmo

        # also it auto-identifies 'format'
        got = read(fp)
        assert got == cosmo

    def test_readwrite_ecsv_renamed_columns(
        self, cosmo_cls, cosmo, read, write, tmp_path, add_cu
    ):
        """Test rename argument to read/write."""
        fp = tmp_path / "test_readwrite_ecsv_rename.ecsv"
        rename = {"name": "cosmo_name"}

        write(fp, format="ascii.ecsv", rename=rename)

        tbl = QTable.read(fp, format="ascii.ecsv")

        assert "name" not in tbl.colnames
        assert "cosmo_name" in tbl.colnames

        # Errors if reading
        with pytest.raises(
            TypeError, match="there are unused parameters {'cosmo_name':"
        ):
            read(fp)

        # Roundtrips
        inv_rename = {v: k for k, v in rename.items()}
        got = read(fp, rename=inv_rename)
        assert got == cosmo

    def test_readwrite_ecsv_subclass_partial_info(
        self, cosmo_cls, cosmo, read, write, tmp_path, add_cu
    ):
        """
        Test writing from an instance and reading from that class.
        This works with missing information.
        """
        fp = tmp_path / "test_read_ecsv_subclass_partial_info.ecsv"

        # test write
        write(fp, format="ascii.ecsv")

        # partial information
        tbl = QTable.read(fp)
        tbl.meta.pop("cosmology", None)
        del tbl["Tcmb0"]
        tbl.write(fp, overwrite=True)

        # read with the same class that wrote fills in the missing info with
        # the default value
        got = cosmo_cls.read(fp, format="ascii.ecsv")
        got2 = read(fp, format="ascii.ecsv", cosmology=cosmo_cls)
        got3 = read(fp, format="ascii.ecsv", cosmology=cosmo_cls.__qualname__)

        assert (got == got2) and (got2 == got3)  # internal consistency

        # not equal, because Tcmb0 is changed, which also changes m_nu
        assert got != cosmo
        assert got.Tcmb0 == cosmo_cls.parameters["Tcmb0"].default
        assert got.clone(name=cosmo.name, Tcmb0=cosmo.Tcmb0, m_nu=cosmo.m_nu) == cosmo
        # but the metadata is the same
        assert got.meta == cosmo.meta

    def test_readwrite_ecsv_mutlirow(self, cosmo, read, write, tmp_path, add_cu):
        """Test if table has multiple rows."""
        fp = tmp_path / "test_readwrite_ecsv_mutlirow.ecsv"

        # Make
        cosmo1 = cosmo.clone(name="row 0")
        cosmo2 = cosmo.clone(name="row 2")
        tbl = vstack(
            [c.to_format("astropy.table") for c in (cosmo1, cosmo, cosmo2)],
            metadata_conflicts="silent",
        )
        tbl.write(fp, format="ascii.ecsv")

        # ------------
        # From Table

        # it will error on a multi-row table
        with pytest.raises(ValueError, match="need to select a specific row"):
            read(fp, format="ascii.ecsv")

        # unless the index argument is provided
        got = read(fp, index=1, format="ascii.ecsv")
        assert got == cosmo

        # the index can be a string
        got = read(fp, index=cosmo.name, format="ascii.ecsv")
        assert got == cosmo

        # it's better if the table already has an index
        # this will be identical to the previous ``got``
        tbl.add_index("name")
        got2 = read(fp, index=cosmo.name, format="ascii.ecsv")
        assert got2 == cosmo


class TestReadWriteECSV(ReadWriteDirectTestBase, ReadWriteECSVTestMixin):
    """
    Directly test ``read/write_ecsv``.
    These are not public API and are discouraged from use, in favor of
    ``Cosmology.read/write(..., format="ascii.ecsv")``, but should be
    tested regardless b/c they are used internally.
    """

    def setup_class(self):
        self.functions = {"read": read_ecsv, "write": write_ecsv}
