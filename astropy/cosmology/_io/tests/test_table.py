# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from astropy.cosmology import Cosmology
from astropy.cosmology._io.table import from_table, to_table
from astropy.cosmology.core import _COSMOLOGY_CLASSES
from astropy.table import QTable, Table, vstack

from .base import ToFromDirectTestBase, ToFromTestMixinBase

###############################################################################


class ToFromTableTestMixin(ToFromTestMixinBase):
    """
    Tests for a Cosmology[To/From]Format with ``format="astropy.table"``.
    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmology`` for an example.
    """

    def test_to_table_bad_index(self, from_format, to_format):
        """Test if argument ``index`` is incorrect"""
        tbl = to_format("astropy.table")

        # single-row table and has a non-0/None index
        with pytest.raises(IndexError, match="index 2 out of range"):
            from_format(tbl, index=2, format="astropy.table")

        # string index where doesn't match
        with pytest.raises(KeyError, match="No matches found for key"):
            from_format(tbl, index="row 0", format="astropy.table")

    # -----------------------

    def test_to_table_failed_cls(self, to_format):
        """Test failed table type."""
        with pytest.raises(TypeError, match="'cls' must be"):
            to_format("astropy.table", cls=list)

    @pytest.mark.parametrize("tbl_cls", [QTable, Table])
    def test_to_table_cls(self, to_format, tbl_cls):
        tbl = to_format("astropy.table", cls=tbl_cls)
        assert isinstance(tbl, tbl_cls)  # test type

    # -----------------------

    @pytest.mark.parametrize("in_meta", [True, False])
    def test_to_table_in_meta(self, cosmo_cls, to_format, in_meta):
        """Test where the cosmology class is placed."""
        tbl = to_format("astropy.table", cosmology_in_meta=in_meta)

        # if it's in metadata, it's not a column. And vice versa.
        if in_meta:
            assert tbl.meta["cosmology"] == cosmo_cls.__qualname__
            assert "cosmology" not in tbl.colnames  # not also a column
        else:
            assert tbl["cosmology"][0] == cosmo_cls.__qualname__
            assert "cosmology" not in tbl.meta

    # -----------------------

    def test_to_table(self, cosmo_cls, cosmo, to_format):
        """Test cosmology -> astropy.table."""
        tbl = to_format("astropy.table")

        # Test properties of Table.
        assert isinstance(tbl, QTable)
        assert tbl.meta["cosmology"] == cosmo_cls.__qualname__
        assert tbl["name"] == cosmo.name
        assert tbl.indices  # indexed

        # Test each Parameter column has expected information.
        for n, P in cosmo_cls.parameters.items():
            col = tbl[n]  # Column

            # Compare the two
            assert col.info.name == P.name
            assert col.info.description == P.__doc__
            assert col.info.meta == (cosmo.meta.get(n) or {})

    # -----------------------

    def test_from_not_table(self, cosmo, from_format):
        """Test not passing a Table to the Table parser."""
        with pytest.raises((TypeError, ValueError)):
            from_format("NOT A TABLE", format="astropy.table")

    def test_tofrom_table_instance(self, cosmo_cls, cosmo, from_format, to_format):
        """Test cosmology -> astropy.table -> cosmology."""
        tbl = to_format("astropy.table")

        # add information
        tbl["mismatching"] = "will error"

        # tests are different if the last argument is a **kwarg
        if tuple(cosmo._init_signature.parameters.values())[-1].kind == 4:
            got = from_format(tbl, format="astropy.table")

            assert got.__class__ is cosmo_cls
            assert got.name == cosmo.name
            assert "mismatching" not in got.meta

            return  # don't continue testing

        # read with mismatching parameters errors
        with pytest.raises(TypeError, match="there are unused parameters"):
            from_format(tbl, format="astropy.table")

        # unless mismatched are moved to meta
        got = from_format(tbl, format="astropy.table", move_to_meta=True)
        assert got == cosmo
        assert got.meta["mismatching"] == "will error"

        # it won't error if everything matches up
        tbl.remove_column("mismatching")
        got = from_format(tbl, format="astropy.table")
        assert got == cosmo

        # and it will also work if the cosmology is a class
        # Note this is not the default output of ``to_format``.
        tbl.meta["cosmology"] = _COSMOLOGY_CLASSES[tbl.meta["cosmology"]]
        got = from_format(tbl, format="astropy.table")
        assert got == cosmo

        # also it auto-identifies 'format'
        got = from_format(tbl)
        assert got == cosmo

    def test_fromformat_table_subclass_partial_info(
        self, cosmo_cls, cosmo, from_format, to_format
    ):
        """
        Test writing from an instance and reading from that class.
        This works with missing information.
        """
        # test to_format
        tbl = to_format("astropy.table")
        assert isinstance(tbl, QTable)

        # partial information
        tbl.meta.pop("cosmology", None)
        del tbl["Tcmb0"]

        # read with the same class that wrote fills in the missing info with
        # the default value
        got = cosmo_cls.from_format(tbl, format="astropy.table")
        got2 = from_format(tbl, format="astropy.table", cosmology=cosmo_cls)
        got3 = from_format(
            tbl, format="astropy.table", cosmology=cosmo_cls.__qualname__
        )

        assert (got == got2) and (got2 == got3)  # internal consistency

        # not equal, because Tcmb0 is changed, which also changes m_nu
        assert got != cosmo
        assert got.Tcmb0 == cosmo_cls._init_signature.parameters["Tcmb0"].default
        assert got.clone(name=cosmo.name, Tcmb0=cosmo.Tcmb0, m_nu=cosmo.m_nu) == cosmo
        # but the metadata is the same
        assert got.meta == cosmo.meta

    @pytest.mark.parametrize("add_index", [True, False])
    def test_tofrom_table_mutlirow(self, cosmo_cls, cosmo, from_format, add_index):
        """Test if table has multiple rows."""
        # ------------
        # To Table

        cosmo1 = cosmo.clone(name="row 0")
        cosmo2 = cosmo.clone(name="row 2")
        tbl = vstack(
            [c.to_format("astropy.table") for c in (cosmo1, cosmo, cosmo2)],
            metadata_conflicts="silent",
        )

        assert isinstance(tbl, QTable)
        assert tbl.meta["cosmology"] == cosmo_cls.__qualname__
        assert tbl[1]["name"] == cosmo.name

        # whether to add an index. `from_format` can work with or without.
        if add_index:
            tbl.add_index("name", unique=True)

        # ------------
        # From Table

        # it will error on a multi-row table
        with pytest.raises(ValueError, match="need to select a specific row"):
            from_format(tbl, format="astropy.table")

        # unless the index argument is provided
        got = from_format(tbl, index=1, format="astropy.table")
        assert got == cosmo

        # the index can be a string
        got = from_format(tbl, index=cosmo.name, format="astropy.table")
        assert got == cosmo

        # when there's more than one cosmology found
        tbls = vstack([tbl, tbl], metadata_conflicts="silent")
        with pytest.raises(ValueError, match="more than one"):
            from_format(tbls, index=cosmo.name, format="astropy.table")

    def test_tofrom_table_rename(self, cosmo, to_format, from_format):
        """Test renaming columns in row."""
        rename = {"name": "cosmo_name"}
        table = to_format("astropy.table", rename=rename)

        assert "name" not in table.colnames
        assert "cosmo_name" in table.colnames

        # Error if just reading
        with pytest.raises(TypeError, match="there are unused parameters"):
            from_format(table)

        # Roundtrip
        inv_rename = {v: k for k, v in rename.items()}
        got = from_format(table, rename=inv_rename)
        assert got == cosmo

    def test_from_table_renamed_index_column(self, cosmo, to_format, from_format):
        """Test reading from a table with a renamed index column."""
        cosmo1 = cosmo.clone(name="row 0")
        cosmo2 = cosmo.clone(name="row 2")
        tbl = vstack(
            [c.to_format("astropy.table") for c in (cosmo1, cosmo, cosmo2)],
            metadata_conflicts="silent",
        )
        tbl.rename_column("name", "cosmo_name")

        inv_rename = {"cosmo_name": "name"}
        newcosmo = from_format(
            tbl, index="row 0", rename=inv_rename, format="astropy.table"
        )
        assert newcosmo == cosmo1

    @pytest.mark.parametrize("format", [True, False, None, "astropy.table"])
    def test_is_equivalent_to_table(self, cosmo, to_format, format):
        """Test :meth:`astropy.cosmology.Cosmology.is_equivalent`.

        This test checks that Cosmology equivalency can be extended to any
        Python object that can be converted to a Cosmology -- in this case
        a |Table|.
        """
        obj = to_format("astropy.table")
        assert not isinstance(obj, Cosmology)

        is_equiv = cosmo.is_equivalent(obj, format=format)
        assert is_equiv is (format is not False)


class TestToFromTable(ToFromDirectTestBase, ToFromTableTestMixin):
    """Directly test ``to/from_table``."""

    def setup_class(self):
        self.functions = {"to": to_table, "from": from_table}
