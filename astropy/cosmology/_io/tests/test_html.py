# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

import astropy.units as u
from astropy.cosmology._io.html import _FORMAT_TABLE, read_html_table, write_html_table
from astropy.table import QTable, Table, vstack
from astropy.utils.compat.optional_deps import HAS_BS4

from .base import ReadWriteDirectTestBase, ReadWriteTestMixinBase

###############################################################################


class ReadWriteHTMLTestMixin(ReadWriteTestMixinBase):
    """
    Tests for a Cosmology[Read/Write] with ``format="ascii.html"``.
    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must dfine a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmology`` for an example.
    """

    @pytest.mark.skipif(not HAS_BS4, reason="requires beautifulsoup4")
    def test_to_html_table_bad_index(self, read, write, tmp_path):
        """Test if argument ``index`` is incorrect"""
        fp = tmp_path / "test_to_html_table_bad_index.html"

        write(fp, format="ascii.html")

        # single-row table and has a non-0/None index
        with pytest.raises(IndexError, match="index 2 out of range"):
            read(fp, index=2, format="ascii.html")

        # string index where doesn't match
        with pytest.raises(KeyError, match="No matches found for key"):
            read(fp, index="row 0", format="ascii.html")

    # -----------------------

    @pytest.mark.skipif(not HAS_BS4, reason="requires beautifulsoup4")
    def test_to_html_table_failed_cls(self, write, tmp_path):
        """Test failed table type."""
        fp = tmp_path / "test_to_html_table_failed_cls.html"

        with pytest.raises(TypeError, match="'cls' must be"):
            write(fp, format="ascii.html", cls=list)

    @pytest.mark.parametrize("tbl_cls", [QTable, Table])
    @pytest.mark.skipif(not HAS_BS4, reason="requires beautifulsoup4")
    def test_to_html_table_cls(self, write, tbl_cls, tmp_path):
        fp = tmp_path / "test_to_html_table_cls.html"
        write(fp, format="ascii.html", cls=tbl_cls)

    # -----------------------

    @pytest.mark.skipif(not HAS_BS4, reason="requires beautifulsoup4")
    def test_readwrite_html_table_instance(
        self, cosmo_cls, cosmo, read, write, tmp_path, add_cu
    ):
        """Test cosmology -> ascii.html -> cosmology."""
        fp = tmp_path / "test_readwrite_html_table_instance.html"

        # ------------
        # To Table

        write(fp, format="ascii.html")

        # some checks on the saved file
        tbl = QTable.read(fp)
        # assert tbl.meta["cosmology"] == cosmo_cls.__qualname__  # metadata read not implemented
        assert tbl["name"] == cosmo.name

        # ------------
        # From Table

        tbl["mismatching"] = "will error"
        tbl.write(fp, format="ascii.html", overwrite=True)

        # tests are different if the last argument is a **kwarg
        if cosmo._init_has_kwargs:
            got = read(fp, format="ascii.html")

            assert got.__class__ is cosmo_cls
            assert got.name == cosmo.name
            # assert "mismatching" not in got.meta # metadata read not implemented

            return  # don't continue testing

        # read with mismatching parameters errors
        with pytest.raises(TypeError, match="there are unused parameters"):
            read(fp, format="ascii.html")

        # unless mismatched are moved to meta
        got = read(fp, format="ascii.html", move_to_meta=True)
        assert got == cosmo
        # assert got.meta["mismatching"] == "will error" # metadata read not implemented

        # it won't error if everything matches up
        tbl.remove_column("mismatching")
        tbl.write(fp, format="ascii.html", overwrite=True)
        got = read(fp, format="ascii.html")
        assert got == cosmo

        # and it will also work if the cosmology is a class
        # Note this is not the default output of ``write``.
        # tbl.meta["cosmology"] = _COSMOLOGY_CLASSES[tbl.meta["cosmology"]] #
        # metadata read not implemented
        got = read(fp, format="ascii.html")
        assert got == cosmo

        got = read(fp)
        assert got == cosmo

    @pytest.mark.skipif(not HAS_BS4, reason="requires beautifulsoup4")
    def test_rename_html_table_columns(self, read, write, tmp_path):
        """Tests renaming columns"""
        fp = tmp_path / "test_rename_html_table_columns.html"

        write(fp, format="ascii.html", latex_names=True)

        tbl = QTable.read(fp)

        # asserts each column name has not been reverted yet
        # For now, Cosmology class and name are stored in first 2 slots
        for column_name in tbl.colnames[2:]:
            assert column_name in _FORMAT_TABLE.values()

        cosmo = read(fp, format="ascii.html")
        converted_tbl = cosmo.to_format("astropy.table")

        # asserts each column name has been reverted
        # cosmology name is still stored in first slot
        for column_name in converted_tbl.colnames[1:]:
            assert column_name in _FORMAT_TABLE.keys()

    @pytest.mark.skipif(not HAS_BS4, reason="requires beautifulsoup4")
    @pytest.mark.parametrize("latex_names", [True, False])
    def test_readwrite_html_subclass_partial_info(
        self, cosmo_cls, cosmo, read, write, latex_names, tmp_path, add_cu
    ):
        """
        Test writing from an instance and reading from that class.
        This works with missing information.
        """
        fp = tmp_path / "test_read_html_subclass_partial_info.html"

        # test write
        write(fp, format="ascii.html", latex_names=latex_names)

        # partial information
        tbl = QTable.read(fp)

        # tbl.meta.pop("cosmology", None) # metadata not implemented
        cname = "$$T_{0}$$" if latex_names else "Tcmb0"
        del tbl[cname]  # format is not converted to original units
        tbl.write(fp, overwrite=True)

        # read with the same class that wrote fills in the missing info with
        # the default value
        got = cosmo_cls.read(fp, format="ascii.html")
        got2 = read(fp, format="ascii.html", cosmology=cosmo_cls)
        got3 = read(fp, format="ascii.html", cosmology=cosmo_cls.__qualname__)

        assert (got == got2) and (got2 == got3)  # internal consistency

        # not equal, because Tcmb0 is changed, which also changes m_nu
        assert got != cosmo
        assert got.Tcmb0 == cosmo_cls.parameters["Tcmb0"].default
        assert got.clone(name=cosmo.name, Tcmb0=cosmo.Tcmb0, m_nu=cosmo.m_nu) == cosmo
        # but the metadata is the same
        # assert got.meta == cosmo.meta # metadata read not implemented

    @pytest.mark.skipif(not HAS_BS4, reason="requires beautifulsoup4")
    def test_readwrite_html_mutlirow(self, cosmo, read, write, tmp_path, add_cu):
        """Test if table has multiple rows."""
        fp = tmp_path / "test_readwrite_html_mutlirow.html"

        # Make
        cosmo1 = cosmo.clone(name="row 0")
        cosmo2 = cosmo.clone(name="row 2")
        table = vstack(
            [c.to_format("astropy.table") for c in (cosmo1, cosmo, cosmo2)],
            metadata_conflicts="silent",
        )

        cosmo_cls = type(cosmo)
        assert cosmo is not None

        for n, col in zip(table.colnames, table.itercols()):
            if n not in cosmo_cls.parameters:
                continue
            param = cosmo_cls.parameters[n]
            if param.unit in (None, u.one):
                continue
            # Replace column with unitless version
            table.replace_column(n, (col << param.unit).value, copy=False)

        table.write(fp, format="ascii.html")

        # ------------
        # From Table

        # it will error on a multi-row table
        with pytest.raises(ValueError, match="need to select a specific row"):
            read(fp, format="ascii.html")

        # unless the index argument is provided
        got = cosmo_cls.read(fp, index=1, format="ascii.html")
        # got = read(fp, index=1, format="ascii.html")
        assert got == cosmo

        # the index can be a string
        got = cosmo_cls.read(fp, index=cosmo.name, format="ascii.html")
        assert got == cosmo

        # it's better if the table already has an index
        # this will be identical to the previous ``got``
        table.add_index("name")
        got2 = cosmo_cls.read(fp, index=cosmo.name, format="ascii.html")
        assert got2 == cosmo


class TestReadWriteHTML(ReadWriteDirectTestBase, ReadWriteHTMLTestMixin):
    """
    Directly test ``read/write_html``.
    These are not public API and are discouraged from use, in favor of
    ``Cosmology.read/write(..., format="ascii.html")``, but should be
    tested regardless b/c they are used internally.
    """

    def setup_class(self):
        self.functions = {"read": read_html_table, "write": write_html_table}

    @pytest.mark.skipif(not HAS_BS4, reason="requires beautifulsoup4")
    def test_rename_direct_html_table_columns(self, read, write, tmp_path):
        """Tests renaming columns"""

        fp = tmp_path / "test_rename_html_table_columns.html"

        write(fp, format="ascii.html", latex_names=True)

        tbl = QTable.read(fp)

        # asserts each column name has not been reverted yet
        for column_name in tbl.colnames[2:]:
            # for now, Cosmology as metadata and name is stored in first 2 slots
            assert column_name in _FORMAT_TABLE.values()

        cosmo = read(fp, format="ascii.html")
        converted_tbl = cosmo.to_format("astropy.table")

        # asserts each column name has been reverted
        for column_name in converted_tbl.colnames[1:]:
            # for now now, metadata is still stored in first slot
            assert column_name in _FORMAT_TABLE.keys()
