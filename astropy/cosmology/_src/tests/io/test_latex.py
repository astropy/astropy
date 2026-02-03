# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from astropy.cosmology._src.io.builtin.latex import (
    _FORMAT_TABLE,
    read_latex,
    write_latex,
)
from astropy.io.registry.base import IORegistryError
from astropy.table import QTable, Table

from .base import ReadWriteDirectTestBase, ReadWriteTestMixinBase


class ReadWriteLATEXTestMixin(ReadWriteTestMixinBase):
    """
    Tests for a Cosmology[Read/Write] with ``format="ascii.latex"``.
    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmology`` for an example.
    """

    def test_to_latex_bad_index(self, read, write, tmp_path):
        """Test if argument ``index`` is incorrect"""

        fp = tmp_path / "test_to_latex_bad_index.tex"

        write(fp, format="ascii.latex")

        # single-row table and has a non-0/None index
        with pytest.raises(IndexError, match="index 2 out of range"):
            read(fp, index=2, format="ascii.latex")

        # string index where doesn't match
        with pytest.raises(KeyError, match="No matches found for key"):
            read(fp, index="row 0", format="ascii.latex")

    def test_read_latex_invalid_path(self, read):
        """Test passing an invalid or non-existent path"""

        # using "blabla" makes it platform independent
        invalid_fp = "blabla.tex"

        with pytest.raises(FileNotFoundError, match="No such file or directory"):
            read(invalid_fp, format="ascii.latex")

    def test_latex_column_mnu(self, read, write, tmp_path):
        """Test for table column m_nu to have a numpy array, essential for proper cosmology conversion"""

        fp = tmp_path / "test_rename_latex_columns.tex"

        write(fp, latex_names=True)

        tbl = QTable.read(fp)

        import numpy as np

        assert isinstance(tbl["$m_{nu}$"].value, np.ndarray)

    # -----------------------

    def test_to_latex_failed_cls(self, write, tmp_path):
        """Test failed table type."""
        fp = tmp_path / "test_to_latex_failed_cls.tex"

        with pytest.raises(TypeError, match="'cls' must be"):
            write(fp, cls=list)

    @pytest.mark.parametrize("tbl_cls", [QTable, Table])
    def test_to_latex_cls(self, write, tbl_cls, tmp_path):
        fp = tmp_path / "test_to_latex_cls.tex"
        write(fp, cls=tbl_cls)

    def test_latex_columns(self, write, tmp_path):
        fp = tmp_path / "test_rename_latex_columns.tex"
        write(fp, latex_names=True)
        tbl = QTable.read(fp)
        # asserts each column name has not been reverted yet
        # For now, Cosmology class and name are stored in first 2 slots
        for column_name in tbl.colnames[2:]:
            assert column_name in _FORMAT_TABLE.values()

    def test_write_latex_invalid_path(self, write):
        """Test passing an invalid path"""
        invalid_fp = ""
        with pytest.raises(FileNotFoundError, match="No such file or directory"):
            write(invalid_fp, format="ascii.latex")

    def test_write_latex_false_overwrite(self, write, tmp_path):
        """Test to write a LaTeX file without overwriting an existing file"""
        # Test that passing an invalid path to write_latex() raises a IOError
        fp = tmp_path / "test_write_latex_false_overwrite.tex"
        write(fp)
        with pytest.raises(OSError, match="overwrite=True"):
            write(fp, overwrite=False)

    def test_write_latex_unsupported_format(self, write, tmp_path):
        """Test for unsupported format"""
        fp = tmp_path / "test_write_latex_unsupported_format.tex"
        invalid_format = "unsupported"
        with pytest.raises((ValueError, IORegistryError)) as exc_info:
            pytest.raises(ValueError, match="format must be 'ascii.latex'")
            pytest.raises(IORegistryError, match="No writer defined for format")
            write(fp, format=invalid_format)


class TestReadWriteLaTex(ReadWriteDirectTestBase, ReadWriteLATEXTestMixin):
    """
    Directly test ``read/write_latex``.
    These are not public API and are discouraged from use, in favor of
    ``Cosmology.read/write(..., format="latex")``, but should be
    tested regardless b/c they are used internally.
    """

    def setup_class(self):
        self.functions = {"write": write_latex, "read": read_latex}

    def test_rename_direct_latex_columns(self, read, write, tmp_path):
        """Tests renaming columns"""
        fp = tmp_path / "test_rename_latex_columns.tex"
        write(fp, latex_names=True)
        tbl = QTable.read(fp)
        # asserts each column name has not been reverted yet
        for column_name in tbl.colnames[2:]:
            # for now, Cosmology as metadata and name is stored in first 2 slots
            assert column_name in _FORMAT_TABLE.values()

        cosmo = read(fp, format="ascii.latex")
        converted_tbl = cosmo.to_format("astropy.table")

        # asserts each column name has been reverted
        for column_name in converted_tbl.colnames[1:]:
            # for now now, metadata is still stored in first slot
            assert column_name in _FORMAT_TABLE.keys()
