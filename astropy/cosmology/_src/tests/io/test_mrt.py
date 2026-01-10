# Licensed under a 3-clause BSD style license - see LICENSE.rst

# THIRD PARTY

import pytest

from astropy.cosmology import Planck18
from astropy.cosmology._src.io.builtin.mrt import (
    read_mrt,
    write_mrt,
)
from astropy.table import QTable, Table

from .base import ReadWriteDirectTestBase, ReadWriteTestMixinBase


class ReadWriteMRTTestMixin(ReadWriteTestMixinBase):
    """
    Tests for a Cosmology[Read/Write] with ``format="mrt"``.
    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmology`` for an example.
    """

    def test_to_mrt_bad_index(self, read, write, tmp_path):
        """Test if argument ``index`` is incorrect"""
        fp = tmp_path / "test_to_mrt_bad_index.mrt"

        write(fp, format="ascii.mrt")

        # single-row table and has a non-0/None index
        with pytest.raises(IndexError, match="index 2 out of range"):
            read(fp, index=2, format="ascii.mrt")

        # string index where doesn't match
        with pytest.raises(KeyError, match="No matches found for key"):
            read(fp, index="row 0", format="ascii.mrt")

    # -----------------------

    def test_to_mrt_failed_cls(self, write, tmp_path):
        """Test failed table type."""
        fp = tmp_path / "test_to_mrt_failed_cls.mrt"

        with pytest.raises(TypeError, match="'cls' must be"):
            write(fp, format="ascii.mrt", cls=list)

    # -----------------------

    @pytest.mark.parametrize("tbl_cls", [QTable, Table])
    def test_to_mrt_cls(self, write, tbl_cls, tmp_path):
        fp = tmp_path / "test_to_mrt_cls.mrt"
        write(fp, format="ascii.mrt", cls=tbl_cls)

    # -----------------------

    def test_readwrite_mrt_instance(self, cosmo_cls, cosmo, read, write, tmp_path):
        """Test cosmology -> ascii.mrt -> cosmology."""
        fp = tmp_path / "test_readwrite_mrt_instance.mrt"

        # ------------
        # To Table

        write(fp, format="ascii.mrt")

        # some checks on the saved file
        tbl = Table.read(fp, format="ascii.mrt")
        assert tbl["name"] == cosmo.name

        # ------------
        # From Table

        tbl["mismatching"] = "will error"
        tbl.write(fp, format="ascii.mrt", overwrite=True)

        # tests are different if the last argument is a **kwarg
        if cosmo._init_has_kwargs:
            got = read(fp, format="ascii.mrt")

            assert got.__class__ is cosmo_cls
            assert got.name == cosmo.name
            assert "mismatching" not in got.meta

            return  # don't continue testing

        # read with mismatching parameters errors
        with pytest.raises(TypeError, match="there are unused parameters"):
            read(fp, format="ascii.mrt")

        # unless mismatched are moved to meta
        got = read(fp, format="ascii.mrt", move_to_meta=True)
        assert got == cosmo
        assert got.meta["mismatching"] == "will error"

        # it won't error if everything matches up
        tbl.remove_column("mismatching")
        tbl.write(fp, format="ascii.mrt", overwrite=True)
        got = read(fp, format="ascii.mrt")
        assert got == cosmo

    # -----------------------

    def test_readwrite_mrt_subclass_partial_info(
        self, cosmo_cls, cosmo, read, write, tmp_path
    ):
        """
        Test writing from an instance and reading from that class.
        This works with missing information.
        """
        fp = tmp_path / "test_read_mrt_subclass_partial_info.mrt"

        # test write
        write(fp, format="ascii.mrt")

        # partial information
        tbl = Table.read(fp, format="ascii.mrt")
        del tbl["Tcmb0"]
        tbl.write(fp, format="ascii.mrt", overwrite=True)

        # read with the same class that wrote fills in the missing info with
        # the default value
        got = cosmo_cls.read(fp, format="ascii.mrt")
        got2 = read(fp, format="ascii.mrt", cosmology=cosmo_cls)

        assert got == got2  # internal consistency

        # not equal, because Tcmb0 is changed, which also changes m_nu
        assert got != cosmo
        assert got.Tcmb0 == cosmo_cls.parameters["Tcmb0"].default
        assert got.clone(name=cosmo.name, Tcmb0=cosmo.Tcmb0, m_nu=cosmo.m_nu) == cosmo


class WriteMRTTestMixin(ReadWriteTestMixinBase):
    def test_to_mrt_failed_cls(self, write, tmp_path):
        """Test failed table type."""
        fp = tmp_path / "test_to_mrt_failed_cls.mrt"

        with pytest.raises(TypeError, match="'cls' must be"):
            write(fp, format="ascii.mrt", cls=list)

    @pytest.mark.parametrize("tbl_cls", [QTable, Table])
    def test_to_mrt_cls(self, write, tbl_cls, tmp_path):
        fp = tmp_path / "test_to_mrt_cls.mrt"
        write(fp, format="ascii.mrt", cls=tbl_cls)

    def test_to_mrt_bad_index(self, read, write, tmp_path):
        """Test if argument ``index`` is incorrect"""
        fp = tmp_path / "test_to_mrt_bad_index.mrt"

        write(fp, format="ascii.mrt")

        # single-row table and has a non-0/None index
        with pytest.raises(IndexError, match="index 2 out of range"):
            read(fp, index=2, format="ascii.mrt")

        # string index where doesn't match
        with pytest.raises(KeyError, match="No matches found for key"):
            read(fp, index="row 0", format="ascii.mrt")

    def test_write_mrt_invalid_path(self, write):
        """Test passing an invalid path"""
        invalid_fp = ""
        with pytest.raises(FileNotFoundError, match="No such file or directory"):
            write(invalid_fp, format="ascii.mrt")

    def test_readwrite_mrt_instance(self, cosmo_cls, cosmo, read, write, tmp_path):
        fp = tmp_path / "test_readwrite_mrt_instance.mrt"
        write(fp, format="ascii.mrt")
        tb1 = Table.read(fp, format="ascii.mrt")
        assert tb1["name"] == cosmo.name


class TestReadWriteMRT(ReadWriteDirectTestBase, WriteMRTTestMixin):
    """
    Directly test ``read/write_mrt``.
    These are not public API and are discouraged from use, in favor of
    ``Cosmology.read/write(..., format="mrt")``, but should be
    tested regardless b/c they are used internally.
    """

    def setup_class(self):
        self.functions = {"read": read_mrt, "write": write_mrt}


def test_write_mrt_invalid_format(tmp_path):
    """Test passing an invalid format"""
    fp = tmp_path / "test_write_mrt_invalid_format.mrt"
    with pytest.raises(ValueError, match="format must be 'ascii.mrt'"):
        write_mrt(Planck18, fp, format="ascii.ecsv")


def test_read_mrt_invalid_format(tmp_path):
    """Test read_mrt with invalid format parameter."""
    fp = tmp_path / "test_read_mrt_invalid_format.mrt"
    Planck18.write(fp, format="ascii.mrt")

    # Test that passing a different format raises ValueError
    with pytest.raises(ValueError, match="format must be 'ascii.mrt'"):
        read_mrt(fp, format="ascii.ecsv")
