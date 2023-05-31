# Licensed under a 3-clause BSD style license - see LICENSE.rst

# THIRD PARTY
import pytest

# LOCAL
from astropy.cosmology.io.mrt import read_mrt, write_mrt
from astropy.table import QTable, Table

from .base import ReadWriteDirectTestBase, ReadWriteTestMixinBase

class ReadWriteMRTTestMixin(ReadWriteTestMixinBase):
    """
    Tests for a Cosmology[Read/Write] with ``format="mrt"``.
    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must dfine a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmology`` for an example.
    """
    def test_to_mrt_bad_index(self, read, write, tmp_path):
        """Test if argument ``index`` is incorrect"""
        fp = tmp_path / "test_to_mrt_bad_index.mrt"

        write(fp, format="mrt")
        
        # single-row table and has a non-0/None index
        with pytest.raises(IndexError, match="index 2 out of range"):
            read(fp, index=2, format="mrt")

        # string index where doesn't match
        with pytest.raises(KeyError, match="No matches found for key"):
            read(fp, index="row 0", format="mrt")
    
    # -----------------------

    def test_to_mrt_failed_cls(self, write, tmp_path):
        """Test failed table type."""
        fp = tmp_path / "test_to_mrt_failed_cls.mrt"
        
        with pytest.raises(TypeError, match="'cls' must be"):
            write(fp, format="mrt", cls=list)
    
    # -----------------------

    @pytest.mark.parametrize("tbl_cls", [QTable, Table])
    def test_to_mrt_cls(self, write, tbl_cls, tmp_path):
        fp = tmp_path / "test_to_mrt_cls.mrt"
        write(fp, format="mrt", cls=tbl_cls)

    # -----------------------

    def test_readwrite_mrt_instance(
        self, cosmo_cls, cosmo, read, write, tmp_path
    ):
        """Test cosmology -> mrt -> cosmology."""
        fp = tmp_path / "test_readwrite_mrt_instance.mrt"

        # ------------
        # To Table

        write(fp, format="mrt")

        # some checks on the saved file
        tbl = Table.read(fp, format="mrt")
        assert tbl["name"] == cosmo.name

        # ------------
        # From Table

        tbl["mismatching"] = "will error"
        tbl.write(fp, format="mrt", overwrite=True)
        
        # tests are different if the last argument is a **kwarg
        if tuple(cosmo._init_signature.parameters.values())[-1].kind == 4:
            got = read(fp, format="mrt")

            assert got.__class__ is cosmo_cls
            assert got.name == cosmo.name
            assert "mismatching" not in got.meta

            return  # don't continue testing

        # read with mismatching parameters errors
        with pytest.raises(TypeError, match="there are unused parameters"):
            read(fp, format="mrt")

        # unless mismatched are moved to meta
        got = read(fp, format="mrt", move_to_meta=True)
        assert got == cosmo
        assert got.meta["mismatching"] == "will error"
        
        # it won't error if everything matches up
        tbl.remove_column("mismatching")
        tbl.write(fp, format="mrt", overwrite=True)
        got = read(fp, format="mrt")
        assert got == cosmo

        # also it auto-identifies 'format'
        got = read(fp)
        assert got == cosmo

    def test_readwrite_mrt_subclass_partial_info(
        self, cosmo_cls, cosmo, read, write, tmp_path
    ):
        """
        Test writing from an instance and reading from that class.
        This works with missing information.
        """
        fp = tmp_path / "test_read_mrt_subclass_partial_info.mrt"

        # test write
        write(fp, format="mrt")

        # partial information
        tbl = Table.read(fp, format="mrt")
        del tbl["Tcmb0"]
        tbl.write(fp, format="mrt", overwrite=True)

        # read with the same class that wrote fills in the missing info with
        # the default value
        got = cosmo_cls.read(fp, format="mrt")
        got2 = read(fp, format="mrt", cosmology=cosmo_cls)

        assert got == got2  # internal consistency

        # not equal, because Tcmb0 is changed, which also changes m_nu
        assert got != cosmo
        assert got.Tcmb0 == cosmo_cls._init_signature.parameters["Tcmb0"].default
        assert got.clone(name=cosmo.name, Tcmb0=cosmo.Tcmb0, m_nu=cosmo.m_nu) == cosmo

class TestReadWriteMRT(ReadWriteDirectTestBase, ReadWriteMRTTestMixin):
    """
    Directly test ``read/write_mrt``.
    These are not public API and are discouraged from use, in favor of
    ``Cosmology.read/write(..., format="mrt")``, but should be
    tested regardless b/c they are used internally.
    """

    def setup_class(self):
        self.functions = {"read": read_mrt, "write": write_mrt}