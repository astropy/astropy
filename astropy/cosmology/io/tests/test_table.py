# Licensed under a 3-clause BSD style license - see LICENSE.rst

# THIRD PARTY
import pytest

# LOCAL
import astropy.units as u
from astropy import cosmology
from astropy.cosmology import Cosmology, realizations, Planck18
from astropy.cosmology.core import _COSMOLOGY_CLASSES, Parameter
from astropy.cosmology.io.table import from_table, to_table
from astropy.table import Table, QTable, vstack
from astropy.cosmology.parameters import available

from .base import IOTestMixinBase, ToFromFormatTestBase

cosmo_instances = [getattr(realizations, name) for name in available]
cosmo_instances.append("TestToFromTable.setup.<locals>.CosmologyWithKwargs")

###############################################################################


class ToFromTableTestMixin(IOTestMixinBase):
    """
    Tests for a Cosmology[To/From]Format with ``format="astropy.table"``.
    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmologyToFromFormat`` or ``TestCosmology`` for examples.
    """

    def test_to_table_bad_index(self, cosmo, to_format, from_format):
        """Test if argument ``index`` is incorrect"""
        tbl = to_format("astropy.table")

        # single-row table and has a non-0/None index
        with pytest.raises(IndexError, match="index 2 out of range"):
            from_format(tbl, index=2, format="astropy.table")

        # string index where doesn't match
        with pytest.raises(IndexError, match="index 0 is out of bounds"):
            from_format(tbl, index="row 0", format="astropy.table")

    # -----------------------

    def test_to_table_failed_cls(self, cosmo, to_format):
        """Test failed table type."""
        with pytest.raises(TypeError, match="'cls' must be"):
            to_format('astropy.table', cls=list)

    @pytest.mark.parametrize("tbl_cls", [QTable, Table])
    def test_to_table_cls(self, cosmo, to_format, tbl_cls):
        tbl = to_format('astropy.table', cls=tbl_cls)
        assert isinstance(tbl, tbl_cls)  # test type

    # -----------------------

    @pytest.mark.parametrize("in_meta", [True, False])
    def test_to_table_in_meta(self, cosmo, in_meta):
        """Test where the cosmology class is placed."""
        tbl = cosmo.to_format('astropy.table', cosmology_in_meta=in_meta)

        # if it's in metadata, it's not a column. And vice versa.
        if in_meta:
            assert tbl.meta["cosmology"] == cosmo.__class__.__qualname__
            assert "cosmology" not in tbl.colnames  # not also a column
        else:
            assert tbl["cosmology"][0] == cosmo.__class__.__qualname__
            assert "cosmology" not in tbl.meta

    # -----------------------

    def test_to_from_table_instance(self, cosmo, to_format, from_format):
        """Test cosmology -> astropy.table -> cosmology."""
        # ------------
        # To Table

        tbl = to_format("astropy.table")
        assert isinstance(tbl, QTable)
        assert tbl.meta["cosmology"] == cosmo.__class__.__qualname__
        assert tbl["name"] == cosmo.name

        # ------------
        # From Table

        tbl["mismatching"] = "will error"

        # tests are different if the last argument is a **kwarg
        if tuple(cosmo._init_signature.parameters.values())[-1].kind == 4:
            got = from_format(tbl, format="astropy.table")

            assert got.__class__ is cosmo.__class__
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

        # and it will also work if the cosmology is a string
        tbl.meta["cosmology"] = _COSMOLOGY_CLASSES[tbl.meta["cosmology"]].__qualname__
        got = from_format(tbl, format="astropy.table")
        assert got == cosmo

        # also it auto-identifies 'format'
        got = from_format(tbl)
        assert got == cosmo

    def test_fromformat_subclass_partial_info_table(self, cosmo):
        """
        Test writing from an instance and reading from that class.
        This works with missing information.
        """
        # test to_format
        tbl = cosmo.to_format("astropy.table")
        assert isinstance(tbl, QTable)

        # partial information
        tbl.meta.pop("cosmology", None)
        del tbl["Tcmb0"]

        # read with the same class that wrote fills in the missing info with
        # the default value
        got = cosmo.__class__.from_format(tbl, format="astropy.table")
        got2 = Cosmology.from_format(tbl, format="astropy.table", cosmology=cosmo.__class__)
        got3 = Cosmology.from_format(tbl, format="astropy.table",
                                     cosmology=cosmo.__class__.__qualname__)

        assert (got == got2) and (got2 == got3)  # internal consistency

        # not equal, because Tcmb0 is changed
        assert got != cosmo
        assert got.Tcmb0 == cosmo.__class__._init_signature.parameters["Tcmb0"].default
        assert got.clone(name=cosmo.name, Tcmb0=cosmo.Tcmb0) == cosmo
        # but the metadata is the same
        assert got.meta == cosmo.meta

    def test_to_from_table_mutlirow(self, cosmo, to_format, from_format):
        """Test if table has multiple rows."""
        # ------------
        # To Table

        cosmo1 = cosmo.clone(name="row 0")
        cosmo2 = cosmo.clone(name="row 2")
        tbl = vstack([c.to_format("astropy.table") for c in (cosmo1, cosmo, cosmo2)],
                     metadata_conflicts='silent')

        assert isinstance(tbl, QTable)
        assert tbl.meta["cosmology"] == cosmo.__class__.__qualname__
        assert tbl[1]["name"] == cosmo.name

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

        # it's better if the table already has an index
        # this will be identical to the previous ``got``
        tbl.add_index("name")
        got2 = from_format(tbl, index=cosmo.name, format="astropy.table")
        assert got2 == cosmo


class TestToFromTable(ToFromFormatTestBase, ToFromTableTestMixin):
    """
    Directly test ``to/from_table``.
    These are not public API and are discouraged from use, in favor of
    ``Cosmology.to/from_format(..., format="astropy.table")``, but should be
    tested regardless b/c 3rd party packages might use these in their Cosmology
    I/O. Also, it's cheap to test.
    """

    def setup_class(self):
        self.io_functions = {"to": to_table, "from": from_table}
