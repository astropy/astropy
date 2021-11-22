# Licensed under a 3-clause BSD style license - see LICENSE.rst

# THIRD PARTY
import pytest

# LOCAL
import astropy.units as u
from astropy import cosmology
from astropy.cosmology import Cosmology, realizations, Planck18
from astropy.cosmology.core import _COSMOLOGY_CLASSES, Parameter
from astropy.cosmology.io.row import from_row, to_row
from astropy.table import Row
from astropy.cosmology.parameters import available

from .base import IOTestMixinBase, IOFormatTestBase

cosmo_instances = [getattr(realizations, name) for name in available]
cosmo_instances.append("TestToFromRow.setup.<locals>.CosmologyWithKwargs")

###############################################################################


class ToFromRowTestMixin(IOTestMixinBase):
    """
    Tests for a Cosmology[To/From]Format with ``format="astropy.row"``.
    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmologyToFromFormat`` or ``TestCosmology`` for examples.
    """

    @pytest.mark.parametrize("in_meta", [True, False])
    def test_to_row_in_meta(self, cosmo_cls, cosmo, in_meta):
        """Test where the cosmology class is placed."""
        row = cosmo.to_format('astropy.row', cosmology_in_meta=in_meta)

        # if it's in metadata, it's not a column. And vice versa.
        if in_meta:
            assert row.meta["cosmology"] == cosmo_cls.__qualname__
            assert "cosmology" not in row.colnames  # not also a column
        else:
            assert row["cosmology"] == cosmo_cls.__qualname__
            assert "cosmology" not in row.meta

    # -----------------------

    def test_from_not_row(self, cosmo, from_format):
        """Test not passing a Row to the Row parser."""
        with pytest.raises(AttributeError):
            from_format("NOT A ROW", format="astropy.row")

    def test_tofrom_row_instance(self, cosmo, to_format, from_format):
        """Test cosmology -> astropy.row -> cosmology."""
        # ------------
        # To Row

        row = to_format("astropy.row")
        assert isinstance(row, Row)
        assert row["cosmology"] == cosmo.__class__.__qualname__
        assert row["name"] == cosmo.name

        # ------------
        # From Row

        row.table["mismatching"] = "will error"

        # tests are different if the last argument is a **kwarg
        if tuple(cosmo._init_signature.parameters.values())[-1].kind == 4:
            got = from_format(row, format="astropy.row")

            assert got.__class__ is cosmo.__class__
            assert got.name == cosmo.name
            assert "mismatching" not in got.meta

            return  # don't continue testing

        # read with mismatching parameters errors
        with pytest.raises(TypeError, match="there are unused parameters"):
            from_format(row, format="astropy.row")

        # unless mismatched are moved to meta
        got = from_format(row, format="astropy.row", move_to_meta=True)
        assert got == cosmo
        assert got.meta["mismatching"] == "will error"

        # it won't error if everything matches up
        row.table.remove_column("mismatching")
        got = from_format(row, format="astropy.row")
        assert got == cosmo

        # and it will also work if the cosmology is a class
        # Note this is not the default output of ``to_format``.
        cosmology = _COSMOLOGY_CLASSES[row["cosmology"]]
        row.table.remove_column("cosmology")
        row.table["cosmology"] = cosmology
        got = from_format(row, format="astropy.row")
        assert got == cosmo

        # also it auto-identifies 'format'
        got = from_format(row)
        assert got == cosmo

    def test_fromformat_row_subclass_partial_info(self, cosmo):
        """
        Test writing from an instance and reading from that class.
        This works with missing information.
        """
        pass  # there are no partial info options


class TestToFromTable(IOFormatTestBase, ToFromRowTestMixin):
    """
    Directly test ``to/from_row``.
    These are not public API and are discouraged from use, in favor of
    ``Cosmology.to/from_format(..., format="astropy.row")``, but should be
    tested regardless b/c 3rd party packages might use these in their Cosmology
    I/O. Also, it's cheap to test.
    """

    def setup_class(self):
        self.functions = {"to": to_row, "from": from_row}
