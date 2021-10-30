# Licensed under a 3-clause BSD style license - see LICENSE.rst

# STDLIB
import copy
from collections import OrderedDict

# THIRD PARTY
import pytest

# LOCAL
import astropy.units as u
from astropy.cosmology import Cosmology, realizations
from astropy.cosmology.core import _COSMOLOGY_CLASSES, Parameter
from astropy.cosmology.io.mapping import from_mapping, to_mapping
from astropy.cosmology.parameters import available
from astropy.table import QTable, vstack

from .base import IOTestMixinBase, IOFormatTestBase

cosmo_instances = [getattr(realizations, name) for name in available]
cosmo_instances.append("TestToFromMapping.setup.<locals>.CosmologyWithKwargs")


###############################################################################


class ToFromMappingTestMixin(IOTestMixinBase):
    """
    Tests for a Cosmology[To/From]Format with ``format="mapping"``.
    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmology`` for an example.
    """

    def test_failed_cls_to_mapping(self, cosmo, to_format):
        """Test incorrect argument ``cls`` in ``to_mapping()``."""
        with pytest.raises(TypeError, match="'cls' must be"):
            to_format('mapping', cls=list)

    @pytest.mark.parametrize("map_cls", [dict, OrderedDict])
    def test_to_mapping_cls(self, cosmo, to_format, map_cls):
        """Test argument ``cls`` in ``to_mapping()``."""
        params = to_format('mapping', cls=map_cls)
        assert isinstance(params, map_cls)  # test type

    def test_tofrom_mapping_instance(self, cosmo, to_format, from_format):
        """Test cosmology -> Mapping -> cosmology."""
        # ------------
        # To Mapping

        params = to_format('mapping')
        assert isinstance(params, dict)  # test type
        assert params["cosmology"] is cosmo.__class__
        assert params["name"] == cosmo.name

        # ------------
        # From Mapping

        params["mismatching"] = "will error"

        # tests are different if the last argument is a **kwarg
        if tuple(cosmo._init_signature.parameters.values())[-1].kind == 4:
            got = from_format(params, format="mapping")

            assert got.name == cosmo.name
            assert "mismatching" not in got.meta

            return  # don't continue testing

        # read with mismatching parameters errors
        with pytest.raises(TypeError, match="there are unused parameters"):
            from_format(params, format="mapping")

        # unless mismatched are moved to meta
        got = from_format(params, format="mapping", move_to_meta=True)
        assert got == cosmo
        assert got.meta["mismatching"] == "will error"

        # it won't error if everything matches up
        params.pop("mismatching")
        got = from_format(params, format="mapping")
        assert got == cosmo
        assert got.meta == cosmo.meta

        # and it will also work if the cosmology is a string
        params["cosmology"] = params["cosmology"].__qualname__
        got = from_format(params, format="mapping")
        assert got == cosmo
        assert got.meta == cosmo.meta

        # also it auto-identifies 'format'
        got = from_format(params)
        assert got == cosmo
        assert got.meta == cosmo.meta

    def test_fromformat_subclass_partial_info_mapping(self, cosmo):
        """
        Test writing from an instance and reading from that class.
        This works with missing information.
        """
        # test to_format
        m = cosmo.to_format("mapping")
        assert isinstance(m, dict)

        # partial information
        m.pop("cosmology", None)
        m.pop("Tcmb0", None)

        # read with the same class that wrote fills in the missing info with
        # the default value
        got = cosmo.__class__.from_format(m, format="mapping")
        got2 = Cosmology.from_format(m, format="mapping", cosmology=cosmo.__class__)
        got3 = Cosmology.from_format(m, format="mapping", cosmology=cosmo.__class__.__qualname__)

        assert (got == got2) and (got2 == got3)  # internal consistency

        # not equal, because Tcmb0 is changed
        assert got != cosmo
        assert got.Tcmb0 == cosmo.__class__._init_signature.parameters["Tcmb0"].default
        assert got.clone(name=cosmo.name, Tcmb0=cosmo.Tcmb0) == cosmo
        # but the metadata is the same
        assert got.meta == cosmo.meta


class TestToFromMapping(IOFormatTestBase, ToFromMappingTestMixin):
    """Directly test ``to/from_mapping``."""

    def setup_class(self):
        self.functions = {"to": to_mapping, "from": from_mapping}
