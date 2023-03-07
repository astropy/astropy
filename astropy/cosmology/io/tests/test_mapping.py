# Licensed under a 3-clause BSD style license - see LICENSE.rst

# STDLIB
import inspect
from collections import OrderedDict

# THIRD PARTY
import numpy as np
import pytest

# LOCAL
from astropy.cosmology import Cosmology
from astropy.cosmology.io.mapping import from_mapping, to_mapping

from .base import ToFromDirectTestBase, ToFromTestMixinBase

###############################################################################


class ToFromMappingTestMixin(ToFromTestMixinBase):
    """Tests for a Cosmology[To/From]Format with ``format="mapping"``.

    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmology`` for an example.
    """

    def test_to_mapping_default(self, cosmo, to_format):
        """Test default usage of Cosmology -> mapping."""
        m = to_format("mapping")
        keys = tuple(m.keys())

        assert isinstance(m, dict)
        # Check equality of all expected items
        assert keys[0] == "cosmology"
        assert m.pop("cosmology") is cosmo.__class__
        assert keys[1] == "name"
        assert m.pop("name") == cosmo.name
        for i, k in enumerate(cosmo.__parameters__, start=2):
            assert keys[i] == k
            assert np.array_equal(m.pop(k), getattr(cosmo, k))
        assert keys[-1] == "meta"
        assert m.pop("meta") == cosmo.meta

        # No unexpected items
        assert not m

    def test_to_mapping_wrong_cls(self, to_format):
        """Test incorrect argument ``cls`` in ``to_mapping()``."""
        with pytest.raises(TypeError, match="'cls' must be"):
            to_format("mapping", cls=list)

    @pytest.mark.parametrize("map_cls", [dict, OrderedDict])
    def test_to_mapping_cls(self, to_format, map_cls):
        """Test argument ``cls`` in ``to_mapping()``."""
        m = to_format("mapping", cls=map_cls)
        assert isinstance(m, map_cls)  # test type

    def test_to_mapping_cosmology_as_str(self, cosmo_cls, to_format):
        """Test argument ``cosmology_as_str`` in ``to_mapping()``."""
        default = to_format("mapping")

        # Cosmology is the class
        m = to_format("mapping", cosmology_as_str=False)
        assert inspect.isclass(m["cosmology"])
        assert cosmo_cls is m["cosmology"]

        assert m == default  # False is the default option

        # Cosmology is a string
        m = to_format("mapping", cosmology_as_str=True)
        assert isinstance(m["cosmology"], str)
        assert m["cosmology"] == cosmo_cls.__qualname__  # Correct class
        assert tuple(m.keys())[0] == "cosmology"  # Stayed at same index

    def test_tofrom_mapping_cosmology_as_str(self, cosmo, to_format, from_format):
        """Test roundtrip with ``cosmology_as_str=True``.

        The test for the default option (`False`) is in ``test_tofrom_mapping_instance``.
        """
        m = to_format("mapping", cosmology_as_str=True)

        got = from_format(m, format="mapping")
        assert got == cosmo
        assert got.meta == cosmo.meta

    def test_to_mapping_move_from_meta(self, to_format):
        """Test argument ``move_from_meta`` in ``to_mapping()``."""
        default = to_format("mapping")

        # Metadata is 'separate' from main mapping
        m = to_format("mapping", move_from_meta=False)
        assert "meta" in m.keys()
        assert not any(k in m for k in m["meta"])  # Not added to main

        assert m == default  # False is the default option

        # Metadata is mixed into main mapping.
        m = to_format("mapping", move_from_meta=True)
        assert "meta" not in m.keys()
        assert all(k in m for k in default["meta"])  # All added to main
        #  The parameters take precedence over the metadata
        assert all(np.array_equal(v, m[k]) for k, v in default.items() if k != "meta")

    def test_tofrom_mapping_move_tofrom_meta(self, cosmo, to_format, from_format):
        """Test roundtrip of ``move_from/to_meta`` in ``to/from_mapping()``."""
        # Metadata is mixed into main mapping.
        m = to_format("mapping", move_from_meta=True)
        # (Just adding something to ensure there's 'metadata')
        m["mismatching"] = "will error"

        # (Tests are different if the last argument is a **kwarg)
        if tuple(cosmo._init_signature.parameters.values())[-1].kind == 4:
            got = from_format(m, format="mapping")

            assert got.name == cosmo.name
            assert "mismatching" not in got.meta

            return  # don't continue testing

        # Reading with mismatching parameters errors...
        with pytest.raises(TypeError, match="there are unused parameters"):
            from_format(m, format="mapping")

        # unless mismatched are moved to meta.
        got = from_format(m, format="mapping", move_to_meta=True)
        assert got == cosmo  # (Doesn't check metadata)
        assert got.meta["mismatching"] == "will error"

    # -----------------------------------------------------

    def test_from_not_mapping(self, cosmo, from_format):
        """Test incorrect map type in ``from_mapping()``."""
        with pytest.raises((TypeError, ValueError)):
            from_format("NOT A MAP", format="mapping")

    def test_from_mapping_default(self, cosmo, to_format, from_format):
        """Test (cosmology -> Mapping) -> cosmology."""
        m = to_format("mapping")

        # Read from exactly as given.
        got = from_format(m, format="mapping")
        assert got == cosmo
        assert got.meta == cosmo.meta

        # Reading auto-identifies 'format'
        got = from_format(m)
        assert got == cosmo
        assert got.meta == cosmo.meta

    def test_fromformat_subclass_partial_info_mapping(self, cosmo):
        """
        Test writing from an instance and reading from that class.
        This works with missing information.
        """
        m = cosmo.to_format("mapping")

        # partial information
        m.pop("cosmology", None)
        m.pop("Tcmb0", None)

        # read with the same class that wrote fills in the missing info with
        # the default value
        got = cosmo.__class__.from_format(m, format="mapping")
        got2 = Cosmology.from_format(m, format="mapping", cosmology=cosmo.__class__)
        got3 = Cosmology.from_format(
            m, format="mapping", cosmology=cosmo.__class__.__qualname__
        )

        assert (got == got2) and (got2 == got3)  # internal consistency

        # not equal, because Tcmb0 is changed, which also changes m_nu
        assert got != cosmo
        assert got.Tcmb0 == cosmo.__class__._init_signature.parameters["Tcmb0"].default
        assert got.clone(name=cosmo.name, Tcmb0=cosmo.Tcmb0, m_nu=cosmo.m_nu) == cosmo
        # but the metadata is the same
        assert got.meta == cosmo.meta

    @pytest.mark.parametrize("format", [True, False, None, "mapping"])
    def test_is_equivalent_to_mapping(self, cosmo, to_format, format):
        """Test :meth:`astropy.cosmology.Cosmology.is_equivalent`.

        This test checks that Cosmology equivalency can be extended to any
        Python object that can be converted to a Cosmology -- in this case
        a mapping.
        """
        obj = to_format("mapping")
        assert not isinstance(obj, Cosmology)

        is_equiv = cosmo.is_equivalent(obj, format=format)
        assert is_equiv is (format is not False)


class TestToFromMapping(ToFromDirectTestBase, ToFromMappingTestMixin):
    """Directly test ``to/from_mapping``."""

    def setup_class(self):
        self.functions = {"to": to_mapping, "from": from_mapping}

    @pytest.mark.skip("N/A")
    def test_fromformat_subclass_partial_info_mapping(self):
        """This test does not apply to the direct functions."""
