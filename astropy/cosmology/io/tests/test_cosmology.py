# Licensed under a 3-clause BSD style license - see LICENSE.rst

# THIRD PARTY
import pytest

# LOCAL
from astropy.cosmology import Cosmology
from astropy.cosmology.io.cosmology import from_cosmology, to_cosmology

from .base import ToFromTestMixinBase, IODirectTestBase


###############################################################################


class ToFromCosmologyTestMixin(ToFromTestMixinBase):
    """
    Tests for a Cosmology[To/From]Format with ``format="astropy.cosmology"``.
    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmology`` for an example.
    """

    def test_to_cosmology_default(self, cosmo, to_format):
        """Test cosmology -> cosmology."""
        newcosmo = to_format("astropy.cosmology")
        assert newcosmo is cosmo

    def test_from_not_cosmology(self, cosmo, from_format):
        """Test incorrect type in ``Cosmology``."""
        with pytest.raises(TypeError):
            from_format("NOT A COSMOLOGY", format="astropy.cosmology")

    def test_from_cosmology_default(self, cosmo, from_format):
        """Test cosmology -> cosmology."""
        newcosmo = from_format(cosmo)
        assert newcosmo is cosmo


class TestToFromCosmology(IODirectTestBase, ToFromCosmologyTestMixin):
    """Directly test ``to/from_cosmology``."""

    def setup_class(self):
        self.functions = {"to": to_cosmology, "from": from_cosmology}
