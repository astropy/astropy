# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from astropy.cosmology._src.io.builtin.cosmology import from_cosmology, to_cosmology

from .base import IODirectTestBase, ToFromTestMixinBase

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

    @pytest.mark.parametrize("format", [True, False, None, "astropy.cosmology"])
    def test_is_equivalent_to_cosmology(self, cosmo, to_format, format):
        """Test :meth:`astropy.cosmology.Cosmology.is_equivalent`.

        This test checks that Cosmology equivalency can be extended to any
        Python object that can be converted to a Cosmology -- in this case
        a Cosmology! Since it's the identity conversion, the cosmology is
        always equivalent to itself, regardless of ``format``.
        """
        obj = to_format("astropy.cosmology")
        assert obj is cosmo

        is_equiv = cosmo.is_equivalent(obj, format=format)
        assert is_equiv is True  # equivalent to self


class TestToFromCosmology(IODirectTestBase, ToFromCosmologyTestMixin):
    """Directly test ``to/from_cosmology``."""

    def setup_class(self):
        self.functions = {"to": to_cosmology, "from": from_cosmology}
