# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from astropy.cosmology import parameters, realizations
from astropy.cosmology.realizations import Planck13, default_cosmology
from astropy.utils.exceptions import AstropyDeprecationWarning

cosmo_realz = [s for s in parameters.available if s != "Planck18_arXiv_v2"]


class Test_default_cosmology(object):
    """Tests for :class:`~astropy.cosmology.realizations.default_cosmology`."""

    @staticmethod
    def test_get_cosmology_from_string():
        """Test method ``get_cosmology_from_string``."""
        cosmo = default_cosmology.get_cosmology_from_string("no_default")
        assert cosmo is None

        cosmo = default_cosmology.get_cosmology_from_string("Planck13")
        assert cosmo is Planck13

        with pytest.raises(ValueError):
            cosmo = default_cosmology.get_cosmology_from_string("fail!")

    @staticmethod
    def test_validate_specific():
        """Test method ``validate`` for specific values."""
        value = default_cosmology.validate(None)
        assert value is realizations.Planck18

        with pytest.warns(AstropyDeprecationWarning) as e:
            value = default_cosmology.validate("Planck18_arXiv_v2")
        assert "use Planck18 instead" in str(e.list[0].message)
        assert value is realizations.Planck18_arXiv_v2

        with pytest.raises(TypeError) as e:
            default_cosmology.validate(TypeError)
        assert "string or Cosmology instance" in str(e.value)

    @pytest.mark.parametrize("name", cosmo_realz)
    def test_validate_str(self, name):
        """Test method ``validate`` for string input."""
        value = default_cosmology.validate(name)
        assert value is getattr(realizations, name)

    @pytest.mark.parametrize("name", cosmo_realz)
    def test_validate_cosmo(self, name):
        """Test method ``validate`` for cosmology instance input."""
        cosmo = getattr(realizations, name)
        value = default_cosmology.validate(cosmo)
        assert value is cosmo
