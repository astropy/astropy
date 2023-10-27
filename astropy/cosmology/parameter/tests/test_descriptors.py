# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.parameter._descriptor`."""

from types import MappingProxyType

import pytest

from astropy.cosmology.parameter import Parameter


class ParametersAttributeTestMixin:
    """Test the descriptor for ``parameters`` on Cosmology classes."""

    @pytest.mark.parametrize("name", ["parameters", "derived_parameters"])
    def test_parameters_from_class(self, cosmo_cls, name):
        """Test descriptor ``parameters`` accessed from the class."""
        # test presence
        assert hasattr(cosmo_cls, name)
        # test Parameter is a MappingProxyType
        parameters = getattr(cosmo_cls, name)
        assert isinstance(parameters, MappingProxyType)
        # Test items
        assert set(parameters.keys()) == set(getattr(cosmo_cls, f"_{name}"))
        assert all(isinstance(p, Parameter) for p in parameters.values())

    @pytest.mark.parametrize("name", ["parameters", "derived_parameters"])
    def test_parameters_from_instance(self, cosmo, name):
        """Test descriptor ``parameters`` accessed from the instance."""
        # test presence
        assert hasattr(cosmo, name)
        # test Parameter is a MappingProxyType
        parameters = getattr(cosmo, name)
        assert isinstance(parameters, MappingProxyType)
        # Test keys
        assert set(parameters) == set(getattr(cosmo, f"_{name}"))

    @pytest.mark.parametrize("name", ["parameters", "derived_parameters"])
    def test_parameters_cannot_set_on_instance(self, cosmo, name):
        """Test descriptor ``parameters`` cannot be set on the instance."""
        with pytest.raises(AttributeError, match=f"cannot set {name!r} of"):
            setattr(cosmo, name, {})
