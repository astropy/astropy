# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.parameter._descriptor`."""

from types import MappingProxyType
from typing import ClassVar

import pytest

from astropy.cosmology._utils import all_cls_vars
from astropy.cosmology.parameter import Parameter
from astropy.cosmology.parameter._descriptors import ParametersAttribute


class Obj:
    a: ClassVar = 1
    b: ClassVar = 2
    c: ClassVar = 3
    _attr_map: ClassVar = ("a", "b", "c")

    attr = ParametersAttribute(attr_name="_attr_map")


class TestParametersAttribute:
    """Test the descriptor ``ParametersAttribute``."""

    @pytest.fixture
    def attr_name(self) -> str:
        return "attr"

    @pytest.fixture
    def obj_class(self) -> type[Obj]:
        return Obj

    @pytest.fixture
    def obj(self, obj_class) -> Obj:
        return obj_class()

    # ===============================================================

    def test_init(self):
        """Test constructing a ParametersAttribute."""
        attr = ParametersAttribute("attr_name")
        assert attr.attr_name == "attr_name"

    def test_get_from_class(self, obj_class):
        """Test the descriptor ``__get__``."""
        assert obj_class.attr == ("a", "b", "c")

    def test_get_from_instance(self, obj):
        """Test the descriptor ``__get__``."""
        assert isinstance(obj.attr, MappingProxyType)
        assert tuple(obj.attr.keys()) == obj._attr_map

    def test_set_from_instance(self, obj):
        """Test the descriptor ``__set__``."""
        with pytest.raises(AttributeError, match="cannot set 'attr' of"):
            obj.attr = {}


##############################################################################


class ParametersAttributeTestMixin:
    """Test the descriptor for ``parameters`` on Cosmology classes.

    This is a mixin class and is mixed into
    :class:`~astropy.cosmology.tests.test_core.CosmologyTest`.
    """

    @pytest.mark.parametrize("name", ["parameters", "_derived_parameters"])
    def test_parameters_from_class(self, cosmo_cls, name):
        """Test descriptor ``parameters`` accessed from the class."""
        descriptor = all_cls_vars(cosmo_cls)[name]
        # test presence
        assert hasattr(cosmo_cls, name)
        # test Parameter is a MappingProxyType
        parameters = getattr(cosmo_cls, name)
        assert isinstance(parameters, MappingProxyType)
        # Test items
        assert set(parameters.keys()) == set(getattr(cosmo_cls, descriptor.attr_name))
        assert all(isinstance(p, Parameter) for p in parameters.values())

    @pytest.mark.parametrize("name", ["parameters", "_derived_parameters"])
    def test_parameters_from_instance(self, cosmo, name):
        """Test descriptor ``parameters`` accessed from the instance."""
        descriptor = all_cls_vars(cosmo)[name]
        # test presence
        assert hasattr(cosmo, name)
        # test Parameter is a MappingProxyType
        parameters = getattr(cosmo, name)
        assert isinstance(parameters, MappingProxyType)
        # Test keys
        assert set(parameters) == set(getattr(cosmo, descriptor.attr_name))

    @pytest.mark.parametrize("name", ["parameters", "_derived_parameters"])
    def test_parameters_cannot_set_on_instance(self, cosmo, name):
        """Test descriptor ``parameters`` cannot be set on the instance."""
        with pytest.raises(AttributeError, match=f"cannot set {name!r} of"):
            setattr(cosmo, name, {})
