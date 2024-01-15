# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.parameter._descriptor`."""

from types import MappingProxyType
from typing import ClassVar

import pytest

from astropy.cosmology._utils import all_cls_vars
from astropy.cosmology.parameter import Parameter
from astropy.cosmology.parameter._descriptors import ParametersAttribute


class Obj:
    """Example class with a ParametersAttribute."""

    # Attributes that will be accessed by ParametersAttribute when called an instance
    # of this class. On a Cosmology these would be the Parameter objects.
    a: ClassVar = 1
    b: ClassVar = 2
    c: ClassVar = 3
    # The class attribute that is accessed by ParametersAttribute when called on the
    # class. On a Cosmology this would be the mapping of Parameter objects.
    # Here it is just the names of the attributes that will be accessed by the
    # ParametersAttribute to better distinguish between the class and instance
    # attributes.
    _attr_map: ClassVar = ("a", "b", "c")

    # The ParametersAttribute descriptor. This will return a mapping of the values of
    # the attributes listed in ``_attr_map`` when called on an instance of this class.
    # When called on the class, it will return ``_attr_map`` itself.
    attr = ParametersAttribute(attr_name="_attr_map")


class TestParametersAttribute:
    """Test the descriptor ``ParametersAttribute``."""

    @pytest.fixture
    def obj_cls(self) -> type[Obj]:
        """The class with the ParametersAttribute."""
        return Obj

    @pytest.fixture
    def obj(self, obj_cls) -> Obj:
        """An instance of the class with the ParametersAttribute."""
        return obj_cls()

    # ===============================================================

    def test_init(self):
        """Test constructing a ParametersAttribute."""
        attr = ParametersAttribute("attr_name")
        assert attr.attr_name == "attr_name"

    def test_get_from_class(self, obj_cls):
        """Test the descriptor ``__get__``."""
        assert obj_cls.attr == ("a", "b", "c")

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
