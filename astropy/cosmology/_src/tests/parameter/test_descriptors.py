# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology._src.parameter.descriptor`."""

from __future__ import annotations

from types import MappingProxyType
from typing import TYPE_CHECKING, ClassVar

import pytest

from astropy.cosmology._src.parameter import Parameter
from astropy.cosmology._src.parameter.descriptors import ParametersAttribute
from astropy.cosmology._src.utils import all_cls_vars

if TYPE_CHECKING:
    from astropy.cosmology._src.core import Cosmology


class Obj:
    """Example class with a ParametersAttribute."""

    # Attributes that will be accessed by ParametersAttribute when called an instance
    # of this class. On a Cosmology these would be the Parameter objects.
    a: ClassVar[int] = 1
    b: ClassVar[int] = 2
    c: ClassVar[int] = 3
    # The class attribute that is accessed by ParametersAttribute when called on the
    # class. On a Cosmology this would be the mapping of Parameter objects.
    # Here it is just the names of the attributes that will be accessed by the
    # ParametersAttribute to better distinguish between the class and instance
    # attributes.
    _attr_map: ClassVar[tuple[str, ...]] = ("a", "b", "c")

    # The ParametersAttribute descriptor. This will return a mapping of the values of
    # the attributes listed in ``_attr_map`` when called on an instance of this class.
    # When called on the class, it will return ``_attr_map`` itself.
    attr = ParametersAttribute(attr_name="_attr_map")


class TestParametersAttribute:
    """Test the descriptor ``ParametersAttribute``."""

    def test_init(self) -> None:
        """Test constructing a ParametersAttribute."""
        # Proper construction
        attr = ParametersAttribute("attr_name")
        assert attr.attr_name == "attr_name"

        # Improper construction
        # There isn't type checking on the attr_name, so this is allowed, but will fail
        # later when the descriptor is used.
        attr = ParametersAttribute(1)  # type: ignore[arg-type]
        assert attr.attr_name == 1

    def test_get_from_class(self) -> None:
        """Test the descriptor ``__get__`` from the class."""
        assert Obj.attr == ("a", "b", "c")

    def test_get_from_instance(self) -> None:
        """Test the descriptor ``__get__``."""
        obj = Obj()  # Construct an instance for the attribute `attr`.
        assert isinstance(obj.attr, MappingProxyType)
        assert tuple(obj.attr.keys()) == obj._attr_map

    def test_set_from_instance(self) -> None:
        """Test the descriptor ``__set__``."""
        obj = Obj()  # Construct an instance for the attribute `attr`.
        with pytest.raises(AttributeError, match="cannot set 'attr' of"):
            obj.attr = {}

    def test_descriptor_attr_name_not_str(self) -> None:
        """Test when ``attr_name`` is not a string and used as a descriptor.

        This is a regression test for #15882.
        """

        class Obj2(Obj):
            attr = ParametersAttribute(attr_name=None)  # type: ignore[arg-type]

        obj = Obj2()
        with pytest.raises(
            TypeError, match=r"attribute name must be string, not 'NoneType'"
        ):
            _ = obj.attr


##############################################################################


class ParametersAttributeTestMixin:
    """Test the descriptor for ``parameters`` on Cosmology classes.

    This is a mixin class and is mixed into
    :class:`~astropy.cosmology._src.tests.test_core.CosmologyTest`.
    """

    @pytest.mark.parametrize("name", ["parameters", "_derived_parameters"])
    def test_parameters_from_class(self, cosmo_cls: type[Cosmology], name: str) -> None:
        """Test descriptor ``parameters`` accessed from the class."""
        # test presence
        assert hasattr(cosmo_cls, name)
        # test Parameter is a MappingProxyType
        parameters = getattr(cosmo_cls, name)
        assert isinstance(parameters, MappingProxyType)
        # Test items
        assert all(isinstance(p, Parameter) for p in parameters.values())
        assert set(parameters) == {
            k
            for k, v in all_cls_vars(cosmo_cls).items()
            if (isinstance(v, Parameter) and (v.derived == ("derived" in name)))
        }

    @pytest.mark.parametrize("name", ["parameters", "_derived_parameters"])
    def test_parameters_from_instance(self, cosmo: Cosmology, name: str) -> None:
        """Test descriptor ``parameters`` accessed from the instance."""
        # test presence
        assert hasattr(cosmo, name)
        # test Parameter is a MappingProxyType
        parameters = getattr(cosmo, name)
        assert isinstance(parameters, MappingProxyType)
        # Test keys
        assert set(parameters) == {
            k
            for k, v in all_cls_vars(cosmo).items()
            if (isinstance(v, Parameter) and (v.derived == ("derived" in name)))
        }

    @pytest.mark.parametrize("name", ["parameters", "_derived_parameters"])
    def test_parameters_cannot_set_on_instance(
        self, cosmo: Cosmology, name: str
    ) -> None:
        """Test descriptor ``parameters`` cannot be set on the instance."""
        with pytest.raises(AttributeError, match=f"cannot assign to field {name!r}"):
            setattr(cosmo, name, {})
