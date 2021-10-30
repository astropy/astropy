# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.parameter`."""

##############################################################################
# IMPORTS

# STDLIB
import inspect

# THIRD PARTY
import pytest

import numpy as np

# LOCAL
import astropy.units as u
from astropy.cosmology import Cosmology
from astropy.cosmology.core import _COSMOLOGY_CLASSES
from astropy.cosmology.parameter import Parameter, _set_to_float, _set_with_unit

##############################################################################
# TESTS
##############################################################################


class ParameterTestMixin:
    """Tests for a :class:`astropy.cosmology.Parameter` on a Cosmology.

    :class:`astropy.cosmology.Parameter` is a descriptor and this test suite
    tests descriptors by class inheritance, so ``ParameterTestMixin`` is mixed
    into ``TestCosmology`` (tests :class:`astropy.cosmology.Cosmology`).
    """

    @pytest.fixture
    def parameter(self, cosmo_cls):
        """Cosmological Parameters"""
        # I wish this would work
        # yield from {getattr(cosmo_cls, n) for n in cosmo_cls.__parameters__}

        # just return one parameter at random
        yield getattr(cosmo_cls, set(cosmo_cls.__parameters__).pop())

    @pytest.fixture
    def all_parameter(self, cosmo_cls):
        """Cosmological All Parameter instances"""
        # I wish this would work
        # yield from {getattr(cosmo_cls, n) for n in cosmo_cls.__all_parameters__}

        # just return one parameter at random
        yield getattr(cosmo_cls, set(cosmo_cls.__all_parameters__).pop())

    # ===============================================================
    # Method Tests

    def test_Parameter_class_attributes(self, all_parameter):
        """Test :class:`astropy.cosmology.Parameter` attributes on class."""
        # _registry_setters
        assert hasattr(all_parameter, "_registry_setters")
        assert isinstance(all_parameter._registry_setters, dict)
        assert all(isinstance(k, str) for k in all_parameter._registry_setters.keys())
        assert all(callable(v) for v in all_parameter._registry_setters.values())

    def test_Parameter_init(self):
        """Test :class:`astropy.cosmology.Parameter` instantiation."""
        # defaults
        parameter = Parameter()
        assert parameter.fget is None
        assert parameter.fset is _set_with_unit
        assert parameter.unit is None
        assert parameter.equivalencies == []
        assert parameter.format_spec == ".3g"
        assert parameter.derived is False
        assert parameter.__wrapped__ is parameter.fget
        assert parameter.__name__ is None

        # setting all kwargs
        parameter = Parameter(fget=lambda x: x, fset="float", doc="DOCSTRING",
                              unit="km", equivalencies=[u.mass_energy()],
                              fmt=".4f", derived=True)
        assert parameter.fget(2) == 2
        assert parameter.fset is _set_to_float
        assert parameter.unit is u.km
        assert parameter.equivalencies == [u.mass_energy()]
        assert parameter.format_spec == ".4f"
        assert parameter.derived is True
        assert parameter.__wrapped__ is parameter.fget
        assert parameter.__name__ == "<lambda>"

    def test_Parameter_instance_attributes(self, all_parameter):
        """Test :class:`astropy.cosmology.Parameter` attributes from init."""
        # property
        assert hasattr(all_parameter, "fget")
        assert all_parameter.fget is None or callable(all_parameter.fget)

        assert hasattr(all_parameter, "fset")
        assert callable(all_parameter.fset)

        assert hasattr(all_parameter, "fdel")
        assert all_parameter.fdel is None

        assert hasattr(all_parameter, "__doc__")

        # Parameter
        assert hasattr(all_parameter, "_unit")
        assert hasattr(all_parameter, "_equivalencies")
        assert hasattr(all_parameter, "_fmt")
        assert hasattr(all_parameter, "_derived")
        assert hasattr(all_parameter, "__wrapped__")
        assert hasattr(all_parameter, "__name__")

        # __set_name__
        assert hasattr(all_parameter, "_attr_name")
        assert hasattr(all_parameter, "_attr_name_private")

    def test_Parameter_fget(self, all_parameter):
        """Test :attr:`astropy.cosmology.Parameter.fget`."""
        assert hasattr(all_parameter, "fget")
        assert callable(all_parameter.fget) or all_parameter.fget is None

    def test_Parameter_fset(self, all_parameter):
        """Test :attr:`astropy.cosmology.Parameter.fset`."""
        assert hasattr(all_parameter, "fset")
        assert callable(all_parameter.fset)

    def test_Parameter_fdel(self, all_parameter):
        """Test :attr:`astropy.cosmology.Parameter.fdel`."""
        assert all_parameter.fdel is None

    def test_Parameter_name(self, all_parameter):
        """Test :attr:`astropy.cosmology.Parameter.name`."""
        assert hasattr(all_parameter, "name")
        assert isinstance(all_parameter.name, str)
        assert all_parameter.name is all_parameter._attr_name

    def test_Parameter_unit(self, all_parameter):
        """Test :attr:`astropy.cosmology.Parameter.unit`."""
        assert hasattr(all_parameter, "unit")
        assert isinstance(all_parameter.unit, (u.UnitBase, type(None)))
        assert all_parameter.unit is all_parameter._unit

    def test_Parameter_equivalencies(self, all_parameter):
        """Test :attr:`astropy.cosmology.Parameter.equivalencies`."""
        assert hasattr(all_parameter, "equivalencies")
        assert isinstance(all_parameter.equivalencies, (list, u.Equivalency))
        assert all_parameter.equivalencies is all_parameter._equivalencies

    def test_Parameter_format_spec(self, all_parameter):
        """Test :attr:`astropy.cosmology.Parameter.format_spec`."""
        assert hasattr(all_parameter, "format_spec")
        assert isinstance(all_parameter.format_spec, str)
        assert all_parameter.format_spec is all_parameter._fmt

    def test_Parameter_derived(self, cosmo_cls, all_parameter):
        """Test :attr:`astropy.cosmology.Parameter.derived`."""
        assert hasattr(all_parameter, "derived")
        assert isinstance(all_parameter.derived, bool)
        assert all_parameter.derived is all_parameter._derived

        # test value
        if all_parameter.name in cosmo_cls.__parameters__:
            assert all_parameter.derived is False
        else:
            assert all_parameter.derived is True

    # -------------------------------------------
    # descriptor methods

    def test_Parameter_descriptor_get(self, cosmo_cls, cosmo, all_parameter):
        """Test :attr:`astropy.cosmology.Parameter.__get__`."""
        # from class
        parameter = getattr(cosmo_cls, all_parameter.name)
        assert isinstance(parameter, Parameter)
        assert parameter is all_parameter

        # from instance
        parameter = getattr(cosmo, all_parameter.name)
        assert np.all(parameter == getattr(cosmo, all_parameter._attr_name_private))

    def test_Parameter_descriptor_set(self, cosmo, all_parameter):
        """Test :attr:`astropy.cosmology.Parameter.__set__`."""
        # test it's already set
        assert hasattr(cosmo, all_parameter._attr_name_private)

        # and raises an error if set again
        with pytest.raises(AttributeError, match="can't set attribute"):
            setattr(cosmo, all_parameter._attr_name, None)

    def test_Parameter_descriptor_delete(self, cosmo, all_parameter):
        """Test :attr:`astropy.cosmology.Parameter.__delete__`."""
        with pytest.raises(AttributeError, match="can't delete attribute"):
            assert delattr(cosmo, all_parameter._attr_name)

    # -------------------------------------------
    # 'property' descriptor overrides

    # def test_Parameter_getter(self, cosmo):
    # def test_Parameter_setter(self, cosmo):
    # def test_Parameter_deleter(self, cosmo):

    # -------------------------------------------
    # set value
    # tested later.

    # -------------------------------------------

    def test_Parameter_repr(self, cosmo_cls, all_parameter):
        """Test Parameter repr."""
        r = repr(getattr(cosmo_cls, all_parameter.name))

        assert all_parameter._attr_name in r
        assert hex(id(all_parameter)) in r

    # ===============================================================
    # Usage Tests

    def test_Parameter_listed(self, cosmo_cls, all_parameter):
        """Test each `astropy.cosmology.Parameter` attached to Cosmology."""
        # just double check that each entry is a Parameter
        assert isinstance(all_parameter, Parameter)

        # the reverse: check that if it is a Parameter, it's listed.
        # note have to check the more inclusive ``__all_parameters__``
        assert all_parameter.name in cosmo_cls.__all_parameters__
        if not all_parameter.derived:
            assert all_parameter.name in cosmo_cls.__parameters__

    def test_parameter_related_attributes_on_Cosmology(self, cosmo_cls):
        """Test `astropy.cosmology.Parameter`-related on Cosmology."""
        # establish has expected attribute
        assert hasattr(cosmo_cls, "__parameters__")
        assert hasattr(cosmo_cls, "__all_parameters__")

    def test_Parameter_not_unique(self, cosmo_cls, clean_registry):
        """Cosmology Parameter not unique to class when subclass defined."""

        # define subclass to show param is same
        class ExampleBase(cosmo_cls):
            param = Parameter()

        class Example(ExampleBase): pass

        assert Example.param is ExampleBase.param
        assert Example.__parameters__ == ExampleBase.__parameters__

    def test_Parameters_reorder_by_signature(self, cosmo_cls, clean_registry):
        """Test parameters are reordered."""

        class Example(cosmo_cls):
            param = Parameter()

            def __init__(self, param, *, name=None, meta=None):
                pass  # never actually initialized

        # param should be 1st, all other parameters next
        Example.__parameters__[0] == "param"
        # Check the other parameters are as expected.
        # only run this test if "param" is not already on the cosmology
        if cosmo_cls.__parameters__[0] != "param":
            assert set(Example.__parameters__[1:]) == set(cosmo_cls.__parameters__)

    def test_make_from_Parameter(self, cosmo_cls, clean_registry):
        """Test the parameter creation process. Uses ``__set__``."""

        class Example(cosmo_cls):
            param = Parameter(unit=u.eV, equivalencies=u.mass_energy())

            def __init__(self, param, *, name=None, meta=None):
                self.param = param

        assert Example(1).param == 1 * u.eV
        assert Example(1 * u.eV).param == 1 * u.eV
        assert Example(1 * u.J).param == (1 * u.J).to(u.eV)
        assert Example(1 * u.kg).param == (1 * u.kg).to(u.eV, u.mass_energy())


# ========================================================================


class TestParameter(ParameterTestMixin):
    """
    Test `astropy.cosmology.Parameter` directly. Adds a lot of specific tests
    that wouldn't be covered by the per-cosmology tests.
    """

    def setup_class(self):
        class Example1(Cosmology):
            param = Parameter(doc="example parameter",
                              unit=u.m, equivalencies=u.mass_energy())

            def __init__(self, param=15):
                self.param = param

        # with setter
        class Example2(Example1):
            def __init__(self, param=15 * u.m):
                self.param = param

            @Example1.param.setter
            def param(self, param, value):
                return value.to(u.km)

        # attributes
        self.classes = {"Example1": Example1, "Example2": Example2}

    def teardown_class(self):
        for cls in self.classes.values():
            _COSMOLOGY_CLASSES.pop(cls.__qualname__)

    @pytest.fixture(params=["Example1", "Example2"])
    def cosmo_cls(self, request):
        """Cosmology class."""
        return self.classes[request.param]

    @pytest.fixture
    def cosmo(self, cosmo_cls):
        """Cosmology instance"""
        return cosmo_cls()

    @pytest.fixture
    def param(self, cosmo_cls):
        """Get Parameter 'param' from cosmology class."""
        return cosmo_cls.param

    # ==============================================================

    def test_Parameter_instance_attributes(self, param):
        """Test :class:`astropy.cosmology.Parameter` attributes from init."""
        super().test_Parameter_instance_attributes(param)

        # property
        assert param.__doc__ == "example parameter"

        # custom from init
        assert param._unit == u.m
        assert param._equivalencies == u.mass_energy()
        assert param._fmt == ".3g"
        assert param._derived == False

        # custom from set_name
        assert param._attr_name == "param"
        assert param._attr_name_private == "_param"
        assert hasattr(param, "__wrapped__")
        assert hasattr(param, "__name__")

    def test_Parameter_fset(self, cosmo, param):
        """Test :attr:`astropy.cosmology.Parameter.fset`."""
        super().test_Parameter_fset(param)

        value = param.fset(cosmo, param, 1000 * u.m)
        assert value == 1 * u.km

    def test_Parameter_name(self, param):
        """Test :attr:`astropy.cosmology.Parameter.name`."""
        super().test_Parameter_name(param)

        assert param.name == "param"

    def test_Parameter_unit(self, param):
        """Test :attr:`astropy.cosmology.Parameter.unit`."""
        super().test_Parameter_unit(param)

        assert param.unit == u.m

    def test_Parameter_equivalencies(self, param):
        """Test :attr:`astropy.cosmology.Parameter.equivalencies`."""
        super().test_Parameter_equivalencies(param)

        assert param.equivalencies == u.mass_energy()

    def test_Parameter_format_spec(self, param):
        """Test :attr:`astropy.cosmology.Parameter.format_spec`."""
        super().test_Parameter_format_spec(param)

        assert param.format_spec == ".3g"

    def test_Parameter_derived(self, cosmo_cls, param):
        """Test :attr:`astropy.cosmology.Parameter.derived`."""
        super().test_Parameter_derived(cosmo_cls, param)

        assert param.derived is False

    # -------------------------------------------
    # descriptor methods

    def test_Parameter_descriptor_get(self, cosmo_cls, cosmo, param):
        """Test :meth:`astropy.cosmology.Parameter.__get__`."""
        super().test_Parameter_descriptor_get(cosmo_cls, cosmo, param)

        # from instance
        value = getattr(cosmo, param.name)
        assert value == 15 * u.m

    # -------------------------------------------
    # property-style methods

    def test_Parameter_getter(self, param):
        """Test :meth:`astropy.cosmology.Parameter.getter`."""
        with pytest.raises(AttributeError, match="can't create custom Parameter getter."):
            param.getter(None)

    def test_Parameter_setter(self, param):
        """Test :meth:`astropy.cosmology.Parameter.setter`."""
        for k in Parameter._registry_setters:
            newparam = param.setter(k)
            assert newparam.fset == newparam._registry_setters[k]

        # error for non-registered str
        with pytest.raises(ValueError, match="`fset` if str"):
            Parameter(fset="NOT REGISTERED")

        # error if wrong type
        with pytest.raises(TypeError, match="`fset` must be a function or"):
            Parameter(fset=object())

    def test_Parameter_deleter(self, param):
        """Test :meth:`astropy.cosmology.Parameter.deleter`."""
        with pytest.raises(AttributeError, match="can't create custom Parameter deleter."):
            param.deleter(None)

    # -------------------------------------------
    # validation

    def test_Parameter_set(self, cosmo, param):
        """Test :meth:`astropy.cosmology.Parameter.set`."""
        value = param.set(cosmo, 1000 * u.m)

        # whether has custom setter
        if param.fset is param._registry_setters["default"]:
            assert value.unit == u.m
            assert value.value == 1000
        else:
            assert value.unit == u.km
            assert value.value == 1

    def test_Parameter_register_setter(self, param):
        """Test :meth:`astropy.cosmology.Parameter.register_setter`."""
        # already registered
        with pytest.raises(KeyError, match="setter 'default' already"):
            param.__class__.register_setter("default", None)

        # setter not None
        try:
            func = lambda x: x
            setter = param.__class__.register_setter("newsetter", func)
            assert setter is func
        finally:
            param.__class__._registry_setters.pop("newsetter", None)

        # used as decorator
        try:
            @param.__class__.register_setter("newsetter")
            def func(cosmology, param, value):
                return value

            assert param.__class__._registry_setters["newsetter"] is func
        finally:
            param.__class__._registry_setters.pop("newsetter", None)

    # ==============================================================

    def test_Parameter_doesnt_change_with_generic_class(self):
        """Descriptors are initialized once and not updated on subclasses."""

        class ExampleBase:
            def __init__(self, param=15):
                self._param = param

            sig = inspect.signature(__init__)
            _init_signature = sig.replace(parameters=list(sig.parameters.values())[1:])

            param = Parameter(doc="example parameter")

        class Example(ExampleBase): pass

        assert Example.param is ExampleBase.param

    def test_Parameter_doesnt_change_with_cosmology(self, cosmo_cls):
        """Cosmology reinitializes all descriptors when a subclass is defined."""

        # define subclass to show param is same
        class Example(cosmo_cls): pass

        assert Example.param is cosmo_cls.param

        # unregister
        _COSMOLOGY_CLASSES.pop(Example.__qualname__)
        assert Example.__qualname__ not in _COSMOLOGY_CLASSES
