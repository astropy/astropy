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
from astropy.cosmology.parameter import Parameter, _validate_to_float, _validate_with_unit

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
        # _registry_validators
        assert hasattr(all_parameter, "_registry_validators")
        assert isinstance(all_parameter._registry_validators, dict)
        assert all(isinstance(k, str) for k in all_parameter._registry_validators.keys())
        assert all(callable(v) for v in all_parameter._registry_validators.values())

    def test_Parameter_init(self):
        """Test :class:`astropy.cosmology.Parameter` instantiation."""
        # defaults
        parameter = Parameter()
        assert parameter.fvalidate is _validate_with_unit
        assert parameter.unit is None
        assert parameter.equivalencies == []
        assert parameter.format_spec == ".3g"
        assert parameter.derived is False

        # setting all kwargs
        parameter = Parameter(fvalidate="float", doc="DOCSTRING",
                              unit="km", equivalencies=[u.mass_energy()],
                              fmt=".4f", derived=True)
        assert parameter.fvalidate is _validate_to_float
        assert parameter.unit is u.km
        assert parameter.equivalencies == [u.mass_energy()]
        assert parameter.format_spec == ".4f"
        assert parameter.derived is True

    def test_Parameter_instance_attributes(self, all_parameter):
        """Test :class:`astropy.cosmology.Parameter` attributes from init."""
        assert hasattr(all_parameter, "fvalidate")
        assert callable(all_parameter.fvalidate)

        assert hasattr(all_parameter, "__doc__")

        # Parameter
        assert hasattr(all_parameter, "_unit")
        assert hasattr(all_parameter, "_equivalencies")
        assert hasattr(all_parameter, "_fmt")
        assert hasattr(all_parameter, "_derived")

        # __set_name__
        assert hasattr(all_parameter, "_attr_name")
        assert hasattr(all_parameter, "_attr_name_private")

    def test_Parameter_fvalidate(self, all_parameter):
        """Test :attr:`astropy.cosmology.Parameter.fvalidate`."""
        assert hasattr(all_parameter, "fvalidate")
        assert callable(all_parameter.fvalidate)

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

    # -------------------------------------------
    # validate value
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

        # with validator
        class Example2(Example1):
            def __init__(self, param=15 * u.m):
                self.param = param

            @Example1.param.validator
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

    def test_Parameter_fvalidate(self, cosmo, param):
        """Test :attr:`astropy.cosmology.Parameter.fvalidate`."""
        super().test_Parameter_fvalidate(param)

        value = param.fvalidate(cosmo, param, 1000 * u.m)
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
    # validation

    def test_Parameter_validator(self, param):
        """Test :meth:`astropy.cosmology.Parameter.validator`."""
        for k in Parameter._registry_validators:
            newparam = param.validator(k)
            assert newparam.fvalidate == newparam._registry_validators[k]

        # error for non-registered str
        with pytest.raises(ValueError, match="`fvalidate`, if str"):
            Parameter(fvalidate="NOT REGISTERED")

        # error if wrong type
        with pytest.raises(TypeError, match="`fvalidate` must be a function or"):
            Parameter(fvalidate=object())

    def test_Parameter_validate(self, cosmo, param):
        """Test :meth:`astropy.cosmology.Parameter.validate`."""
        value = param.validate(cosmo, 1000 * u.m)

        # whether has custom validator
        if param.fvalidate is param._registry_validators["default"]:
            assert value.unit == u.m
            assert value.value == 1000
        else:
            assert value.unit == u.km
            assert value.value == 1

    def test_Parameter_register_validator(self, param):
        """Test :meth:`astropy.cosmology.Parameter.register_validator`."""
        # already registered
        with pytest.raises(KeyError, match="validator 'default' already"):
            param.__class__.register_validator("default", None)

        # validator not None
        try:
            func = lambda x: x
            validator = param.__class__.register_validator("newvalidator", func)
            assert validator is func
        finally:
            param.__class__._registry_validators.pop("newvalidator", None)

        # used as decorator
        try:
            @param.__class__.register_validator("newvalidator")
            def func(cosmology, param, value):
                return value

            assert param.__class__._registry_validators["newvalidator"] is func
        finally:
            param.__class__._registry_validators.pop("newvalidator", None)

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
