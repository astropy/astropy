# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.parameter`."""

from collections.abc import Callable

import numpy as np
import pytest

import astropy.units as u
from astropy.cosmology import Cosmology, Parameter
from astropy.cosmology.core import _COSMOLOGY_CLASSES, dataclass_decorator
from astropy.cosmology.parameter._converter import (
    _REGISTRY_FVALIDATORS,
    _validate_with_unit,
)
from astropy.cosmology.parameter._core import MISSING

##############################################################################


def test_registry_validators():
    """Test :class:`astropy.cosmology.Parameter` attributes on class."""
    # _registry_validators
    assert isinstance(_REGISTRY_FVALIDATORS, dict)
    assert all(isinstance(k, str) for k in _REGISTRY_FVALIDATORS.keys())
    assert all(callable(v) for v in _REGISTRY_FVALIDATORS.values())


class Test_Parameter:
    """Test :class:`astropy.cosmology.Parameter` not on a cosmology."""

    @pytest.mark.parametrize(
        "kwargs",
        [
            {},
            dict(
                default=1.0,
                fvalidate="float",
                doc="DOCSTRING",
                unit="km",
                equivalencies=[u.mass_energy()],
                derived=True,
            ),
        ],
    )
    def test_Parameter_init(self, kwargs):
        """Test :class:`astropy.cosmology.Parameter` instantiation."""
        unit = kwargs.get("unit")

        param = Parameter(**kwargs)
        assert param.default == kwargs.get("default", MISSING)
        assert param.fvalidate is _REGISTRY_FVALIDATORS.get(
            kwargs.get("fvalidate"), _validate_with_unit
        )
        assert param.doc == kwargs.get("doc")
        assert param.unit is (u.Unit(unit) if unit is not None else None)
        assert param.equivalencies == kwargs.get("equivalencies", [])
        assert param.derived is kwargs.get("derived", False)
        assert param.name is None

    def test_Parameter_default(self):
        """Test :attr:`astropy.cosmology.Parameter.default`."""
        parameter = Parameter()
        assert parameter.default is MISSING
        assert repr(parameter.default) == "<MISSING>"


class ParameterTestMixin:
    """Tests for a :class:`astropy.cosmology.Parameter` on a Cosmology.

    :class:`astropy.cosmology.Parameter` is a descriptor and this test suite
    tests descriptors by class inheritance, so ``ParameterTestMixin`` is mixed
    into ``TestCosmology`` (tests :class:`astropy.cosmology.Cosmology`).
    """

    @pytest.fixture
    def parameter(self, cosmo_cls):
        """Cosmological Parameters"""
        yield from cosmo_cls.parameters.values()

    @pytest.fixture
    def all_parameter(self, cosmo_cls):
        """Cosmological All Parameter instances"""
        # just return one parameter at random
        n = set(cosmo_cls._parameters_all).pop()
        try:
            yield cosmo_cls.parameters[n]
        except KeyError:
            yield cosmo_cls._derived_parameters[n]

    # ===============================================================
    # Method Tests

    def test_Parameter_instance_attributes(self, all_parameter):
        """Test :class:`astropy.cosmology.Parameter` attributes from init."""
        assert hasattr(all_parameter, "fvalidate")
        assert callable(all_parameter.fvalidate)

        assert hasattr(all_parameter, "__doc__")

        # Parameter
        assert hasattr(all_parameter, "_unit")
        assert hasattr(all_parameter, "equivalencies")
        assert hasattr(all_parameter, "derived")

        # __set_name__
        assert hasattr(all_parameter, "name")
        assert hasattr(all_parameter, "_attr_name")

    def test_Parameter_fvalidate(self, all_parameter):
        """Test :attr:`astropy.cosmology.Parameter.fvalidate`."""
        assert hasattr(all_parameter, "fvalidate")
        assert callable(all_parameter.fvalidate)
        assert hasattr(all_parameter, "_fvalidate_in")
        assert isinstance(all_parameter._fvalidate_in, (str, Callable))

    def test_Parameter_name(self, all_parameter):
        """Test :attr:`astropy.cosmology.Parameter.name`."""
        assert hasattr(all_parameter, "name")
        assert isinstance(all_parameter.name, str)
        assert all_parameter._attr_name == f"_{all_parameter.name}"

    def test_Parameter_unit(self, all_parameter):
        """Test :attr:`astropy.cosmology.Parameter.unit`."""
        assert hasattr(all_parameter, "unit")
        assert isinstance(all_parameter.unit, (u.UnitBase, type(None)))
        assert all_parameter.unit is all_parameter._unit

    def test_Parameter_equivalencies(self, all_parameter):
        """Test :attr:`astropy.cosmology.Parameter.equivalencies`."""
        assert hasattr(all_parameter, "equivalencies")
        assert isinstance(all_parameter.equivalencies, (list, u.Equivalency))

    def test_Parameter_derived(self, cosmo_cls, all_parameter):
        """Test :attr:`astropy.cosmology.Parameter.derived`."""
        assert hasattr(all_parameter, "derived")
        assert isinstance(all_parameter.derived, bool)

        # test value
        assert all_parameter.derived is (all_parameter.name not in cosmo_cls.parameters)

    def test_Parameter_default(self, cosmo_cls, all_parameter):
        """Test :attr:`astropy.cosmology.Parameter.default`."""
        assert hasattr(all_parameter, "default")
        assert all_parameter.default is MISSING or isinstance(
            all_parameter.default, (type(None), int, float, u.Quantity)
        )

    # -------------------------------------------
    # descriptor methods

    def test_Parameter_descriptor_get(self, cosmo_cls, cosmo, all_parameter):
        """Test :attr:`astropy.cosmology.Parameter.__get__`."""
        # from class
        np.testing.assert_array_equal(
            getattr(cosmo_cls, all_parameter.name).default, all_parameter.default
        )

        # from instance
        parameter = getattr(cosmo, all_parameter.name)
        assert np.all(parameter == getattr(cosmo, all_parameter._attr_name))

    def test_Parameter_descriptor_set(self, cosmo, all_parameter):
        """Test :attr:`astropy.cosmology.Parameter.__set__`."""
        # test it's already set
        assert hasattr(cosmo, all_parameter._attr_name)

    # -------------------------------------------
    # validate value
    # tested later.

    # ===============================================================
    # Usage Tests

    def test_Parameter_listed(self, cosmo_cls, all_parameter):
        """Test each `astropy.cosmology.Parameter` attached to Cosmology."""
        # just double check that each entry is a Parameter
        assert isinstance(all_parameter, Parameter)

        # the reverse: check that if it is a Parameter, it's listed.
        if all_parameter.derived:
            assert all_parameter.name in cosmo_cls._derived_parameters
        else:
            assert all_parameter.name in cosmo_cls.parameters


# ========================================================================


class TestParameter(ParameterTestMixin):
    """
    Test `astropy.cosmology.Parameter` directly. Adds a lot of specific tests
    that wouldn't be covered by the per-cosmology tests.
    """

    def setup_class(self):
        theparam = Parameter(
            default=15,
            doc="Description of example parameter.",
            unit=u.m,
            equivalencies=u.mass_energy(),
        )

        @dataclass_decorator
        class Example1(Cosmology):
            param: Parameter = theparam.clone()

            @property
            def is_flat(self):
                return super().is_flat()

        # with validator
        @dataclass_decorator
        class Example2(Example1):
            param: Parameter = theparam.clone(default=15 * u.m)

            @param.validator
            def param(self, param, value):
                return value.to(u.km)

        # attributes
        self.classes = {"Example1": Example1, "Example2": Example2}

    def teardown_class(self):
        for cls in self.classes.values():
            _COSMOLOGY_CLASSES.pop(cls.__qualname__, None)

    @pytest.fixture(scope="class", params=["Example1", "Example2"])
    def cosmo_cls(self, request):
        """Cosmology class."""
        return self.classes[request.param]

    @pytest.fixture(scope="class")
    def cosmo(self, cosmo_cls):
        """Cosmology instance"""
        return cosmo_cls()

    @pytest.fixture(scope="class")
    def param(self, cosmo_cls):
        """Get Parameter 'param' from cosmology class."""
        return cosmo_cls.parameters["param"]

    @pytest.fixture(scope="class")
    def param_cls(self, param):
        """Get Parameter class from cosmology class."""
        return type(param)

    # ==============================================================

    def test_Parameter_instance_attributes(self, param):
        """Test :class:`astropy.cosmology.Parameter` attributes from init."""
        super().test_Parameter_instance_attributes(param)

        # property
        assert param.__doc__ == "Description of example parameter."

        # custom from init
        assert param.unit == u.m
        assert param.equivalencies == u.mass_energy()
        assert param.derived == np.False_

        # custom from set_name
        assert param.name == "param"
        assert param._attr_name == "_param"

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
        for k in _REGISTRY_FVALIDATORS:
            newparam = param.validator(k)
            assert newparam.fvalidate == _REGISTRY_FVALIDATORS[k]

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
        if param.fvalidate is _REGISTRY_FVALIDATORS["default"]:
            assert value.unit == u.m
            assert value.value == 1000
        else:
            assert value.unit == u.km
            assert value.value == 1

    def test_Parameter_register_validator(self, param_cls):
        """Test :meth:`astropy.cosmology.Parameter.register_validator`."""
        # already registered
        with pytest.raises(KeyError, match="validator 'default' already"):
            param_cls.register_validator("default", None)

        # validator not None
        def notnonefunc(x):
            return x

        try:
            validator = param_cls.register_validator("newvalidator", notnonefunc)
            assert validator is notnonefunc
        finally:
            _REGISTRY_FVALIDATORS.pop("newvalidator", None)

        # used as decorator
        try:

            @param_cls.register_validator("newvalidator")
            def func(cosmology, param, value):
                return value

            assert _REGISTRY_FVALIDATORS["newvalidator"] is func
        finally:
            _REGISTRY_FVALIDATORS.pop("newvalidator", None)

    # -------------------------------------------

    def test_Parameter_clone(self, param):
        """Test :meth:`astropy.cosmology.Parameter.clone`."""
        # this implicitly relies on `__eq__` testing properly. Which is tested.

        # basic test that nothing changes
        assert param.clone() == param
        assert param.clone() is not param  # but it's not a 'singleton'

        # passing kwargs will change stuff
        newparam = param.clone(unit="km/(yr sr)")
        assert newparam.unit == u.km / u.yr / u.sr
        assert param.unit != u.km / u.yr / u.sr  # original is unchanged

        # expected failure for not-an-argument
        with pytest.raises(TypeError):
            param.clone(not_a_valid_parameter=True)

    # -------------------------------------------

    def test_Parameter_equality(self):
        """Test Parameter equality.

        Determined from the processed initialization args (including defaults).
        """
        p1 = Parameter(unit="km / (s Mpc)")
        p2 = Parameter(unit="km / (s Mpc)")
        assert p1 == p2

        # not equal parameters
        p3 = Parameter(unit="km / s")
        assert p3 != p1

        # misc
        assert p1 != 2  # show doesn't error

    # -------------------------------------------

    def test_Parameter_repr(self, cosmo_cls, param):
        """Test Parameter repr."""
        r = repr(param)

        assert "Parameter(" in r
        for subs in (
            "derived=False",
            'unit=Unit("m")',
            'equivalencies=[(Unit("kg"), Unit("J")',
            "doc='Description of example parameter.'",
        ):
            assert subs in r, subs

        # `fvalidate` is a little tricker b/c one of them is custom!
        if param.fvalidate in _REGISTRY_FVALIDATORS.values():  # not custom
            assert "fvalidate='default'" in r
        else:
            assert "fvalidate=<" in r  # Some function, don't care about details.

    def test_Parameter_repr_roundtrip(self, param):
        """Test ``eval(repr(Parameter))`` can round trip to ``Parameter``."""
        P = Parameter(doc="A description of this parameter.", derived=True)
        NP = eval(repr(P))  # Evaluate string representation back into a param.

        assert P == NP

    # ========================================================================

    def test_make_from_Parameter(self, cosmo_cls, clean_registry):
        """Test the parameter creation process. Uses ``__set__``."""

        @dataclass_decorator
        class Example(cosmo_cls):
            param: Parameter = Parameter(unit=u.eV, equivalencies=u.mass_energy())

            @property
            def is_flat(self):
                return super().is_flat()

        assert Example(1).param == 1 * u.eV
        assert Example(1 * u.eV).param == 1 * u.eV
        assert Example(1 * u.J).param == (1 * u.J).to(u.eV)
        assert Example(1 * u.kg).param == (1 * u.kg).to(u.eV, u.mass_energy())
