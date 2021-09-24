# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.core`."""

##############################################################################
# IMPORTS

import abc
import inspect
from types import MappingProxyType

import pytest

import numpy as np

import astropy.units as u
from astropy.cosmology import Cosmology, core
from astropy.cosmology.core import _COSMOLOGY_CLASSES, Parameter

from .test_connect import ReadWriteTestMixin, ToFromFormatTestMixin

##############################################################################
# TESTS
##############################################################################


class TestParameter:
    """Test `astropy.cosmology.Parameter`."""

    def setup_class(self):
        class Example1(Cosmology):
            param = Parameter(doc="example parameter",
                              unit=u.m, equivalencies=u.mass_energy())

            def __init__(self, param=15):
                self._param = param

        class Example2(Example1):
            def __init__(self, param=15 * u.m):
                self._param = param.to_value(u.km)

            @Example1.param.getter
            def param(self):
                return self._param << u.km

        self.classes = {"Example1": Example1, "Example2": Example2}

    def teardown_class(self):
        for cls in self.classes.values():
            _COSMOLOGY_CLASSES.pop(cls.__qualname__)

    @pytest.fixture(params=["Example1", "Example2"])
    def cosmo_cls(self, request):
        yield self.classes[request.param]

    @pytest.fixture
    def cosmo(self, cosmo_cls):
        yield cosmo_cls()

    @pytest.fixture
    def parameter(self, cosmo_cls):
        yield cosmo_cls.param

    # ==============================================================

    def test_has_expected_attributes(self, parameter):
        # property
        assert hasattr(parameter, "fget")  # None or callable
        assert parameter.__doc__ == "example parameter"

        # custom from init
        assert parameter._fmt == ".3g"
        assert parameter._unit == u.m
        assert hasattr(parameter, "_equivalencies")

        # custom from set_name
        assert parameter._attr_name == "param"
        assert parameter._attr_name_private == "_param"

    def test_name(self, parameter):
        """Test :attr:`astropy.cosmology.Parameter.name`."""
        assert parameter.name == "param"

    def test_unit(self, parameter):
        """Test :attr:`astropy.cosmology.Parameter.unit`."""
        assert parameter.unit is parameter._unit
        assert parameter.unit == u.m

    def test_equivalencies(self, parameter):
        """Test :attr:`astropy.cosmology.Parameter.equivalencies`."""
        assert parameter.equivalencies == u.mass_energy()

    def test_format_spec(self, parameter):
        """Test :attr:`astropy.cosmology.Parameter.format_spec`."""
        # see test_format for more in-depth tests
        assert parameter.format_spec is parameter._fmt
        assert parameter._fmt == ".3g"

    # -------------------------------------------
    # descriptor methods

    def test_get(self, cosmo_cls):
        """Test :meth:`astropy.cosmology.Parameter.__get__`."""
        # from the class
        parameter = cosmo_cls.param
        assert isinstance(parameter, Parameter)

        # from instance
        value = cosmo_cls().param
        if parameter.fget is None:
            assert value == 15
        else:
            assert value == 15 * u.m

    def test_set(self, cosmo):
        """Test :meth:`astropy.cosmology.Parameter.__set__`."""
        with pytest.raises(AttributeError, match="can't set attribute"):
            cosmo.param = 2

    def test_delete(self, cosmo):
        """Test :meth:`astropy.cosmology.Parameter.__delete__`."""
        with pytest.raises(AttributeError, match="can't delete attribute"):
            del cosmo.param

    # -------------------------------------------
    # property-style methods

    def test_getter_method(self, parameter):
        """Test :meth:`astropy.cosmology.Parameter.getter`."""
        newparam = parameter.getter("NOT NONE")
        assert newparam.fget == "NOT NONE"

    # -------------------------------------------
    # misc

    def test_repr(self, cosmo_cls):
        r = repr(cosmo_cls.param)
        assert f"Parameter 'param'" in r
        assert str(hex(id(cosmo_cls.param))) in r

    # ==============================================================

    def test_parameter_doesnt_change_with_generic_class(self):
        """Descriptors are initialized once and not updated on subclasses."""

        class ExampleBase:
            def __init__(self, param=15):
                self._param = param

            sig = inspect.signature(__init__)
            _init_signature = sig.replace(parameters=list(sig.parameters.values())[1:])

            param = Parameter(doc="example parameter")

        class Example(ExampleBase): pass

        assert Example.param is ExampleBase.param

    def test_parameter_doesnt_change_with_cosmology(self, cosmo_cls):
        """Cosmology reinitializes all descriptors when a subclass is defined."""

        # define subclass to show param is same
        class Example(cosmo_cls): pass

        assert Example.param is cosmo_cls.param

        # unregister
        _COSMOLOGY_CLASSES.pop(Example.__qualname__)
        assert Example.__qualname__ not in _COSMOLOGY_CLASSES


# ========================================================================


class ParameterTestMixin:
    """Tests for a :class:`astropy.cosmology.Parameter` on a Cosmology."""

    def test_Parameters(self, cosmo_cls):
        """Test `astropy.cosmology.Parameter` attached to Cosmology."""
        # first establish has expected attribute
        assert hasattr(cosmo_cls, "__parameters__")

        # check that each entry is a Parameter
        for n in cosmo_cls.__parameters__:
            assert hasattr(cosmo_cls, n)  # checks on the instance
            assert isinstance(getattr(cosmo_cls, n), Parameter)

        # the reverse: check that if it is a Parameter, it's listed.
        for n in dir(cosmo_cls):
            if isinstance(getattr(cosmo_cls, n), Parameter):
                assert n in cosmo_cls.__parameters__

    def test_Parameter_not_unique(self, cosmo_cls, clean_registry):
        """Cosmology reinitializes Parameter when a class is defined."""

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
        assert set(Example.__parameters__[1:]) == set(cosmo_cls.__parameters__)

    def test_make_from_Parameter(self, cosmo_cls, clean_registry):
        """Test the parameter creation process."""

        class Example(cosmo_cls):
            param = Parameter(unit=u.eV, equivalencies=u.mass_energy())

            def __init__(self, param, *, name=None, meta=None):
                cls = self.__class__
                with u.add_enabled_equivalencies(cls.param.equivalencies):
                    self._param = param << cls.param.unit

        assert Example(1).param == 1 * u.eV
        assert Example(1 * u.eV).param == 1 * u.eV
        assert Example(1 * u.J).param == (1 * u.J).to(u.eV)
        assert Example(1 * u.kg).param == (1 * u.kg).to(u.eV, u.mass_energy())


class TestCosmology(ParameterTestMixin, ReadWriteTestMixin, ToFromFormatTestMixin,
                    metaclass=abc.ABCMeta):
    """Test :class:`astropy.cosmology.Cosmology`."""

    def setup_class(self):
        """
        Setup for testing.
        Cosmology should not be instantiated, so tests are done on a subclass.
        """
        class SubCosmology(Cosmology):

            H0 = Parameter(unit=u.km / u.s / u.Mpc)
            Tcmb0 = Parameter(unit=u.K)

            def __init__(self, H0, Tcmb0=0*u.K, name=None, meta=None):
                super().__init__(name=name, meta=meta)
                self._H0 = H0
                self._Tcmb0 = Tcmb0

        self.cls = SubCosmology
        self.cls_args = (70 * (u.km / u.s / u.Mpc), 2.7 * u.K)
        self.cls_kwargs = dict(name=self.__class__.__name__, meta={"a": "b"})

    def teardown_class(self):
        _COSMOLOGY_CLASSES.pop("TestCosmology.setup_class.<locals>.SubCosmology", None)

    @pytest.fixture
    def cosmo_cls(self):
        return self.cls

    @pytest.fixture
    def cosmo(self):
        """The cosmology instance with which to test."""
        return self.cls(*self.cls_args, **self.cls_kwargs)

    # ===============================================================
    # Method & Attribute Tests

    def _cosmo_test_init_attr(self, cosmo):
        """Helper function for testing ``__init__``."""
        assert hasattr(cosmo, "_name")
        assert cosmo._name is None or isinstance(cosmo._name, str)

        assert hasattr(cosmo, "meta")
        assert isinstance(cosmo.meta, dict)

    def test_init(self, cosmo_cls):
        """Test initialization."""
        cosmo1 = cosmo_cls(70, 2.7)
        self._cosmo_test_init_attr(cosmo1)

        # but only "name" and "meta" are used
        cosmo2 = cosmo_cls(70, name="test", meta={"m": 1})
        self._cosmo_test_init_attr(cosmo2)
        assert cosmo2.name == "test"
        assert cosmo2.meta["m"] == 1

    # ---------------------------------------------------------------

    def test_is_equivalent(self, cosmo):
        """
        Test :meth:`astropy.cosmology.Cosmology.is_equivalent`.
        """
        # to self
        assert cosmo.is_equivalent(cosmo)

        # same class, different instance
        newclone = cosmo.clone(name="cloned")
        assert cosmo.is_equivalent(newclone)
        assert newclone.is_equivalent(cosmo)

        # different class
        assert not cosmo.is_equivalent(2)


class CosmologySubclassTest(TestCosmology):
    """
    Test subclasses of :class:`astropy.cosmology.Cosmology`.
    This is broken away from ``TestCosmology``, because |Cosmology| is/will be
    an ABC and subclasses must override some methods.
    """

    @abc.abstractmethod
    def setup_class(self):
        """Setup for testing."""
        pass


# -----------------------------------------------------------------------------


class FlatCosmologyMixinTest:
    """Test :class:`astropy.cosmology.core.FlatCosmologyMixin`."""

    def test_is_equivalent(self, cosmo):
        """Test :meth:`astropy.cosmology.core.FlatCosmologyMixin.is_equivalent`.

        normally this would pass up via super(), but ``__equiv__`` is meant
        to be overridden, so we skip super().
        e.g. FlatFLRWMixinTest -> FlatCosmologyMixinTest -> TestCosmology
        vs   FlatFLRWMixinTest -> FlatCosmologyMixinTest -> TestFLRW -> TestCosmology
        """
        CosmologySubclassTest.test_is_equivalent(self, cosmo)
