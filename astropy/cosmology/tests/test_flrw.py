# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.flrw`."""

##############################################################################
# IMPORTS

# STDLIB
import abc
import copy

# THIRD PARTY
import pytest

import numpy as np

# LOCAL
import astropy.cosmology.units as cu
import astropy.units as u
import astropy.constants as const
from astropy.cosmology import (FLRW, FlatLambdaCDM, Flatw0waCDM, FlatwCDM,
                               LambdaCDM, Planck18, w0waCDM, w0wzCDM, wCDM, wpwaCDM)
from astropy.cosmology.core import _COSMOLOGY_CLASSES
from astropy.cosmology.flrw import ellipkinc, hyp2f1, quad
from astropy.cosmology.parameter import Parameter
from astropy.utils.compat.optional_deps import HAS_SCIPY

from .test_core import CosmologySubclassTest as CosmologyTest
from .test_core import FlatCosmologyMixinTest, ParameterTestMixin

##############################################################################
# TESTS
##############################################################################


@pytest.mark.skipif(HAS_SCIPY, reason="scipy is installed")
def test_optional_deps_functions():
    """Test stand-in functions when optional dependencies not installed."""
    with pytest.raises(ModuleNotFoundError, match="No module named 'scipy.integrate'"):
        quad()

    with pytest.raises(ModuleNotFoundError, match="No module named 'scipy.special'"):
        ellipkinc()

    with pytest.raises(ModuleNotFoundError, match="No module named 'scipy.special'"):
        hyp2f1()


##############################################################################


class ParameterH0TestMixin(ParameterTestMixin):

    def test_H0(self, cosmo_cls, cosmo):
        """Test Parameter ``H0``."""
        unit = u.Unit("km/(s Mpc)")

        # on the class
        assert isinstance(cosmo_cls.H0, Parameter)
        assert "Hubble constant" in cosmo_cls.H0.__doc__
        assert cosmo_cls.H0.unit == unit

        # validation
        assert cosmo_cls.H0.validate(cosmo, 1) == 1 * unit
        assert cosmo_cls.H0.validate(cosmo, 10 * unit) == 10 * unit
        with pytest.raises(ValueError, match="H0 is a non-scalar quantity"):
            cosmo_cls.H0.validate(cosmo, [1, 2])

        # on the instance
        assert cosmo.H0 is cosmo._H0
        assert cosmo.H0 == self._cls_args["H0"]
        assert isinstance(cosmo.H0, u.Quantity) and cosmo.H0.unit == unit

    def test_init_H0(self, cosmo_cls):
        """Test initialization for values of ``H0``."""
        ba = cosmo_cls._init_signature.bind(*self.cls_args, **self.cls_kwargs)

        ba.arguments["H0"] = u.Quantity([70, 100], u.km / u.s / u.Mpc)
        with pytest.raises(ValueError, match="H0 is a non-scalar quantity"):
            cosmo_cls(*ba.args, **ba.kwargs)


class ParameterOm0TestMixin(ParameterTestMixin):

    def test_Om0(self, cosmo_cls, cosmo):
        """Test Parameter ``Om0``."""
        # on the class
        assert isinstance(cosmo_cls.Om0, Parameter)
        assert "Omega matter" in cosmo_cls.Om0.__doc__
        assert cosmo_cls.Tcmb0.format_spec == "0.4g"

        # validation
        assert cosmo_cls.Om0.validate(cosmo, 1) == 1
        assert cosmo_cls.Om0.validate(cosmo, 10 * u.one) == 10
        with pytest.raises(ValueError, match="Om0 cannot be negative"):
            cosmo_cls.Om0.validate(cosmo, -1)

        # on the instance
        assert cosmo.Om0 is cosmo._Om0
        assert cosmo.Om0 == self._cls_args["Om0"]
        assert isinstance(cosmo.Om0, float)

    def test_init_Om0(self, cosmo_cls):
        """Test initialization for values of ``Om0``."""
        ba = cosmo_cls._init_signature.bind(*self.cls_args, **self.cls_kwargs)

        ba.arguments["Om0"] = -0.27
        with pytest.raises(ValueError, match="Om0 cannot be negative."):
            cosmo_cls(*ba.args, **ba.kwargs)


class ParameterOde0TestMixin(ParameterTestMixin):

    def test_Ode0(self, cosmo_cls, cosmo):
        """Test Parameter ``Ode0``."""
        # on the class
        assert isinstance(cosmo_cls.Ode0, Parameter)
        assert "Omega dark energy" in cosmo_cls.Ode0.__doc__

        # validation
        assert cosmo_cls.Ode0.validate(cosmo, 1.1) == 1.1
        assert cosmo_cls.Ode0.validate(cosmo, 10 * u.one) == 10.0
        with pytest.raises(TypeError, match="only dimensionless"):
            cosmo_cls.Ode0.validate(cosmo, 10 * u.km)

        # on the instance
        # if Ode0 is a parameter, test its value
        assert cosmo.Ode0 is cosmo._Ode0
        assert cosmo.Ode0 == self._cls_args["Ode0"]
        assert isinstance(cosmo.Ode0, float)

    def test_init_Ode0(self, cosmo_cls):
        """Test initialization for values of ``Ode0``."""
        ba = cosmo_cls._init_signature.bind(*self.cls_args, **self.cls_kwargs)

        ba.arguments["Ode0"] = 10 * u.km
        with pytest.raises(TypeError, match="only dimensionless"):
            cosmo_cls(*ba.args, **ba.kwargs)


class ParameterTcmb0TestMixin(ParameterTestMixin):

    def test_Tcmb0(self, cosmo_cls, cosmo):
        """Test Parameter ``Tcmb0``."""
        # on the class
        assert isinstance(cosmo_cls.Tcmb0, Parameter)
        assert "Temperature of the CMB" in cosmo_cls.Tcmb0.__doc__
        assert cosmo_cls.Tcmb0.unit == u.K

        # validation
        assert cosmo_cls.Tcmb0.validate(cosmo, 1) == 1 * u.K
        assert cosmo_cls.Tcmb0.validate(cosmo, 10 * u.K) == 10 * u.K
        with pytest.raises(ValueError, match="Tcmb0 is a non-scalar quantity"):
            cosmo_cls.Tcmb0.validate(cosmo, [1, 2])

        # on the instance
        assert cosmo.Tcmb0 is cosmo._Tcmb0
        assert cosmo.Tcmb0 == self.cls_kwargs["Tcmb0"]
        assert isinstance(cosmo.Tcmb0, u.Quantity) and cosmo.Tcmb0.unit == u.K

    def test_init_Tcmb0(self, cosmo_cls):
        """Test initialization for values of ``Tcmb0``."""
        ba = cosmo_cls._init_signature.bind(*self.cls_args, **self.cls_kwargs)

        ba.arguments["Tcmb0"] = u.Quantity([0.0, 2], u.K)
        with pytest.raises(ValueError, match="Tcmb0 is a non-scalar quantity"):
            cosmo_cls(*ba.args, **ba.kwargs)


class ParameterNeffTestMixin(ParameterTestMixin):

    def test_Neff(self, cosmo_cls, cosmo):
        """Test Parameter ``Neff``."""
        # on the class
        assert isinstance(cosmo_cls.Neff, Parameter)
        assert "Number of effective neutrino species" in cosmo_cls.Neff.__doc__

        # validation
        assert cosmo_cls.Neff.validate(cosmo, 1) == 1
        assert cosmo_cls.Neff.validate(cosmo, 10 * u.one) == 10
        with pytest.raises(ValueError, match="Neff cannot be negative"):
            cosmo_cls.Neff.validate(cosmo, -1)

        # on the instance
        assert cosmo.Neff is cosmo._Neff
        assert cosmo.Neff == self.cls_kwargs.get("Neff", 3.04)
        assert isinstance(cosmo.Neff, float)

    def test_init_Neff(self, cosmo_cls):
        """Test initialization for values of ``Neff``."""
        ba = cosmo_cls._init_signature.bind(*self.cls_args, **self.cls_kwargs)

        ba.arguments["Neff"] = -1
        with pytest.raises(ValueError):
            cosmo_cls(*ba.args, **ba.kwargs)


class Parameterm_nuTestMixin(ParameterTestMixin):

    def test_m_nu(self, cosmo_cls, cosmo):
        """Test Parameter ``m_nu``."""
        # on the class
        assert isinstance(cosmo_cls.m_nu, Parameter)
        assert "Mass of neutrino species" in cosmo_cls.m_nu.__doc__
        assert cosmo_cls.m_nu.unit == u.eV
        assert cosmo_cls.m_nu.equivalencies == u.mass_energy()
        assert cosmo_cls.m_nu.format_spec == ""

        # on the instance
        # assert cosmo.m_nu is cosmo._m_nu
        assert u.allclose(cosmo.m_nu, [0.0, 0.0, 0.0] * u.eV)

        # set differently depending on the other inputs
        if cosmo.Tnu0.value == 0:
            assert cosmo.m_nu is None
        elif not cosmo._massivenu:  # only massless
            assert u.allclose(cosmo.m_nu, 0 * u.eV)
        elif self._nmasslessnu == 0:  # only massive
            assert cosmo.m_nu == cosmo._massivenu_mass
        else:  # a mix -- the most complicated case
            assert u.allclose(cosmo.m_nu[:self._nmasslessnu], 0 * u.eV)
            assert u.allclose(cosmo.m_nu[self._nmasslessnu], cosmo._massivenu_mass)

    def test_init_m_nu(self, cosmo_cls):
        """Test initialization for values of ``m_nu``."""
        ba = cosmo_cls._init_signature.bind(*self.cls_args, **self.cls_kwargs)

        # negative m_nu
        tba = copy.deepcopy(ba)
        tba.arguments["m_nu"] = u.Quantity([-0.3, 0.2, 0.1], u.eV)
        with pytest.raises(ValueError):
            cosmo_cls(*tba.args, **tba.kwargs)

        # mismatch with Neff and Tcmb0
        tba = copy.deepcopy(ba)
        tba.arguments["Tcmb0"] = 3
        tba.arguments["Neff"] = 2
        tba.arguments["m_nu"] = u.Quantity([0.15, 0.2, 0.1], u.eV)
        with pytest.raises(ValueError):
            cosmo_cls(*tba.args, **tba.kwargs)

        # wrong number of neutrinos
        tba = copy.deepcopy(ba)
        tba.arguments["Tcmb0"] = 3
        tba.arguments["m_nu"] = u.Quantity([-0.3, 0.2], u.eV)  # 2, expecting 3
        with pytest.raises(ValueError):
            cosmo_cls(*tba.args, **tba.kwargs)

        # TODO! transfer tests for initializing neutrinos


class ParameterOb0TestMixin(ParameterTestMixin):

    def test_Ob0(self, cosmo_cls, cosmo):
        """Test Parameter ``Ob0``."""
        # on the class
        assert isinstance(cosmo_cls.Ob0, Parameter)
        assert "Omega baryon;" in cosmo_cls.Ob0.__doc__

        # validation
        assert cosmo_cls.Ob0.validate(cosmo, None) is None
        assert cosmo_cls.Ob0.validate(cosmo, 0.1) == 0.1
        assert cosmo_cls.Ob0.validate(cosmo, 0.1 * u.one) == 0.1
        with pytest.raises(ValueError, match="Ob0 cannot be negative"):
            cosmo_cls.Ob0.validate(cosmo, -1)
        with pytest.raises(ValueError, match="baryonic density can not be larger"):
            cosmo_cls.Ob0.validate(cosmo, cosmo.Om0 + 1)

        # on the instance
        assert cosmo.Ob0 is cosmo._Ob0
        assert cosmo.Ob0 is None

    def test_init_Ob0(self, cosmo_cls):
        """Test initialization for values of ``Ob0``."""
        ba = cosmo_cls._init_signature.bind(*self.cls_args, **self.cls_kwargs)

        # negative Ob0
        tba = copy.deepcopy(ba)
        tba.arguments["Ob0"] = -0.04
        with pytest.raises(ValueError, match="Ob0 cannot be negative"):
            cosmo_cls(*tba.args, **tba.kwargs)

        # Ob0 > Om0
        tba.arguments["Ob0"] = tba.arguments["Om0"] + 0.1
        with pytest.raises(ValueError, match="baryonic density can not be larger"):
            cosmo_cls(*tba.args, **tba.kwargs)

        # no baryons specified
        tba = copy.deepcopy(ba)
        tba.arguments.pop("Ob0", None)
        cosmo = cosmo_cls(*tba.args, **tba.kwargs)
        with pytest.raises(ValueError):
            cosmo.Ob(1)

        # also means DM fraction is un
        with pytest.raises(ValueError):
            cosmo.Odm(1)


class TestFLRW(CosmologyTest,
               ParameterH0TestMixin, ParameterOm0TestMixin, ParameterOde0TestMixin,
               ParameterTcmb0TestMixin, ParameterNeffTestMixin, Parameterm_nuTestMixin,
               ParameterOb0TestMixin):
    """Test :class:`astropy.cosmology.FLRW`."""

    def setup_class(self):
        """
        Setup for testing.
        FLRW is abstract, so tests are done on a subclass.
        """

        class SubFLRW(FLRW):
            def w(self, z):
                return super().w(z)

        self.cls = SubFLRW
        # H0, Om0, Ode0
        self._cls_args = dict(H0=70 * u.km / u.s / u.Mpc, Om0=0.27 * u.one, Ode0=0.689 * u.one)
        self.cls_kwargs = dict(Tcmb0=3.0 * u.K, name=self.__class__.__name__, meta={"a": "b"})

    def teardown_class(self):
        super().teardown_class(self)
        _COSMOLOGY_CLASSES.pop("TestFLRW.setup_class.<locals>.SubFLRW", None)

    @pytest.fixture
    def nonflatcosmo(self):
        """A non-flat cosmology used in equivalence tests."""
        return LambdaCDM(70, 0.4, 0.8)

    # ===============================================================
    # Method & Attribute Tests

    def test_init(self, cosmo_cls):
        """Test initialization."""
        super().test_init(cosmo_cls)

        # TODO! tests for initializing calculated values, e.g. `h`
        # TODO! transfer tests for initializing neutrinos

    # ---------------------------------------------------------------
    # from Cosmology

    def test_clone_change_param(self, cosmo):
        """Test method ``.clone()`` changing a(many) Parameter(s)."""
        super().test_clone_change_param(cosmo)

        # don't change any values
        kwargs = cosmo._init_arguments.copy()
        kwargs.pop("name", None)  # make sure not setting name
        c = cosmo.clone(**kwargs)
        assert c.__class__ == cosmo.__class__
        assert c.name == cosmo.name + " (modified)"
        assert c.is_equivalent(cosmo)

        # change ``H0``
        # Note that H0 affects Ode0 because it changes Ogamma0
        c = cosmo.clone(H0=100)
        assert c.__class__ == cosmo.__class__
        assert c.name == cosmo.name + " (modified)"
        assert c.H0.value == 100
        for n in ("Om0", "Ode0", "Tcmb0", "Neff", "m_nu", "Ok0", "Ob0"):
            v = getattr(c, n)
            if v is None:
                assert v is getattr(cosmo, n)
                continue
            assert u.allclose(v, getattr(cosmo, n), atol=1e-4 * getattr(v, "unit", 1))
        assert not u.allclose(c.Ogamma0, cosmo.Ogamma0)
        assert not u.allclose(c.Onu0, cosmo.Onu0)

        # change multiple things
        c = cosmo.clone(name="new name", H0=100, Tcmb0=2.8, meta=dict(zz="tops"))
        assert c.__class__ == cosmo.__class__
        assert c.name == "new name"
        assert c.H0.value == 100
        assert c.Tcmb0.value == 2.8
        assert c.meta == {**cosmo.meta, **dict(zz="tops")}
        for n in ("Om0", "Ode0", "Neff", "m_nu", "Ok0", "Ob0"):
            v = getattr(c, n)
            if v is None:
                assert v is getattr(cosmo, n)
                continue
            assert u.allclose(v, getattr(cosmo, n), atol=1e-4 * getattr(v, "unit", 1))
        assert not u.allclose(c.Ogamma0, cosmo.Ogamma0)
        assert not u.allclose(c.Onu0, cosmo.Onu0)
        assert not u.allclose(c.Tcmb0.value, cosmo.Tcmb0.value)

    def test_is_equivalent(self, cosmo):
        """Test :meth:`astropy.cosmology.FLRW.is_equivalent`."""
        super().test_is_equivalent(cosmo)  # pass to CosmologySubclassTest

        # test against a FlatFLRWMixin
        # case (3) in FLRW.is_equivalent
        if isinstance(cosmo, FlatLambdaCDM):
            assert cosmo.is_equivalent(Planck18)
            assert Planck18.is_equivalent(cosmo)
        else:
            assert not cosmo.is_equivalent(Planck18)
            assert not Planck18.is_equivalent(cosmo)

    # ---------------------------------------------------------------

    def test_efunc_vs_invefunc(self, cosmo):
        """
        Test that efunc and inv_efunc give inverse values.
        Here they just fail b/c no ``w(z)`` or no scipy.
        """
        exception = NotImplementedError if HAS_SCIPY else ModuleNotFoundError

        with pytest.raises(exception):
            cosmo.efunc(0.5)

        with pytest.raises(exception):
            cosmo.inv_efunc(0.5)


class FLRWSubclassTest(TestFLRW):
    """
    Test subclasses of :class:`astropy.cosmology.FLRW`.
    This is broken away from ``TestFLRW``, because ``FLRW`` is an ABC and
    subclasses must override some methods.
    """

    @abc.abstractmethod
    def setup_class(self):
        """Setup for testing."""
        pass

    # ===============================================================

    def test_efunc_vs_invefunc(self, cosmo):
        """Test that ``efunc`` and ``inv_efunc`` give inverse values.

        Note that the test doesn't need scipy because it doesn't need to call
        ``de_density_scale``.
        """
        # super().test_efunc_vs_invefunc(cosmo)  # NOT b/c abstract `w(z)`
        z0 = 0.5
        z = np.array([0.5, 1.0, 2.0, 5.0])

        assert np.allclose(cosmo.efunc(z0), 1.0 / cosmo.inv_efunc(z0))
        assert np.allclose(cosmo.efunc(z), 1.0 / cosmo.inv_efunc(z))


# -----------------------------------------------------------------------------

class ParameterFlatOde0TestMixin(ParameterOde0TestMixin):

    def test_Ode0(self, cosmo_cls, cosmo):
        """Test no-longer-Parameter ``Ode0``."""
        # on the class
        assert isinstance(cosmo_cls.Ode0, Parameter)
        assert cosmo_cls.Ode0.derived == True
        assert "Omega dark energy" in cosmo_cls.Ode0.__doc__

        # on the instance
        assert cosmo.Ode0 is cosmo._Ode0
        assert cosmo.Ode0 == 1.0 - (cosmo.Om0 + cosmo.Ogamma0 + cosmo.Onu0)

    def test_init_Ode0(self, cosmo_cls):
        """Test initialization for values of ``Ode0``."""
        ba = cosmo_cls._init_signature.bind(*self.cls_args, **self.cls_kwargs)

        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.Ode0 == 1.0 - (cosmo.Om0 + cosmo.Ogamma0 + cosmo.Onu0)

        # Ode0 is not in the signature
        with pytest.raises(TypeError, match="Ode0"):
            cosmo_cls(*ba.args, **ba.kwargs, Ode0=1)


class FlatFLRWMixinTest(FlatCosmologyMixinTest, ParameterFlatOde0TestMixin):
    """Tests for :class:`astropy.cosmology.FlatFLRWMixin`."""

    # ---------------------------------------------------------------
    # class-level

    def test_init_subclass(self, cosmo_cls):
        """Test initializing subclass, mostly that can't have Ode0 in init."""
        super().test_init_subclass(cosmo_cls)

        with pytest.raises(TypeError, match="subclasses of"):

            class HASOde0SubClass(cosmo_cls):
                def __init__(self, Ode0):
                    pass

            _COSMOLOGY_CLASSES.pop(HASOde0SubClass.__qualname__, None)

    # ---------------------------------------------------------------
    # instance-level

    def test_init(self, cosmo_cls):
        super().test_init(cosmo_cls)

        cosmo = cosmo_cls(*self.cls_args, **self.cls_kwargs)
        assert cosmo._Ode0 == 1.0 - cosmo._Om0 - cosmo._Ogamma0 - cosmo._Onu0
        assert cosmo._Ok0 == 0.0

    def test_is_equivalent(self, cosmo, nonflatcosmo):
        """Test :meth:`astropy.cosmology.FLRW.is_equivalent`."""
        super().test_is_equivalent(cosmo)  # pass to TestFLRW

        # against non-flat Cosmology
        assert not cosmo.is_equivalent(nonflatcosmo)
        assert not nonflatcosmo.is_equivalent(cosmo)

        # non-flat version of class
        nonflat_cosmo_cls = cosmo.__class__.mro()[3]
        # keys check in `test_is_equivalent_nonflat_class_different_params`

        # non-flat
        nonflat = nonflat_cosmo_cls(*self.cls_args, Ode0=0.9, **self.cls_kwargs)
        assert not nonflat.is_equivalent(cosmo)
        assert not cosmo.is_equivalent(nonflat)

        # flat, but not FlatFLRWMixin
        flat = nonflat_cosmo_cls(*self.cls_args,
                                   Ode0=1.0 - cosmo.Om0 - cosmo.Ogamma0 - cosmo.Onu0,
                                   **self.cls_kwargs)
        flat._Ok0 = 0.0
        assert flat.is_equivalent(cosmo)
        assert cosmo.is_equivalent(flat)

    def test_repr(self, cosmo_cls, cosmo):
        """
        Test method ``.__repr__()``. Skip non-flat superclass test.
        e.g. `TestFlatLambdaCDDM` -> `FlatFLRWMixinTest`
        vs   `TestFlatLambdaCDDM` -> `TestLambdaCDDM` -> `FlatFLRWMixinTest`
        """
        FLRWSubclassTest.test_repr(self, cosmo_cls, cosmo)

        # test eliminated Ode0 from parameters
        assert "Ode0" not in repr(cosmo)


# -----------------------------------------------------------------------------


class TestLambdaCDM(FLRWSubclassTest):
    """Test :class:`astropy.cosmology.LambdaCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = LambdaCDM
        self._cls_args = dict(H0=70 * (u.km / u.s / u.Mpc), Om0=0.27, Ode0=0.73)
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name=self.__class__.__name__, meta={"a": "b"})

    def test_repr(self, cosmo_cls, cosmo):
        """Test method ``.__repr__()``."""
        super().test_repr(cosmo_cls, cosmo)

        expected = ("LambdaCDM(name=\"ABCMeta\", H0=70 km / (Mpc s), Om0=0.27,"
                    " Ode0=0.73, Tcmb0=3 K, Neff=3.04, m_nu=[0. 0. 0.] eV,"
                    " Ob0=None)")
        assert repr(cosmo) == expected


# -----------------------------------------------------------------------------


class TestFlatLambdaCDM(FlatFLRWMixinTest, TestLambdaCDM):
    """Test :class:`astropy.cosmology.FlatLambdaCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = FlatLambdaCDM
        self._cls_args = dict(H0=70 * (u.km / u.s / u.Mpc), Om0=0.27)
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name=self.__class__.__name__, meta={"a": "b"})

    def test_repr(self, cosmo_cls, cosmo):
        """Test method ``.__repr__()``."""
        super().test_repr(cosmo_cls, cosmo)

        expected = ("FlatLambdaCDM(name=\"ABCMeta\", H0=70 km / (Mpc s),"
                    " Om0=0.27, Tcmb0=3 K, Neff=3.04, m_nu=[0. 0. 0.] eV,"
                    " Ob0=None)")
        assert repr(cosmo) == expected


# -----------------------------------------------------------------------------


class Parameterw0TestMixin(ParameterTestMixin):

    def test_w0(self, cosmo_cls, cosmo):
        """Test Parameter ``w0``."""
        # on the class
        assert isinstance(cosmo_cls.w0, Parameter)
        assert "Dark energy equation of state" in cosmo_cls.w0.__doc__
        assert cosmo_cls.w0.unit is None

        # on the instance
        assert cosmo.w0 is cosmo._w0
        assert cosmo.w0 == -1.0

    def test_init_w0(self, cosmo_cls):
        """Test initialization for values of ``w0``."""
        ba = cosmo_cls._init_signature.bind(*self.cls_args, **self.cls_kwargs)

        ba.arguments["w0"] = 10 * u.km
        with pytest.raises(TypeError):
            cosmo_cls(*ba.args, **ba.kwargs)


class TestwCDM(FLRWSubclassTest, Parameterw0TestMixin):
    """Test :class:`astropy.cosmology.wCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = wCDM
        self._cls_args = dict(H0=70 * (u.km / u.s / u.Mpc), Om0=0.27, Ode0=0.73)
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name=self.__class__.__name__, meta={"a": "b"})

    def test_clone_change_param(self, cosmo):
        """Test method ``.clone()`` changing a(many) Parameter(s)."""
        super().test_clone_change_param(cosmo)

        # `w` params
        c = cosmo.clone(w0=0.1)
        assert c.w0 == 0.1
        for n in ("H0", "Om0", "Ode0", "Tcmb0", "Neff", "m_nu", "Ok0", "Ob0"):
            v = getattr(c, n)
            if v is None:
                assert v is getattr(cosmo, n)
                continue
            assert u.allclose(v, getattr(cosmo, n), atol=1e-4 * getattr(v, "unit", 1))

    def test_repr(self, cosmo_cls, cosmo):
        """Test method ``.__repr__()``."""
        super().test_repr(cosmo_cls, cosmo)

        expected = ("wCDM(name=\"ABCMeta\", H0=70 km / (Mpc s), Om0=0.27,"
                    " Ode0=0.73, w0=-1, Tcmb0=3 K, Neff=3.04,"
                    " m_nu=[0. 0. 0.] eV, Ob0=None)")
        assert repr(cosmo) == expected


# -----------------------------------------------------------------------------


class TestFlatwCDM(FlatFLRWMixinTest, TestwCDM):
    """Test :class:`astropy.cosmology.FlatwCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = FlatwCDM
        self._cls_args = dict(H0=70 * (u.km / u.s / u.Mpc), Om0=0.27)
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name=self.__class__.__name__, meta={"a": "b"})

    def test_repr(self, cosmo_cls, cosmo):
        """Test method ``.__repr__()``."""
        super().test_repr(cosmo_cls, cosmo)

        expected = ("FlatwCDM(name=\"ABCMeta\", H0=70 km / (Mpc s), Om0=0.27,"
                    " w0=-1, Tcmb0=3 K, Neff=3.04, m_nu=[0. 0. 0.] eV,"
                    " Ob0=None)")
        assert repr(cosmo) == expected


# -----------------------------------------------------------------------------


class ParameterwaTestMixin(ParameterTestMixin):

    def test_wa(self, cosmo_cls, cosmo):
        """Test Parameter ``wa``."""
        # on the class
        assert isinstance(cosmo_cls.wa, Parameter)
        assert "Negative derivative" in cosmo_cls.wa.__doc__
        assert cosmo_cls.wa.unit is None

        # on the instance
        assert cosmo.wa is cosmo._wa
        assert cosmo.wa == 0.0

    def test_init_wa(self, cosmo_cls):
        """Test initialization for values of ``wa``."""
        ba = cosmo_cls._init_signature.bind(*self.cls_args, **self.cls_kwargs)

        ba.arguments["wa"] = 10 * u.km
        with pytest.raises(TypeError):
            cosmo_cls(*ba.args, **ba.kwargs)


class Testw0waCDM(FLRWSubclassTest, Parameterw0TestMixin, ParameterwaTestMixin):
    """Test :class:`astropy.cosmology.w0waCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = w0waCDM
        self._cls_args = dict(H0=70 * (u.km / u.s / u.Mpc), Om0=0.27, Ode0=0.73)
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name=self.__class__.__name__, meta={"a": "b"})

    def test_clone_change_param(self, cosmo):
        """Test method ``.clone()`` changing a(many) Parameter(s)."""
        super().test_clone_change_param(cosmo)

        # `w` params
        c = cosmo.clone(w0=0.1, wa=0.2)
        assert c.w0 == 0.1
        assert c.wa == 0.2
        for n in ("H0", "Om0", "Ode0", "Tcmb0", "Neff", "m_nu", "Ok0", "Ob0"):
            v = getattr(c, n)
            if v is None:
                assert v is getattr(cosmo, n)
                continue
            assert u.allclose(v, getattr(cosmo, n), atol=1e-4 * getattr(v, "unit", 1))

    def test_repr(self, cosmo_cls, cosmo):
        """Test method ``.__repr__()``."""
        super().test_repr(cosmo_cls, cosmo)

        expected = ("w0waCDM(name=\"ABCMeta\", H0=70 km / (Mpc s), Om0=0.27,"
                    " Ode0=0.73, w0=-1, wa=0, Tcmb0=3 K, Neff=3.04,"
                    " m_nu=[0. 0. 0.] eV, Ob0=None)")
        assert repr(cosmo) == expected


# -----------------------------------------------------------------------------


class TestFlatw0waCDM(FlatFLRWMixinTest, Testw0waCDM):
    """Test :class:`astropy.cosmology.Flatw0waCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = Flatw0waCDM
        self._cls_args = dict(H0=70 * (u.km / u.s / u.Mpc), Om0=0.27)
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name=self.__class__.__name__, meta={"a": "b"})

    def test_repr(self, cosmo_cls, cosmo):
        """Test method ``.__repr__()``."""
        super().test_repr(cosmo_cls, cosmo)

        expected = ("Flatw0waCDM(name=\"ABCMeta\", H0=70 km / (Mpc s),"
                    " Om0=0.27, w0=-1, wa=0, Tcmb0=3 K, Neff=3.04,"
                    " m_nu=[0. 0. 0.] eV, Ob0=None)")
        assert repr(cosmo) == expected


# -----------------------------------------------------------------------------


class ParameterwpTestMixin(ParameterTestMixin):

    def test_wp(self, cosmo_cls, cosmo):
        """Test Parameter ``wp``."""
        # on the class
        assert isinstance(cosmo_cls.wp, Parameter)
        assert "at the pivot" in cosmo_cls.wp.__doc__
        assert cosmo_cls.wp.unit is None

        # on the instance
        assert cosmo.wp is cosmo._wp
        assert cosmo.wp == -1.0

    def test_init_wp(self, cosmo_cls):
        """Test initialization for values of ``wp``."""
        ba = cosmo_cls._init_signature.bind(*self.cls_args, **self.cls_kwargs)

        ba.arguments["wp"] = 10 * u.km
        with pytest.raises(TypeError):
            cosmo_cls(*ba.args, **ba.kwargs)


class ParameterzpTestMixin(ParameterTestMixin):

    def test_zp(self, cosmo_cls, cosmo):
        """Test Parameter ``zp``."""
        # on the class
        assert isinstance(cosmo_cls.zp, Parameter)
        assert "pivot redshift" in cosmo_cls.zp.__doc__
        assert cosmo_cls.zp.unit == cu.redshift

        # on the instance
        assert cosmo.zp is cosmo._zp
        assert cosmo.zp == 0 * cu.redshift

    def test_init_zp(self, cosmo_cls):
        """Test initialization for values of ``zp``."""
        ba = cosmo_cls._init_signature.bind(*self.cls_args, **self.cls_kwargs)

        ba.arguments["zp"] = 10 * u.km
        with pytest.raises(u.UnitConversionError):
            cosmo_cls(*ba.args, **ba.kwargs)


class TestwpwaCDM(FLRWSubclassTest,
                  ParameterwpTestMixin, ParameterwaTestMixin, ParameterzpTestMixin):
    """Test :class:`astropy.cosmology.wpwaCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = wpwaCDM
        self._cls_args = dict(H0=70 * (u.km / u.s / u.Mpc), Om0=0.27, Ode0=0.73)
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name=self.__class__.__name__, meta={"a": "b"})

    def test_clone_change_param(self, cosmo):
        """Test method ``.clone()`` changing a(many) Parameter(s)."""
        super().test_clone_change_param(cosmo)

        # `w` params
        c = cosmo.clone(wp=0.1, wa=0.2, zp=14)
        assert c.wp == 0.1
        assert c.wa == 0.2
        assert c.zp == 14
        for n in ("H0", "Om0", "Ode0", "Tcmb0", "Neff", "m_nu", "Ok0", "Ob0"):
            v = getattr(c, n)
            if v is None:
                assert v is getattr(cosmo, n)
                continue
            assert u.allclose(v, getattr(cosmo, n), atol=1e-4 * getattr(v, "unit", 1))

    def test_repr(self, cosmo_cls, cosmo):
        """Test method ``.__repr__()``."""
        super().test_repr(cosmo_cls, cosmo)

        expected = ("wpwaCDM(name=\"ABCMeta\", H0=70 km / (Mpc s), Om0=0.27,"
                    " Ode0=0.73, wp=-1, wa=0, zp=0 redshift, Tcmb0=3 K,"
                    " Neff=3.04, m_nu=[0. 0. 0.] eV, Ob0=None)")
        assert repr(cosmo) == expected


# -----------------------------------------------------------------------------


class ParameterwzTestMixin(ParameterTestMixin):

    def test_wz(self, cosmo_cls, cosmo):
        """Test Parameter ``wz``."""
        # on the class
        assert isinstance(cosmo_cls.wz, Parameter)
        assert "Derivative of the dark energy" in cosmo_cls.wz.__doc__
        assert cosmo_cls.wz.unit is None

        # on the instance
        assert cosmo.wz is cosmo._wz
        assert cosmo.wz == 0.0

    def test_init_wz(self, cosmo_cls):
        """Test initialization for values of ``wz``."""
        ba = cosmo_cls._init_signature.bind(*self.cls_args, **self.cls_kwargs)

        ba.arguments["wz"] = 10 * u.km
        with pytest.raises(TypeError):
            cosmo_cls(*ba.args, **ba.kwargs)


class Testw0wzCDM(FLRWSubclassTest, Parameterw0TestMixin, ParameterwzTestMixin):
    """Test :class:`astropy.cosmology.w0wzCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = w0wzCDM
        self._cls_args = dict(H0=70 * (u.km / u.s / u.Mpc), Om0=0.27, Ode0=0.73)
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name=self.__class__.__name__, meta={"a": "b"})

    def test_clone_change_param(self, cosmo):
        """Test method ``.clone()`` changing a(many) Parameter(s)."""
        super().test_clone_change_param(cosmo)

        # `w` params
        c = cosmo.clone(w0=0.1, wz=0.2)
        assert c.w0 == 0.1
        assert c.wz == 0.2
        for n in ("H0", "Om0", "Ode0", "Tcmb0", "Neff", "m_nu", "Ok0", "Ob0"):
            v = getattr(c, n)
            if v is None:
                assert v is getattr(cosmo, n)
                continue
            assert u.allclose(v, getattr(cosmo, n), atol=1e-4 * getattr(v, "unit", 1))

    def test_repr(self, cosmo_cls, cosmo):
        """Test method ``.__repr__()``."""
        super().test_repr(cosmo_cls, cosmo)

        expected = ("w0wzCDM(name=\"ABCMeta\", H0=70 km / (Mpc s), Om0=0.27,"
                    " Ode0=0.73, w0=-1, wz=0, Tcmb0=3 K, Neff=3.04,"
                    " m_nu=[0. 0. 0.] eV, Ob0=None)")
        assert repr(cosmo) == expected
