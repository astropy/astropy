# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.flrw.base`."""

##############################################################################
# IMPORTS

# STDLIB
import abc
import copy

# THIRD PARTY
import numpy as np
import pytest

import astropy.constants as const

# LOCAL
import astropy.units as u
from astropy.cosmology import FLRW, FlatLambdaCDM, LambdaCDM, Planck18
from astropy.cosmology.core import _COSMOLOGY_CLASSES
from astropy.cosmology.flrw.base import _a_B_c2, _critdens_const, _H0units_to_invs, quad
from astropy.cosmology.parameter import Parameter
from astropy.cosmology.tests.helper import get_redshift_methods
from astropy.cosmology.tests.test_core import (
    CosmologyTest,
    FlatCosmologyMixinTest,
    ParameterTestMixin,
    invalid_zs,
    valid_zs,
)
from astropy.utils.compat.optional_deps import HAS_SCIPY

##############################################################################
# SETUP / TEARDOWN


class SubFLRW(FLRW):
    def w(self, z):
        return super().w(z)


##############################################################################
# TESTS
##############################################################################


@pytest.mark.skipif(HAS_SCIPY, reason="scipy is installed")
def test_optional_deps_functions():
    """Test stand-in functions when optional dependencies not installed."""
    with pytest.raises(ModuleNotFoundError, match="No module named 'scipy.integrate'"):
        quad()


##############################################################################


class ParameterH0TestMixin(ParameterTestMixin):
    """Tests for `astropy.cosmology.Parameter` H0 on a Cosmology.

    H0 is a descriptor, which are tested by mixin, here with ``TestFLRW``.
    These tests expect dicts ``_cls_args`` and ``cls_kwargs`` which give the
    args and kwargs for the cosmology class, respectively. See ``TestFLRW``.
    """

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

    def test_init_H0(self, cosmo_cls, ba):
        """Test initialization for values of ``H0``."""
        # test that it works with units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.H0 == ba.arguments["H0"]

        # also without units
        ba.arguments["H0"] = ba.arguments["H0"].value  # strip units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.H0.value == ba.arguments["H0"]

        # fails for non-scalar
        ba.arguments["H0"] = u.Quantity([70, 100], u.km / u.s / u.Mpc)
        with pytest.raises(ValueError, match="H0 is a non-scalar quantity"):
            cosmo_cls(*ba.args, **ba.kwargs)


class ParameterOm0TestMixin(ParameterTestMixin):
    """Tests for `astropy.cosmology.Parameter` Om0 on a Cosmology.

    Om0 is a descriptor, which are tested by mixin, here with ``TestFLRW``.
    These tests expect dicts ``_cls_args`` and ``cls_kwargs`` which give the
    args and kwargs for the cosmology class, respectively. See ``TestFLRW``.
    """

    def test_Om0(self, cosmo_cls, cosmo):
        """Test Parameter ``Om0``."""
        # on the class
        assert isinstance(cosmo_cls.Om0, Parameter)
        assert "Omega matter" in cosmo_cls.Om0.__doc__

        # validation
        assert cosmo_cls.Om0.validate(cosmo, 1) == 1
        assert cosmo_cls.Om0.validate(cosmo, 10 * u.one) == 10
        with pytest.raises(ValueError, match="Om0 cannot be negative"):
            cosmo_cls.Om0.validate(cosmo, -1)

        # on the instance
        assert cosmo.Om0 is cosmo._Om0
        assert cosmo.Om0 == self._cls_args["Om0"]
        assert isinstance(cosmo.Om0, float)

    def test_init_Om0(self, cosmo_cls, ba):
        """Test initialization for values of ``Om0``."""
        # test that it works with units
        ba.arguments["Om0"] = ba.arguments["Om0"] << u.one  # ensure units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.Om0 == ba.arguments["Om0"]

        # also without units
        ba.arguments["Om0"] = ba.arguments["Om0"].value  # strip units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.Om0 == ba.arguments["Om0"]

        # fails for negative numbers
        ba.arguments["Om0"] = -0.27
        with pytest.raises(ValueError, match="Om0 cannot be negative."):
            cosmo_cls(*ba.args, **ba.kwargs)


class ParameterOde0TestMixin(ParameterTestMixin):
    """Tests for `astropy.cosmology.Parameter` Ode0 on a Cosmology.

    Ode0 is a descriptor, which are tested by mixin, here with ``TestFLRW``.
    These tests expect dicts ``_cls_args`` and ``cls_kwargs`` which give the
    args and kwargs for the cosmology class, respectively. See ``TestFLRW``.
    """

    def test_Parameter_Ode0(self, cosmo_cls):
        """Test Parameter ``Ode0`` on the class."""
        assert isinstance(cosmo_cls.Ode0, Parameter)
        assert "Omega dark energy" in cosmo_cls.Ode0.__doc__

    def test_Parameter_Ode0_validation(self, cosmo_cls, cosmo):
        """Test Parameter ``Ode0`` validation."""
        assert cosmo_cls.Ode0.validate(cosmo, 1.1) == 1.1
        assert cosmo_cls.Ode0.validate(cosmo, 10 * u.one) == 10.0
        with pytest.raises(TypeError, match="only dimensionless"):
            cosmo_cls.Ode0.validate(cosmo, 10 * u.km)

    def test_Ode0(self, cosmo):
        """Test Parameter ``Ode0`` validation."""
        # if Ode0 is a parameter, test its value
        assert cosmo.Ode0 is cosmo._Ode0
        assert cosmo.Ode0 == self._cls_args["Ode0"]
        assert isinstance(cosmo.Ode0, float)

    def test_init_Ode0(self, cosmo_cls, ba):
        """Test initialization for values of ``Ode0``."""
        # test that it works with units
        ba.arguments["Ode0"] = ba.arguments["Ode0"] << u.one  # ensure units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.Ode0 == ba.arguments["Ode0"]

        # also without units
        ba.arguments["Ode0"] = ba.arguments["Ode0"].value  # strip units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.Ode0 == ba.arguments["Ode0"]

        # Setting param to 0 respects that. Note this test uses ``Ode()``.
        ba.arguments["Ode0"] = 0.0
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert u.allclose(cosmo.Ode([0, 1, 2, 3]), [0, 0, 0, 0])
        assert u.allclose(cosmo.Ode(1), 0)

        # Must be dimensionless or have no units. Errors otherwise.
        ba.arguments["Ode0"] = 10 * u.km
        with pytest.raises(TypeError, match="only dimensionless"):
            cosmo_cls(*ba.args, **ba.kwargs)


class ParameterTcmb0TestMixin(ParameterTestMixin):
    """Tests for `astropy.cosmology.Parameter` Tcmb0 on a Cosmology.

    Tcmb0 is a descriptor, which are tested by mixin, here with ``TestFLRW``.
    These tests expect dicts ``_cls_args`` and ``cls_kwargs`` which give the
    args and kwargs for the cosmology class, respectively. See ``TestFLRW``.
    """

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

    def test_init_Tcmb0(self, cosmo_cls, ba):
        """Test initialization for values of ``Tcmb0``."""
        # test that it works with units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.Tcmb0 == ba.arguments["Tcmb0"]

        # also without units
        ba.arguments["Tcmb0"] = ba.arguments["Tcmb0"].value  # strip units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.Tcmb0.value == ba.arguments["Tcmb0"]

        # must be a scalar
        ba.arguments["Tcmb0"] = u.Quantity([0.0, 2], u.K)
        with pytest.raises(ValueError, match="Tcmb0 is a non-scalar quantity"):
            cosmo_cls(*ba.args, **ba.kwargs)


class ParameterNeffTestMixin(ParameterTestMixin):
    """Tests for `astropy.cosmology.Parameter` Neff on a Cosmology.

    Neff is a descriptor, which are tested by mixin, here with ``TestFLRW``.
    These tests expect dicts ``_cls_args`` and ``cls_kwargs`` which give the
    args and kwargs for the cosmology class, respectively. See ``TestFLRW``.
    """

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

    def test_init_Neff(self, cosmo_cls, ba):
        """Test initialization for values of ``Neff``."""
        # test that it works with units
        ba.arguments["Neff"] = ba.arguments["Neff"] << u.one  # ensure units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.Neff == ba.arguments["Neff"]

        # also without units
        ba.arguments["Neff"] = ba.arguments["Neff"].value  # strip units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.Neff == ba.arguments["Neff"]

        ba.arguments["Neff"] = -1
        with pytest.raises(ValueError):
            cosmo_cls(*ba.args, **ba.kwargs)


class Parameterm_nuTestMixin(ParameterTestMixin):
    """Tests for `astropy.cosmology.Parameter` m_nu on a Cosmology.

    m_nu is a descriptor, which are tested by mixin, here with ``TestFLRW``.
    These tests expect dicts ``_cls_args`` and ``cls_kwargs`` which give the
    args and kwargs for the cosmology class, respectively. See ``TestFLRW``.
    """

    def test_m_nu(self, cosmo_cls, cosmo):
        """Test Parameter ``m_nu``."""
        # on the class
        assert isinstance(cosmo_cls.m_nu, Parameter)
        assert "Mass of neutrino species" in cosmo_cls.m_nu.__doc__
        assert cosmo_cls.m_nu.unit == u.eV
        assert cosmo_cls.m_nu.equivalencies == u.mass_energy()

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
            assert u.allclose(cosmo.m_nu[: self._nmasslessnu], 0 * u.eV)
            assert u.allclose(cosmo.m_nu[self._nmasslessnu], cosmo._massivenu_mass)

    def test_init_m_nu(self, cosmo_cls, ba):
        """Test initialization for values of ``m_nu``.

        Note this requires the class to have a property ``has_massive_nu``.
        """
        # Test that it works when m_nu has units.
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert np.all(cosmo.m_nu == ba.arguments["m_nu"])  # (& checks len, unit)
        assert not cosmo.has_massive_nu
        assert cosmo.m_nu.unit == u.eV  # explicitly check unit once.

        # And it works when m_nu doesn't have units.
        ba.arguments["m_nu"] = ba.arguments["m_nu"].value  # strip units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert np.all(cosmo.m_nu.value == ba.arguments["m_nu"])
        assert not cosmo.has_massive_nu

        # A negative m_nu raises an exception.
        tba = copy.copy(ba)
        tba.arguments["m_nu"] = u.Quantity([-0.3, 0.2, 0.1], u.eV)
        with pytest.raises(ValueError, match="invalid"):
            cosmo_cls(*tba.args, **tba.kwargs)

    def test_init_m_nu_and_Neff(self, cosmo_cls, ba):
        """Test initialization for values of ``m_nu`` and ``Neff``.

        Note this test requires ``Neff`` as constructor input, and a property
        ``has_massive_nu``.
        """
        # Mismatch with Neff = wrong number of neutrinos
        tba = copy.copy(ba)
        tba.arguments["Neff"] = 4.05
        tba.arguments["m_nu"] = u.Quantity([0.15, 0.2, 0.1], u.eV)
        with pytest.raises(ValueError, match="unexpected number of neutrino"):
            cosmo_cls(*tba.args, **tba.kwargs)

        # No neutrinos, but Neff
        tba.arguments["m_nu"] = 0
        cosmo = cosmo_cls(*tba.args, **tba.kwargs)
        assert not cosmo.has_massive_nu
        assert len(cosmo.m_nu) == 4
        assert cosmo.m_nu.unit == u.eV
        assert u.allclose(cosmo.m_nu, 0 * u.eV)
        # TODO! move this test when create ``test_nu_relative_density``
        assert u.allclose(
            cosmo.nu_relative_density(1.0), 0.22710731766 * 4.05, rtol=1e-6
        )

        # All massive neutrinos case, len from Neff
        tba.arguments["m_nu"] = 0.1 * u.eV
        cosmo = cosmo_cls(*tba.args, **tba.kwargs)
        assert cosmo.has_massive_nu
        assert len(cosmo.m_nu) == 4
        assert cosmo.m_nu.unit == u.eV
        assert u.allclose(cosmo.m_nu, [0.1, 0.1, 0.1, 0.1] * u.eV)

    def test_init_m_nu_override_by_Tcmb0(self, cosmo_cls, ba):
        """Test initialization for values of ``m_nu``.

        Note this test requires ``Tcmb0`` as constructor input, and a property
        ``has_massive_nu``.
        """
        # If Neff = 0, m_nu is None.
        tba = copy.copy(ba)
        tba.arguments["Neff"] = 0
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.m_nu is None
        assert not cosmo.has_massive_nu

        # If Tcmb0 = 0, m_nu is None
        tba = copy.copy(ba)
        tba.arguments["Tcmb0"] = 0
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.m_nu is None
        assert not cosmo.has_massive_nu


class ParameterOb0TestMixin(ParameterTestMixin):
    """Tests for `astropy.cosmology.Parameter` Ob0 on a Cosmology.

    Ob0 is a descriptor, which are tested by mixin, here with ``TestFLRW``.
    These tests expect dicts ``_cls_args`` and ``cls_kwargs`` which give the
    args and kwargs for the cosmology class, respectively. See ``TestFLRW``.
    """

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
        assert cosmo.Ob0 == 0.03

    def test_init_Ob0(self, cosmo_cls, ba):
        """Test initialization for values of ``Ob0``."""
        # test that it works with units
        assert isinstance(ba.arguments["Ob0"], u.Quantity)
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.Ob0 == ba.arguments["Ob0"]

        # also without units
        ba.arguments["Ob0"] = ba.arguments["Ob0"].value  # strip units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.Ob0 == ba.arguments["Ob0"]

        # Setting param to 0 respects that.  Note this test uses ``Ob()``.
        ba.arguments["Ob0"] = 0.0
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.Ob0 == 0.0
        if not self.abstract_w:
            assert u.allclose(cosmo.Ob(1), 0)
            assert u.allclose(cosmo.Ob([0, 1, 2, 3]), [0, 0, 0, 0])

        # Negative Ob0 errors
        tba = copy.copy(ba)
        tba.arguments["Ob0"] = -0.04
        with pytest.raises(ValueError, match="Ob0 cannot be negative"):
            cosmo_cls(*tba.args, **tba.kwargs)

        # Ob0 > Om0 errors
        tba.arguments["Ob0"] = tba.arguments["Om0"] + 0.1
        with pytest.raises(ValueError, match="baryonic density can not be larger"):
            cosmo_cls(*tba.args, **tba.kwargs)

        # No baryons specified means baryon-specific methods fail.
        tba = copy.copy(ba)
        tba.arguments.pop("Ob0", None)
        cosmo = cosmo_cls(*tba.args, **tba.kwargs)
        with pytest.raises(ValueError):
            cosmo.Ob(1)

        # also means DM fraction is undefined
        with pytest.raises(ValueError):
            cosmo.Odm(1)

        # The default value is None
        assert cosmo_cls._init_signature.parameters["Ob0"].default is None


class FLRWTest(
    CosmologyTest,
    ParameterH0TestMixin,
    ParameterOm0TestMixin,
    ParameterOde0TestMixin,
    ParameterTcmb0TestMixin,
    ParameterNeffTestMixin,
    Parameterm_nuTestMixin,
    ParameterOb0TestMixin,
):
    abstract_w = False

    @abc.abstractmethod
    def setup_class(self):
        """Setup for testing."""
        super().setup_class(self)

        # Default cosmology args and kwargs
        self._cls_args = dict(
            H0=70 * u.km / u.s / u.Mpc, Om0=0.27 * u.one, Ode0=0.73 * u.one
        )
        self.cls_kwargs = dict(
            Tcmb0=3.0 * u.K,
            Ob0=0.03 * u.one,
            name=self.__class__.__name__,
            meta={"a": "b"},
        )

    @pytest.fixture(scope="class")
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

    def test_init_Tcmb0_zeroing(self, cosmo_cls, ba):
        """Test if setting Tcmb0 parameter to 0 influences other parameters.

        TODO: consider moving this test to ``FLRWTest``
        """
        ba.arguments["Tcmb0"] = 0.0
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)

        assert cosmo.Ogamma0 == 0.0
        assert cosmo.Onu0 == 0.0

        if not self.abstract_w:
            assert u.allclose(cosmo.Ogamma(1.5), [0, 0, 0, 0])
            assert u.allclose(cosmo.Ogamma([0, 1, 2, 3]), [0, 0, 0, 0])
            assert u.allclose(cosmo.Onu(1.5), [0, 0, 0, 0])
            assert u.allclose(cosmo.Onu([0, 1, 2, 3]), [0, 0, 0, 0])

    # ---------------------------------------------------------------
    # Properties

    def test_Odm0(self, cosmo_cls, cosmo):
        """Test property ``Odm0``."""
        # on the class
        assert isinstance(cosmo_cls.Odm0, property)
        assert cosmo_cls.Odm0.fset is None  # immutable

        # on the instance
        assert cosmo.Odm0 is cosmo._Odm0
        # Odm0 can be None, if Ob0 is None. Otherwise DM = matter - baryons.
        if cosmo.Ob0 is None:
            assert cosmo.Odm0 is None
        else:
            assert np.allclose(cosmo.Odm0, cosmo.Om0 - cosmo.Ob0)

    def test_Ok0(self, cosmo_cls, cosmo):
        """Test property ``Ok0``."""
        # on the class
        assert isinstance(cosmo_cls.Ok0, property)
        assert cosmo_cls.Ok0.fset is None  # immutable

        # on the instance
        assert cosmo.Ok0 is cosmo._Ok0
        assert np.allclose(
            cosmo.Ok0, 1.0 - (cosmo.Om0 + cosmo.Ode0 + cosmo.Ogamma0 + cosmo.Onu0)
        )

    def test_is_flat(self, cosmo_cls, cosmo):
        """Test property ``is_flat``."""
        # on the class
        assert isinstance(cosmo_cls.is_flat, property)
        assert cosmo_cls.is_flat.fset is None  # immutable

        # on the instance
        assert isinstance(cosmo.is_flat, bool)
        assert cosmo.is_flat is bool((cosmo.Ok0 == 0.0) and (cosmo.Otot0 == 1.0))

    def test_Tnu0(self, cosmo_cls, cosmo):
        """Test property ``Tnu0``."""
        # on the class
        assert isinstance(cosmo_cls.Tnu0, property)
        assert cosmo_cls.Tnu0.fset is None  # immutable

        # on the instance
        assert cosmo.Tnu0 is cosmo._Tnu0
        assert cosmo.Tnu0.unit == u.K
        assert u.allclose(cosmo.Tnu0, 0.7137658555036082 * cosmo.Tcmb0, rtol=1e-5)

    def test_has_massive_nu(self, cosmo_cls, cosmo):
        """Test property ``has_massive_nu``."""
        # on the class
        assert isinstance(cosmo_cls.has_massive_nu, property)
        assert cosmo_cls.has_massive_nu.fset is None  # immutable

        # on the instance
        if cosmo.Tnu0 == 0:
            assert cosmo.has_massive_nu is False
        else:
            assert cosmo.has_massive_nu is cosmo._massivenu

    def test_h(self, cosmo_cls, cosmo):
        """Test property ``h``."""
        # on the class
        assert isinstance(cosmo_cls.h, property)
        assert cosmo_cls.h.fset is None  # immutable

        # on the instance
        assert cosmo.h is cosmo._h
        assert np.allclose(cosmo.h, cosmo.H0.value / 100.0)

    def test_hubble_time(self, cosmo_cls, cosmo):
        """Test property ``hubble_time``."""
        # on the class
        assert isinstance(cosmo_cls.hubble_time, property)
        assert cosmo_cls.hubble_time.fset is None  # immutable

        # on the instance
        assert cosmo.hubble_time is cosmo._hubble_time
        assert u.allclose(cosmo.hubble_time, (1 / cosmo.H0) << u.Gyr)

    def test_hubble_distance(self, cosmo_cls, cosmo):
        """Test property ``hubble_distance``."""
        # on the class
        assert isinstance(cosmo_cls.hubble_distance, property)
        assert cosmo_cls.hubble_distance.fset is None  # immutable

        # on the instance
        assert cosmo.hubble_distance is cosmo._hubble_distance
        assert cosmo.hubble_distance == (const.c / cosmo._H0).to(u.Mpc)

    def test_critical_density0(self, cosmo_cls, cosmo):
        """Test property ``critical_density0``."""
        # on the class
        assert isinstance(cosmo_cls.critical_density0, property)
        assert cosmo_cls.critical_density0.fset is None  # immutable

        # on the instance
        assert cosmo.critical_density0 is cosmo._critical_density0
        assert cosmo.critical_density0.unit == u.g / u.cm**3

        cd0value = _critdens_const * (cosmo.H0.value * _H0units_to_invs) ** 2
        assert cosmo.critical_density0.value == cd0value

    def test_Ogamma0(self, cosmo_cls, cosmo):
        """Test property ``Ogamma0``."""
        # on the class
        assert isinstance(cosmo_cls.Ogamma0, property)
        assert cosmo_cls.Ogamma0.fset is None  # immutable

        # on the instance
        assert cosmo.Ogamma0 is cosmo._Ogamma0
        # Ogamma cor \propto T^4/rhocrit
        expect = _a_B_c2 * cosmo.Tcmb0.value**4 / cosmo.critical_density0.value
        assert np.allclose(cosmo.Ogamma0, expect)
        # check absolute equality to 0 if Tcmb0 is 0
        if cosmo.Tcmb0 == 0:
            assert cosmo.Ogamma0 == 0

    def test_Onu0(self, cosmo_cls, cosmo):
        """Test property ``Onu0``."""
        # on the class
        assert isinstance(cosmo_cls.Onu0, property)
        assert cosmo_cls.Onu0.fset is None  # immutable

        # on the instance
        assert cosmo.Onu0 is cosmo._Onu0
        # neutrino temperature <= photon temperature since the neutrinos
        # decouple first.
        if cosmo.has_massive_nu:  # Tcmb0 > 0 & has massive
            # check the expected formula
            assert cosmo.Onu0 == cosmo.Ogamma0 * cosmo.nu_relative_density(0)
            # a sanity check on on the ratio of neutrinos to photons
            # technically it could be 1, but not for any of the tested cases.
            assert cosmo.nu_relative_density(0) <= 1
        elif cosmo.Tcmb0 == 0:
            assert cosmo.Onu0 == 0
        else:
            # check the expected formula
            assert cosmo.Onu0 == 0.22710731766 * cosmo._Neff * cosmo.Ogamma0
            # and check compatibility with nu_relative_density
            assert np.allclose(
                cosmo.nu_relative_density(0), 0.22710731766 * cosmo._Neff
            )

    def test_Otot0(self, cosmo):
        """Test :attr:`astropy.cosmology.FLRW.Otot0`."""
        assert (
            cosmo.Otot0
            == cosmo.Om0 + cosmo.Ogamma0 + cosmo.Onu0 + cosmo.Ode0 + cosmo.Ok0
        )

    # ---------------------------------------------------------------
    # Methods

    _FLRW_redshift_methods = get_redshift_methods(
        FLRW, include_private=True, include_z2=False
    )

    @pytest.mark.skipif(not HAS_SCIPY, reason="scipy is not installed")
    @pytest.mark.parametrize("z, exc", invalid_zs)
    @pytest.mark.parametrize("method", _FLRW_redshift_methods)
    def test_redshift_method_bad_input(self, cosmo, method, z, exc):
        """Test all the redshift methods for bad input."""
        with pytest.raises(exc):
            getattr(cosmo, method)(z)

    @pytest.mark.parametrize("z", valid_zs)
    @abc.abstractmethod
    def test_w(self, cosmo, z):
        """Test :meth:`astropy.cosmology.FLRW.w`.

        Since ``w`` is abstract, each test class needs to define further tests.
        """
        # super().test_w(cosmo, z)  # NOT b/c abstract `w(z)`
        w = cosmo.w(z)
        assert np.shape(w) == np.shape(z)  # test same shape
        assert u.Quantity(w).unit == u.one  # test no units or dimensionless

    # -------------------------------------------

    @pytest.mark.parametrize("z", valid_zs)
    def test_Otot(self, cosmo, z):
        """Test :meth:`astropy.cosmology.FLRW.Otot`."""
        # super().test_Otot(cosmo)  # NOT b/c abstract `w(z)`
        assert np.allclose(
            cosmo.Otot(z),
            cosmo.Om(z) + cosmo.Ogamma(z) + cosmo.Onu(z) + cosmo.Ode(z) + cosmo.Ok(z),
        )

    # ---------------------------------------------------------------

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

    # ---------------------------------------------------------------
    # from Cosmology

    def test_clone_change_param(self, cosmo):
        """Test method ``.clone()`` changing a(many) Parameter(s)."""
        super().test_clone_change_param(cosmo)

        # don't change any values
        kwargs = cosmo._init_arguments.copy()
        kwargs.pop("name", None)  # make sure not setting name
        kwargs.pop("meta", None)  # make sure not setting name
        c = cosmo.clone(**kwargs)
        assert c.__class__ == cosmo.__class__
        assert c == cosmo

        # change ``H0``
        # Note that H0 affects Ode0 because it changes Ogamma0
        c = cosmo.clone(H0=100)
        assert c.__class__ == cosmo.__class__
        assert c.name == cosmo.name + " (modified)"
        assert c.H0.value == 100
        for n in set(cosmo.__parameters__) - {"H0"}:
            v = getattr(c, n)
            if v is None:
                assert v is getattr(cosmo, n)
            else:
                assert u.allclose(
                    v, getattr(cosmo, n), atol=1e-4 * getattr(v, "unit", 1)
                )
        assert not u.allclose(c.Ogamma0, cosmo.Ogamma0)
        assert not u.allclose(c.Onu0, cosmo.Onu0)

        # change multiple things
        c = cosmo.clone(name="new name", H0=100, Tcmb0=2.8, meta=dict(zz="tops"))
        assert c.__class__ == cosmo.__class__
        assert c.name == "new name"
        assert c.H0.value == 100
        assert c.Tcmb0.value == 2.8
        assert c.meta == {**cosmo.meta, **dict(zz="tops")}
        for n in set(cosmo.__parameters__) - {"H0", "Tcmb0"}:
            v = getattr(c, n)
            if v is None:
                assert v is getattr(cosmo, n)
            else:
                assert u.allclose(
                    v, getattr(cosmo, n), atol=1e-4 * getattr(v, "unit", 1)
                )
        assert not u.allclose(c.Ogamma0, cosmo.Ogamma0)
        assert not u.allclose(c.Onu0, cosmo.Onu0)
        assert not u.allclose(c.Tcmb0.value, cosmo.Tcmb0.value)

    def test_is_equivalent(self, cosmo):
        """Test :meth:`astropy.cosmology.FLRW.is_equivalent`."""
        super().test_is_equivalent(cosmo)  # pass to CosmologyTest

        # test against a FlatFLRWMixin
        # case (3) in FLRW.is_equivalent
        if isinstance(cosmo, FlatLambdaCDM):
            assert cosmo.is_equivalent(Planck18)
            assert Planck18.is_equivalent(cosmo)
        else:
            assert not cosmo.is_equivalent(Planck18)
            assert not Planck18.is_equivalent(cosmo)

    # ===============================================================
    # Usage Tests

    # TODO: this test should be subsumed by other tests
    @pytest.mark.parametrize("method", ("Om", "Ode", "w", "de_density_scale"))
    def test_distance_broadcast(self, cosmo, method):
        """Test distance methods broadcast z correctly."""
        g = getattr(cosmo, method)
        z = np.linspace(0.1, 1, 6)
        z2d = z.reshape(2, 3)
        z3d = z.reshape(3, 2, 1)

        value_flat = g(z)
        assert value_flat.shape == z.shape

        value_2d = g(z2d)
        assert value_2d.shape == z2d.shape

        value_3d = g(z3d)
        assert value_3d.shape == z3d.shape
        assert u.allclose(value_flat, value_2d.flatten())
        assert u.allclose(value_flat, value_3d.flatten())

    @pytest.mark.skipif(not HAS_SCIPY, reason="scipy required for this test.")
    def test_comoving_distance_example(self, cosmo_cls, args, kwargs, expected):
        """Test :meth:`astropy.cosmology.LambdaCDM.comoving_distance`.

        These do not come from external codes -- they are just internal checks to make
        sure nothing changes if we muck with the distance calculators.
        """
        z = np.array([1.0, 2.0, 3.0, 4.0])

        cosmo = cosmo_cls(*args, **kwargs)
        assert u.allclose(cosmo.comoving_distance(z), expected, rtol=1e-4)


class TestFLRW(FLRWTest):
    """Test :class:`astropy.cosmology.FLRW`."""

    abstract_w = True

    def setup_class(self):
        """
        Setup for testing.
        FLRW is abstract, so tests are done on a subclass.
        """
        super().setup_class(self)

        # make sure SubCosmology is known
        _COSMOLOGY_CLASSES["SubFLRW"] = SubFLRW

        self.cls = SubFLRW

    def teardown_class(self):
        super().teardown_class(self)
        _COSMOLOGY_CLASSES.pop("SubFLRW", None)

    # ===============================================================
    # Method & Attribute Tests

    # ---------------------------------------------------------------
    # Methods

    def test_w(self, cosmo):
        """Test abstract :meth:`astropy.cosmology.FLRW.w`."""
        with pytest.raises(NotImplementedError, match="not implemented"):
            cosmo.w(1)

    def test_Otot(self, cosmo):
        """Test :meth:`astropy.cosmology.FLRW.Otot`."""
        exception = NotImplementedError if HAS_SCIPY else ModuleNotFoundError
        with pytest.raises(exception):
            assert cosmo.Otot(1)

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

    _FLRW_redshift_methods = get_redshift_methods(
        FLRW, include_private=True, include_z2=False
    ) - {"w"}

    @pytest.mark.skipif(not HAS_SCIPY, reason="scipy is not installed")
    @pytest.mark.parametrize("z, exc", invalid_zs)
    @pytest.mark.parametrize("method", _FLRW_redshift_methods)
    def test_redshift_method_bad_input(self, cosmo, method, z, exc):
        """Test all the redshift methods for bad input."""
        with pytest.raises(exc):
            getattr(cosmo, method)(z)

    # ===============================================================
    # Usage Tests

    @pytest.mark.skipif(not HAS_SCIPY, reason="scipy required for this test.")
    @pytest.mark.parametrize("method", ("Om", "Ode", "w", "de_density_scale"))
    def test_distance_broadcast(self, cosmo, method):
        with pytest.raises(NotImplementedError):
            super().test_distance_broadcast(cosmo, method)

    @pytest.mark.skipif(not HAS_SCIPY, reason="scipy required for this test.")
    @pytest.mark.parametrize(
        ("args", "kwargs", "expected"),
        [((70, 0.27, 0.73), {"Tcmb0": 3.0, "Ob0": 0.03}, None)],
    )
    def test_comoving_distance_example(self, cosmo_cls, args, kwargs, expected):
        with pytest.raises(NotImplementedError):
            super().test_comoving_distance_example(cosmo_cls, args, kwargs, expected)


# -----------------------------------------------------------------------------


class ParameterFlatOde0TestMixin(ParameterOde0TestMixin):
    """Tests for `astropy.cosmology.Parameter` Ode0 on a flat Cosmology.

    This will augment or override some tests in ``ParameterOde0TestMixin``.

    Ode0 is a descriptor, which are tested by mixin, here with ``TestFLRW``.
    These tests expect dicts ``_cls_args`` and ``cls_kwargs`` which give the
    args and kwargs for the cosmology class, respectively. See ``TestFLRW``.
    """

    def test_Parameter_Ode0(self, cosmo_cls):
        """Test Parameter ``Ode0`` on the class."""
        super().test_Parameter_Ode0(cosmo_cls)
        assert cosmo_cls.Ode0.derived in (True, np.True_)

    def test_Ode0(self, cosmo):
        """Test no-longer-Parameter ``Ode0``."""
        assert cosmo.Ode0 is cosmo._Ode0
        assert cosmo.Ode0 == 1.0 - (cosmo.Om0 + cosmo.Ogamma0 + cosmo.Onu0)

    def test_init_Ode0(self, cosmo_cls, ba):
        """Test initialization for values of ``Ode0``."""
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.Ode0 == 1.0 - (cosmo.Om0 + cosmo.Ogamma0 + cosmo.Onu0 + cosmo.Ok0)

        # Ode0 is not in the signature
        with pytest.raises(TypeError, match="Ode0"):
            cosmo_cls(*ba.args, **ba.kwargs, Ode0=1)


class FlatFLRWMixinTest(FlatCosmologyMixinTest, ParameterFlatOde0TestMixin):
    """Tests for :class:`astropy.cosmology.FlatFLRWMixin` subclasses.

    E.g to use this class::

        class TestFlatSomeFLRW(FlatFLRWMixinTest, TestSomeFLRW):
            ...
    """

    def setup_class(self):
        """Setup for testing.

        Set up as for regular FLRW test class, but remove dark energy component
        since flat cosmologies are forbidden Ode0 as an argument,
        see ``test_init_subclass``.
        """
        super().setup_class(self)
        self._cls_args.pop("Ode0")

    # ===============================================================
    # Method & Attribute Tests

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
        assert cosmo._Ok0 == 0.0
        assert cosmo._Ode0 == 1.0 - (
            cosmo._Om0 + cosmo._Ogamma0 + cosmo._Onu0 + cosmo._Ok0
        )

    def test_Ok0(self, cosmo_cls, cosmo):
        """Test property ``Ok0``."""
        super().test_Ok0(cosmo_cls, cosmo)

        # for flat cosmologies, Ok0 is not *close* to 0, it *is* 0
        assert cosmo.Ok0 == 0.0

    def test_Otot0(self, cosmo):
        """Test :attr:`astropy.cosmology.FLRW.Otot0`. Should always be 1."""
        super().test_Otot0(cosmo)

        # for flat cosmologies, Otot0 is not *close* to 1, it *is* 1
        assert cosmo.Otot0 == 1.0

    @pytest.mark.parametrize("z", valid_zs)
    def test_Otot(self, cosmo, z):
        """Test :meth:`astropy.cosmology.FLRW.Otot`. Should always be 1."""
        super().test_Otot(cosmo, z)

        # for flat cosmologies, Otot is 1, within precision.
        assert u.allclose(cosmo.Otot(z), 1.0)

    @pytest.mark.skipif(not HAS_SCIPY, reason="scipy is not installed")
    @pytest.mark.parametrize("z, exc", invalid_zs)
    @pytest.mark.parametrize("method", FLRWTest._FLRW_redshift_methods - {"Otot"})
    def test_redshift_method_bad_input(self, cosmo, method, z, exc):
        """Test all the redshift methods for bad input."""
        super().test_redshift_method_bad_input(cosmo, method, z, exc)

    # ---------------------------------------------------------------

    def test_clone_to_nonflat_change_param(self, cosmo):
        """Test method ``.clone()`` changing a(many) Parameter(s)."""
        super().test_clone_to_nonflat_change_param(cosmo)

        # change Ode0, without non-flat
        with pytest.raises(TypeError):
            cosmo.clone(Ode0=1)

        # change to non-flat
        nc = cosmo.clone(to_nonflat=True, Ode0=cosmo.Ode0)
        assert isinstance(nc, cosmo.__nonflatclass__)
        assert nc == cosmo.nonflat

        nc = cosmo.clone(to_nonflat=True, Ode0=1)
        assert nc.Ode0 == 1.0
        assert nc.name == cosmo.name + " (modified)"

    # ---------------------------------------------------------------

    def test_is_equivalent(self, cosmo, nonflatcosmo):
        """Test :meth:`astropy.cosmology.FLRW.is_equivalent`."""
        super().test_is_equivalent(cosmo)  # pass to TestFLRW

        # against non-flat Cosmology
        assert not cosmo.is_equivalent(nonflatcosmo)
        assert not nonflatcosmo.is_equivalent(cosmo)

        # non-flat version of class
        nonflat_cosmo_cls = cosmo.__nonflatclass__
        # keys check in `test_is_equivalent_nonflat_class_different_params`

        # non-flat
        nonflat = nonflat_cosmo_cls(*self.cls_args, Ode0=0.9, **self.cls_kwargs)
        assert not nonflat.is_equivalent(cosmo)
        assert not cosmo.is_equivalent(nonflat)

        # flat, but not FlatFLRWMixin
        flat = nonflat_cosmo_cls(
            *self.cls_args,
            Ode0=1.0 - cosmo.Om0 - cosmo.Ogamma0 - cosmo.Onu0,
            **self.cls_kwargs
        )
        flat._Ok0 = 0.0
        assert flat.is_equivalent(cosmo)
        assert cosmo.is_equivalent(flat)

    def test_repr(self, cosmo_cls, cosmo):
        """
        Test method ``.__repr__()``. Skip non-flat superclass test.
        e.g. `TestFlatLambdaCDDM` -> `FlatFLRWMixinTest`
        vs   `TestFlatLambdaCDDM` -> `TestLambdaCDDM` -> `FlatFLRWMixinTest`
        """
        FLRWTest.test_repr(self, cosmo_cls, cosmo)

        # test eliminated Ode0 from parameters
        assert "Ode0" not in repr(cosmo)
