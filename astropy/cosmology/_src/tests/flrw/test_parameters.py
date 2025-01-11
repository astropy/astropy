# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Parameter test mixin classes."""

from __future__ import annotations

import copy
from typing import TYPE_CHECKING

import numpy as np
import pytest

import astropy.units as u
from astropy.cosmology import FlatFLRWMixin, Parameter
from astropy.cosmology._src.parameter import MISSING
from astropy.cosmology._src.tests.test_core import ParameterTestMixin
from astropy.tests.helper import assert_quantity_allclose

if TYPE_CHECKING:
    from inspect import BoundArguments

    from astropy.cosmology import Cosmology


class ParameterH0TestMixin(ParameterTestMixin):
    """Tests for `astropy.cosmology.Parameter` H0 on a Cosmology.

    H0 is a descriptor, which are tested by mixin, here with ``TestFLRW``.
    These tests expect dicts ``_cls_args`` and ``cls_kwargs`` which give the
    args and kwargs for the cosmology class, respectively. See ``TestFLRW``.
    """

    def test_H0(self, cosmo_cls: type[Cosmology], cosmo: Cosmology):
        """Test Parameter ``H0``."""
        unit = u.Unit("km/(s Mpc)")

        # on the class
        H0 = cosmo_cls.parameters["H0"]
        assert isinstance(H0, Parameter)
        assert "Hubble constant" in H0.__doc__
        assert H0.unit == unit
        assert H0.default is MISSING

        # validation
        assert H0.validate(cosmo, 1) == 1 * unit
        assert H0.validate(cosmo, 10 * unit) == 10 * unit
        with pytest.raises(ValueError, match="H0 is a non-scalar quantity"):
            H0.validate(cosmo, [1, 2])

        # on the instance
        assert cosmo.H0 is cosmo.__dict__["H0"]
        assert cosmo.H0 == self._cls_args["H0"]
        assert isinstance(cosmo.H0, u.Quantity) and cosmo.H0.unit == unit

    def test_init_H0(self, cosmo_cls: type[Cosmology], ba: BoundArguments):
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


# =============================================================================


class ParameterOm0TestMixin(ParameterTestMixin):
    """Tests for `astropy.cosmology.Parameter` Om0 on a Cosmology.

    Om0 is a descriptor, which are tested by mixin, here with ``TestFLRW``.
    These tests expect dicts ``_cls_args`` and ``cls_kwargs`` which give the
    args and kwargs for the cosmology class, respectively. See ``TestFLRW``.
    """

    def test_Om0(self, cosmo_cls: type[Cosmology], cosmo: Cosmology):
        """Test Parameter ``Om0``."""
        # on the class
        Om0 = cosmo_cls.parameters["Om0"]
        assert isinstance(Om0, Parameter)
        assert "Omega matter" in Om0.__doc__
        assert Om0.default is MISSING

        # validation
        assert Om0.validate(cosmo, 1) == 1
        assert Om0.validate(cosmo, 10 * u.one) == 10
        with pytest.raises(ValueError, match="Om0 cannot be negative"):
            Om0.validate(cosmo, -1)

        # on the instance
        assert cosmo.Om0 is cosmo.__dict__["Om0"]
        assert cosmo.Om0 == self._cls_args["Om0"]
        assert isinstance(cosmo.Om0, float)

    def test_init_Om0(self, cosmo_cls: type[Cosmology], ba: BoundArguments):
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


# =============================================================================


class ParameterOde0TestMixin(ParameterTestMixin):
    """Tests for `astropy.cosmology.Parameter` Ode0 on a Cosmology.

    Ode0 is a descriptor, which are tested by mixin, here with ``TestFLRW``.
    These tests expect dicts ``_cls_args`` and ``cls_kwargs`` which give the
    args and kwargs for the cosmology class, respectively. See ``TestFLRW``.
    """

    def test_Parameter_Ode0(self, cosmo_cls: type[Cosmology]):
        """Test Parameter ``Ode0`` on the class."""
        Ode0 = cosmo_cls.parameters.get(
            "Ode0", cosmo_cls._derived_parameters.get("Ode0")
        )
        assert isinstance(Ode0, Parameter)
        assert "Omega dark energy" in Ode0.__doc__
        if issubclass(cosmo_cls, FlatFLRWMixin):
            assert Ode0.default == 0
        else:
            assert Ode0.default is MISSING

    def test_Parameter_Ode0_validation(
        self, cosmo_cls: type[Cosmology], cosmo: Cosmology
    ):
        """Test Parameter ``Ode0`` validation."""
        Ode0 = cosmo_cls.parameters.get(
            "Ode0", cosmo_cls._derived_parameters.get("Ode0")
        )
        assert Ode0.validate(cosmo, 1.1) == 1.1
        assert Ode0.validate(cosmo, 10 * u.one) == 10.0
        with pytest.raises(TypeError, match="only dimensionless"):
            Ode0.validate(cosmo, 10 * u.km)

    def test_Ode0(self, cosmo: Cosmology):
        """Test Parameter ``Ode0`` validation."""
        # if Ode0 is a parameter, test its value
        assert cosmo.Ode0 is cosmo.__dict__["Ode0"]
        assert cosmo.Ode0 == self._cls_args["Ode0"]
        assert isinstance(cosmo.Ode0, float)

    def test_init_Ode0(self, cosmo_cls: type[Cosmology], ba: BoundArguments):
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
        assert_quantity_allclose(cosmo.Ode([0, 1, 2, 3]), [0, 0, 0, 0])
        assert_quantity_allclose(cosmo.Ode(1), 0)

        # Must be dimensionless or have no units. Errors otherwise.
        ba.arguments["Ode0"] = 10 * u.km
        with pytest.raises(TypeError, match="only dimensionless"):
            cosmo_cls(*ba.args, **ba.kwargs)


# =============================================================================


class ParameterTcmb0TestMixin(ParameterTestMixin):
    """Tests for `astropy.cosmology.Parameter` Tcmb0 on a Cosmology.

    Tcmb0 is a descriptor, which are tested by mixin, here with ``TestFLRW``.
    These tests expect dicts ``_cls_args`` and ``cls_kwargs`` which give the
    args and kwargs for the cosmology class, respectively. See ``TestFLRW``.
    """

    def test_Tcmb0(self, cosmo_cls: type[Cosmology], cosmo: Cosmology):
        """Test Parameter ``Tcmb0``."""
        # on the class
        Tcmb0 = cosmo_cls.parameters["Tcmb0"]
        assert isinstance(Tcmb0, Parameter)
        assert "Temperature of the CMB" in Tcmb0.__doc__
        assert Tcmb0.unit == u.K
        assert Tcmb0.default == 0.0 * u.K

        # validation
        assert Tcmb0.validate(cosmo, 1) == 1 * u.K
        assert Tcmb0.validate(cosmo, 10 * u.K) == 10 * u.K
        with pytest.raises(ValueError, match="Tcmb0 is a non-scalar quantity"):
            Tcmb0.validate(cosmo, [1, 2])

        # on the instance
        assert cosmo.Tcmb0 is cosmo.__dict__["Tcmb0"]
        assert cosmo.Tcmb0 == self.cls_kwargs["Tcmb0"]
        assert isinstance(cosmo.Tcmb0, u.Quantity) and cosmo.Tcmb0.unit == u.K

    def test_init_Tcmb0(self, cosmo_cls: type[Cosmology], ba: BoundArguments):
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


# =============================================================================


class ParameterNeffTestMixin(ParameterTestMixin):
    """Tests for `astropy.cosmology.Parameter` Neff on a Cosmology.

    Neff is a descriptor, which are tested by mixin, here with ``TestFLRW``.
    These tests expect dicts ``_cls_args`` and ``cls_kwargs`` which give the
    args and kwargs for the cosmology class, respectively. See ``TestFLRW``.
    """

    def test_Neff(self, cosmo_cls: type[Cosmology], cosmo: Cosmology):
        """Test Parameter ``Neff``."""
        # on the class
        Neff = cosmo_cls.parameters["Neff"]
        assert isinstance(Neff, Parameter)
        assert "Number of effective neutrino species" in Neff.__doc__
        assert Neff.default == 3.04

        # validation
        assert Neff.validate(cosmo, 1) == 1
        assert Neff.validate(cosmo, 10 * u.one) == 10
        with pytest.raises(ValueError, match="Neff cannot be negative"):
            Neff.validate(cosmo, -1)

        # on the instance
        assert cosmo.Neff is cosmo.__dict__["Neff"]
        assert cosmo.Neff == self.cls_kwargs.get("Neff", 3.04)
        assert isinstance(cosmo.Neff, float)

    def test_init_Neff(self, cosmo_cls: type[Cosmology], ba: BoundArguments):
        """Test initialization for values of ``Neff``."""
        # test that it works with units
        ba.arguments["Neff"] = (
            cosmo_cls.parameters["Neff"].default << u.one
        )  # ensure units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.Neff == ba.arguments["Neff"]

        # also without units
        ba.arguments["Neff"] = ba.arguments["Neff"].value  # strip units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.Neff == ba.arguments["Neff"]

        ba.arguments["Neff"] = -1
        with pytest.raises(ValueError):
            cosmo_cls(*ba.args, **ba.kwargs)


# =============================================================================


class Parameterm_nuTestMixin(ParameterTestMixin):
    """Tests for `astropy.cosmology.Parameter` m_nu on a Cosmology.

    m_nu is a descriptor, which are tested by mixin, here with ``TestFLRW``.
    These tests expect dicts ``_cls_args`` and ``cls_kwargs`` which give the
    args and kwargs for the cosmology class, respectively. See ``TestFLRW``.
    """

    def test_m_nu(self, cosmo_cls: type[Cosmology], cosmo: Cosmology):
        """Test Parameter ``m_nu``."""
        # on the class
        m_nu = cosmo_cls.parameters["m_nu"]
        assert isinstance(m_nu, Parameter)
        assert "Mass of neutrino species" in m_nu.__doc__
        assert m_nu.unit == u.eV
        assert m_nu.equivalencies == u.mass_energy()
        assert m_nu.default == 0.0 * u.eV

        # on the instance
        # assert cosmo.m_nu is cosmo._m_nu
        assert_quantity_allclose(cosmo.m_nu, [0.0, 0.0, 0.0] * u.eV)

        # set differently depending on the other inputs
        if cosmo.Tnu0.value == 0:
            assert cosmo.m_nu is None
        elif not cosmo._massivenu:  # only massless
            assert_quantity_allclose(cosmo.m_nu, 0 * u.eV)
        elif self._nmasslessnu == 0:  # only massive
            assert cosmo.m_nu == cosmo._massivenu_mass
        else:  # a mix -- the most complicated case
            assert_quantity_allclose(cosmo.m_nu[: self._nmasslessnu], 0 * u.eV)
            assert_quantity_allclose(
                cosmo.m_nu[self._nmasslessnu], cosmo._massivenu_mass
            )

    def test_init_m_nu(self, cosmo_cls: type[Cosmology], ba: BoundArguments):
        """Test initialization for values of ``m_nu``.

        Note this requires the class to have a property ``has_massive_nu``.
        """
        # Test that it works when m_nu has units.
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        np.testing.assert_array_equal(
            cosmo.m_nu, cosmo_cls.parameters["m_nu"].default
        )  # (& checks len, unit)
        assert not cosmo.has_massive_nu
        assert cosmo.m_nu.unit == u.eV  # explicitly check unit once.

        # And it works when m_nu doesn't have units.
        ba.arguments["m_nu"] = cosmo_cls.parameters["m_nu"].default.value  # strip units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        np.testing.assert_array_equal(cosmo.m_nu.value, ba.arguments["m_nu"])
        assert not cosmo.has_massive_nu

        # A negative m_nu raises an exception.
        tba = copy.copy(ba)
        tba.arguments["m_nu"] = u.Quantity([-0.3, 0.2, 0.1], u.eV)
        with pytest.raises(ValueError, match="invalid"):
            cosmo_cls(*tba.args, **tba.kwargs)

    def test_init_m_nu_and_Neff(self, cosmo_cls: type[Cosmology], ba: BoundArguments):
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
        assert_quantity_allclose(cosmo.m_nu, 0 * u.eV)
        # TODO! move this test when create ``test_nu_relative_density``
        assert_quantity_allclose(
            cosmo.nu_relative_density(1.0), 0.22710731766 * 4.05, rtol=1e-6
        )

        # All massive neutrinos case, len from Neff
        tba.arguments["m_nu"] = 0.1 * u.eV
        cosmo = cosmo_cls(*tba.args, **tba.kwargs)
        assert cosmo.has_massive_nu
        assert len(cosmo.m_nu) == 4
        assert cosmo.m_nu.unit == u.eV
        assert_quantity_allclose(cosmo.m_nu, [0.1, 0.1, 0.1, 0.1] * u.eV)

    def test_init_m_nu_override_by_Tcmb0(
        self, cosmo_cls: type[Cosmology], ba: BoundArguments
    ):
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


# =============================================================================


class ParameterOb0TestMixin(ParameterTestMixin):
    """Tests for `astropy.cosmology.Parameter` Ob0 on a Cosmology.

    Ob0 is a descriptor, which are tested by mixin, here with ``TestFLRW``.
    These tests expect dicts ``_cls_args`` and ``cls_kwargs`` which give the
    args and kwargs for the cosmology class, respectively. See ``TestFLRW``.
    """

    def test_Ob0(self, cosmo_cls: type[Cosmology], cosmo: Cosmology):
        """Test Parameter ``Ob0``."""
        # on the class
        Ob0 = cosmo_cls.parameters["Ob0"]
        assert isinstance(Ob0, Parameter)
        assert "Omega baryon;" in Ob0.__doc__
        assert Ob0.default == 0.0

        # validation
        with pytest.warns(DeprecationWarning, match="Ob0=None is deprecated"):
            assert Ob0.validate(cosmo, None) == 0.0
        assert Ob0.validate(cosmo, 0.1) == 0.1
        assert Ob0.validate(cosmo, 0.1 * u.one) == 0.1
        with pytest.raises(ValueError, match="Ob0 cannot be negative"):
            Ob0.validate(cosmo, -1)
        with pytest.raises(ValueError, match="baryonic density can not be larger"):
            Ob0.validate(cosmo, cosmo.Om0 + 1)

        # on the instance
        assert cosmo.Ob0 is cosmo.__dict__["Ob0"]
        assert cosmo.Ob0 == 0.03

    def test_init_Ob0(self, cosmo_cls: type[Cosmology], ba: BoundArguments):
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
            assert_quantity_allclose(cosmo.Ob(1), 0)
            assert_quantity_allclose(cosmo.Ob([0, 1, 2, 3]), [0, 0, 0, 0])

        # Negative Ob0 errors
        tba = copy.copy(ba)
        tba.arguments["Ob0"] = -0.04
        with pytest.raises(ValueError, match="Ob0 cannot be negative"):
            cosmo_cls(*tba.args, **tba.kwargs)

        # Ob0 > Om0 errors
        tba.arguments["Ob0"] = tba.arguments["Om0"] + 0.1
        with pytest.raises(ValueError, match="baryonic density can not be larger"):
            cosmo_cls(*tba.args, **tba.kwargs)

        # No baryons specified means baryons = 0 and gives a warning.
        tba = copy.copy(ba)
        tba.arguments["Ob0"] = None
        with pytest.warns(DeprecationWarning, match="Ob0=None is deprecated"):
            cosmo = cosmo_cls(*tba.args, **tba.kwargs)

        # In FLRW `Ob(z)` requires `w(z)`.
        if not self.abstract_w:
            assert cosmo.Ob(1) == 0

        # The default value is None
        assert cosmo_cls.parameters["Ob0"].default == 0.0


# =============================================================================


class ParameterFlatOde0TestMixin(ParameterOde0TestMixin):
    """Tests for `astropy.cosmology.Parameter` Ode0 on a flat Cosmology.

    This will augment or override some tests in ``ParameterOde0TestMixin``.

    Ode0 is a descriptor, which are tested by mixin, here with ``TestFLRW``.
    These tests expect dicts ``_cls_args`` and ``cls_kwargs`` which give the
    args and kwargs for the cosmology class, respectively. See ``TestFLRW``.
    """

    def test_Parameter_Ode0(self, cosmo_cls: type[Cosmology]):
        """Test Parameter ``Ode0`` on the class."""
        super().test_Parameter_Ode0(cosmo_cls)
        Ode0 = cosmo_cls.parameters.get("Ode0", cosmo_cls._derived_parameters["Ode0"])
        assert Ode0.derived in (True, np.True_)

    def test_Ode0(self, cosmo: Cosmology):
        """Test no-longer-Parameter ``Ode0``."""
        assert cosmo.Ode0 is cosmo.__dict__["Ode0"]
        assert cosmo.Ode0 == 1.0 - (cosmo.Om0 + cosmo.Ogamma0 + cosmo.Onu0)

    def test_init_Ode0(self, cosmo_cls: type[Cosmology], ba: BoundArguments):
        """Test initialization for values of ``Ode0``."""
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.Ode0 == 1.0 - (cosmo.Om0 + cosmo.Ogamma0 + cosmo.Onu0 + cosmo.Ok0)

        # Ode0 is not in the signature
        with pytest.raises(TypeError, match="Ode0"):
            cosmo_cls(*ba.args, **ba.kwargs, Ode0=1)
