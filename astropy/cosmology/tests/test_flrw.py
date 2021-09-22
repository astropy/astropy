# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.flrw`."""

##############################################################################
# IMPORTS

import pytest

import astropy.units as u
from astropy.cosmology import (FLRW, FlatLambdaCDM, Flatw0waCDM, FlatwCDM,
                               LambdaCDM, w0waCDM, w0wzCDM, wCDM, wpwaCDM)

from .test_core import FlatCosmologyMixinTest
from .test_core import TestCosmology as CosmologyTest

##############################################################################
# TESTS
##############################################################################


class TestFLRW(CosmologyTest):
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
        self.cls_args = (70 * u.km / u.s / u.Mpc, 0.27 * u.one, 0.689 * u.one)
        self.cls_kwargs = dict(Tcmb0=3.0 * u.K, name="test", meta={"a": "b"})

    def cleanup_class(self):
        _COSMOLOGY_CLASSES.pop("TestFLRW.setup_class.<locals>.SubFLRW")

    # ===============================================================
    # Method & Attribute Tests

    def test_init(self, cosmo_cls):
        # Cosmology accepts any args, kwargs
        cosmo1 = cosmo_cls(*self.cls_args, **self.cls_kwargs)
        self._cosmo_test_init_attr(cosmo1)


# -----------------------------------------------------------------------------


class FlatFLRWMixinTest(FlatCosmologyMixinTest):

    def test_init(self, cosmo_cls):
        super().test_init(cosmo_cls)

        cosmo = cosmo_cls(*self.cls_args, **self.cls_kwargs)
        assert cosmo._Ode0 == 1.0 - cosmo._Om0 - cosmo._Ogamma0 - cosmo._Onu0
        assert cosmo._Ok0 == 0.0


# -----------------------------------------------------------------------------


class TestLambdaCDM(TestFLRW):
    """Test :class:`astropy.cosmology.LambdaCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = LambdaCDM
        self.cls_args = (70 * (u.km / u.s / u.Mpc), 0.27, 0.73)  # H0, Om0, Ode0
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name="test", meta={"a": "b"})


# -----------------------------------------------------------------------------


class TestFlatLambdaCDM(FlatFLRWMixinTest, TestLambdaCDM):
    """Test :class:`astropy.cosmology.FlatLambdaCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = FlatLambdaCDM
        self.cls_args = (70 * (u.km / u.s / u.Mpc), 0.27)  # H0, Om0, Ode0
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name="test", meta={"a": "b"})


# -----------------------------------------------------------------------------


class TestwCDM(TestFLRW):
    """Test :class:`astropy.cosmology.wCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = wCDM
        self.cls_args = (70 * (u.km / u.s / u.Mpc), 0.27, 0.73)  # H0, Om0, Ode0
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name="test", meta={"a": "b"})


# -----------------------------------------------------------------------------


class TestFlatwCDM(FlatFLRWMixinTest, TestwCDM):
    """Test :class:`astropy.cosmology.FlatwCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = FlatwCDM
        self.cls_args = (70 * (u.km / u.s / u.Mpc), 0.27)  # H0, Om0
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name="test", meta={"a": "b"})


# -----------------------------------------------------------------------------


class Testw0waCDM(TestFLRW):
    """Test :class:`astropy.cosmology.w0waCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = w0waCDM
        self.cls_args = (70 * (u.km / u.s / u.Mpc), 0.27, 0.73)  # H0, Om0, Ode0
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name="test", meta={"a": "b"})


# -----------------------------------------------------------------------------


class TestFlatw0waCDM(FlatFLRWMixinTest, Testw0waCDM):
    """Test :class:`astropy.cosmology.Flatw0waCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = Flatw0waCDM
        self.cls_args = (70 * (u.km / u.s / u.Mpc), 0.27)  # H0, Om0
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name="test", meta={"a": "b"})


# -----------------------------------------------------------------------------


class TestwpwaCDM(TestFLRW):
    """Test :class:`astropy.cosmology.wpwaCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = wpwaCDM
        self.cls_args = (70 * (u.km / u.s / u.Mpc), 0.27, 0.73)  # H0, Om0, Ode0
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name="test", meta={"a": "b"})


# -----------------------------------------------------------------------------


class Testw0wzCDM(TestFLRW):
    """Test :class:`astropy.cosmology.w0wzCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = w0wzCDM
        self.cls_args = (70 * (u.km / u.s / u.Mpc), 0.27, 0.73)  # H0, Om0, Ode0
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name="test", meta={"a": "b"})
