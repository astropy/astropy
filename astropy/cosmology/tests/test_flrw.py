# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.flrw`."""

##############################################################################
# IMPORTS

import abc

import pytest

import numpy as np

import astropy.units as u
from astropy.cosmology import (FLRW, FlatLambdaCDM, Flatw0waCDM, FlatwCDM,
                               LambdaCDM, Planck18, w0waCDM, w0wzCDM, wCDM, wpwaCDM)
from astropy.cosmology.core import _COSMOLOGY_CLASSES
from astropy.cosmology.flrw import ellipkinc, hyp2f1, quad
from astropy.utils.compat.optional_deps import HAS_SCIPY

from .test_core import CosmologySubclassTest as CosmologyTest
from .test_core import FlatCosmologyMixinTest

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
        # Cosmology accepts any args, kwargs
        cosmo1 = cosmo_cls(*self.cls_args, **self.cls_kwargs)
        self._cosmo_test_init_attr(cosmo1)

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


class FlatFLRWMixinTest(FlatCosmologyMixinTest):
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


# -----------------------------------------------------------------------------


class TestLambdaCDM(FLRWSubclassTest):
    """Test :class:`astropy.cosmology.LambdaCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = LambdaCDM
        self.cls_args = (70 * (u.km / u.s / u.Mpc), 0.27, 0.73)  # H0, Om0, Ode0
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name=self.__class__.__name__, meta={"a": "b"})


# -----------------------------------------------------------------------------


class TestFlatLambdaCDM(FlatFLRWMixinTest, TestLambdaCDM):
    """Test :class:`astropy.cosmology.FlatLambdaCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = FlatLambdaCDM
        self.cls_args = (70 * (u.km / u.s / u.Mpc), 0.27)  # H0, Om0, Ode0
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name=self.__class__.__name__, meta={"a": "b"})


# -----------------------------------------------------------------------------


class TestwCDM(FLRWSubclassTest):
    """Test :class:`astropy.cosmology.wCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = wCDM
        self.cls_args = (70 * (u.km / u.s / u.Mpc), 0.27, 0.73)  # H0, Om0, Ode0
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name=self.__class__.__name__, meta={"a": "b"})


# -----------------------------------------------------------------------------


class TestFlatwCDM(FlatFLRWMixinTest, TestwCDM):
    """Test :class:`astropy.cosmology.FlatwCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = FlatwCDM
        self.cls_args = (70 * (u.km / u.s / u.Mpc), 0.27)  # H0, Om0
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name=self.__class__.__name__, meta={"a": "b"})


# -----------------------------------------------------------------------------


class Testw0waCDM(FLRWSubclassTest):
    """Test :class:`astropy.cosmology.w0waCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = w0waCDM
        self.cls_args = (70 * (u.km / u.s / u.Mpc), 0.27, 0.73)  # H0, Om0, Ode0
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name=self.__class__.__name__, meta={"a": "b"})


# -----------------------------------------------------------------------------


class TestFlatw0waCDM(FlatFLRWMixinTest, Testw0waCDM):
    """Test :class:`astropy.cosmology.Flatw0waCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = Flatw0waCDM
        self.cls_args = (70 * (u.km / u.s / u.Mpc), 0.27)  # H0, Om0
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name=self.__class__.__name__, meta={"a": "b"})


# -----------------------------------------------------------------------------


class TestwpwaCDM(FLRWSubclassTest):
    """Test :class:`astropy.cosmology.wpwaCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = wpwaCDM
        self.cls_args = (70 * (u.km / u.s / u.Mpc), 0.27, 0.73)  # H0, Om0, Ode0
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name=self.__class__.__name__, meta={"a": "b"})


# -----------------------------------------------------------------------------


class Testw0wzCDM(FLRWSubclassTest):
    """Test :class:`astropy.cosmology.w0wzCDM`."""

    def setup_class(self):
        """Setup for testing."""
        self.cls = w0wzCDM
        self.cls_args = (70 * (u.km / u.s / u.Mpc), 0.27, 0.73)  # H0, Om0, Ode0
        self.cls_kwargs = dict(Tcmb0=3 * u.K, name=self.__class__.__name__, meta={"a": "b"})
