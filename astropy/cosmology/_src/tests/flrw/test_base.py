# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.flrw.base`.

This module sets up the tests for subclasses of :class:`astropy.cosmology.FLRW`. The
tests for the specific abstract class :class:`astropy.cosmology.FLRW` are in
``test_flrw``.

"""

import abc
from functools import cached_property

import numpy as np
import pytest

import astropy.constants as const
import astropy.units as u
from astropy.cosmology import FLRW, FlatLambdaCDM, LambdaCDM, Planck18
from astropy.cosmology._src.core import _COSMOLOGY_CLASSES, dataclass_decorator
from astropy.cosmology._src.flrw.base import _a_B_c2, quad
from astropy.cosmology._src.tests.helper import get_redshift_methods
from astropy.cosmology._src.tests.test_core import (
    CosmologyTest,
    FlatCosmologyMixinTest,
    invalid_zs,
    valid_zs,
)
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.compat.optional_deps import HAS_PANDAS, HAS_SCIPY

from .conftest import filter_keys_from_items
from .test_parameters import (
    ParameterFlatOde0TestMixin,
    ParameterH0TestMixin,
    Parameterm_nuTestMixin,
    ParameterNeffTestMixin,
    ParameterOb0TestMixin,
    ParameterOde0TestMixin,
    ParameterOm0TestMixin,
    ParameterTcmb0TestMixin,
)

##############################################################################
# TESTS
##############################################################################


@pytest.mark.skipif(HAS_SCIPY, reason="scipy is installed")
def test_optional_deps_functions():
    """Test stand-in functions when optional dependencies not installed."""
    with pytest.raises(ModuleNotFoundError, match="No module named 'scipy.integrate'"):
        quad()


##############################################################################


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
        """Test ``cached_property`` ``Odm0``."""
        # on the class
        assert isinstance(cosmo_cls.Odm0, cached_property)

        # on the instance
        assert np.allclose(cosmo.Odm0, cosmo.Om0 - cosmo.Ob0)

    def test_Ok0(self, cosmo_cls, cosmo):
        """Test ``cached_property`` ``Ok0``."""
        # on the class
        assert isinstance(cosmo_cls.Ok0, cached_property)

        # on the instance
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
        """Test ``cached_property`` ``Tnu0``."""
        # on the class
        assert isinstance(cosmo_cls.Tnu0, cached_property)

        # on the instance
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
        """Test ``cached_property`` ``h``."""
        # on the class
        assert isinstance(cosmo_cls.h, cached_property)

        # on the instance
        assert np.allclose(cosmo.h, cosmo.H0.value / 100.0)

    def test_hubble_time(self, cosmo_cls, cosmo):
        """Test ``cached_property`` ``hubble_time``."""
        # on the class
        assert isinstance(cosmo_cls.hubble_time, cached_property)

        # on the instance
        assert u.allclose(cosmo.hubble_time, (1 / cosmo.H0) << u.Gyr)

    def test_hubble_distance(self, cosmo_cls, cosmo):
        """Test ``cached_property`` ``hubble_distance``."""
        # on the class
        assert isinstance(cosmo_cls.hubble_distance, cached_property)

        # on the instance
        assert cosmo.hubble_distance == (const.c / cosmo.H0).to(u.Mpc)

    def test_critical_density0(self, cosmo_cls, cosmo):
        """Test ``cached_property`` ``critical_density0``."""
        # on the class
        assert isinstance(cosmo_cls.critical_density0, cached_property)

        # on the instance
        assert cosmo.critical_density0.unit == u.g / u.cm**3
        assert u.allclose(  # sanity check
            cosmo.critical_density0, 3 * cosmo.H0**2 / (8 * np.pi * const.G)
        )

    def test_Ogamma0(self, cosmo_cls, cosmo):
        """Test ``cached_property`` ``Ogamma0``."""
        # on the class
        assert isinstance(cosmo_cls.Ogamma0, cached_property)

        # on the instance
        # Ogamma cor \propto T^4/rhocrit
        expect = _a_B_c2 * cosmo.Tcmb0.value**4 / cosmo.critical_density0.value
        assert np.allclose(cosmo.Ogamma0, expect)
        # check absolute equality to 0 if Tcmb0 is 0
        if cosmo.Tcmb0 == 0:
            assert cosmo.Ogamma0 == 0

    def test_Onu0(self, cosmo_cls, cosmo):
        """Test ``cached_property`` ``Onu0``."""
        # on the class
        assert isinstance(cosmo_cls.Onu0, cached_property)

        # on the instance
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
            assert cosmo.Onu0 == 0.22710731766 * cosmo.__dict__["Neff"] * cosmo.Ogamma0
            # and check compatibility with nu_relative_density
            assert np.allclose(
                cosmo.nu_relative_density(0), 0.22710731766 * cosmo.__dict__["Neff"]
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
    @pytest.mark.parametrize("method", sorted(_FLRW_redshift_methods))
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

    def test_scale_factor0(self, cosmo):
        """Test :meth:`astropy.cosmology.FLRW.scale_factor`."""
        assert isinstance(cosmo.scale_factor0, u.Quantity)
        assert cosmo.scale_factor0.unit == u.one
        assert cosmo.scale_factor0 == 1
        assert np.allclose(cosmo.scale_factor0, cosmo.scale_factor(0))

    @pytest.mark.parametrize("z", valid_zs)
    def test_scale_factor(self, cosmo, z):
        """Test :meth:`astropy.cosmology.FLRW.scale_factor`."""
        assert np.allclose(cosmo.scale_factor(z), 1 / (1 + np.array(z)))

    # -------------------------------------------

    @pytest.mark.skipif(not HAS_SCIPY, reason="scipy required for this test.")
    def test_comoving_distance_1arg_equal_to_2arg(self, cosmo):
        """Test :meth:`astropy.cosmology.FLRW.comoving_distance`."""
        # Special case of z1 = 0
        z = np.linspace(0, 1, 10)
        assert u.allclose(cosmo.comoving_distance(z), cosmo.comoving_distance(0, z))

        # General case of z1, z2
        z1 = z
        z2 = z + 1
        assert u.allclose(
            cosmo.comoving_distance(z2) - cosmo.comoving_distance(z1),
            cosmo.comoving_distance(z1, z2),
        )

    @pytest.mark.skipif(
        not (HAS_PANDAS and HAS_SCIPY), reason="requires pandas and scipy"
    )
    def test_luminosity_distance_pandas(self, cosmo):
        """Test :meth:`astropy.cosmology.FLRW.luminosity_distance`.

        Regression test for https://github.com/astropy/astropy/issues/15576.
        """
        import pandas as pd

        z = pd.Series([0.1, 0.2, 0.3])
        d = cosmo.luminosity_distance(z)

        assert isinstance(d, u.Quantity)
        assert d.unit == u.Mpc
        np.testing.assert_array_equal(d, cosmo.luminosity_distance(np.array(z)))

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
        kwargs = dict(cosmo.parameters)
        c = cosmo.clone(**kwargs)
        assert c.__class__ == cosmo.__class__
        assert c == cosmo

        # change ``H0``
        # Note that H0 affects Ode0 because it changes Ogamma0
        c = cosmo.clone(H0=100)
        assert c.__class__ == cosmo.__class__
        assert c.name == cosmo.name + " (modified)"
        assert c.H0.value == 100
        for n, v in filter_keys_from_items(c.parameters, ("H0",)):
            v_expect = getattr(cosmo, n)
            assert_quantity_allclose(v, v_expect, atol=1e-4 * getattr(v, "unit", 1))
        assert not u.allclose(c.Ogamma0, cosmo.Ogamma0)
        assert not u.allclose(c.Onu0, cosmo.Onu0)

        # change multiple things
        c = cosmo.clone(name="new name", H0=100, Tcmb0=2.8, meta=dict(zz="tops"))
        assert c.__class__ == cosmo.__class__
        assert c.name == "new name"
        assert c.H0.value == 100
        assert c.Tcmb0.value == 2.8
        assert c.meta == {**cosmo.meta, **dict(zz="tops")}
        for n, v in filter_keys_from_items(c.parameters, ("H0", "Tcmb0")):
            v_expect = getattr(cosmo, n)
            assert_quantity_allclose(v, v_expect, atol=1e-4 * getattr(v, "unit", 1))
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


# ==============================================================================


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

            @dataclass_decorator
            class HASOde0SubClass(cosmo_cls):
                def __init__(self, Ode0):
                    pass

            _COSMOLOGY_CLASSES.pop(HASOde0SubClass.__qualname__, None)

    # ---------------------------------------------------------------
    # instance-level

    def test_init(self, cosmo_cls):
        super().test_init(cosmo_cls)

        cosmo = cosmo_cls(*self.cls_args, **self.cls_kwargs)
        assert cosmo.Ok0 == 0.0
        assert cosmo.Ode0 == 1.0 - (cosmo.Om0 + cosmo.Ogamma0 + cosmo.Onu0 + cosmo.Ok0)

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
    @pytest.mark.parametrize(
        "method", sorted(FLRWTest._FLRW_redshift_methods - {"Otot"})
    )
    def test_redshift_method_bad_input(self, cosmo, method, z, exc):
        """Test all the redshift methods for bad input."""
        super().test_redshift_method_bad_input(cosmo, method, z, exc)

    # ---------------------------------------------------------------

    def test_clone_to_nonflat_change_param(self, cosmo):
        """Test method ``.clone()`` changing a(many) Parameter(s)."""
        super().test_clone_to_nonflat_change_param(cosmo)

        # change Ode0, without non-flat
        msg = "Cannot set 'Ode0' in clone unless 'to_nonflat=True'. "
        with pytest.raises(ValueError, match=msg):
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

        # Flat, but not FlatFLRWMixin
        # This will require forcing flatness by overriding attribute values.
        # Since Cosmology is frozen, the easiest way is via __dict__.
        flat = nonflat_cosmo_cls(
            *self.cls_args,
            Ode0=1.0 - cosmo.Om0 - cosmo.Ogamma0 - cosmo.Onu0,
            **self.cls_kwargs,
        )
        flat.__dict__["Ok0"] = 0.0  # manually forcing flatness by setting `Ok0`.
        assert flat.is_equivalent(cosmo)
        assert cosmo.is_equivalent(flat)

    def test_repr(self, cosmo_cls, cosmo):
        """
        Test method ``.__repr__()``. Skip non-flat superclass test.
        e.g. `TestFlatLambdaCDDM` -> `FlatFLRWMixinTest`
        vs   `TestFlatLambdaCDDM` -> `TestLambdaCDDM` -> `FlatFLRWMixinTest`
        """
        # test eliminated Ode0 from parameters
        assert "Ode0" not in repr(cosmo)
