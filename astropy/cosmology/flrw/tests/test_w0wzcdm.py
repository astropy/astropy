# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.flrw.w0wzcdm`."""

# THIRD PARTY
import numpy as np
import pytest

# LOCAL
import astropy.units as u
from astropy.cosmology import Flatw0wzCDM, w0wzCDM
from astropy.cosmology.parameter import Parameter
from astropy.cosmology.tests.test_core import ParameterTestMixin, make_valid_zs
from astropy.utils.compat.optional_deps import HAS_SCIPY

from .test_base import FlatFLRWMixinTest, FLRWTest
from .test_w0cdm import Parameterw0TestMixin

##############################################################################
# PARAMETERS

COMOVING_DISTANCE_EXAMPLE_KWARGS = {"w0": -0.9, "wz": 0.1, "Tcmb0": 0.0}

valid_zs = make_valid_zs(max_z=400)[-1]


##############################################################################
# TESTS
##############################################################################


class ParameterwzTestMixin(ParameterTestMixin):
    """Tests for `astropy.cosmology.Parameter` wz on a Cosmology.

    wz is a descriptor, which are tested by mixin, here with ``TestFLRW``.
    These tests expect dicts ``_cls_args`` and ``cls_kwargs`` which give the
    args and kwargs for the cosmology class, respectively. See ``TestFLRW``.
    """

    def test_wz(self, cosmo_cls, cosmo):
        """Test Parameter ``wz``."""
        # on the class
        assert isinstance(cosmo_cls.wz, Parameter)
        assert "Derivative of the dark energy" in cosmo_cls.wz.__doc__
        assert cosmo_cls.wz.unit is None

        # on the instance
        assert cosmo.wz is cosmo._wz
        assert cosmo.wz == self.cls_kwargs["wz"]

    def test_init_wz(self, cosmo_cls, ba):
        """Test initialization for values of ``wz``."""
        # test that it works with units
        ba.arguments["wz"] = ba.arguments["wz"] << u.one  # ensure units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.wz == ba.arguments["wz"]

        # also without units
        ba.arguments["wz"] = ba.arguments["wz"].value  # strip units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.wz == ba.arguments["wz"]

        # must be dimensionless
        ba.arguments["wz"] = 10 * u.km
        with pytest.raises(TypeError):
            cosmo_cls(*ba.args, **ba.kwargs)


class Testw0wzCDM(FLRWTest, Parameterw0TestMixin, ParameterwzTestMixin):
    """Test :class:`astropy.cosmology.w0wzCDM`."""

    def setup_class(self):
        """Setup for testing."""
        super().setup_class(self)
        self.cls = w0wzCDM
        self.cls_kwargs.update(w0=-1, wz=0.5)

    # ===============================================================
    # Method & Attribute Tests

    def test_clone_change_param(self, cosmo):
        """Test method ``.clone()`` changing a(many) Parameter(s)."""
        super().test_clone_change_param(cosmo)

        # `w` params
        c = cosmo.clone(w0=0.1, wz=0.2)
        assert c.w0 == 0.1
        assert c.wz == 0.2
        for n in set(cosmo.__parameters__) - {"w0", "wz"}:
            v = getattr(c, n)
            if v is None:
                assert v is getattr(cosmo, n)
            else:
                assert u.allclose(
                    v, getattr(cosmo, n), atol=1e-4 * getattr(v, "unit", 1)
                )

    # @pytest.mark.parametrize("z", valid_zs)  # TODO! recompute comparisons below
    def test_w(self, cosmo):
        """Test :meth:`astropy.cosmology.w0wzCDM.w`."""
        # super().test_w(cosmo, z)

        assert u.allclose(cosmo.w(1.0), -0.5)
        assert u.allclose(
            cosmo.w([0.0, 0.5, 1.0, 1.5, 2.3]), [-1.0, -0.75, -0.5, -0.25, 0.15]
        )

    def test_repr(self, cosmo_cls, cosmo):
        """Test method ``.__repr__()``."""
        super().test_repr(cosmo_cls, cosmo)

        expected = (
            'w0wzCDM(name="ABCMeta", H0=70.0 km / (Mpc s), Om0=0.27,'
            " Ode0=0.73, w0=-1.0, wz=0.5, Tcmb0=3.0 K, Neff=3.04,"
            " m_nu=[0. 0. 0.] eV, Ob0=0.03)"
        )
        assert repr(cosmo) == expected

    # ---------------------------------------------------------------

    @pytest.mark.parametrize("z", valid_zs)
    def test_Otot(self, cosmo, z):
        """Test :meth:`astropy.cosmology.w0wzCDM.Otot`.

        This is tested in the base class, but we need to override it here because
        this class is quite unstable.
        """
        super().test_Otot(cosmo, z)

    def test_Otot_overflow(self, cosmo):
        """Test :meth:`astropy.cosmology.w0wzCDM.Otot` for overflow."""
        with pytest.warns(RuntimeWarning, match="overflow encountered in exp"):
            cosmo.Otot(1e3)

    # ===============================================================
    # Usage Tests

    @pytest.mark.skipif(not HAS_SCIPY, reason="scipy is not installed")
    @pytest.mark.parametrize(
        ("args", "kwargs", "expected"),
        [
            (  # no relativistic species
                (75.0, 0.3, 0.6),
                {},
                [3051.68786716, 4756.17714818, 5822.38084257, 6562.70873734] * u.Mpc,
            ),
            (  # massless neutrinos
                (75.0, 0.25, 0.5),
                {"Tcmb0": 3.0, "Neff": 3, "m_nu": 0 * u.eV},
                [2997.8115653, 4686.45599916, 5764.54388557, 6524.17408738] * u.Mpc,
            ),
            (  # massive neutrinos
                (75.0, 0.25, 0.5),
                {"Tcmb0": 3.0, "Neff": 4, "m_nu": 5 * u.eV},
                [2676.73467639, 3940.57967585, 4686.90810278, 5191.54178243] * u.Mpc,
            ),
        ],
    )
    def test_comoving_distance_example(self, cosmo_cls, args, kwargs, expected):
        """Test :meth:`astropy.cosmology.LambdaCDM.comoving_distance`.

        These do not come from external codes -- they are just internal checks to make
        sure nothing changes if we muck with the distance calculators.
        """
        super().test_comoving_distance_example(
            cosmo_cls, args, {**COMOVING_DISTANCE_EXAMPLE_KWARGS, **kwargs}, expected
        )


class TestFlatw0wzCDM(FlatFLRWMixinTest, Testw0wzCDM):
    """Test :class:`astropy.cosmology.Flatw0wzCDM`."""

    def setup_class(self):
        """Setup for testing."""
        super().setup_class(self)
        self.cls = Flatw0wzCDM

    def test_repr(self, cosmo_cls, cosmo):
        """Test method ``.__repr__()``."""
        super().test_repr(cosmo_cls, cosmo)

        expected = (
            'Flatw0wzCDM(name="ABCMeta", H0=70.0 km / (Mpc s),'
            " Om0=0.27, w0=-1.0, wz=0.5, Tcmb0=3.0 K,"
            " Neff=3.04, m_nu=[0. 0. 0.] eV, Ob0=0.03)"
        )
        assert repr(cosmo) == expected

    # ---------------------------------------------------------------

    @pytest.mark.parametrize("z", valid_zs)
    def test_Otot(self, cosmo, z):
        """Test :meth:`astropy.cosmology.Flatw0wzCDM.Otot`.

        This is tested in the base class, but we need to override it here because
        this class is quite unstable.
        """
        super().test_Otot(cosmo, z)

    def test_Otot_overflow(self, cosmo):
        """Test :meth:`astropy.cosmology.Flatw0wzCDM.Otot` for NOT overflowing."""
        cosmo.Otot(1e5)

    # ---------------------------------------------------------------

    @pytest.mark.skipif(not HAS_SCIPY, reason="scipy is not installed")
    @pytest.mark.parametrize(
        ("args", "kwargs", "expected"),
        [
            (  # no relativistic species
                (75.0, 0.3),
                {},
                [3156.41804372, 4951.19475878, 6064.40591021, 6831.18710042] * u.Mpc,
            ),
            (  # massless neutrinos
                (75.0, 0.25),
                {"Tcmb0": 3.0, "Neff": 3, "m_nu": 0 * u.eV},
                [3268.38450997, 5205.96494068, 6419.75447923, 7257.77819438] * u.Mpc,
            ),
            (  # massive neutrinos
                (75.0, 0.25),
                {"Tcmb0": 3.0, "Neff": 4, "m_nu": 5 * u.eV},
                [2536.77159626, 3721.76294016, 4432.3526772, 4917.90352107] * u.Mpc,
            ),
        ],
    )
    def test_comoving_distance_example(self, cosmo_cls, args, kwargs, expected):
        """Test :meth:`astropy.cosmology.LambdaCDM.comoving_distance`.

        These do not come from external codes -- they are just internal checks to make
        sure nothing changes if we muck with the distance calculators.
        """
        super().test_comoving_distance_example(
            cosmo_cls, args, {**COMOVING_DISTANCE_EXAMPLE_KWARGS, **kwargs}, expected
        )


##############################################################################
# Miscellaneous
# TODO: these should be better integrated into the new test framework


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_de_densityscale():
    cosmo = w0wzCDM(H0=70, Om0=0.3, Ode0=0.50, w0=-1, wz=0.5)

    z = np.array([0.1, 0.2, 0.5, 1.5, 2.5])
    assert u.allclose(
        cosmo.de_density_scale(z),
        [1.00705953, 1.02687239, 1.15234885, 2.40022841, 6.49384982],
        rtol=1e-4,
    )

    assert u.allclose(cosmo.de_density_scale(3), cosmo.de_density_scale(3.0), rtol=1e-7)
    assert u.allclose(
        cosmo.de_density_scale([1, 2, 3]),
        cosmo.de_density_scale([1.0, 2.0, 3.0]),
        rtol=1e-7,
    )

    # Flat tests
    cosmo = w0wzCDM(H0=70, Om0=0.3, Ode0=0.7, w0=-1, wz=0.5)
    flatcosmo = Flatw0wzCDM(H0=70, Om0=0.3, w0=-1, wz=0.5)

    assert u.allclose(
        cosmo.de_density_scale(z), flatcosmo.de_density_scale(z), rtol=1e-4
    )
