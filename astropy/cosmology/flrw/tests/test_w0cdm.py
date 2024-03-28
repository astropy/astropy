# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.flrw.w0cdm`."""

import numpy as np
import pytest

import astropy.units as u
from astropy.cosmology import FlatwCDM, wCDM
from astropy.cosmology.parameter import Parameter
from astropy.cosmology.tests.test_core import ParameterTestMixin, valid_zs
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.compat.optional_deps import HAS_SCIPY

from .conftest import filter_keys_from_items
from .test_base import FlatFLRWMixinTest, FLRWTest

##############################################################################
# TESTS
##############################################################################


class Parameterw0TestMixin(ParameterTestMixin):
    """Tests for `astropy.cosmology.Parameter` w0 on a Cosmology.

    w0 is a descriptor, which are tested by mixin, here with ``TestFLRW``.
    These tests expect dicts ``_cls_args`` and ``cls_kwargs`` which give the
    args and kwargs for the cosmology class, respectively. See ``TestFLRW``.
    """

    def test_w0(self, cosmo_cls, cosmo):
        """Test Parameter ``w0``."""
        # on the class
        w0 = cosmo_cls.parameters["w0"]
        assert isinstance(w0, Parameter)
        assert "Dark energy equation of state" in w0.__doc__
        assert w0.unit is None
        assert w0.default == -1.0

        # on the instance
        assert cosmo.w0 is cosmo._w0
        assert cosmo.w0 == self.cls_kwargs["w0"]

    def test_init_w0(self, cosmo_cls, ba):
        """Test initialization for values of ``w0``."""
        # test that it works with units
        ba.arguments["w0"] = ba.arguments["w0"] << u.one  # ensure units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.w0 == ba.arguments["w0"]

        # also without units
        ba.arguments["w0"] = ba.arguments["w0"].value  # strip units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.w0 == ba.arguments["w0"]

        # must be dimensionless
        ba.arguments["w0"] = 10 * u.km
        with pytest.raises(TypeError):
            cosmo_cls(*ba.args, **ba.kwargs)


class TestwCDM(FLRWTest, Parameterw0TestMixin):
    """Test :class:`astropy.cosmology.wCDM`."""

    def setup_class(self):
        """Setup for testing."""
        super().setup_class(self)

        self.cls = wCDM
        self.cls_kwargs.update(w0=-0.5)

    # ===============================================================
    # Method & Attribute Tests

    def test_clone_change_param(self, cosmo):
        """Test method ``.clone()`` changing a(many) Parameter(s)."""
        super().test_clone_change_param(cosmo)

        # `w` params
        c = cosmo.clone(w0=0.1)
        assert c.w0 == 0.1
        for n, v in filter_keys_from_items(c.parameters, ("w0",)):
            v_expect = getattr(cosmo, n)
            if v is None:
                assert v is v_expect
            else:
                assert_quantity_allclose(v, v_expect, atol=1e-4 * getattr(v, "unit", 1))

    @pytest.mark.parametrize("z", valid_zs)
    def test_w(self, cosmo, z):
        """Test :meth:`astropy.cosmology.wCDM.w`."""
        super().test_w(cosmo, z)

        w = cosmo.w(z)
        assert u.allclose(w, self.cls_kwargs["w0"])

    def test_repr(self, cosmo_cls, cosmo):
        """Test method ``.__repr__()``."""
        assert repr(cosmo) == (
            "wCDM(name='ABCMeta', H0=<Quantity 70. km / (Mpc s)>, Om0=0.27, "
            "Ode0=0.73, Tcmb0=<Quantity 3. K>, Neff=3.04, "
            "m_nu=<Quantity [0., 0., 0.] eV>, Ob0=0.03, w0=-0.5)"
        )

    # ===============================================================
    # Usage Tests

    @pytest.mark.skipif(not HAS_SCIPY, reason="scipy is not installed")
    @pytest.mark.parametrize(
        ("args", "kwargs", "expected"),
        [
            (  # no relativistic species
                (75.0, 0.25, 0.4),
                {"w0": -0.9, "Tcmb0": 0.0},
                [2849.6163356, 4428.71661565, 5450.97862778, 6179.37072324] * u.Mpc,
            ),
            (  # massless neutrinos
                (75.0, 0.25, 0.4),
                {"w0": -1.1, "Tcmb0": 3.0, "Neff": 3, "m_nu": u.Quantity(0.0, u.eV)},
                [2904.35580229, 4511.11471267, 5543.43643353, 6275.9206788] * u.Mpc,
            ),
            (  # massive neutrinos
                (75.0, 0.25, 0.4),
                {"w0": -0.9, "Tcmb0": 3.0, "Neff": 3, "m_nu": u.Quantity(10.0, u.eV)},
                [2473.32522734, 3581.54519631, 4232.41674426, 4671.83818117] * u.Mpc,
            ),
        ],
    )
    def test_comoving_distance_example(self, cosmo_cls, args, kwargs, expected):
        """Test :meth:`astropy.cosmology.LambdaCDM.comoving_distance`.

        These do not come from external codes -- they are just internal checks to make
        sure nothing changes if we muck with the distance calculators.
        """
        super().test_comoving_distance_example(cosmo_cls, args, kwargs, expected)


# -----------------------------------------------------------------------------


class TestFlatwCDM(FlatFLRWMixinTest, TestwCDM):
    """Test :class:`astropy.cosmology.FlatwCDM`."""

    def setup_class(self):
        """Setup for testing."""
        super().setup_class(self)
        self.cls = FlatwCDM
        self.cls_kwargs.update(w0=-0.5)

    def test_repr(self, cosmo_cls, cosmo):
        """Test method ``.__repr__()``."""
        super().test_repr(cosmo_cls, cosmo)

        assert repr(cosmo) == (
            "FlatwCDM(name='ABCMeta', H0=<Quantity 70. km / (Mpc s)>, Om0=0.27, "
            "Tcmb0=<Quantity 3. K>, Neff=3.04, m_nu=<Quantity [0., 0., 0.] eV>, "
            "Ob0=0.03, w0=-0.5)"
        )

    # ===============================================================
    # Usage Tests

    @pytest.mark.skipif(not HAS_SCIPY, reason="scipy is not installed")
    @pytest.mark.parametrize(
        ("args", "kwargs", "expected"),
        [
            (  # no relativistic species
                (75.0, 0.25),
                {"w0": -1.05, "Tcmb0": 0.0},
                [3216.8296894, 5117.2097601, 6317.05995437, 7149.68648536] * u.Mpc,
            ),
            (  # massless neutrinos
                (75.0, 0.25),
                {"w0": -0.95, "Tcmb0": 3.0, "Neff": 3, "m_nu": u.Quantity(0.0, u.eV)},
                [3143.56537758, 5000.32196494, 6184.11444601, 7009.80166062] * u.Mpc,
            ),
            (  # massive neutrinos
                (75.0, 0.25),
                {"w0": -0.9, "Tcmb0": 3.0, "Neff": 3, "m_nu": u.Quantity(10.0, u.eV)},
                [2337.76035371, 3372.1971387, 3988.71362289, 4409.40817174] * u.Mpc,
            ),
        ],
    )
    def test_comoving_distance_example(self, cosmo_cls, args, kwargs, expected):
        """Test :meth:`astropy.cosmology.LambdaCDM.comoving_distance`.

        These do not come from external codes -- they are just internal checks to make
        sure nothing changes if we muck with the distance calculators.
        """
        super().test_comoving_distance_example(cosmo_cls, args, kwargs, expected)


##############################################################################
# Miscellaneous
# TODO: these should be better integrated into the new test framework


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_de_densityscale():
    cosmo = wCDM(H0=70, Om0=0.3, Ode0=0.60, w0=-0.5)

    z = np.array([0.1, 0.2, 0.5, 1.5, 2.5])
    assert u.allclose(
        cosmo.de_density_scale(z),
        [1.15369, 1.31453, 1.83712, 3.95285, 6.5479],
        rtol=1e-4,
    )

    assert u.allclose(cosmo.de_density_scale(3), cosmo.de_density_scale(3.0), rtol=1e-7)
    assert u.allclose(
        cosmo.de_density_scale([1, 2, 3]),
        cosmo.de_density_scale([1.0, 2.0, 3.0]),
        rtol=1e-7,
    )
