# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.flrw.wpwazpcdm`."""

##############################################################################
# IMPORTS

# STDLIB

# THIRD PARTY
import pytest

# LOCAL
import astropy.cosmology.units as cu
import astropy.units as u
from astropy.cosmology import wpwaCDM
from astropy.cosmology.parameter import Parameter
from astropy.cosmology.tests.test_core import ParameterTestMixin

from .test_base import FLRWSubclassTest
from .test_w0wacdm import ParameterwaTestMixin

##############################################################################
# TESTS
##############################################################################


class ParameterwpTestMixin(ParameterTestMixin):
    """Tests for `astropy.cosmology.Parameter` wp on a Cosmology.

    wp is a descriptor, which are tested by mixin, here with ``TestFLRW``.
    These tests expect dicts ``_cls_args`` and ``cls_kwargs`` which give the
    args and kwargs for the cosmology class, respectively. See ``TestFLRW``.
    """

    def test_wp(self, cosmo_cls, cosmo):
        """Test Parameter ``wp``."""
        # on the class
        assert isinstance(cosmo_cls.wp, Parameter)
        assert "at the pivot" in cosmo_cls.wp.__doc__
        assert cosmo_cls.wp.unit is None

        # on the instance
        assert cosmo.wp is cosmo._wp
        assert cosmo.wp == self.cls_kwargs["wp"]

    def test_init_wp(self, cosmo_cls, ba):
        """Test initialization for values of ``wp``."""
        # test that it works with units
        ba.arguments["wp"] = ba.arguments["wp"] << u.one  # ensure units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.wp == ba.arguments["wp"]

        # also without units
        ba.arguments["wp"] = ba.arguments["wp"].value  # strip units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.wp == ba.arguments["wp"]

        # must be dimensionless
        ba.arguments["wp"] = 10 * u.km
        with pytest.raises(TypeError):
            cosmo_cls(*ba.args, **ba.kwargs)


class ParameterzpTestMixin(ParameterTestMixin):
    """Tests for `astropy.cosmology.Parameter` zp on a Cosmology.

    zp is a descriptor, which are tested by mixin, here with ``TestFLRW``.
    These tests expect dicts ``_cls_args`` and ``cls_kwargs`` which give the
    args and kwargs for the cosmology class, respectively. See ``TestFLRW``.
    """

    def test_zp(self, cosmo_cls, cosmo):
        """Test Parameter ``zp``."""
        # on the class
        assert isinstance(cosmo_cls.zp, Parameter)
        assert "pivot redshift" in cosmo_cls.zp.__doc__
        assert cosmo_cls.zp.unit == cu.redshift

        # on the instance
        assert cosmo.zp is cosmo._zp
        assert cosmo.zp == self.cls_kwargs["zp"] << cu.redshift

    def test_init_zp(self, cosmo_cls, ba):
        """Test initialization for values of ``zp``."""
        # test that it works with units
        ba.arguments["zp"] = ba.arguments["zp"] << u.one  # ensure units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.zp == ba.arguments["zp"]

        # also without units
        ba.arguments["zp"] = ba.arguments["zp"].value  # strip units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.zp.value == ba.arguments["zp"]

        # must be dimensionless
        ba.arguments["zp"] = 10 * u.km
        with pytest.raises(u.UnitConversionError):
            cosmo_cls(*ba.args, **ba.kwargs)


class TestwpwaCDM(FLRWSubclassTest,
                  ParameterwpTestMixin, ParameterwaTestMixin, ParameterzpTestMixin):
    """Test :class:`astropy.cosmology.wpwaCDM`."""

    def setup_class(self):
        """Setup for testing."""
        super().setup_class(self)
        self.cls = wpwaCDM
        self.cls_kwargs.update(wp=-0.9, wa=0.2, zp=0.5)

    # ===============================================================
    # Method & Attribute Tests

    def test_clone_change_param(self, cosmo):
        """Test method ``.clone()`` changing a(many) Parameter(s)."""
        super().test_clone_change_param(cosmo)

        # `w` params
        c = cosmo.clone(wp=0.1, wa=0.2, zp=14)
        assert c.wp == 0.1
        assert c.wa == 0.2
        assert c.zp == 14
        for n in (set(cosmo.__parameters__) - {"wp", "wa", "zp"}):
            v = getattr(c, n)
            if v is None:
                assert v is getattr(cosmo, n)
            else:
                assert u.allclose(v, getattr(cosmo, n), atol=1e-4 * getattr(v, "unit", 1))

    # @pytest.mark.parametrize("z", valid_zs)  # TODO! recompute comparisons below
    def test_w(self, cosmo):
        """Test :meth:`astropy.cosmology.wpwaCDM.w`."""
        # super().test_w(cosmo, z)

        assert u.allclose(cosmo.w(0.5), -0.9)
        assert u.allclose(cosmo.w([0.1, 0.2, 0.5, 1.5, 2.5, 11.5]),
                          [-0.94848485, -0.93333333, -0.9, -0.84666667,
                           -0.82380952, -0.78266667])

    def test_repr(self, cosmo_cls, cosmo):
        """Test method ``.__repr__()``."""
        super().test_repr(cosmo_cls, cosmo)

        expected = ("wpwaCDM(name=\"ABCMeta\", H0=70.0 km / (Mpc s), Om0=0.27,"
                    " Ode0=0.73, wp=-0.9, wa=0.2, zp=0.5 redshift, Tcmb0=3.0 K,"
                    " Neff=3.04, m_nu=[0. 0. 0.] eV, Ob0=0.03)")
        assert repr(cosmo) == expected
