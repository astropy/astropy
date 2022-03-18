# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.flrw.w0wzcdm`."""

##############################################################################
# IMPORTS

# STDLIB

# THIRD PARTY
import pytest

# LOCAL
import astropy.units as u
from astropy.cosmology import w0wzCDM
from astropy.cosmology.parameter import Parameter
from astropy.cosmology.tests.test_core import ParameterTestMixin

from .test_base import FLRWSubclassTest
from .test_w0cdm import Parameterw0TestMixin

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


class Testw0wzCDM(FLRWSubclassTest,
                  Parameterw0TestMixin, ParameterwzTestMixin):
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
        for n in (set(cosmo.__parameters__) - {"w0", "wz"}):
            v = getattr(c, n)
            if v is None:
                assert v is getattr(cosmo, n)
            else:
                assert u.allclose(v, getattr(cosmo, n), atol=1e-4 * getattr(v, "unit", 1))

    # @pytest.mark.parametrize("z", valid_zs)  # TODO! recompute comparisons below
    def test_w(self, cosmo):
        """Test :meth:`astropy.cosmology.w0wzCDM.w`."""
        # super().test_w(cosmo, z)

        assert u.allclose(cosmo.w(1.0), -0.5)
        assert u.allclose(cosmo.w([0.0, 0.5, 1.0, 1.5, 2.3]),
                          [-1.0, -0.75, -0.5, -0.25, 0.15])

    def test_repr(self, cosmo_cls, cosmo):
        """Test method ``.__repr__()``."""
        super().test_repr(cosmo_cls, cosmo)

        expected = ("w0wzCDM(name=\"ABCMeta\", H0=70.0 km / (Mpc s), Om0=0.27,"
                    " Ode0=0.73, w0=-1.0, wz=0.5, Tcmb0=3.0 K, Neff=3.04,"
                    " m_nu=[0. 0. 0.] eV, Ob0=0.03)")
        assert repr(cosmo) == expected
