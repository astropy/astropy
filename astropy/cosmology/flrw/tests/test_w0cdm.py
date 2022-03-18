# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.flrw.w0cdm`."""

##############################################################################
# IMPORTS

# STDLIB

# THIRD PARTY
import pytest

# LOCAL
import astropy.units as u
from astropy.cosmology import FlatwCDM, wCDM
from astropy.cosmology.parameter import Parameter
from astropy.cosmology.tests.test_core import ParameterTestMixin, valid_zs

from .test_base import FlatFLRWMixinTest, FLRWSubclassTest

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
        assert isinstance(cosmo_cls.w0, Parameter)
        assert "Dark energy equation of state" in cosmo_cls.w0.__doc__
        assert cosmo_cls.w0.unit is None

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


class TestwCDM(FLRWSubclassTest, Parameterw0TestMixin):
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
        for n in (set(cosmo.__parameters__) - {"w0"}):
            v = getattr(c, n)
            if v is None:
                assert v is getattr(cosmo, n)
            else:
                assert u.allclose(v, getattr(cosmo, n), atol=1e-4 * getattr(v, "unit", 1))

    @pytest.mark.parametrize("z", valid_zs)
    def test_w(self, cosmo, z):
        """Test :meth:`astropy.cosmology.wCDM.w`."""
        super().test_w(cosmo, z)

        w = cosmo.w(z)
        assert u.allclose(w, self.cls_kwargs["w0"])

    def test_repr(self, cosmo_cls, cosmo):
        """Test method ``.__repr__()``."""
        super().test_repr(cosmo_cls, cosmo)

        expected = ("wCDM(name=\"ABCMeta\", H0=70.0 km / (Mpc s), Om0=0.27,"
                    " Ode0=0.73, w0=-0.5, Tcmb0=3.0 K, Neff=3.04,"
                    " m_nu=[0. 0. 0.] eV, Ob0=0.03)")
        assert repr(cosmo) == expected


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

        expected = ("FlatwCDM(name=\"ABCMeta\", H0=70.0 km / (Mpc s), Om0=0.27,"
                    " w0=-0.5, Tcmb0=3.0 K, Neff=3.04, m_nu=[0. 0. 0.] eV,"
                    " Ob0=0.03)")
        assert repr(cosmo) == expected
