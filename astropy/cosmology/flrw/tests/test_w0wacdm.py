# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.flrw.w0wacdm`."""

##############################################################################
# IMPORTS

# STDLIB

# THIRD PARTY
import pytest

# LOCAL
import astropy.units as u
from astropy.cosmology import Flatw0waCDM, w0waCDM
from astropy.cosmology.parameter import Parameter
from astropy.cosmology.tests.test_core import ParameterTestMixin

from .test_base import FlatFLRWMixinTest, FLRWSubclassTest
from .test_w0cdm import Parameterw0TestMixin

##############################################################################
# TESTS
##############################################################################


class ParameterwaTestMixin(ParameterTestMixin):
    """Tests for `astropy.cosmology.Parameter` wa on a Cosmology.

    wa is a descriptor, which are tested by mixin, here with ``TestFLRW``.
    These tests expect dicts ``_cls_args`` and ``cls_kwargs`` which give the
    args and kwargs for the cosmology class, respectively. See ``TestFLRW``.
    """

    def test_wa(self, cosmo_cls, cosmo):
        """Test Parameter ``wa``."""
        # on the class
        assert isinstance(cosmo_cls.wa, Parameter)
        assert "Negative derivative" in cosmo_cls.wa.__doc__
        assert cosmo_cls.wa.unit is None

        # on the instance
        assert cosmo.wa is cosmo._wa
        assert cosmo.wa == self.cls_kwargs["wa"]

    def test_init_wa(self, cosmo_cls, ba):
        """Test initialization for values of ``wa``."""
        # test that it works with units
        ba.arguments["wa"] = ba.arguments["wa"] << u.one  # ensure units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.wa == ba.arguments["wa"]

        # also without units
        ba.arguments["wa"] = ba.arguments["wa"].value  # strip units
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        assert cosmo.wa == ba.arguments["wa"]

        # must be dimensionless
        ba.arguments["wa"] = 10 * u.km
        with pytest.raises(TypeError):
            cosmo_cls(*ba.args, **ba.kwargs)


class Testw0waCDM(FLRWSubclassTest, Parameterw0TestMixin, ParameterwaTestMixin):
    """Test :class:`astropy.cosmology.w0waCDM`."""

    def setup_class(self):
        """Setup for testing."""
        super().setup_class(self)
        self.cls = w0waCDM
        self.cls_kwargs.update(w0=-1, wa=-0.5)

    # ===============================================================
    # Method & Attribute Tests

    def test_clone_change_param(self, cosmo):
        """Test method ``.clone()`` changing a(many) Parameter(s)."""
        super().test_clone_change_param(cosmo)

        # `w` params
        c = cosmo.clone(w0=0.1, wa=0.2)
        assert c.w0 == 0.1
        assert c.wa == 0.2
        for n in (set(cosmo.__parameters__) - {"w0", "wa"}):
            v = getattr(c, n)
            if v is None:
                assert v is getattr(cosmo, n)
            else:
                assert u.allclose(v, getattr(cosmo, n), atol=1e-4 * getattr(v, "unit", 1))

    # @pytest.mark.parametrize("z", valid_zs)  # TODO! recompute comparisons below
    def test_w(self, cosmo):
        """Test :meth:`astropy.cosmology.w0waCDM.w`."""
        # super().test_w(cosmo, z)

        assert u.allclose(cosmo.w(1.0), -1.25)
        assert u.allclose(cosmo.w([0.0, 0.5, 1.0, 1.5, 2.3]),
                          [-1, -1.16666667, -1.25, -1.3, -1.34848485])

    def test_repr(self, cosmo_cls, cosmo):
        """Test method ``.__repr__()``."""
        super().test_repr(cosmo_cls, cosmo)

        expected = ("w0waCDM(name=\"ABCMeta\", H0=70.0 km / (Mpc s), Om0=0.27,"
                    " Ode0=0.73, w0=-1.0, wa=-0.5, Tcmb0=3.0 K, Neff=3.04,"
                    " m_nu=[0. 0. 0.] eV, Ob0=0.03)")
        assert repr(cosmo) == expected


# -----------------------------------------------------------------------------


class TestFlatw0waCDM(FlatFLRWMixinTest, Testw0waCDM):
    """Test :class:`astropy.cosmology.Flatw0waCDM`."""

    def setup_class(self):
        """Setup for testing."""
        super().setup_class(self)
        self.cls = Flatw0waCDM
        self.cls_kwargs.update(w0=-1, wa=-0.5)

    def test_repr(self, cosmo_cls, cosmo):
        """Test method ``.__repr__()``."""
        super().test_repr(cosmo_cls, cosmo)

        expected = ("Flatw0waCDM(name=\"ABCMeta\", H0=70.0 km / (Mpc s),"
                    " Om0=0.27, w0=-1.0, wa=-0.5, Tcmb0=3.0 K, Neff=3.04,"
                    " m_nu=[0. 0. 0.] eV, Ob0=0.03)")
        assert repr(cosmo) == expected
