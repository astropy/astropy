# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.flrw.lambdacdm`."""

##############################################################################
# IMPORTS

import numpy as np

# THIRD PARTY
import pytest

# LOCAL
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM, LambdaCDM
from astropy.cosmology.flrw.lambdacdm import ellipkinc, hyp2f1
from astropy.cosmology.tests.helper import get_redshift_methods
from astropy.cosmology.tests.test_core import invalid_zs, valid_zs
from astropy.utils.compat.optional_deps import HAS_SCIPY

from .test_base import FlatFLRWMixinTest, FLRWTest

##############################################################################
# TESTS
##############################################################################


@pytest.mark.skipif(HAS_SCIPY, reason="scipy is installed")
def test_optional_deps_functions():
    """Test stand-in functions when optional dependencies not installed."""
    with pytest.raises(ModuleNotFoundError, match="No module named 'scipy.special'"):
        ellipkinc()

    with pytest.raises(ModuleNotFoundError, match="No module named 'scipy.special'"):
        hyp2f1()


##############################################################################


class TestLambdaCDM(FLRWTest):
    """Test :class:`astropy.cosmology.LambdaCDM`."""

    def setup_class(self):
        """Setup for testing."""
        super().setup_class(self)
        self.cls = LambdaCDM

    # ===============================================================
    # Method & Attribute Tests

    _FLRW_redshift_methods = get_redshift_methods(
        LambdaCDM, include_private=True, include_z2=False
    ) - {"_dS_age"}
    # `_dS_age` is removed because it doesn't strictly rely on the value of `z`,
    # so any input that doesn't trip up ``np.shape`` is "valid"

    @pytest.mark.skipif(not HAS_SCIPY, reason="scipy is not installed")
    @pytest.mark.parametrize("z, exc", invalid_zs)
    @pytest.mark.parametrize("method", _FLRW_redshift_methods)
    def test_redshift_method_bad_input(self, cosmo, method, z, exc):
        """Test all the redshift methods for bad input."""
        super().test_redshift_method_bad_input(cosmo, method, z, exc)

    @pytest.mark.parametrize("z", valid_zs)
    def test_w(self, cosmo, z):
        """Test :meth:`astropy.cosmology.LambdaCDM.w`."""
        super().test_w(cosmo, z)

        w = cosmo.w(z)
        assert u.allclose(w, -1.0)

    def test_repr(self, cosmo_cls, cosmo):
        """Test method ``.__repr__()``."""
        super().test_repr(cosmo_cls, cosmo)

        expected = (
            'LambdaCDM(name="ABCMeta", H0=70.0 km / (Mpc s), Om0=0.27,'
            " Ode0=0.73, Tcmb0=3.0 K, Neff=3.04, m_nu=[0. 0. 0.] eV,"
            " Ob0=0.03)"
        )
        assert repr(cosmo) == expected


# -----------------------------------------------------------------------------


class TestFlatLambdaCDM(FlatFLRWMixinTest, TestLambdaCDM):
    """Test :class:`astropy.cosmology.FlatLambdaCDM`."""

    def setup_class(self):
        """Setup for testing."""
        super().setup_class(self)
        self.cls = FlatLambdaCDM

    @pytest.mark.skipif(not HAS_SCIPY, reason="scipy is not installed")
    @pytest.mark.parametrize("z, exc", invalid_zs)
    @pytest.mark.parametrize("method", TestLambdaCDM._FLRW_redshift_methods - {"Otot"})
    def test_redshift_method_bad_input(self, cosmo, method, z, exc):
        """Test all the redshift methods for bad input."""
        super().test_redshift_method_bad_input(cosmo, method, z, exc)

    # ===============================================================
    # Method & Attribute Tests

    def test_repr(self, cosmo_cls, cosmo):
        """Test method ``.__repr__()``."""
        super().test_repr(cosmo_cls, cosmo)

        expected = (
            'FlatLambdaCDM(name="ABCMeta", H0=70.0 km / (Mpc s),'
            " Om0=0.27, Tcmb0=3.0 K, Neff=3.04, m_nu=[0. 0. 0.] eV,"
            " Ob0=0.03)"
        )
        assert repr(cosmo) == expected


##############################################################################
# Comparison to Other Codes


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy.")
def test_flat_z1():
    """Test a flat cosmology at z=1 against several other on-line calculators.

    Test values were taken from the following web cosmology calculators on
    2012-02-11:

    Wright: http://www.astro.ucla.edu/~wright/CosmoCalc.html
            (https://ui.adsabs.harvard.edu/abs/2006PASP..118.1711W)
    Kempner: http://www.kempner.net/cosmic.php
    iCosmos: http://www.icosmos.co.uk/index.html
    """
    cosmo = FlatLambdaCDM(H0=70, Om0=0.27, Tcmb0=0.0)

    # The order of values below is Wright, Kempner, iCosmos'
    assert u.allclose(
        cosmo.comoving_distance(1), [3364.5, 3364.8, 3364.7988] * u.Mpc, rtol=1e-4
    )
    assert u.allclose(
        cosmo.angular_diameter_distance(1),
        [1682.3, 1682.4, 1682.3994] * u.Mpc,
        rtol=1e-4,
    )
    assert u.allclose(
        cosmo.luminosity_distance(1), [6729.2, 6729.6, 6729.5976] * u.Mpc, rtol=1e-4
    )
    assert u.allclose(
        cosmo.lookback_time(1), [7.841, 7.84178, 7.843] * u.Gyr, rtol=1e-3
    )
    assert u.allclose(
        cosmo.lookback_distance(1), [2404.0, 2404.24, 2404.4] * u.Mpc, rtol=1e-3
    )


##############################################################################
# Regression Tests


SPECIALIZED_COMOVING_DISTANCE_COSMOLOGIES = [
    FlatLambdaCDM(H0=70, Om0=0.0, Tcmb0=0.0),  # de Sitter
    FlatLambdaCDM(H0=70, Om0=1.0, Tcmb0=0.0),  # Einstein - de Sitter
    FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=0.0),  # Hypergeometric
    LambdaCDM(H0=70, Om0=0.3, Ode0=0.6, Tcmb0=0.0),  # Elliptic
]

ITERABLE_REDSHIFTS = [
    (0, 1, 2, 3, 4),  # tuple
    [0, 1, 2, 3, 4],  # list
    np.array([0, 1, 2, 3, 4]),  # array
]


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
@pytest.mark.parametrize("cosmo", SPECIALIZED_COMOVING_DISTANCE_COSMOLOGIES)
@pytest.mark.parametrize("z", ITERABLE_REDSHIFTS)
def test_comoving_distance_iterable_argument(cosmo, z):
    """
    Regression test for #10980
    Test that specialized comoving distance methods handle iterable arguments.
    """

    assert u.allclose(
        cosmo.comoving_distance(z), cosmo._integral_comoving_distance_z1z2(0.0, z)
    )


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
@pytest.mark.parametrize("cosmo", SPECIALIZED_COMOVING_DISTANCE_COSMOLOGIES)
def test_comoving_distance_broadcast(cosmo):
    """
    Regression test for #10980
    Test that specialized comoving distance methods broadcast array arguments.
    """

    z1 = np.zeros((2, 5))
    z2 = np.ones((3, 1, 5))
    z3 = np.ones((7, 5))
    output_shape = np.broadcast(z1, z2).shape

    # Check compatible array arguments return an array with the correct shape
    assert cosmo._comoving_distance_z1z2(z1, z2).shape == output_shape

    # Check incompatible array arguments raise an error
    with pytest.raises(ValueError, match="z1 and z2 have different shapes"):
        cosmo._comoving_distance_z1z2(z1, z3)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_elliptic_comoving_distance_z1z2():
    """Regression test for #8388."""
    cosmo = LambdaCDM(70.0, 2.3, 0.05, Tcmb0=0)
    z = 0.2
    assert u.allclose(
        cosmo.comoving_distance(z), cosmo._integral_comoving_distance_z1z2(0.0, z)
    )
    assert u.allclose(
        cosmo._elliptic_comoving_distance_z1z2(0.0, z),
        cosmo._integral_comoving_distance_z1z2(0.0, z),
    )


##############################################################################
# Miscellaneous
# TODO: these should be better integrated into the new test framework


def test_xtfuncs():
    """Test of absorption and lookback integrand"""
    cosmo = LambdaCDM(70, 0.3, 0.5, Tcmb0=2.725)
    z = np.array([2.0, 3.2])
    assert u.allclose(cosmo.lookback_time_integrand(3), 0.052218976654969378, rtol=1e-4)
    assert u.allclose(
        cosmo.lookback_time_integrand(z), [0.10333179, 0.04644541], rtol=1e-4
    )
    assert u.allclose(cosmo.abs_distance_integrand(3), 3.3420145059180402, rtol=1e-4)
    assert u.allclose(
        cosmo.abs_distance_integrand(z), [2.7899584, 3.44104758], rtol=1e-4
    )


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_matter():
    # Test non-relativistic matter evolution
    tcos = FlatLambdaCDM(70.0, 0.3, Ob0=0.045)

    assert u.allclose(tcos.Om(0), 0.3)
    assert u.allclose(tcos.Ob(0), 0.045)

    z = np.array([0.0, 0.5, 1.0, 2.0])
    assert u.allclose(tcos.Om(z), [0.3, 0.59124088, 0.77419355, 0.92045455], rtol=1e-4)
    assert u.allclose(
        tcos.Ob(z), [0.045, 0.08868613, 0.11612903, 0.13806818], rtol=1e-4
    )
    assert u.allclose(
        tcos.Odm(z), [0.255, 0.50255474, 0.65806452, 0.78238636], rtol=1e-4
    )
    # Consistency of dark and baryonic matter evolution with all
    # non-relativistic matter
    assert u.allclose(tcos.Ob(z) + tcos.Odm(z), tcos.Om(z))
