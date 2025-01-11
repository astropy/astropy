# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.flrw.lambdacdm`."""

import pathlib

import numpy as np
import pytest

import astropy.constants as const
import astropy.cosmology.units as cu
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM, LambdaCDM
from astropy.cosmology._src.flrw.lambdacdm import ellipkinc, hyp2f1
from astropy.cosmology._src.tests.helper import get_redshift_methods
from astropy.cosmology._src.tests.test_core import invalid_zs, valid_zs
from astropy.table import QTable
from astropy.utils.compat.optional_deps import HAS_SCIPY
from astropy.utils.exceptions import AstropyUserWarning

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
    @pytest.mark.parametrize("method", sorted(_FLRW_redshift_methods))
    def test_redshift_method_bad_input(self, cosmo, method, z, exc):
        """Test all the redshift methods for bad input."""
        super().test_redshift_method_bad_input(cosmo, method, z, exc)

    @pytest.mark.parametrize("z", valid_zs)
    def test_w(self, cosmo, z):
        """Test :meth:`astropy.cosmology.LambdaCDM.w`."""
        super().test_w(cosmo, z)

        w = cosmo.w(z)
        assert u.allclose(w, -1.0)

    def test_repr(self, cosmo):
        """Test method ``.__repr__()``."""
        assert repr(cosmo) == (
            "LambdaCDM(name='ABCMeta', H0=<Quantity 70. km / (Mpc s)>, Om0=0.27,"
            " Ode0=0.73, Tcmb0=<Quantity 3. K>, Neff=3.04,"
            " m_nu=<Quantity [0., 0., 0.] eV>, Ob0=0.03)"
        )

    def test_str(self, cosmo):
        """Test method ``.__str__()``."""
        assert str(cosmo) == (
            'LambdaCDM(name="ABCMeta", H0=70.0 km / (Mpc s), Om0=0.27, Ode0=0.73,'
            " Tcmb0=3.0 K, Neff=3.04, m_nu=[0. 0. 0.] eV, Ob0=0.03)"
        )

    @pytest.mark.skipif(not HAS_SCIPY, reason="scipy is not installed")
    @pytest.mark.parametrize(
        ("args", "kwargs", "expected"),
        [
            (  # no relativistic species
                (75.0, 0.25, 0.5),
                {"Tcmb0": 0.0},
                [2953.93001902, 4616.7134253, 5685.07765971, 6440.80611897] * u.Mpc,
            ),
            (  # massless neutrinos
                (75.0, 0.25, 0.6),
                {"Tcmb0": 3.0, "Neff": 3, "m_nu": u.Quantity(0.0, u.eV)},
                [3037.12620424, 4776.86236327, 5889.55164479, 6671.85418235] * u.Mpc,
            ),
            (  # massive neutrinos
                (75.0, 0.3, 0.4),
                {"Tcmb0": 3.0, "Neff": 3, "m_nu": u.Quantity(10.0, u.eV)},
                [2471.80626824, 3567.1902565, 4207.15995626, 4638.20476018] * u.Mpc,
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


class TestFlatLambdaCDM(FlatFLRWMixinTest, TestLambdaCDM):
    """Test :class:`astropy.cosmology.FlatLambdaCDM`."""

    def setup_class(self):
        """Setup for testing."""
        super().setup_class(self)
        self.cls = FlatLambdaCDM

    @pytest.mark.skipif(not HAS_SCIPY, reason="scipy is not installed")
    @pytest.mark.parametrize("z, exc", invalid_zs)
    @pytest.mark.parametrize(
        "method", sorted(TestLambdaCDM._FLRW_redshift_methods - {"Otot"})
    )
    def test_redshift_method_bad_input(self, cosmo, method, z, exc):
        """Test all the redshift methods for bad input."""
        super().test_redshift_method_bad_input(cosmo, method, z, exc)

    # ===============================================================
    # Method & Attribute Tests

    def test_repr(self, cosmo):
        """Test method ``.__repr__()``."""
        assert repr(cosmo) == (
            "FlatLambdaCDM(name='ABCMeta', H0=<Quantity 70. km / (Mpc s)>, Om0=0.27,"
            " Tcmb0=<Quantity 3. K>, Neff=3.04, m_nu=<Quantity [0., 0., 0.] eV>,"
            " Ob0=0.03)"
        )

    def test_str(self, cosmo):
        """Test method ``.__str__()``."""
        assert str(cosmo) == (
            'FlatLambdaCDM(name="ABCMeta", H0=70.0 km / (Mpc s), Om0=0.27, '
            "Tcmb0=3.0 K, Neff=3.04, m_nu=[0. 0. 0.] eV, Ob0=0.03)"
        )

    # ===============================================================
    # Usage Tests

    @pytest.mark.skipif(not HAS_SCIPY, reason="scipy is not installed")
    @pytest.mark.parametrize(
        ("args", "kwargs", "expected"),
        [
            (  # no relativistic species
                (75.0, 0.25),
                {"Tcmb0": 0.0},
                [3180.83488552, 5060.82054204, 6253.6721173, 7083.5374303] * u.Mpc,
            ),
            (  # massless neutrinos
                (75.0, 0.25),
                {"Tcmb0": 3.0, "Neff": 3, "m_nu": u.Quantity(0.0, u.eV)},
                [3180.42662867, 5059.60529655, 6251.62766102, 7080.71698117] * u.Mpc,
            ),
            (  # massive neutrinos
                (75.0, 0.25),
                {"Tcmb0": 3.0, "Neff": 3, "m_nu": u.Quantity(10.0, u.eV)},
                [2337.54183142, 3371.91131264, 3988.40711188, 4409.09346922] * u.Mpc,
            ),
            (  # work the scalar nu density functions
                (75.0, 0.25),
                {"Tcmb0": 3.0, "m_nu": u.Quantity([10.0, 0, 0], u.eV)},
                [2777.71589173, 4186.91111666, 5046.0300719, 5636.10397302] * u.Mpc,
            ),
            (  # work the scalar nu density functions
                (75.0, 0.25),
                {"Tcmb0": 3.0, "m_nu": u.Quantity([10.0, 5, 0], u.eV)},
                [2636.48149391, 3913.14102091, 4684.59108974, 5213.07557084] * u.Mpc,
            ),
            (  # work the scalar nu density functions
                (75.0, 0.25),
                {"Tcmb0": 3.0, "m_nu": u.Quantity([4.0, 5, 9], u.eV)},
                [2563.5093049, 3776.63362071, 4506.83448243, 5006.50158829] * u.Mpc,
            ),
            (  # work the scalar nu density functions
                (75.0, 0.25),
                {"Tcmb0": 3.0, "Neff": 4.2, "m_nu": u.Quantity([1.0, 4.0, 5, 9], u.eV)},
                [2525.58017482, 3706.87633298, 4416.58398847, 4901.96669755] * u.Mpc,
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


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_ogamma():
    """Tests the effects of changing the temperature of the CMB"""

    # Tested against Ned Wright's advanced cosmology calculator,
    # Sep 7 2012.  The accuracy of our comparison is limited by
    # how many digits it outputs, which limits our test to about
    # 0.2% accuracy.  The NWACC does not allow one
    # to change the number of nuetrino species, fixing that at 3.
    # Also, inspection of the NWACC code shows it uses inaccurate
    # constants at the 0.2% level (specifically, a_B),
    # so we shouldn't expect to match it that well. The integral is
    # also done rather crudely.  Therefore, we should not expect
    # the NWACC to be accurate to better than about 0.5%, which is
    # unfortunate, but reflects a problem with it rather than this code.
    # More accurate tests below using Mathematica
    z = np.array([1.0, 10.0, 500.0, 1000.0])

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=0, Neff=3)
    assert u.allclose(
        cosmo.angular_diameter_distance(z),
        [1651.9, 858.2, 26.855, 13.642] * u.Mpc,
        rtol=5e-4,
    )

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725, Neff=3)
    assert u.allclose(
        cosmo.angular_diameter_distance(z),
        [1651.8, 857.9, 26.767, 13.582] * u.Mpc,
        rtol=5e-4,
    )

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=4.0, Neff=3)
    assert u.allclose(
        cosmo.angular_diameter_distance(z),
        [1651.4, 856.6, 26.489, 13.405] * u.Mpc,
        rtol=5e-4,
    )

    # Next compare with doing the integral numerically in Mathematica,
    # which allows more precision in the test.  It is at least as
    # good as 0.01%, possibly better
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=0, Neff=3.04)
    assert u.allclose(
        cosmo.angular_diameter_distance(z),
        [1651.91, 858.205, 26.8586, 13.6469] * u.Mpc,
        rtol=1e-5,
    )

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725, Neff=3.04)
    assert u.allclose(
        cosmo.angular_diameter_distance(z),
        [1651.76, 857.817, 26.7688, 13.5841] * u.Mpc,
        rtol=1e-5,
    )

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=4.0, Neff=3.04)
    assert u.allclose(
        cosmo.angular_diameter_distance(z),
        [1651.21, 856.411, 26.4845, 13.4028] * u.Mpc,
        rtol=1e-5,
    )

    # Just to be really sure, we also do a version where the integral
    # is analytic, which is a Ode = 0 flat universe.  In this case
    # Integrate(1/E(x),{x,0,z}) = 2 ( sqrt((1+Or z)/(1+z)) - 1 )/(Or - 1)
    # Recall that c/H0 * Integrate(1/E) is FLRW.comoving_distance.
    Ogamma0h2 = 4 * 5.670373e-8 / 299792458.0**3 * 2.725**4 / 1.87837e-26
    Onu0h2 = Ogamma0h2 * 7.0 / 8.0 * (4.0 / 11.0) ** (4.0 / 3.0) * 3.04
    Or0 = (Ogamma0h2 + Onu0h2) / 0.7**2
    Om0 = 1.0 - Or0
    hubdis = (299792.458 / 70.0) * u.Mpc
    cosmo = FlatLambdaCDM(H0=70, Om0=Om0, Tcmb0=2.725, Neff=3.04)
    targvals = 2.0 * hubdis * (np.sqrt((1.0 + Or0 * z) / (1.0 + z)) - 1.0) / (Or0 - 1.0)
    assert u.allclose(cosmo.comoving_distance(z), targvals, rtol=1e-5)

    # And integers for z
    assert u.allclose(cosmo.comoving_distance(z.astype(int)), targvals, rtol=1e-5)

    # Try Tcmb0 = 4
    Or0 *= (4.0 / 2.725) ** 4
    Om0 = 1.0 - Or0
    cosmo = FlatLambdaCDM(H0=70, Om0=Om0, Tcmb0=4.0, Neff=3.04)
    targvals = 2.0 * hubdis * (np.sqrt((1.0 + Or0 * z) / (1.0 + z)) - 1.0) / (Or0 - 1.0)
    assert u.allclose(cosmo.comoving_distance(z), targvals, rtol=1e-5)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
@pytest.mark.parametrize(
    "file_name", ["cosmo_flat.ecsv", "cosmo_open.ecsv", "cosmo_closed.ecsv"]
)
def test_flat_open_closed_icosmo(file_name):
    """Test against the tabulated values generated from icosmo.org
    with three example cosmologies (flat, open and closed).
    """
    with u.add_enabled_units(cu):
        tbl = QTable.read(pathlib.Path(__file__).parent / "data" / file_name)
    cosmo = LambdaCDM(
        H0=100 * tbl.meta["h"], Om0=tbl.meta["Om"], Ode0=tbl.meta["Ol"], Tcmb0=0.0
    )
    assert u.allclose(cosmo.comoving_transverse_distance(tbl["redshift"]), tbl["dm"])
    assert u.allclose(cosmo.angular_diameter_distance(tbl["redshift"]), tbl["da"])
    assert u.allclose(cosmo.luminosity_distance(tbl["redshift"]), tbl["dl"])


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_comoving_transverse_distance_z1z2():
    tcos = FlatLambdaCDM(100, 0.3, Tcmb0=0.0)

    with pytest.raises(ValueError):  # test diff size z1, z2 fail
        tcos._comoving_transverse_distance_z1z2((1, 2), (3, 4, 5))

    # Tests that should actually work, target values computed with
    # http://www.astro.multivax.de:8000/phillip/angsiz_prog/README.HTML
    # Kayser, Helbig, and Schramm (Astron.Astrophys. 318 (1997) 680-686)
    assert u.allclose(
        tcos._comoving_transverse_distance_z1z2(1, 2), 1313.2232194828466 * u.Mpc
    )

    # In a flat universe comoving distance and comoving transverse
    # distance are identical
    z1 = 0, 0, 2, 0.5, 1
    z2 = 2, 1, 1, 2.5, 1.1

    assert u.allclose(
        tcos._comoving_distance_z1z2(z1, z2),
        tcos._comoving_transverse_distance_z1z2(z1, z2),
    )

    # Test Flat Universe with Omega_M > 1.  Rarely used, but perfectly valid.
    tcos = FlatLambdaCDM(100, 1.5, Tcmb0=0.0)
    results = (
        2202.72682564,
        1559.51679971,
        -643.21002593,
        1408.36365679,
        85.09286258,
    ) * u.Mpc

    assert u.allclose(tcos._comoving_transverse_distance_z1z2(z1, z2), results)

    # In a flat universe comoving distance and comoving transverse
    # distance are identical
    z1 = 0, 0, 2, 0.5, 1
    z2 = 2, 1, 1, 2.5, 1.1

    assert u.allclose(
        tcos._comoving_distance_z1z2(z1, z2),
        tcos._comoving_transverse_distance_z1z2(z1, z2),
    )
    # Test non-flat cases to avoid simply testing
    # comoving_distance_z1z2. Test array, array case.
    tcos = LambdaCDM(100, 0.3, 0.5, Tcmb0=0.0)
    results = (
        3535.931375645655,
        2226.430046551708,
        -1208.6817970036532,
        2595.567367601969,
        151.36592003406884,
    ) * u.Mpc

    assert u.allclose(tcos._comoving_transverse_distance_z1z2(z1, z2), results)

    # Test positive curvature with scalar, array combination.
    tcos = LambdaCDM(100, 1.0, 0.2, Tcmb0=0.0)
    z1 = 0.1
    z2 = 0, 0.1, 0.2, 0.5, 1.1, 2
    results = (
        -281.31602666724865,
        0.0,
        248.58093707820436,
        843.9331377460543,
        1618.6104987686672,
        2287.5626543279927,
    ) * u.Mpc

    assert u.allclose(tcos._comoving_transverse_distance_z1z2(z1, z2), results)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_angular_diameter_distance_z1z2():
    tcos = FlatLambdaCDM(70.4, 0.272, Tcmb0=0.0)

    with pytest.raises(ValueError):  # test diff size z1, z2 fail
        tcos.angular_diameter_distance_z1z2([1, 2], [3, 4, 5])

    # Tests that should actually work, target values computed with
    # http://www.astro.multivax.de:8000/phillip/angsiz_prog/README.HTML
    # Kayser, Helbig, and Schramm (Astron.Astrophys. 318 (1997) 680-686)
    assert u.allclose(
        tcos.angular_diameter_distance_z1z2(1, 2), 646.22968662822018 * u.Mpc
    )

    z1 = 2  # Separate test for z2<z1, returns negative value with warning
    z2 = 1
    results = -969.34452994 * u.Mpc
    with pytest.warns(AstropyUserWarning, match="less than first redshift"):
        assert u.allclose(tcos.angular_diameter_distance_z1z2(z1, z2), results)

    z1 = 0, 0, 0.5, 1
    z2 = 2, 1, 2.5, 1.1
    results = (
        1760.0628637762106,
        1670.7497657219858,
        1159.0970895962193,
        115.72768186186921,
    ) * u.Mpc

    assert u.allclose(tcos.angular_diameter_distance_z1z2(z1, z2), results)

    z1 = 0.1
    z2 = 0.1, 0.2, 0.5, 1.1, 2
    results = (0.0, 332.09893173, 986.35635069, 1508.37010062, 1621.07937976) * u.Mpc
    assert u.allclose(tcos.angular_diameter_distance_z1z2(0.1, z2), results)

    # Non-flat (positive Ok0) test
    tcos = LambdaCDM(H0=70.4, Om0=0.2, Ode0=0.5, Tcmb0=0.0)
    assert u.allclose(
        tcos.angular_diameter_distance_z1z2(1, 2), 620.1175337852428 * u.Mpc
    )
    # Non-flat (negative Ok0) test
    tcos = LambdaCDM(H0=100, Om0=2, Ode0=1, Tcmb0=0.0)
    assert u.allclose(
        tcos.angular_diameter_distance_z1z2(1, 2), 228.42914659246014 * u.Mpc
    )


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_massivenu_density():
    # Testing neutrino density calculation

    # Simple test cosmology, where we compare rho_nu and rho_gamma
    # against the exact formula (eq 24/25 of Komatsu et al. 2011)
    # computed using Mathematica.  The approximation we use for f(y)
    # is only good to ~ 0.5% (with some redshift dependence), so that's
    # what we test to.
    ztest = np.array([0.0, 1.0, 2.0, 10.0, 1000.0])
    nuprefac = 7.0 / 8.0 * (4.0 / 11.0) ** (4.0 / 3.0)
    #  First try 3 massive neutrinos, all 100 eV -- note this is a universe
    #  seriously dominated by neutrinos!
    tcos = FlatLambdaCDM(75.0, 0.25, Tcmb0=3.0, Neff=3, m_nu=u.Quantity(100.0, u.eV))
    assert tcos.has_massive_nu
    assert tcos.Neff == 3
    nurel_exp = (
        nuprefac * tcos.Neff * np.array([171969, 85984.5, 57323, 15633.5, 171.801])
    )
    assert u.allclose(tcos.nu_relative_density(ztest), nurel_exp, rtol=5e-3)
    assert u.allclose(tcos.efunc([0.0, 1.0]), [1.0, 7.46144727668], rtol=5e-3)

    # Next, slightly less massive
    tcos = FlatLambdaCDM(75.0, 0.25, Tcmb0=3.0, Neff=3, m_nu=u.Quantity(0.25, u.eV))
    nurel_exp = (
        nuprefac * tcos.Neff * np.array([429.924, 214.964, 143.312, 39.1005, 1.11086])
    )
    assert u.allclose(tcos.nu_relative_density(ztest), nurel_exp, rtol=5e-3)

    # For this one also test Onu directly
    onu_exp = np.array([0.01890217, 0.05244681, 0.0638236, 0.06999286, 0.1344951])
    assert u.allclose(tcos.Onu(ztest), onu_exp, rtol=5e-3)

    # And fairly light
    tcos = FlatLambdaCDM(80.0, 0.30, Tcmb0=3.0, Neff=3, m_nu=u.Quantity(0.01, u.eV))

    nurel_exp = (
        nuprefac * tcos.Neff * np.array([17.2347, 8.67345, 5.84348, 1.90671, 1.00021])
    )
    assert u.allclose(tcos.nu_relative_density(ztest), nurel_exp, rtol=5e-3)
    onu_exp = np.array([0.00066599, 0.00172677, 0.0020732, 0.00268404, 0.0978313])
    assert u.allclose(tcos.Onu(ztest), onu_exp, rtol=5e-3)
    assert u.allclose(tcos.efunc([1.0, 2.0]), [1.76225893, 2.97022048], rtol=1e-4)
    assert u.allclose(tcos.inv_efunc([1.0, 2.0]), [0.5674535, 0.33667534], rtol=1e-4)

    # Now a mixture of neutrino masses, with non-integer Neff
    tcos = FlatLambdaCDM(
        80.0, 0.30, Tcmb0=3.0, Neff=3.04, m_nu=u.Quantity([0.0, 0.01, 0.25], u.eV)
    )
    nurel_exp = (
        nuprefac
        * tcos.Neff
        * np.array([149.386233, 74.87915, 50.0518, 14.002403, 1.03702333])
    )
    assert u.allclose(tcos.nu_relative_density(ztest), nurel_exp, rtol=5e-3)
    onu_exp = np.array([0.00584959, 0.01493142, 0.01772291, 0.01963451, 0.10227728])
    assert u.allclose(tcos.Onu(ztest), onu_exp, rtol=5e-3)

    # Integer redshifts
    ztest = ztest.astype(int)
    assert u.allclose(tcos.nu_relative_density(ztest), nurel_exp, rtol=5e-3)
    assert u.allclose(tcos.Onu(ztest), onu_exp, rtol=5e-3)


##############################################################################
# Miscellaneous
# TODO: these should be better integrated into the new test framework


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_units():
    """Test if the right units are being returned"""
    cosmo = FlatLambdaCDM(H0=70, Om0=0.27, Tcmb0=2.0)

    assert cosmo.comoving_distance(1.0).unit == u.Mpc
    assert cosmo._comoving_distance_z1z2(1.0, 2.0).unit == u.Mpc
    assert cosmo.comoving_transverse_distance(1.0).unit == u.Mpc
    assert cosmo._comoving_transverse_distance_z1z2(1.0, 2.0).unit == u.Mpc
    assert cosmo.angular_diameter_distance(1.0).unit == u.Mpc
    assert cosmo.angular_diameter_distance_z1z2(1.0, 2.0).unit == u.Mpc
    assert cosmo.luminosity_distance(1.0).unit == u.Mpc
    assert cosmo.lookback_time(1.0).unit == u.Gyr
    assert cosmo.lookback_distance(1.0).unit == u.Mpc
    assert cosmo.H(1.0).unit == u.km / u.Mpc / u.s
    assert cosmo.Tcmb(1.0).unit == u.K
    assert cosmo.Tcmb([0.0, 1.0]).unit == u.K
    assert cosmo.Tnu(1.0).unit == u.K
    assert cosmo.Tnu([0.0, 1.0]).unit == u.K
    assert cosmo.arcsec_per_kpc_comoving(1.0).unit == u.arcsec / u.kpc
    assert cosmo.arcsec_per_kpc_proper(1.0).unit == u.arcsec / u.kpc
    assert cosmo.kpc_comoving_per_arcmin(1.0).unit == u.kpc / u.arcmin
    assert cosmo.kpc_proper_per_arcmin(1.0).unit == u.kpc / u.arcmin
    assert cosmo.critical_density(1.0).unit == u.g / u.cm**3
    assert cosmo.comoving_volume(1.0).unit == u.Mpc**3
    assert cosmo.age(1.0).unit == u.Gyr
    assert cosmo.distmod(1.0).unit == u.mag


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


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_ocurv():
    # Test Ok evolution
    # Flat, boring case
    tcos = FlatLambdaCDM(70.0, 0.3)

    assert u.allclose(tcos.Ok0, 0.0)
    assert u.allclose(tcos.Ok(0), 0.0)
    z = np.array([0.0, 0.5, 1.0, 2.0])
    assert u.allclose(tcos.Ok(z), [0.0, 0.0, 0.0, 0.0], rtol=1e-6)

    # Not flat
    tcos = LambdaCDM(70.0, 0.3, 0.5, Tcmb0=u.Quantity(0.0, u.K))
    assert u.allclose(tcos.Ok0, 0.2)
    assert u.allclose(tcos.Ok(0), 0.2)
    assert u.allclose(tcos.Ok(z), [0.2, 0.22929936, 0.21621622, 0.17307692], rtol=1e-4)

    # Test the sum; note that Ogamma/Onu are 0
    assert u.allclose(
        tcos.Ok(z) + tcos.Om(z) + tcos.Ode(z), [1.0, 1.0, 1.0, 1.0], rtol=1e-5
    )


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_ode():
    # Test Ode evolution, turn off neutrinos, cmb
    tcos = FlatLambdaCDM(70.0, 0.3, Tcmb0=0)

    assert u.allclose(tcos.Ode0, 0.7)
    assert u.allclose(tcos.Ode(0), 0.7)
    z = np.array([0.0, 0.5, 1.0, 2.0])
    assert u.allclose(tcos.Ode(z), [0.7, 0.408759, 0.2258065, 0.07954545], rtol=1e-5)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_tcmb():
    cosmo = FlatLambdaCDM(70.4, 0.272, Tcmb0=2.5)

    assert u.allclose(cosmo.Tcmb0, 2.5 * u.K)
    assert u.allclose(cosmo.Tcmb(2), 7.5 * u.K)
    z = [0.0, 1.0, 2.0, 3.0, 9.0]
    assert u.allclose(cosmo.Tcmb(z), [2.5, 5.0, 7.5, 10.0, 25.0] * u.K, rtol=1e-6)
    # Make sure it's the same for integers
    z = [0, 1, 2, 3, 9]
    assert u.allclose(cosmo.Tcmb(z), [2.5, 5.0, 7.5, 10.0, 25.0] * u.K, rtol=1e-6)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_tnu():
    cosmo = FlatLambdaCDM(70.4, 0.272, Tcmb0=3.0)

    assert u.allclose(cosmo.Tnu0, 2.1412975665108247 * u.K, rtol=1e-6)
    assert u.allclose(cosmo.Tnu(2), 6.423892699532474 * u.K, rtol=1e-6)
    z = [0.0, 1.0, 2.0, 3.0]
    expected = [2.14129757, 4.28259513, 6.4238927, 8.56519027] * u.K
    assert u.allclose(cosmo.Tnu(z), expected, rtol=1e-6)

    # Test for integers
    z = [0, 1, 2, 3]
    assert u.allclose(cosmo.Tnu(z), expected, rtol=1e-6)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_kpc_methods():
    cosmo = FlatLambdaCDM(70.4, 0.272, Tcmb0=0.0)

    assert u.allclose(cosmo.arcsec_per_kpc_comoving(3), 0.0317179167 * u.arcsec / u.kpc)
    assert u.allclose(cosmo.arcsec_per_kpc_proper(3), 0.1268716668 * u.arcsec / u.kpc)
    assert u.allclose(cosmo.kpc_comoving_per_arcmin(3), 1891.6753126 * u.kpc / u.arcmin)
    assert u.allclose(cosmo.kpc_proper_per_arcmin(3), 472.918828 * u.kpc / u.arcmin)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_comoving_volume():
    c_flat = LambdaCDM(H0=70, Om0=0.27, Ode0=0.73, Tcmb0=0.0)
    c_open = LambdaCDM(H0=70, Om0=0.27, Ode0=0.0, Tcmb0=0.0)
    c_closed = LambdaCDM(H0=70, Om0=2, Ode0=0.0, Tcmb0=0.0)

    # test against ned wright's calculator (cubic Gpc)
    redshifts = np.array([0.5, 1, 2, 3, 5, 9])
    wright_flat = (
        np.array([29.123, 159.529, 630.427, 1178.531, 2181.485, 3654.802]) * u.Gpc**3
    )
    wright_open = (
        np.array([20.501, 99.019, 380.278, 747.049, 1558.363, 3123.814]) * u.Gpc**3
    )
    wright_closed = (
        np.array([12.619, 44.708, 114.904, 173.709, 258.82, 358.992]) * u.Gpc**3
    )
    # The wright calculator isn't very accurate, so we use a rather
    # modest precision
    assert u.allclose(c_flat.comoving_volume(redshifts), wright_flat, rtol=1e-2)
    assert u.allclose(c_open.comoving_volume(redshifts), wright_open, rtol=1e-2)
    assert u.allclose(c_closed.comoving_volume(redshifts), wright_closed, rtol=1e-2)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_differential_comoving_volume():
    from scipy.integrate import quad

    c_flat = LambdaCDM(H0=70, Om0=0.27, Ode0=0.73, Tcmb0=0.0)
    c_open = LambdaCDM(H0=70, Om0=0.27, Ode0=0.0, Tcmb0=0.0)
    c_closed = LambdaCDM(H0=70, Om0=2, Ode0=0.0, Tcmb0=0.0)

    # test that integration of differential_comoving_volume()
    #  yields same as comoving_volume()
    redshifts = np.array([0.5, 1, 2, 3, 5, 9])
    wright_flat = (
        np.array([29.123, 159.529, 630.427, 1178.531, 2181.485, 3654.802]) * u.Gpc**3
    )
    wright_open = (
        np.array([20.501, 99.019, 380.278, 747.049, 1558.363, 3123.814]) * u.Gpc**3
    )
    wright_closed = (
        np.array([12.619, 44.708, 114.904, 173.709, 258.82, 358.992]) * u.Gpc**3
    )

    # The wright calculator isn't very accurate, so we use a rather
    # modest precision.
    def ftemp(x):
        return c_flat.differential_comoving_volume(x).value

    def otemp(x):
        return c_open.differential_comoving_volume(x).value

    def ctemp(x):
        return c_closed.differential_comoving_volume(x).value

    # Multiply by solid_angle (4 * pi)
    assert u.allclose(
        np.array([4.0 * np.pi * quad(ftemp, 0, redshift)[0] for redshift in redshifts])
        * u.Mpc**3,
        wright_flat,
        rtol=1e-2,
    )
    assert u.allclose(
        np.array([4.0 * np.pi * quad(otemp, 0, redshift)[0] for redshift in redshifts])
        * u.Mpc**3,
        wright_open,
        rtol=1e-2,
    )
    assert u.allclose(
        np.array([4.0 * np.pi * quad(ctemp, 0, redshift)[0] for redshift in redshifts])
        * u.Mpc**3,
        wright_closed,
        rtol=1e-2,
    )


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_age():
    # WMAP7 but with Omega_relativisitic = 0
    tcos = FlatLambdaCDM(70.4, 0.272, Tcmb0=0.0)

    assert u.allclose(tcos.hubble_time, 13.889094057856937 * u.Gyr)
    assert u.allclose(tcos.age(4), 1.5823603508870991 * u.Gyr)
    assert u.allclose(tcos.age([1.0, 5.0]), [5.97113193, 1.20553129] * u.Gyr)
    assert u.allclose(tcos.age([1, 5]), [5.97113193, 1.20553129] * u.Gyr)

    # Add relativistic species
    tcos = FlatLambdaCDM(70.4, 0.272, Tcmb0=3.0)
    assert u.allclose(tcos.age(4), 1.5773003779230699 * u.Gyr)
    assert u.allclose(tcos.age([1, 5]), [5.96344942, 1.20093077] * u.Gyr)

    # And massive neutrinos
    tcos = FlatLambdaCDM(70.4, 0.272, Tcmb0=3.0, m_nu=0.1 * u.eV)
    assert u.allclose(tcos.age(4), 1.5546485439853412 * u.Gyr)
    assert u.allclose(tcos.age([1, 5]), [5.88448152, 1.18383759] * u.Gyr)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_distmod():
    # WMAP7 but with Omega_relativisitic = 0
    tcos = FlatLambdaCDM(70.4, 0.272, Tcmb0=0.0)

    assert u.allclose(tcos.hubble_distance, 4258.415596590909 * u.Mpc)
    assert u.allclose(tcos.distmod([1, 5]), [44.124857, 48.40167258] * u.mag)
    assert u.allclose(tcos.distmod([1.0, 5.0]), [44.124857, 48.40167258] * u.mag)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_neg_distmod():
    # Cosmology with negative luminosity distances (perfectly okay,
    #  if obscure)
    tcos = LambdaCDM(70, 0.2, 1.3, Tcmb0=0)

    assert u.allclose(
        tcos.luminosity_distance([50, 100]), [16612.44047622, -46890.79092244] * u.Mpc
    )
    assert u.allclose(tcos.distmod([50, 100]), [46.102167189, 48.355437790944] * u.mag)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_critical_density():
    from astropy.constants import codata2014

    # WMAP7 but with Omega_relativistic = 0
    # These tests will fail if astropy.const starts returning non-mks
    #  units by default; see the comment at the top of core.py.
    # critical_density0 is inversely proportional to G.
    tcos = FlatLambdaCDM(70.4, 0.272, Tcmb0=0.0)

    fac = (const.G / codata2014.G).to(u.dimensionless_unscaled).value
    assert u.allclose(
        tcos.critical_density0 * fac, 9.309668456020899e-30 * (u.g / u.cm**3)
    )
    assert u.allclose(tcos.critical_density0, tcos.critical_density(0))
    assert u.allclose(
        tcos.critical_density([1, 5]) * fac,
        [2.70352772e-29, 5.53739080e-28] * (u.g / u.cm**3),
    )
    assert u.allclose(
        tcos.critical_density([1.0, 5.0]) * fac,
        [2.70352772e-29, 5.53739080e-28] * (u.g / u.cm**3),
    )


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_integral():
    # Test integer vs. floating point inputs
    cosmo = LambdaCDM(H0=73.2, Om0=0.3, Ode0=0.50)

    assert u.allclose(
        cosmo.comoving_distance(3), cosmo.comoving_distance(3.0), rtol=1e-7
    )
    assert u.allclose(
        cosmo.comoving_distance([1, 2, 3, 5]),
        cosmo.comoving_distance([1.0, 2.0, 3.0, 5.0]),
        rtol=1e-7,
    )
    assert u.allclose(cosmo.efunc(6), cosmo.efunc(6.0), rtol=1e-7)
    assert u.allclose(cosmo.efunc([1, 2, 6]), cosmo.efunc([1.0, 2.0, 6.0]), rtol=1e-7)
    assert u.allclose(
        cosmo.inv_efunc([1, 2, 6]), cosmo.inv_efunc([1.0, 2.0, 6.0]), rtol=1e-7
    )


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_de_densityscale():
    cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.70)

    z = np.array([0.1, 0.2, 0.5, 1.5, 2.5])
    assert u.allclose(cosmo.de_density_scale(z), [1.0, 1.0, 1.0, 1.0, 1.0])

    # Integer check
    assert u.allclose(cosmo.de_density_scale(3), cosmo.de_density_scale(3.0), rtol=1e-7)
    assert u.allclose(
        cosmo.de_density_scale([1, 2, 3]),
        cosmo.de_density_scale([1.0, 2.0, 3.0]),
        rtol=1e-7,
    )


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_comoving_distance_z1z2():
    tcos = LambdaCDM(100, 0.3, 0.8, Tcmb0=0.0)

    with pytest.raises(ValueError):  # test diff size z1, z2 fail
        tcos._comoving_distance_z1z2((1, 2), (3, 4, 5))

    # Comoving distances are invertible
    assert u.allclose(
        tcos._comoving_distance_z1z2(1, 2), -tcos._comoving_distance_z1z2(2, 1)
    )

    z1 = 0, 0, 2, 0.5, 1
    z2 = 2, 1, 1, 2.5, 1.1
    results = (
        3767.90579253,
        2386.25591391,
        -1381.64987862,
        2893.11776663,
        174.1524683,
    ) * u.Mpc

    assert u.allclose(tcos._comoving_distance_z1z2(z1, z2), results)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_age_in_special_cosmologies():
    """Check that age in de Sitter and Einstein-de Sitter Universes work.

    Some analytic solutions fail at these critical points.
    """
    c_dS = FlatLambdaCDM(100, 0, Tcmb0=0)
    assert u.allclose(c_dS.age(0), np.inf * u.Gyr)
    assert u.allclose(c_dS.age(1), np.inf * u.Gyr)
    assert u.allclose(c_dS.lookback_time(0), 0 * u.Gyr)
    assert u.allclose(c_dS.lookback_time(1), 6.777539216261741 * u.Gyr)

    c_EdS = FlatLambdaCDM(100, 1, Tcmb0=0)
    assert u.allclose(c_EdS.age(0), 6.518614811154189 * u.Gyr)
    assert u.allclose(c_EdS.age(1), 2.3046783684542738 * u.Gyr)
    assert u.allclose(c_EdS.lookback_time(0), 0 * u.Gyr)
    assert u.allclose(c_EdS.lookback_time(1), 4.213936442699092 * u.Gyr)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_distance_in_special_cosmologies():
    """Check that de Sitter and Einstein-de Sitter Universes both work.

    Some analytic solutions fail at these critical points.
    """
    c_dS = FlatLambdaCDM(100, 0, Tcmb0=0)
    assert u.allclose(c_dS.comoving_distance(0), 0 * u.Mpc)
    assert u.allclose(c_dS.comoving_distance(1), 2997.92458 * u.Mpc)

    c_EdS = FlatLambdaCDM(100, 1, Tcmb0=0)
    assert u.allclose(c_EdS.comoving_distance(0), 0 * u.Mpc)
    assert u.allclose(c_EdS.comoving_distance(1), 1756.1435599923348 * u.Mpc)

    c_dS = LambdaCDM(100, 0, 1, Tcmb0=0)
    assert u.allclose(c_dS.comoving_distance(0), 0 * u.Mpc)
    assert u.allclose(c_dS.comoving_distance(1), 2997.92458 * u.Mpc)

    c_EdS = LambdaCDM(100, 1, 0, Tcmb0=0)
    assert u.allclose(c_EdS.comoving_distance(0), 0 * u.Mpc)
    assert u.allclose(c_EdS.comoving_distance(1), 1756.1435599923348 * u.Mpc)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_absorption_distance():
    tcos = FlatLambdaCDM(70.4, 0.272, Tcmb0=0.0)

    assert u.allclose(tcos.absorption_distance([1, 3]), [1.72576635, 7.98685853])
    assert u.allclose(tcos.absorption_distance([1.0, 3.0]), [1.72576635, 7.98685853])
    assert u.allclose(tcos.absorption_distance(3), 7.98685853)
    assert u.allclose(tcos.absorption_distance(3.0), 7.98685853)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_distance_broadcast():
    """Test array shape broadcasting for functions with single
    redshift inputs"""

    cosmo = FlatLambdaCDM(H0=70, Om0=0.27, m_nu=u.Quantity([0.0, 0.1, 0.011], u.eV))
    z = np.linspace(0.1, 1, 6)
    z_reshape2d = z.reshape(2, 3)
    z_reshape3d = z.reshape(3, 2, 1)
    # Things with units
    methods = [
        "comoving_distance",
        "luminosity_distance",
        "comoving_transverse_distance",
        "angular_diameter_distance",
        "distmod",
        "lookback_time",
        "age",
        "comoving_volume",
        "differential_comoving_volume",
        "kpc_comoving_per_arcmin",
    ]
    for method in methods:
        g = getattr(cosmo, method)
        value_flat = g(z)
        assert value_flat.shape == z.shape
        value_2d = g(z_reshape2d)
        assert value_2d.shape == z_reshape2d.shape
        value_3d = g(z_reshape3d)
        assert value_3d.shape == z_reshape3d.shape
        assert value_flat.unit == value_2d.unit
        assert value_flat.unit == value_3d.unit
        assert u.allclose(value_flat, value_2d.flatten())
        assert u.allclose(value_flat, value_3d.flatten())

    # Also test unitless ones
    methods = [
        "absorption_distance",
        "Om",
        "Ode",
        "Ok",
        "H",
        "w",
        "de_density_scale",
        "Onu",
        "Ogamma",
        "nu_relative_density",
    ]
    for method in methods:
        g = getattr(cosmo, method)
        value_flat = g(z)
        assert value_flat.shape == z.shape
        value_2d = g(z_reshape2d)
        assert value_2d.shape == z_reshape2d.shape
        value_3d = g(z_reshape3d)
        assert value_3d.shape == z_reshape3d.shape
        assert u.allclose(value_flat, value_2d.flatten())
        assert u.allclose(value_flat, value_3d.flatten())
