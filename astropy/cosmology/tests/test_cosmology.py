# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Stand-alone overall systems tests for :mod:`astropy.cosmology`."""

import numpy as np
import pytest

import astropy.units as u
from astropy.cosmology import flrw
from astropy.units import allclose
from astropy.utils.compat.optional_deps import HAS_SCIPY
from astropy.utils.exceptions import AstropyUserWarning

###############################################################################
# TODO! sort and refactor following tests.
# overall systems tests stay here, specific tests go to new test suite.


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_units():
    """Test if the right units are being returned"""

    cosmo = flrw.FlatLambdaCDM(H0=70, Om0=0.27, Tcmb0=2.0)
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


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_distance_broadcast():
    """Test array shape broadcasting for functions with single
    redshift inputs"""

    cosmo = flrw.FlatLambdaCDM(
        H0=70, Om0=0.27, m_nu=u.Quantity([0.0, 0.1, 0.011], u.eV)
    )
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
        assert allclose(value_flat, value_2d.flatten())
        assert allclose(value_flat, value_3d.flatten())

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
        assert allclose(value_flat, value_2d.flatten())
        assert allclose(value_flat, value_3d.flatten())

    # Test some dark energy models
    methods = ["Om", "Ode", "w", "de_density_scale"]
    for tcosmo in [
        flrw.LambdaCDM(H0=70, Om0=0.27, Ode0=0.5),
        flrw.wCDM(H0=70, Om0=0.27, Ode0=0.5, w0=-1.2),
        flrw.w0waCDM(H0=70, Om0=0.27, Ode0=0.5, w0=-1.2, wa=-0.2),
        flrw.wpwaCDM(H0=70, Om0=0.27, Ode0=0.5, wp=-1.2, wa=-0.2, zp=0.9),
        flrw.w0wzCDM(H0=70, Om0=0.27, Ode0=0.5, w0=-1.2, wz=0.1),
    ]:
        for method in methods:
            g = getattr(cosmo, method)
            value_flat = g(z)
            assert value_flat.shape == z.shape
            value_2d = g(z_reshape2d)
            assert value_2d.shape == z_reshape2d.shape
            value_3d = g(z_reshape3d)
            assert value_3d.shape == z_reshape3d.shape
            assert allclose(value_flat, value_2d.flatten())
            assert allclose(value_flat, value_3d.flatten())


# This class is to test whether the routines work correctly
# if one only overloads w(z)
class test_cos_sub(flrw.FLRW):
    def __init__(self):
        super().__init__(70.0, 0.27, 0.73, Tcmb0=0.0, name="test_cos")
        self._w0 = -0.9

    def w(self, z):
        return self._w0 * np.ones_like(z)


# Similar, but with neutrinos
class test_cos_subnu(flrw.FLRW):
    def __init__(self):
        super().__init__(
            70.0, 0.27, 0.73, Tcmb0=3.0, m_nu=0.1 * u.eV, name="test_cos_nu"
        )
        self._w0 = -0.8

    def w(self, z):
        return self._w0 * np.ones_like(z)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_de_subclass():
    # This is the comparison object
    z = [0.2, 0.4, 0.6, 0.9]
    cosmo = flrw.wCDM(H0=70, Om0=0.27, Ode0=0.73, w0=-0.9, Tcmb0=0.0)
    # Values taken from Ned Wrights advanced cosmo calculator, Aug 17 2012
    assert allclose(
        cosmo.luminosity_distance(z), [975.5, 2158.2, 3507.3, 5773.1] * u.Mpc, rtol=1e-3
    )
    # Now try the subclass that only gives w(z)
    cosmo = test_cos_sub()
    assert allclose(
        cosmo.luminosity_distance(z), [975.5, 2158.2, 3507.3, 5773.1] * u.Mpc, rtol=1e-3
    )
    # Test efunc
    assert allclose(cosmo.efunc(1.0), 1.7489240754, rtol=1e-5)
    assert allclose(cosmo.efunc([0.5, 1.0]), [1.31744953, 1.7489240754], rtol=1e-5)
    assert allclose(cosmo.inv_efunc([0.5, 1.0]), [0.75904236, 0.57178011], rtol=1e-5)
    # Test de_density_scale
    assert allclose(cosmo.de_density_scale(1.0), 1.23114444, rtol=1e-4)
    assert allclose(
        cosmo.de_density_scale([0.5, 1.0]), [1.12934694, 1.23114444], rtol=1e-4
    )

    # Add neutrinos for efunc, inv_efunc


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_efunc_vs_invefunc_flrw():
    """Test that efunc and inv_efunc give inverse values"""
    z0 = 0.5
    z = np.array([0.5, 1.0, 2.0, 5.0])

    # FLRW is abstract, so requires test_cos_sub defined earlier
    # This requires scipy, unlike the built-ins, because it
    # calls de_density_scale, which has an integral in it
    cosmo = test_cos_sub()
    assert allclose(cosmo.efunc(z0), 1.0 / cosmo.inv_efunc(z0))
    assert allclose(cosmo.efunc(z), 1.0 / cosmo.inv_efunc(z))
    # Add neutrinos
    cosmo = test_cos_subnu()
    assert allclose(cosmo.efunc(z0), 1.0 / cosmo.inv_efunc(z0))
    assert allclose(cosmo.efunc(z), 1.0 / cosmo.inv_efunc(z))


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_integral():
    # Test integer vs. floating point inputs
    cosmo = flrw.LambdaCDM(H0=73.2, Om0=0.3, Ode0=0.50)
    assert allclose(cosmo.comoving_distance(3), cosmo.comoving_distance(3.0), rtol=1e-7)
    assert allclose(
        cosmo.comoving_distance([1, 2, 3, 5]),
        cosmo.comoving_distance([1.0, 2.0, 3.0, 5.0]),
        rtol=1e-7,
    )
    assert allclose(cosmo.efunc(6), cosmo.efunc(6.0), rtol=1e-7)
    assert allclose(cosmo.efunc([1, 2, 6]), cosmo.efunc([1.0, 2.0, 6.0]), rtol=1e-7)
    assert allclose(
        cosmo.inv_efunc([1, 2, 6]), cosmo.inv_efunc([1.0, 2.0, 6.0]), rtol=1e-7
    )


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_de_densityscale():
    cosmo = flrw.LambdaCDM(H0=70, Om0=0.3, Ode0=0.70)
    z = np.array([0.1, 0.2, 0.5, 1.5, 2.5])
    assert allclose(cosmo.de_density_scale(z), [1.0, 1.0, 1.0, 1.0, 1.0])
    # Integer check
    assert allclose(cosmo.de_density_scale(3), cosmo.de_density_scale(3.0), rtol=1e-7)
    assert allclose(
        cosmo.de_density_scale([1, 2, 3]),
        cosmo.de_density_scale([1.0, 2.0, 3.0]),
        rtol=1e-7,
    )

    cosmo = flrw.wCDM(H0=70, Om0=0.3, Ode0=0.60, w0=-0.5)
    assert allclose(
        cosmo.de_density_scale(z),
        [1.15369, 1.31453, 1.83712, 3.95285, 6.5479],
        rtol=1e-4,
    )
    assert allclose(cosmo.de_density_scale(3), cosmo.de_density_scale(3.0), rtol=1e-7)
    assert allclose(
        cosmo.de_density_scale([1, 2, 3]),
        cosmo.de_density_scale([1.0, 2.0, 3.0]),
        rtol=1e-7,
    )

    cosmo = flrw.w0wzCDM(H0=70, Om0=0.3, Ode0=0.50, w0=-1, wz=0.5)
    assert allclose(
        cosmo.de_density_scale(z),
        [0.746048, 0.5635595, 0.25712378, 0.026664129, 0.0035916468],
        rtol=1e-4,
    )
    assert allclose(cosmo.de_density_scale(3), cosmo.de_density_scale(3.0), rtol=1e-7)
    assert allclose(
        cosmo.de_density_scale([1, 2, 3]),
        cosmo.de_density_scale([1.0, 2.0, 3.0]),
        rtol=1e-7,
    )

    cosmo = flrw.w0waCDM(H0=70, Om0=0.3, Ode0=0.70, w0=-1, wa=-0.5)
    assert allclose(
        cosmo.de_density_scale(z),
        [0.9934201, 0.9767912, 0.897450, 0.622236, 0.4458753],
        rtol=1e-4,
    )
    assert allclose(cosmo.de_density_scale(3), cosmo.de_density_scale(3.0), rtol=1e-7)
    assert allclose(
        cosmo.de_density_scale([1, 2, 3]),
        cosmo.de_density_scale([1.0, 2.0, 3.0]),
        rtol=1e-7,
    )

    cosmo = flrw.wpwaCDM(H0=70, Om0=0.3, Ode0=0.70, wp=-0.9, wa=0.2, zp=0.5)
    assert allclose(
        cosmo.de_density_scale(z),
        [1.012246048, 1.0280102, 1.087439, 1.324988, 1.565746],
        rtol=1e-4,
    )
    assert allclose(cosmo.de_density_scale(3), cosmo.de_density_scale(3.0), rtol=1e-7)
    assert allclose(
        cosmo.de_density_scale([1, 2, 3]),
        cosmo.de_density_scale([1.0, 2.0, 3.0]),
        rtol=1e-7,
    )


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_comoving_distance_z1z2():
    tcos = flrw.LambdaCDM(100, 0.3, 0.8, Tcmb0=0.0)
    with pytest.raises(ValueError):  # test diff size z1, z2 fail
        tcos._comoving_distance_z1z2((1, 2), (3, 4, 5))
    # Comoving distances are invertible
    assert allclose(
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

    assert allclose(tcos._comoving_distance_z1z2(z1, z2), results)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_age_in_special_cosmologies():
    """Check that age in de Sitter and Einstein-de Sitter Universes work.

    Some analytic solutions fail at these critical points.
    """
    c_dS = flrw.FlatLambdaCDM(100, 0, Tcmb0=0)
    assert allclose(c_dS.age(z=0), np.inf * u.Gyr)
    assert allclose(c_dS.age(z=1), np.inf * u.Gyr)
    assert allclose(c_dS.lookback_time(z=0), 0 * u.Gyr)
    assert allclose(c_dS.lookback_time(z=1), 6.777539216261741 * u.Gyr)

    c_EdS = flrw.FlatLambdaCDM(100, 1, Tcmb0=0)
    assert allclose(c_EdS.age(z=0), 6.518614811154189 * u.Gyr)
    assert allclose(c_EdS.age(z=1), 2.3046783684542738 * u.Gyr)
    assert allclose(c_EdS.lookback_time(z=0), 0 * u.Gyr)
    assert allclose(c_EdS.lookback_time(z=1), 4.213936442699092 * u.Gyr)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_distance_in_special_cosmologies():
    """Check that de Sitter and Einstein-de Sitter Universes both work.

    Some analytic solutions fail at these critical points.
    """
    c_dS = flrw.FlatLambdaCDM(100, 0, Tcmb0=0)
    assert allclose(c_dS.comoving_distance(z=0), 0 * u.Mpc)
    assert allclose(c_dS.comoving_distance(z=1), 2997.92458 * u.Mpc)

    c_EdS = flrw.FlatLambdaCDM(100, 1, Tcmb0=0)
    assert allclose(c_EdS.comoving_distance(z=0), 0 * u.Mpc)
    assert allclose(c_EdS.comoving_distance(z=1), 1756.1435599923348 * u.Mpc)

    c_dS = flrw.LambdaCDM(100, 0, 1, Tcmb0=0)
    assert allclose(c_dS.comoving_distance(z=0), 0 * u.Mpc)
    assert allclose(c_dS.comoving_distance(z=1), 2997.92458 * u.Mpc)

    c_EdS = flrw.LambdaCDM(100, 1, 0, Tcmb0=0)
    assert allclose(c_EdS.comoving_distance(z=0), 0 * u.Mpc)
    assert allclose(c_EdS.comoving_distance(z=1), 1756.1435599923348 * u.Mpc)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_comoving_transverse_distance_z1z2():
    tcos = flrw.FlatLambdaCDM(100, 0.3, Tcmb0=0.0)
    with pytest.raises(ValueError):  # test diff size z1, z2 fail
        tcos._comoving_transverse_distance_z1z2((1, 2), (3, 4, 5))
    # Tests that should actually work, target values computed with
    # http://www.astro.multivax.de:8000/phillip/angsiz_prog/README.HTML
    # Kayser, Helbig, and Schramm (Astron.Astrophys. 318 (1997) 680-686)
    assert allclose(
        tcos._comoving_transverse_distance_z1z2(1, 2), 1313.2232194828466 * u.Mpc
    )

    # In a flat universe comoving distance and comoving transverse
    # distance are identical
    z1 = 0, 0, 2, 0.5, 1
    z2 = 2, 1, 1, 2.5, 1.1

    assert allclose(
        tcos._comoving_distance_z1z2(z1, z2),
        tcos._comoving_transverse_distance_z1z2(z1, z2),
    )

    # Test Flat Universe with Omega_M > 1.  Rarely used, but perfectly valid.
    tcos = flrw.FlatLambdaCDM(100, 1.5, Tcmb0=0.0)
    results = (
        2202.72682564,
        1559.51679971,
        -643.21002593,
        1408.36365679,
        85.09286258,
    ) * u.Mpc

    assert allclose(tcos._comoving_transverse_distance_z1z2(z1, z2), results)

    # In a flat universe comoving distance and comoving transverse
    # distance are identical
    z1 = 0, 0, 2, 0.5, 1
    z2 = 2, 1, 1, 2.5, 1.1

    assert allclose(
        tcos._comoving_distance_z1z2(z1, z2),
        tcos._comoving_transverse_distance_z1z2(z1, z2),
    )
    # Test non-flat cases to avoid simply testing
    # comoving_distance_z1z2. Test array, array case.
    tcos = flrw.LambdaCDM(100, 0.3, 0.5, Tcmb0=0.0)
    results = (
        3535.931375645655,
        2226.430046551708,
        -1208.6817970036532,
        2595.567367601969,
        151.36592003406884,
    ) * u.Mpc

    assert allclose(tcos._comoving_transverse_distance_z1z2(z1, z2), results)

    # Test positive curvature with scalar, array combination.
    tcos = flrw.LambdaCDM(100, 1.0, 0.2, Tcmb0=0.0)
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

    assert allclose(tcos._comoving_transverse_distance_z1z2(z1, z2), results)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_angular_diameter_distance_z1z2():
    tcos = flrw.FlatLambdaCDM(70.4, 0.272, Tcmb0=0.0)
    with pytest.raises(ValueError):  # test diff size z1, z2 fail
        tcos.angular_diameter_distance_z1z2([1, 2], [3, 4, 5])
    # Tests that should actually work
    assert allclose(
        tcos.angular_diameter_distance_z1z2(1, 2), 646.22968662822018 * u.Mpc
    )

    z1 = 2  # Separate test for z2<z1, returns negative value with warning
    z2 = 1
    results = -969.34452994 * u.Mpc
    with pytest.warns(AstropyUserWarning, match="less than first redshift"):
        assert allclose(tcos.angular_diameter_distance_z1z2(z1, z2), results)

    z1 = 0, 0, 0.5, 1
    z2 = 2, 1, 2.5, 1.1
    results = (
        1760.0628637762106,
        1670.7497657219858,
        1159.0970895962193,
        115.72768186186921,
    ) * u.Mpc

    assert allclose(tcos.angular_diameter_distance_z1z2(z1, z2), results)

    z1 = 0.1
    z2 = 0.1, 0.2, 0.5, 1.1, 2
    results = (0.0, 332.09893173, 986.35635069, 1508.37010062, 1621.07937976) * u.Mpc
    assert allclose(tcos.angular_diameter_distance_z1z2(0.1, z2), results)

    # Non-flat (positive Ok0) test
    tcos = flrw.LambdaCDM(H0=70.4, Om0=0.2, Ode0=0.5, Tcmb0=0.0)
    assert allclose(
        tcos.angular_diameter_distance_z1z2(1, 2), 620.1175337852428 * u.Mpc
    )
    # Non-flat (negative Ok0) test
    tcos = flrw.LambdaCDM(H0=100, Om0=2, Ode0=1, Tcmb0=0.0)
    assert allclose(
        tcos.angular_diameter_distance_z1z2(1, 2), 228.42914659246014 * u.Mpc
    )


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_absorption_distance():
    tcos = flrw.FlatLambdaCDM(70.4, 0.272, Tcmb0=0.0)
    assert allclose(tcos.absorption_distance([1, 3]), [1.72576635, 7.98685853])
    assert allclose(tcos.absorption_distance([1.0, 3.0]), [1.72576635, 7.98685853])
    assert allclose(tcos.absorption_distance(3), 7.98685853)
    assert allclose(tcos.absorption_distance(3.0), 7.98685853)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_distances():
    # Test distance calculations for various special case
    #  scenarios (no relativistic species, normal, massive neutrinos)
    # These do not come from external codes -- they are just internal
    #  checks to make sure nothing changes if we muck with the distance
    #  calculators

    z = np.array([1.0, 2.0, 3.0, 4.0])

    # The pattern here is: no relativistic species, the relativistic
    # species with massless neutrinos, then massive neutrinos
    cos = flrw.LambdaCDM(75.0, 0.25, 0.5, Tcmb0=0.0)
    assert allclose(
        cos.comoving_distance(z),
        [2953.93001902, 4616.7134253, 5685.07765971, 6440.80611897] * u.Mpc,
        rtol=1e-4,
    )
    cos = flrw.LambdaCDM(75.0, 0.25, 0.6, Tcmb0=3.0, Neff=3, m_nu=u.Quantity(0.0, u.eV))
    assert allclose(
        cos.comoving_distance(z),
        [3037.12620424, 4776.86236327, 5889.55164479, 6671.85418235] * u.Mpc,
        rtol=1e-4,
    )
    cos = flrw.LambdaCDM(75.0, 0.3, 0.4, Tcmb0=3.0, Neff=3, m_nu=u.Quantity(10.0, u.eV))
    assert allclose(
        cos.comoving_distance(z),
        [2471.80626824, 3567.1902565, 4207.15995626, 4638.20476018] * u.Mpc,
        rtol=1e-4,
    )
    # Flat
    cos = flrw.FlatLambdaCDM(75.0, 0.25, Tcmb0=0.0)
    assert allclose(
        cos.comoving_distance(z),
        [3180.83488552, 5060.82054204, 6253.6721173, 7083.5374303] * u.Mpc,
        rtol=1e-4,
    )
    cos = flrw.FlatLambdaCDM(75.0, 0.25, Tcmb0=3.0, Neff=3, m_nu=u.Quantity(0.0, u.eV))
    assert allclose(
        cos.comoving_distance(z),
        [3180.42662867, 5059.60529655, 6251.62766102, 7080.71698117] * u.Mpc,
        rtol=1e-4,
    )
    cos = flrw.FlatLambdaCDM(75.0, 0.25, Tcmb0=3.0, Neff=3, m_nu=u.Quantity(10.0, u.eV))
    assert allclose(
        cos.comoving_distance(z),
        [2337.54183142, 3371.91131264, 3988.40711188, 4409.09346922] * u.Mpc,
        rtol=1e-4,
    )
    # Add w
    cos = flrw.FlatwCDM(75.0, 0.25, w0=-1.05, Tcmb0=0.0)
    assert allclose(
        cos.comoving_distance(z),
        [3216.8296894, 5117.2097601, 6317.05995437, 7149.68648536] * u.Mpc,
        rtol=1e-4,
    )
    cos = flrw.FlatwCDM(
        75.0, 0.25, w0=-0.95, Tcmb0=3.0, Neff=3, m_nu=u.Quantity(0.0, u.eV)
    )
    assert allclose(
        cos.comoving_distance(z),
        [3143.56537758, 5000.32196494, 6184.11444601, 7009.80166062] * u.Mpc,
        rtol=1e-4,
    )
    cos = flrw.FlatwCDM(
        75.0, 0.25, w0=-0.9, Tcmb0=3.0, Neff=3, m_nu=u.Quantity(10.0, u.eV)
    )
    assert allclose(
        cos.comoving_distance(z),
        [2337.76035371, 3372.1971387, 3988.71362289, 4409.40817174] * u.Mpc,
        rtol=1e-4,
    )
    # Non-flat w
    cos = flrw.wCDM(75.0, 0.25, 0.4, w0=-0.9, Tcmb0=0.0)
    assert allclose(
        cos.comoving_distance(z),
        [2849.6163356, 4428.71661565, 5450.97862778, 6179.37072324] * u.Mpc,
        rtol=1e-4,
    )
    cos = flrw.wCDM(
        75.0, 0.25, 0.4, w0=-1.1, Tcmb0=3.0, Neff=3, m_nu=u.Quantity(0.0, u.eV)
    )
    assert allclose(
        cos.comoving_distance(z),
        [2904.35580229, 4511.11471267, 5543.43643353, 6275.9206788] * u.Mpc,
        rtol=1e-4,
    )
    cos = flrw.wCDM(
        75.0, 0.25, 0.4, w0=-0.9, Tcmb0=3.0, Neff=3, m_nu=u.Quantity(10.0, u.eV)
    )
    assert allclose(
        cos.comoving_distance(z),
        [2473.32522734, 3581.54519631, 4232.41674426, 4671.83818117] * u.Mpc,
        rtol=1e-4,
    )
    # w0wa
    cos = flrw.w0waCDM(75.0, 0.3, 0.6, w0=-0.9, wa=0.1, Tcmb0=0.0)
    assert allclose(
        cos.comoving_distance(z),
        [2937.7807638, 4572.59950903, 5611.52821924, 6339.8549956] * u.Mpc,
        rtol=1e-4,
    )
    cos = flrw.w0waCDM(
        75.0, 0.25, 0.5, w0=-0.9, wa=0.1, Tcmb0=3.0, Neff=3, m_nu=u.Quantity(0.0, u.eV)
    )
    assert allclose(
        cos.comoving_distance(z),
        [2907.34722624, 4539.01723198, 5593.51611281, 6342.3228444] * u.Mpc,
        rtol=1e-4,
    )
    cos = flrw.w0waCDM(
        75.0, 0.25, 0.5, w0=-0.9, wa=0.1, Tcmb0=3.0, Neff=3, m_nu=u.Quantity(10.0, u.eV)
    )
    assert allclose(
        cos.comoving_distance(z),
        [2507.18336722, 3633.33231695, 4292.44746919, 4736.35404638] * u.Mpc,
        rtol=1e-4,
    )
    # Flatw0wa
    cos = flrw.Flatw0waCDM(75.0, 0.25, w0=-0.95, wa=0.15, Tcmb0=0.0)
    assert allclose(
        cos.comoving_distance(z),
        [3123.29892781, 4956.15204302, 6128.15563818, 6948.26480378] * u.Mpc,
        rtol=1e-4,
    )
    cos = flrw.Flatw0waCDM(
        75.0, 0.25, w0=-0.95, wa=0.15, Tcmb0=3.0, Neff=3, m_nu=u.Quantity(0.0, u.eV)
    )
    assert allclose(
        cos.comoving_distance(z),
        [3122.92671907, 4955.03768936, 6126.25719576, 6945.61856513] * u.Mpc,
        rtol=1e-4,
    )
    cos = flrw.Flatw0waCDM(
        75.0, 0.25, w0=-0.95, wa=0.15, Tcmb0=3.0, Neff=3, m_nu=u.Quantity(10.0, u.eV)
    )
    assert allclose(
        cos.comoving_distance(z),
        [2337.70072701, 3372.13719963, 3988.6571093, 4409.35399673] * u.Mpc,
        rtol=1e-4,
    )
    # wpwa
    cos = flrw.wpwaCDM(75.0, 0.3, 0.6, wp=-0.9, zp=0.5, wa=0.1, Tcmb0=0.0)
    assert allclose(
        cos.comoving_distance(z),
        [2954.68975298, 4599.83254834, 5643.04013201, 6373.36147627] * u.Mpc,
        rtol=1e-4,
    )
    cos = flrw.wpwaCDM(
        75.0,
        0.25,
        0.5,
        wp=-0.9,
        zp=0.4,
        wa=0.1,
        Tcmb0=3.0,
        Neff=3,
        m_nu=u.Quantity(0.0, u.eV),
    )
    assert allclose(
        cos.comoving_distance(z),
        [2919.00656215, 4558.0218123, 5615.73412391, 6366.10224229] * u.Mpc,
        rtol=1e-4,
    )
    cos = flrw.wpwaCDM(
        75.0,
        0.25,
        0.5,
        wp=-0.9,
        zp=1.0,
        wa=0.1,
        Tcmb0=3.0,
        Neff=4,
        m_nu=u.Quantity(5.0, u.eV),
    )
    assert allclose(
        cos.comoving_distance(z),
        [2629.48489827, 3874.13392319, 4614.31562397, 5116.51184842] * u.Mpc,
        rtol=1e-4,
    )

    # w0wz
    cos = flrw.w0wzCDM(75.0, 0.3, 0.6, w0=-0.9, wz=0.1, Tcmb0=0.0)
    assert allclose(
        cos.comoving_distance(z),
        [3051.68786716, 4756.17714818, 5822.38084257, 6562.70873734] * u.Mpc,
        rtol=1e-4,
    )
    cos = flrw.w0wzCDM(
        75.0, 0.25, 0.5, w0=-0.9, wz=0.1, Tcmb0=3.0, Neff=3, m_nu=u.Quantity(0.0, u.eV)
    )
    assert allclose(
        cos.comoving_distance(z),
        [2997.8115653, 4686.45599916, 5764.54388557, 6524.17408738] * u.Mpc,
        rtol=1e-4,
    )
    cos = flrw.w0wzCDM(
        75.0, 0.25, 0.5, w0=-0.9, wz=0.1, Tcmb0=3.0, Neff=4, m_nu=u.Quantity(5.0, u.eV)
    )
    assert allclose(
        cos.comoving_distance(z),
        [2676.73467639, 3940.57967585, 4686.90810278, 5191.54178243] * u.Mpc,
        rtol=1e-4,
    )

    # Also test different numbers of massive neutrinos
    # for FlatLambdaCDM to give the scalar nu density functions a
    # work out
    cos = flrw.FlatLambdaCDM(75.0, 0.25, Tcmb0=3.0, m_nu=u.Quantity([10.0, 0, 0], u.eV))
    assert allclose(
        cos.comoving_distance(z),
        [2777.71589173, 4186.91111666, 5046.0300719, 5636.10397302] * u.Mpc,
        rtol=1e-4,
    )
    cos = flrw.FlatLambdaCDM(75.0, 0.25, Tcmb0=3.0, m_nu=u.Quantity([10.0, 5, 0], u.eV))
    assert allclose(
        cos.comoving_distance(z),
        [2636.48149391, 3913.14102091, 4684.59108974, 5213.07557084] * u.Mpc,
        rtol=1e-4,
    )
    cos = flrw.FlatLambdaCDM(75.0, 0.25, Tcmb0=3.0, m_nu=u.Quantity([4.0, 5, 9], u.eV))
    assert allclose(
        cos.comoving_distance(z),
        [2563.5093049, 3776.63362071, 4506.83448243, 5006.50158829] * u.Mpc,
        rtol=1e-4,
    )
    cos = flrw.FlatLambdaCDM(
        75.0, 0.25, Tcmb0=3.0, Neff=4.2, m_nu=u.Quantity([1.0, 4.0, 5, 9], u.eV)
    )
    assert allclose(
        cos.comoving_distance(z),
        [2525.58017482, 3706.87633298, 4416.58398847, 4901.96669755] * u.Mpc,
        rtol=1e-4,
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
    tcos = flrw.FlatLambdaCDM(
        75.0, 0.25, Tcmb0=3.0, Neff=3, m_nu=u.Quantity(100.0, u.eV)
    )
    assert tcos.has_massive_nu
    assert tcos.Neff == 3
    nurel_exp = (
        nuprefac * tcos.Neff * np.array([171969, 85984.5, 57323, 15633.5, 171.801])
    )
    assert allclose(tcos.nu_relative_density(ztest), nurel_exp, rtol=5e-3)
    assert allclose(tcos.efunc([0.0, 1.0]), [1.0, 7.46144727668], rtol=5e-3)

    # Next, slightly less massive
    tcos = flrw.FlatLambdaCDM(
        75.0, 0.25, Tcmb0=3.0, Neff=3, m_nu=u.Quantity(0.25, u.eV)
    )
    nurel_exp = (
        nuprefac * tcos.Neff * np.array([429.924, 214.964, 143.312, 39.1005, 1.11086])
    )
    assert allclose(tcos.nu_relative_density(ztest), nurel_exp, rtol=5e-3)

    # For this one also test Onu directly
    onu_exp = np.array([0.01890217, 0.05244681, 0.0638236, 0.06999286, 0.1344951])
    assert allclose(tcos.Onu(ztest), onu_exp, rtol=5e-3)

    # And fairly light
    tcos = flrw.FlatLambdaCDM(
        80.0, 0.30, Tcmb0=3.0, Neff=3, m_nu=u.Quantity(0.01, u.eV)
    )

    nurel_exp = (
        nuprefac * tcos.Neff * np.array([17.2347, 8.67345, 5.84348, 1.90671, 1.00021])
    )
    assert allclose(tcos.nu_relative_density(ztest), nurel_exp, rtol=5e-3)
    onu_exp = np.array([0.00066599, 0.00172677, 0.0020732, 0.00268404, 0.0978313])
    assert allclose(tcos.Onu(ztest), onu_exp, rtol=5e-3)
    assert allclose(tcos.efunc([1.0, 2.0]), [1.76225893, 2.97022048], rtol=1e-4)
    assert allclose(tcos.inv_efunc([1.0, 2.0]), [0.5674535, 0.33667534], rtol=1e-4)

    # Now a mixture of neutrino masses, with non-integer Neff
    tcos = flrw.FlatLambdaCDM(
        80.0, 0.30, Tcmb0=3.0, Neff=3.04, m_nu=u.Quantity([0.0, 0.01, 0.25], u.eV)
    )
    nurel_exp = (
        nuprefac
        * tcos.Neff
        * np.array([149.386233, 74.87915, 50.0518, 14.002403, 1.03702333])
    )
    assert allclose(tcos.nu_relative_density(ztest), nurel_exp, rtol=5e-3)
    onu_exp = np.array([0.00584959, 0.01493142, 0.01772291, 0.01963451, 0.10227728])
    assert allclose(tcos.Onu(ztest), onu_exp, rtol=5e-3)

    # Integer redshifts
    ztest = ztest.astype(int)
    assert allclose(tcos.nu_relative_density(ztest), nurel_exp, rtol=5e-3)
    assert allclose(tcos.Onu(ztest), onu_exp, rtol=5e-3)
