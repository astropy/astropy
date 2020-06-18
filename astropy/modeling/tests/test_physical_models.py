# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for physical functions."""
# pylint: disable=no-member, invalid-name
import pytest
import numpy as np

from astropy.modeling.physical_models import BlackBody, NFW
from astropy.modeling.fitting import LevMarLSQFitter

from astropy.tests.helper import assert_quantity_allclose
from astropy import units as u
from astropy.utils.exceptions import AstropyUserWarning
from astropy import cosmology

try:
    from scipy import optimize, integrate  # noqa

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

__doctest_skip__ = ["*"]


# BlackBody tests


@pytest.mark.parametrize("temperature", (3000 * u.K, 2726.85 * u.deg_C))
def test_blackbody_evaluate(temperature):

    b = BlackBody(temperature=temperature, scale=1.0)

    assert_quantity_allclose(b(1.4 * u.micron), 486787299458.15656 * u.MJy / u.sr)
    assert_quantity_allclose(b(214.13747 * u.THz), 486787299458.15656 * u.MJy / u.sr)


def test_blackbody_weins_law():
    b = BlackBody(293.0 * u.K)
    assert_quantity_allclose(b.lambda_max, 9.890006672986939 * u.micron)
    assert_quantity_allclose(b.nu_max, 17.22525080856469 * u.THz)


def test_blackbody_sefanboltzman_law():
    b = BlackBody(293.0 * u.K)
    assert_quantity_allclose(b.bolometric_flux, 133.02471751812573 * u.W / (u.m * u.m))


def test_blackbody_return_units():
    # return of evaluate has no units when temperature has no units
    b = BlackBody(1000.0 * u.K, scale=1.0)
    assert not isinstance(b.evaluate(1.0 * u.micron, 1000.0, 1.0), u.Quantity)

    # return has "standard" units when scale has no units
    b = BlackBody(1000.0 * u.K, scale=1.0)
    assert isinstance(b(1.0 * u.micron), u.Quantity)
    assert b(1.0 * u.micron).unit == u.erg / (u.cm ** 2 * u.s * u.Hz * u.sr)

    # return has scale units when scale has units
    b = BlackBody(1000.0 * u.K, scale=1.0 * u.MJy / u.sr)
    assert isinstance(b(1.0 * u.micron), u.Quantity)
    assert b(1.0 * u.micron).unit == u.MJy / u.sr


@pytest.mark.skipif("not HAS_SCIPY")
def test_blackbody_fit():

    fitter = LevMarLSQFitter()

    b = BlackBody(3000 * u.K, scale=5e-17 * u.Jy / u.sr)

    wav = np.array([0.5, 5, 10]) * u.micron
    fnu = np.array([1, 10, 5]) * u.Jy / u.sr

    b_fit = fitter(b, wav, fnu, maxiter=1000)

    assert_quantity_allclose(b_fit.temperature, 2840.7438355865065 * u.K)
    assert_quantity_allclose(b_fit.scale, 5.803783292762381e-17 * u.Jy / u.sr)


def test_blackbody_overflow():
    """Test Planck function with overflow."""
    photlam = u.photon / (u.cm ** 2 * u.s * u.AA)
    wave = [0.0, 1000.0, 100000.0, 1e55]  # Angstrom
    temp = 10000.0  # Kelvin
    bb = BlackBody(temperature=temp * u.K, scale=1.0)
    with pytest.warns(
            AstropyUserWarning,
            match=r'Input contains invalid wavelength/frequency value\(s\)'):
        with np.errstate(all="ignore"):
            bb_lam = bb(wave) * u.sr
    flux = bb_lam.to(photlam, u.spectral_density(wave * u.AA)) / u.sr

    # First element is NaN, last element is very small, others normal
    assert np.isnan(flux[0])
    with np.errstate(all="ignore"):
        assert np.log10(flux[-1].value) < -134
    np.testing.assert_allclose(
        flux.value[1:-1], [0.00046368, 0.04636773], rtol=1e-3
    )  # 0.1% accuracy in PHOTLAM/sr
    with np.errstate(all="ignore"):
        flux = bb(1.0 * u.AA)
    assert flux.value == 0


def test_blackbody_exceptions_and_warnings():
    """Test exceptions."""

    # Negative temperature
    with pytest.raises(ValueError) as exc:
        bb = BlackBody(-100 * u.K)
        bb(1.0 * u.micron)
    assert exc.value.args[0] == "Temperature should be positive: [-100.] K"

    bb = BlackBody(5000 * u.K)

    # Zero wavelength given for conversion to Hz
    with pytest.warns(AstropyUserWarning, match='invalid') as w:
        bb(0 * u.AA)
    assert len(w) == 3  # 2 of these are RuntimeWarning from zero divide

    # Negative wavelength given for conversion to Hz
    with pytest.warns(AstropyUserWarning, match='invalid') as w:
        bb(-1.0 * u.AA)
    assert len(w) == 1

    # Test that a non surface brightness converatable scale unit
    with pytest.raises(ValueError) as exc:
        bb = BlackBody(5000 * u.K, scale=1.0 * u.Jy)
        bb(1.0 * u.micron)
    assert exc.value.args[0] == "scale units not surface brightness: Jy"


def test_blackbody_array_temperature():
    """Regression test to make sure that the temperature can be an array."""
    multibb = BlackBody([100, 200, 300] * u.K)
    flux = multibb(1.2 * u.mm)
    np.testing.assert_allclose(
        flux.value, [1.804908e-12, 3.721328e-12, 5.638513e-12], rtol=1e-5
    )

    flux = multibb([2, 4, 6] * u.mm)
    np.testing.assert_allclose(
        flux.value, [6.657915e-13, 3.420677e-13, 2.291897e-13], rtol=1e-5
    )

    multibb = BlackBody(np.ones(4) * u.K)
    flux = multibb(np.ones((3, 4)) * u.mm)
    assert flux.shape == (3, 4)


@pytest.mark.parametrize("mass", (2.0000000000000E15 * u.M_sun, 3.976819741e+45 * u.kg))
def test_NFW_evaluate(mass):
    """Evalution, density, and radii validation of NFW model."""
    # Test parameters
    concentration = 8.5
    redshift = 0.63
    cosmo = cosmology.Planck15

    # Parsec tests
    # 200c Overdensity
    massfactor = ("critical", 200)

    n200c = NFW(mass=mass, concentration=concentration, redshift=redshift, cosmo=cosmo,
                massfactor=massfactor)
    assert_quantity_allclose(n200c(3.0 * u.Mpc), (3.709693508e+12 * (u.solMass / u.Mpc ** 3),
                                                  7.376391187e+42 * (u.kg / u.Mpc ** 3)))
    assert_quantity_allclose(n200c.rho_scale, (7800150779863018.0 * (u.solMass / u.Mpc ** 3)))
    assert_quantity_allclose(n200c.r_s, (0.24684627641195428 * u.Mpc))
    assert_quantity_allclose(n200c.r_virial, (2.0981933495016114 * u.Mpc))

    # 200m Overdensity
    massfactor = ("mean", 200)

    n200m = NFW(mass=mass, concentration=concentration, redshift=redshift, cosmo=cosmo,
                massfactor=massfactor)
    assert_quantity_allclose(n200m(3.0 * u.Mpc), (3.626093406e+12 * (u.solMass / u.Mpc**3),
                                                  7.210159921e+42 * (u.kg / u.Mpc**3)))
    assert_quantity_allclose(n200m.rho_scale, (5118547639858115.0 * (u.solMass / u.Mpc ** 3)))
    assert_quantity_allclose(n200m.r_s, (0.2840612517326848 * u.Mpc))
    assert_quantity_allclose(n200m.r_virial, (2.414520639727821 * u.Mpc))

    # Virial mass
    massfactor = ("virial")

    nvir = NFW(mass=mass, concentration=concentration, redshift=redshift, cosmo=cosmo,
               massfactor=massfactor)
    assert_quantity_allclose(nvir(3.0 * u.Mpc), (3.646475546e+12 * (u.solMass / u.Mpc**3),
                                                 7.250687967e+42 * (u.kg / u.Mpc**3)))
    assert_quantity_allclose(nvir.rho_scale, (5649367524651067.0 * (u.solMass / u.Mpc ** 3)))
    assert_quantity_allclose(nvir.r_s, (0.2748701862303786 * u.Mpc))
    assert_quantity_allclose(nvir.r_virial, (2.3363965829582183 * u.Mpc))

    # kpc tests
    # 200c Overdensity
    massfactor = ("critical", 200)

    n200c = NFW(mass=mass, concentration=concentration, redshift=redshift, cosmo=cosmo,
                massfactor=massfactor)
    assert_quantity_allclose(n200c(3141 * u.kpc), (3254.373619264334 * (u.solMass / u.kpc ** 3),
                                                   6.471028627484543e+33 * (u.kg / u.kpc ** 3)))
    assert_quantity_allclose(n200c.rho_scale, (7800150.779863021 * (u.solMass / u.kpc ** 3)))
    assert_quantity_allclose(n200c.r_s, (246.84627641195425 * u.kpc))
    assert_quantity_allclose(n200c.r_virial, (2098.193349501611 * u.kpc))

    # 200m Overdensity
    massfactor = ("mean", 200)

    n200m = NFW(mass=mass, concentration=concentration, redshift=redshift, cosmo=cosmo,
                massfactor=massfactor)
    assert_quantity_allclose(n200m(3141 * u.kpc), (3184.0370866188623 * (u.solMass / u.kpc**3),
                                                   6.33117077170161e+33 * (u.kg / u.kpc**3)))
    assert_quantity_allclose(n200m.rho_scale, (5118547.639858116 * (u.solMass / u.kpc ** 3)))
    assert_quantity_allclose(n200m.r_s, (284.0612517326848 * u.kpc))
    assert_quantity_allclose(n200m.r_virial, (2414.5206397278207 * u.kpc))

    # Virial mass
    massfactor = ("virial")

    nvir = NFW(mass=mass, concentration=concentration, redshift=redshift, cosmo=cosmo,
               massfactor=massfactor)
    assert_quantity_allclose(nvir(3141 * u.kpc), (3201.1946851294997 * (u.solMass / u.kpc**3),
                                                  6.365287109937637e+33 * (u.kg / u.kpc**3)))
    assert_quantity_allclose(nvir.rho_scale, (5649367.5246510655 * (u.solMass / u.kpc ** 3)))
    assert_quantity_allclose(nvir.r_s, (274.87018623037864 * u.kpc))
    assert_quantity_allclose(nvir.r_virial, (2336.3965829582185 * u.kpc))

    # Meter tests
    # 200c Overdensity
    massfactor = ("critical", 200)

    n200c = NFW(mass=mass, concentration=concentration, redshift=redshift, cosmo=cosmo,
                massfactor=massfactor)
    assert_quantity_allclose(n200c(4.2e+23 * u.m), (1.527649658673012e-57 * (u.solMass / u.m ** 3),
                                                    3.0375936602739256e-27 * (u.kg / u.m ** 3)))
    assert_quantity_allclose(n200c.rho_scale, (2.654919529637763e-52 * (u.solMass / u.m ** 3)))
    assert_quantity_allclose(n200c.r_s, (7.616880211930209e+21 * u.m))
    assert_quantity_allclose(n200c.r_virial, (6.474348180140678e+22 * u.m))

    # 200m Overdensity
    massfactor = ("mean", 200)

    n200m = NFW(mass=mass, concentration=concentration, redshift=redshift, cosmo=cosmo,
                massfactor=massfactor)
    assert_quantity_allclose(n200m(4.2e+23 * u.m), (1.5194778058079436e-57 * (u.solMass / u.m ** 3),
                                                    3.0213446673751314e-27 * (u.kg / u.m ** 3)))
    assert_quantity_allclose(n200m.rho_scale, (1.742188385322371e-52 * (u.solMass / u.m ** 3)))
    assert_quantity_allclose(n200m.r_s, (8.76521436235054e+21 * u.m))
    assert_quantity_allclose(n200m.r_virial, (7.450432207997959e+22 * u.m))

    # Virial mass
    massfactor = ("virial")

    nvir = NFW(mass=mass, concentration=concentration, redshift=redshift, cosmo=cosmo,
               massfactor=massfactor)
    assert_quantity_allclose(nvir(4.2e+23 * u.m), (1.5214899184117633e-57 * (u.solMass / u.m ** 3),
                                                   3.0253455719375224e-27 * (u.kg / u.m ** 3)))
    assert_quantity_allclose(nvir.rho_scale, (1.922862338766335e-52 * (u.solMass / u.m ** 3)))
    assert_quantity_allclose(nvir.r_s, (8.481607714647913e+21 * u.m))
    assert_quantity_allclose(nvir.r_virial, (7.209366557450727e+22 * u.m))

    # Verify string input of overdensity type
    # 200c Overdensity
    massfactor = "200c"

    n200c = NFW(mass=mass, concentration=concentration, redshift=redshift, cosmo=cosmo,
                massfactor=massfactor)
    assert_quantity_allclose(n200c(3.0 * u.Mpc), (3.709693508e+12 * (u.solMass / u.Mpc ** 3),
                                                  7.376391187e+42 * (u.kg / u.Mpc ** 3)))

    # 200m Overdensity
    massfactor = "200m"

    n200m = NFW(mass=mass, concentration=concentration, redshift=redshift, cosmo=cosmo,
                massfactor=massfactor)
    assert_quantity_allclose(n200m(3.0 * u.Mpc), (3.626093406e+12 * (u.solMass / u.Mpc**3),
                                                  7.210159921e+42 * (u.kg / u.Mpc**3)))

    # Virial mass
    massfactor = "virial"

    nvir = NFW(mass=mass, concentration=concentration, redshift=redshift,  cosmo=cosmo,
               massfactor=massfactor)
    assert_quantity_allclose(nvir(3.0 * u.Mpc), (3.646475546e+12 * (u.solMass / u.Mpc**3),
                                                 7.250687967e+42 * (u.kg / u.Mpc**3)))


@pytest.mark.skipif("not HAS_SCIPY")
def test_NFW_fit():
    """Test linear fitting of NFW model."""
    # Fixed parameters
    redshift = 0.63
    cosmo = cosmology.Planck15

    # Radial set
    r = np.array([1.00e+01, 1.00e+02, 2.00e+02, 2.50e+02, 3.00e+02, 4.00e+02, 5.00e+02,
                  7.50e+02, 1.00e+03, 1.50e+03, 2.50e+03, 6.50e+03, 1.15e+04]) * u.kpc

    # 200c Overdensity
    massfactor = ("critical", 200)

    density_r = np.array([1.77842761e+08, 9.75233623e+06, 2.93789626e+06, 1.90107238e+06,
                          1.30776878e+06, 7.01004140e+05, 4.20678479e+05, 1.57421880e+05,
                          7.54669701e+04, 2.56319769e+04, 6.21976562e+03, 3.96522424e+02,
                          7.39336808e+01]) * (u.solMass / u.kpc ** 3)

    fitter = LevMarLSQFitter()

    n200c = NFW(mass=1.8E15 * u.M_sun, concentration=7.0, redshift=redshift, cosmo=cosmo,
                massfactor=massfactor)
    n200c.redshift.fixed = True

    n_fit = fitter(n200c, r, density_r, maxiter=1000)

    assert_quantity_allclose(n_fit.mass, 2.0000000000000E15 * u.M_sun)
    assert_quantity_allclose(n_fit.concentration, 8.5)

    # 200m Overdensity
    massfactor = ("mean", 200)

    density_r = np.array([1.35677282e+08, 7.95392979e+06, 2.50352599e+06, 1.64535870e+06,
                          1.14642248e+06, 6.26805453e+05, 3.81691731e+05, 1.46294819e+05,
                          7.11559560e+04, 2.45737796e+04, 6.05459585e+03, 3.92183991e+02,
                          7.34674416e+01]) * (u.solMass / u.kpc ** 3)

    fitter = LevMarLSQFitter()

    n200m = NFW(mass=1.8E15 * u.M_sun, concentration=7.0, redshift=redshift, cosmo=cosmo,
                massfactor=massfactor)
    n200m.redshift.fixed = True

    n_fit = fitter(n200m, r, density_r, maxiter=1000)

    assert_quantity_allclose(n_fit.mass, 2.0000000000000E15 * u.M_sun)
    assert_quantity_allclose(n_fit.concentration, 8.5)

    # Virial mass
    massfactor = ("virial", 200)

    density_r = np.array([1.44573515e+08, 8.34873998e+06, 2.60137484e+06, 1.70348738e+06,
                          1.18337370e+06, 6.43994654e+05, 3.90800249e+05, 1.48930537e+05,
                          7.21856397e+04, 2.48289464e+04, 6.09477095e+03, 3.93248818e+02,
                          7.35821787e+01]) * (u.solMass / u.kpc ** 3)

    fitter = LevMarLSQFitter()

    nvir = NFW(mass=1.8E15 * u.M_sun, concentration=7.0, redshift=redshift, cosmo=cosmo,
               massfactor=massfactor)
    nvir.redshift.fixed = True

    n_fit = fitter(nvir, r, density_r, maxiter=1000)

    assert_quantity_allclose(n_fit.mass, 2.0000000000000E15 * u.M_sun)
    assert_quantity_allclose(n_fit.concentration, 8.5)


def test_NFW_circular_velocity():
    """Test circular velocity and radial validation of NFW model."""
    # Test parameters
    mass = 2.0000000000000E15 * u.M_sun
    concentration = 8.5
    redshift = 0.63
    cosmo = cosmology.Planck15

    r_r = np.array([0.01, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.5, 6.5, 11.5]) * u.Mpc

    # 200c Overdensity tests
    massfactor = ("critical", 200)

    n200c = NFW(mass=mass, concentration=concentration, redshift=redshift, cosmo=cosmo,
                massfactor=massfactor)
    circ_v_200c = np.array([702.45487454, 1812.4138346, 2150.50929296, 2231.5802568, 2283.96950242,
                            2338.45989696, 2355.78876772, 2332.41766543, 2276.89433811,
                            2154.53909153, 1950.07947819, 1512.37442943,
                            1260.94034541]) * (u.km / u.s)
    assert_quantity_allclose(n200c.circular_velocity(r_r), circ_v_200c)
    assert_quantity_allclose(n200c.r_max, (0.5338248204429641 * u.Mpc))
    assert_quantity_allclose(n200c.v_max, (2356.7204380904027 * (u.km / u.s)))

    # 200m Overdensity tests
    massfactor = ("mean", 200)
    mass = 1.0e14 * u.M_sun
    concentration = 12.3
    redshift = 1.5

    n200m = NFW(mass=mass, concentration=concentration, redshift=redshift, cosmo=cosmo,
                massfactor=massfactor)
    circ_v_200m = np.array([670.18236647, 1088.9843324, 1046.82334367, 1016.88890732, 987.97273478,
                            936.00207134, 891.80115232, 806.63307977, 744.91002191, 659.33401039,
                            557.82823549, 395.9735786, 318.29863006]) * (u.km / u.s)
    assert_quantity_allclose(n200m.circular_velocity(r_r), circ_v_200m)
    assert_quantity_allclose(n200m.r_max, (0.10196917920081808 * u.Mpc))
    assert_quantity_allclose(n200m.v_max, (1089.0224395818727 * (u.km / u.s)))

    # Virial Overdensity tests
    massfactor = ("virial")
    mass = 1.2e+45 * u.kg
    concentration = 2.4
    redshift = 0.34

    r_r = np.array([3.08567758e+20, 3.08567758e+21, 6.17135516e+21, 7.71419395e+21,
                    9.25703274e+21, 1.23427103e+22, 1.54283879e+22, 2.31425819e+22,
                    3.08567758e+22, 4.62851637e+22, 7.71419395e+22, 2.00569043e+23,
                    3.54852922e+23]) * u.m

    nvir = NFW(mass=mass, concentration=concentration, redshift=redshift, cosmo=cosmo,
               massfactor=massfactor)
    circ_v_vir = np.array([205.87461783, 604.65091823, 793.9190629, 857.52516521, 908.90280843,
                           986.53582718, 1041.69089845, 1124.19719446, 1164.58270747, 1191.33193561,
                           1174.02934755, 1023.69360527, 895.52206321]) * (u.km / u.s)
    assert_quantity_allclose(nvir.circular_velocity(r_r), circ_v_vir)
    assert_quantity_allclose(nvir.r_max, (1.6484542328623448 * u.Mpc))
    assert_quantity_allclose(nvir.v_max, (1192.3130989914962 * (u.km / u.s)))


def test_NFW_exceptions_and_warnings_and_misc():
    """Test NFW exceptions."""

    # Arbitrary Test parameters
    mass = 2.0000000000000E15 * u.M_sun
    concentration = 8.5
    redshift = 0.63
    cosmo = cosmology.Planck15
    massfactor = ("critical", 200)

    r_r = np.array([1.00e+01, 1.00e+02, 2.00e+02, 2.50e+02, 3.00e+02, 4.00e+02, 5.00e+02,
                    7.50e+02, 1.00e+03, 1.50e+03, 2.50e+03, 6.50e+03, 1.15e+04]) * u.kpc

    # Massfactor exception tests
    with pytest.raises(ValueError) as exc:
        NFW(mass=mass, concentration=concentration, redshift=redshift, cosmo=cosmo,
            massfactor=("not", "virial"))
    assert exc.value.args[0] == "Massfactor 'not' not one of 'critical', 'mean', or 'virial'"
    with pytest.raises(ValueError) as exc:
        NFW(mass=mass, concentration=concentration, redshift=redshift, cosmo=cosmo,
            massfactor="not virial")
    assert exc.value.args[0] == "Massfactor not virial string not of the form '#m', '#c', " \
                                "or 'virial'"
    with pytest.raises(TypeError) as exc:
        NFW(mass=mass, concentration=concentration, redshift=redshift, cosmo=cosmo,
            massfactor=200)
    assert exc.value.args[0] == "Massfactor 200 not a tuple or string"

    # Verify unitless mass
    # Density test
    n200c = NFW(mass=mass.value, concentration=concentration, redshift=redshift, cosmo=cosmo,
                massfactor=massfactor)
    assert_quantity_allclose(n200c(3000.0), (3.709693508e+12 * (u.solMass / u.Mpc ** 3),
                                             7.376391187e+42 * (u.kg / u.Mpc ** 3)))

    # Circular velocity test with unitless mass
    circ_v_200c = np.array([702.45487454, 1812.4138346, 2150.50929296, 2231.5802568, 2283.96950242,
                            2338.45989696, 2355.78876772, 2332.41766543, 2276.89433811,
                            2154.53909153, 1950.07947819, 1512.37442943,
                            1260.94034541]) * (u.km / u.s)
    assert_quantity_allclose(n200c.circular_velocity(r_r), circ_v_200c)

    # Test Default Cosmology
    ncos = NFW(mass=mass, concentration=concentration, redshift=redshift)
    assert_quantity_allclose(ncos.A_NFW(concentration), 1.356554956501232)
