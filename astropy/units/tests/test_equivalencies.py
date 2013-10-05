# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Separate tests specifically for equivalencies
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np
from numpy.testing.utils import assert_allclose

from ...tests.helper import pytest

from ...extern.six.moves import zip
from ... import units as u

functions = [u.doppler_optical, u.doppler_radio, u.doppler_relativistic]


@pytest.mark.parametrize(('function'), functions)
def test_doppler_frequency_0(function):
    rest = 105.01 * u.GHz
    velo0 = rest.to(u.km/u.s, equivalencies=function(rest))
    assert velo0.value == 0


@pytest.mark.parametrize(('function'), functions)
def test_doppler_wavelength_0(function):
    rest = 105.01 * u.GHz
    q1 = 0.00285489437196 * u.m
    velo0 = q1.to(u.km/u.s, equivalencies=function(rest))
    np.testing.assert_almost_equal(velo0.value, 0, decimal=6)


@pytest.mark.parametrize(('function'), functions)
def test_doppler_energy_0(function):
    rest = 105.01 * u.GHz
    q1 = 0.000434286445543 * u.eV
    velo0 = q1.to(u.km/u.s, equivalencies=function(rest))
    np.testing.assert_almost_equal(velo0.value, 0, decimal=6)


@pytest.mark.parametrize(('function'), functions)
def test_doppler_frequency_circle(function):
    rest = 105.01 * u.GHz
    shifted = 105.03 * u.GHz
    velo = shifted.to(u.km/u.s, equivalencies=function(rest))
    freq = velo.to(u.GHz, equivalencies=function(rest))
    np.testing.assert_almost_equal(freq.value, shifted.value, decimal=7)


@pytest.mark.parametrize(('function'), functions)
def test_doppler_wavelength_circle(function):
    rest = 105.01 * u.nm
    shifted = 105.03 * u.nm
    velo = shifted.to(u.km / u.s, equivalencies=function(rest))
    wav = velo.to(u.nm, equivalencies=function(rest))
    np.testing.assert_almost_equal(wav.value, shifted.value, decimal=7)


@pytest.mark.parametrize(('function'), functions)
def test_doppler_energy_circle(function):
    rest = 1.0501 * u.eV
    shifted = 1.0503 * u.eV
    velo = shifted.to(u.km / u.s, equivalencies=function(rest))
    en = velo.to(u.eV, equivalencies=function(rest))
    np.testing.assert_almost_equal(en.value, shifted.value, decimal=7)


values_ghz = (999.899940784289,999.8999307714406,999.8999357778647)
@pytest.mark.parametrize(('function', 'value'),
                         list(zip(functions, values_ghz)))
def test_30kms(function, value):
    rest = 1000 * u.GHz
    velo = 30 * u.km/u.s
    shifted = velo.to(u.GHz, equivalencies=function(rest))
    np.testing.assert_almost_equal(shifted.value, value, decimal=7)

def test_massenergy():
    # The relative tolerance of these tests is set by the uncertainties
    # in the charge of the electron, which is known to about
    # 3e-9 (relative tolerance).  Therefore, we limit the
    # precision of the tests to 1e-7 to be safe.  The masses are
    # (loosely) known to ~ 5e-8 rel tolerance, so we couldn't test to
    # 1e-7 if we used the values from astropy.constants; that is,
    # they might change by more than 1e-7 in some future update, so instead
    # they are hardwired here.

    # Electron, proton, neutron, muon, 1g
    mass_eV = u.Quantity([510.998928e3, 938.272046e6, 939.565378e6,
                          105.6583715e6, 5.60958884539e32], u.eV)
    mass_g = u.Quantity([9.10938291e-28, 1.672621777e-24, 1.674927351e-24,
                         1.88353147e-25, 1], u.g)
    # Test both ways
    assert np.allclose(mass_eV.to(u.g, equivalencies=u.mass_energy()).value,
                       mass_g.value, rtol=1e-7)
    assert np.allclose(mass_g.to(u.eV, equivalencies=u.mass_energy()).value,
                       mass_eV.value, rtol=1e-7)

    # Basic tests of 'derived' equivalencies
    # Surface density
    sdens_eV = u.Quantity(5.60958884539e32, u.eV / u.m**2)
    sdens_g = u.Quantity(1e-4, u.g / u.cm**2)
    assert np.allclose(sdens_eV.to(u.g / u.cm**2,
                                   equivalencies=u.mass_energy()).value,
                       sdens_g.value, rtol=1e-7)
    assert np.allclose(sdens_g.to(u.eV / u.m**2,
                                  equivalencies=u.mass_energy()).value,
                       sdens_eV.value, rtol=1e-7)

    # Density
    dens_eV = u.Quantity(5.60958884539e32, u.eV / u.m**3)
    dens_g = u.Quantity(1e-6, u.g / u.cm**3)
    assert np.allclose(dens_eV.to(u.g / u.cm**3,
                                   equivalencies=u.mass_energy()).value,
                       dens_g.value, rtol=1e-7)
    assert np.allclose(dens_g.to(u.eV / u.m**3,
                                  equivalencies=u.mass_energy()).value,
                       dens_eV.value, rtol=1e-7)

    # Power
    pow_eV = u.Quantity(5.60958884539e32, u.eV / u.s)
    pow_g = u.Quantity(1, u.g / u.s)
    assert np.allclose(pow_eV.to(u.g / u.s,
                                   equivalencies=u.mass_energy()).value,
                       pow_g.value, rtol=1e-7)
    assert np.allclose(pow_g.to(u.eV / u.s,
                                  equivalencies=u.mass_energy()).value,
                       pow_eV.value, rtol=1e-7)

def test_is_equivalent():

    assert u.m.is_equivalent(u.inch)
    assert not (u.Hz.is_equivalent(u.J))
    assert u.Hz.is_equivalent(u.J, u.spectral())
    assert u.J.is_equivalent(u.Hz, u.spectral())
    assert u.pc.is_equivalent(u.arcsecond, u.parallax())
    assert u.arcminute.is_equivalent(u.au, u.parallax())

    # Pass a tuple for multiple possibilities
    assert u.cm.is_equivalent((u.m, u.s, u.kg))
    assert u.ms.is_equivalent((u.m, u.s, u.kg))
    assert u.g.is_equivalent((u.m, u.s, u.kg))

def test_parallax():
    a = u.arcsecond.to(u.pc, 10, u.parallax())
    assert_allclose(a, 0.10)
    b = u.pc.to(u.arcsecond, a, u.parallax())
    assert_allclose(b, 10)

    a = u.arcminute.to(u.au, 1, u.parallax())
    assert_allclose(a, 3437.7467916)
    b = u.au.to(u.arcminute, a, u.parallax())
    assert_allclose(b, 1)


def test_parallax2():
    a = u.arcsecond.to(u.pc, [0.1, 2.5], u.parallax())
    assert_allclose(a, [10, 0.4])


def test_spectral():
    a = u.AA.to(u.Hz, 1, u.spectral())
    assert_allclose(a, 2.9979245799999995e+18)
    b = u.Hz.to(u.AA, a, u.spectral())
    assert_allclose(b, 1)

    a = u.AA.to(u.MHz, 1, u.spectral())
    assert_allclose(a, 2.9979245799999995e+12)
    b = u.MHz.to(u.AA, a, u.spectral())
    assert_allclose(b, 1)

    a = u.m.to(u.Hz, 1, u.spectral())
    assert_allclose(a, 2.9979245799999995e+8)
    b = u.Hz.to(u.m, a, u.spectral())
    assert_allclose(b, 1)


def test_spectral2():
    a = u.nm.to(u.J, 500, u.spectral())
    assert_allclose(a, 3.972891366538605e-19)
    b = u.J.to(u.nm, a, u.spectral())
    assert_allclose(b, 500)

    a = u.AA.to(u.Hz, 1, u.spectral())
    b = u.Hz.to(u.J, a, u.spectral())
    c = u.AA.to(u.J, 1, u.spectral())
    assert_allclose(b, c)


def test_spectral3():
    a = u.nm.to(u.Hz, [1000, 2000], u.spectral())
    assert_allclose(a, [2.99792458e+14, 1.49896229e+14])


def test_spectraldensity():

    a = u.AA.to(u.Jy, 1, u.spectral_density(u.eV, 2.2))
    assert_allclose(a, 1059416252057.8357, rtol=1e-4)

    b = u.Jy.to(u.AA, a, u.spectral_density(u.eV, 2.2))
    assert_allclose(b, 1)


def test_spectraldensity2():
    flambda = u.erg / u.angstrom / u.cm ** 2 / u.s
    fnu = u.erg / u.Hz / u.cm ** 2 / u.s

    a = flambda.to(fnu, 1, u.spectral_density(u.AA, 3500))
    assert_allclose(a, 4.086160166177361e-12)


def test_spectraldensity3():

    # Define F_nu in Jy
    f_nu = u.Jy

    # Convert to ergs / cm^2 / s / Hz
    assert_allclose(f_nu.to(u.erg / u.cm ** 2 / u.s / u.Hz, 1.), 1.e-23, 10)

    # Convert to ergs / cm^2 / s at 10 Ghz
    assert_allclose(f_nu.to(u.erg / u.cm ** 2 / u.s, 1.,
                    equivalencies=u.spectral_density(u.GHz, 10)), 1.e-13, 10)

    # Convert to ergs / cm^2 / s / micron at 1 Ghz
    assert_allclose(f_nu.to(u.erg / u.cm ** 2 / u.s / u.micron, 1.,
                    equivalencies=u.spectral_density(u.Hz, 1.e9)),
                    3.335640951981521e-20, 10)

    # Define F_lambda in ergs / cm^2 / s / micron
    f_lambda = u.erg / u.cm ** 2 / u.s / u.micron

    # Convert to Jy at 1 Ghz
    assert_allclose(f_lambda.to(u.Jy, 1.,
                    equivalencies=u.spectral_density(u.Hz, 1.e9)),
                    1. / 3.335640951981521e-20, 10)

    # Convert to ergs / cm^2 / s at 10 microns
    assert_allclose(f_lambda.to(u.erg / u.cm ** 2 / u.s, 1.,
                    equivalencies=u.spectral_density(u.micron, 10.)), 10., 10)


def test_equivalent_units():
    units = u.g.find_equivalent_units()
    units_set = set(units)
    match = set(
        [u.M_e, u.M_p, u.g, u.kg, u.lb, u.oz,
         u.solMass, u.t, u.ton, u.u])
    assert units_set == match

    r = repr(units)
    assert r.count('\n') == len(units) + 2


def test_equivalent_units2():
    units = set(u.Hz.find_equivalent_units(u.spectral()))
    match = set(
        [u.AU, u.Angstrom, u.BTU, u.Hz, u.J, u.Ry, u.cal, u.cm, u.eV,
         u.erg, u.ft, u.inch, u.kcal, u.lyr, u.m, u.mi, u.micron,
         u.pc, u.solRad, u.yd, u.Bq, u.Ci, u.nmi])
    assert units == match


def test_trivial_equivalency():
    assert u.m.to(u.kg, equivalencies=[(u.m, u.kg)]) == 1.0


def test_invalid_equivalency():
    with pytest.raises(ValueError):
        u.m.to(u.kg, equivalencies=[(u.m,)])

    with pytest.raises(ValueError):
        u.m.to(u.kg, equivalencies=[(u.m, 5.0)])


def test_irrelevant_equivalency():
    with pytest.raises(u.UnitsException):
        u.m.to(u.kg, equivalencies=[(u.m, u.l)])

def test_brightness_temperature():
    omega_B = np.pi*(50*u.arcsec)**2
    nu = u.GHz * 5
    tb = 7.05258885885 * u.K
    np.testing.assert_almost_equal(tb.value,
                                   (1*u.Jy).to(u.K,
                                               equivalencies=u.brightness_temperature(omega_B, nu)).value)
    np.testing.assert_almost_equal(1.0,
                                   tb.to(u.Jy,
                                         equivalencies=u.brightness_temperature(omega_B, nu)).value)
