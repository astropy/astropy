# Licensed under a 3-clause BSD style license - see LICENSE.rst

from StringIO import StringIO

import numpy as np

from .. import core, funcs
from ...tests.helper import pytest
from ... import units as u

try:
    import scipy  # pylint: disable=W0611
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True


def setup_function(function):
    # Make sure that tests don't affect default cosmology
    core.set_current('no_default')


def teardown_function(function):
    # Make sure that tests don't affect default cosmology
    core.set_current('no_default')


def test_basic():
    cosmo = core.FlatLambdaCDM(H0=70, Om0=0.27, Tcmb0=2.0, Neff=3.04)
    assert np.allclose(cosmo.Om0, 0.27)
    assert np.allclose(cosmo.Ode0, 0.729975, rtol=1e-4)
    assert np.allclose(cosmo.Ogamma0, 1.463285e-5, rtol=1e-4)
    assert np.allclose(cosmo.Onu0, 1.01026e-5, rtol=1e-4)
    assert np.allclose(cosmo.Ok0, 0.0)
    assert np.allclose(cosmo.Om0 + cosmo.Ode0 + cosmo.Ogamma0 + cosmo.Onu0,
                       1.0, rtol=1e-6)
    assert np.allclose(cosmo.Om(1) + cosmo.Ode(1) + cosmo.Ogamma(1) +
                       cosmo.Onu(1), 1.0, rtol=1e-6)
    assert np.allclose(cosmo.Tcmb0.value, 2.0)
    assert np.allclose(cosmo.Tnu0.value, 1.4275317, rtol=1e-5)
    assert np.allclose(cosmo.Neff, 3.04)
    assert np.allclose(cosmo.h, 0.7)
    assert np.allclose(cosmo.H0.value, 70.0)

    # Make sure setting them as quantities gives the same results
    H0 = u.Quantity(70, u.km / (u.s * u.Mpc))
    T = u.Quantity(2.0, u.K)
    cosmo = core.FlatLambdaCDM(H0=H0, Om0=0.27, Tcmb0=T, Neff=3.04)
    assert np.allclose(cosmo.Om0, 0.27)
    assert np.allclose(cosmo.Ode0, 0.729975, rtol=1e-4)
    assert np.allclose(cosmo.Ogamma0, 1.463285e-5, rtol=1e-4)
    assert np.allclose(cosmo.Onu0, 1.01026e-5, rtol=1e-4)
    assert np.allclose(cosmo.Ok0, 0.0)
    assert np.allclose(cosmo.Om0 + cosmo.Ode0 + cosmo.Ogamma0 + cosmo.Onu0,
                       1.0, rtol=1e-6)
    assert np.allclose(cosmo.Om(1) + cosmo.Ode(1) + cosmo.Ogamma(1) +
                       cosmo.Onu(1), 1.0, rtol=1e-6)
    assert np.allclose(cosmo.Tcmb0.value, 2.0)
    assert np.allclose(cosmo.Tnu0.value, 1.4275317, rtol=1e-5)
    assert np.allclose(cosmo.Neff, 3.04)
    assert np.allclose(cosmo.h, 0.7)
    assert np.allclose(cosmo.H0.value, 70.0)


@pytest.mark.skipif('not HAS_SCIPY')
def test_units():
    """ Test if the right units are being returned"""

    cosmo = core.FlatLambdaCDM(H0=70, Om0=0.27, Tcmb0=2.0)
    assert cosmo.comoving_distance(1.0).unit == u.Mpc
    assert cosmo.angular_diameter_distance(1.0).unit == u.Mpc
    assert cosmo.angular_diameter_distance_z1z2(1.0, 2.0).unit == u.Mpc
    assert cosmo.comoving_distance(1.0).unit == u.Mpc
    assert cosmo.luminosity_distance(1.0).unit == u.Mpc
    assert cosmo.lookback_time(1.0).unit == u.Gyr
    assert cosmo.H0.unit == u.km / u.Mpc / u.s
    assert cosmo.H(1.0).unit == u.km / u.Mpc / u.s
    assert cosmo.Tcmb0.unit == u.K
    assert cosmo.Tcmb(1.0).unit == u.K
    assert cosmo.Tcmb([0.0, 1.0]).unit == u.K
    assert cosmo.Tnu0.unit == u.K
    assert cosmo.Tnu(1.0).unit == u.K
    assert cosmo.Tnu([0.0, 1.0]).unit == u.K
    assert cosmo.arcsec_per_kpc_comoving(1.0).unit == u.arcsec / u.kpc
    assert cosmo.arcsec_per_kpc_proper(1.0).unit == u.arcsec / u.kpc
    assert cosmo.kpc_comoving_per_arcmin(1.0).unit == u.kpc / u.arcmin
    assert cosmo.kpc_proper_per_arcmin(1.0).unit == u.kpc / u.arcmin
    assert cosmo.critical_density(1.0).unit == u.g / u.cm ** 3
    assert cosmo.comoving_volume(1.0).unit == u.Mpc ** 3
    assert cosmo.age(1.0).unit == u.Gyr
    assert cosmo.distmod(1.0).unit == u.mag

def test_repr():
    """ Test string representation of built in classes"""
    cosmo = core.LambdaCDM(70, 0.3, 0.5)
    expected = 'LambdaCDM(H0=70 km / (Mpc s), Om0=0.3, Ode0=0.5, Tcmb0=2.725 K, Neff=3.04, m_nu=[ 0.  0.  0.] eV)'
    assert "{0:s}".format(str(cosmo)) == expected

    cosmo = core.LambdaCDM(70, 0.3, 0.5, m_nu=u.Quantity(0.01, u.eV))
    expected = 'LambdaCDM(H0=70 km / (Mpc s), Om0=0.3, Ode0=0.5, Tcmb0=2.725 K, Neff=3.04, m_nu=[ 0.01  0.01  0.01] eV)'
    assert "{0:s}".format(str(cosmo)) == expected

    cosmo = core.FlatLambdaCDM(50.0, 0.27)
    expected = 'FlatLambdaCDM(H0=50 km / (Mpc s), Om0=0.27, Tcmb0=2.725 K, Neff=3.04, m_nu=[ 0.  0.  0.] eV)'
    assert "{0:s}".format(str(cosmo)) == expected

    cosmo = core.wCDM(60.0, 0.27, 0.6, w0=-0.8, name='test1')
    expected = 'wCDM(name="test1", H0=60 km / (Mpc s), Om0=0.27, Ode0=0.6, w0=-0.8, Tcmb0=2.725 K, Neff=3.04, m_nu=[ 0.  0.  0.] eV)'
    assert "{0:s}".format(str(cosmo)) == expected

    cosmo = core.FlatwCDM(65.0, 0.27, w0=-0.6, name='test2')
    expected = 'FlatwCDM(name="test2", H0=65 km / (Mpc s), Om0=0.27, w0=-0.6, Tcmb0=2.725 K, Neff=3.04, m_nu=[ 0.  0.  0.] eV)'
    assert "{0:s}".format(str(cosmo)) == expected

    cosmo = core.w0waCDM(60.0, 0.25, 0.4, w0=-0.6, wa=0.1, name='test3')
    expected = 'w0waCDM(name="test3", H0=60 km / (Mpc s), Om0=0.25, Ode0=0.4, w0=-0.6, wa=0.1, Tcmb0=2.725 K, Neff=3.04, m_nu=[ 0.  0.  0.] eV)'
    assert "{0:s}".format(str(cosmo)) == expected

    cosmo = core.Flatw0waCDM(55.0, 0.35, w0=-0.9, wa=-0.2, name='test4')
    expected = 'Flatw0waCDM(name="test4", H0=55 km / (Mpc s), Om0=0.35, w0=-0.9, Tcmb0=2.725 K, Neff=3.04, m_nu=[ 0.  0.  0.] eV)'
    assert "{0:s}".format(str(cosmo)) == expected

    cosmo = core.wpwaCDM(50.0, 0.3, 0.3, wp=-0.9, wa=-0.2, zp=0.3, name='test5')
    expected = 'wpwaCDM(name="test5", H0=50 km / (Mpc s), Om0=0.3, Ode0=0.3, wp=-0.9, wa=-0.2, zp=0.3, Tcmb0=2.725 K, Neff=3.04, m_nu=[ 0.  0.  0.] eV)'
    assert "{0:s}".format(str(cosmo)) == expected

    cosmo = core.w0wzCDM(55.0, 0.4, 0.8, w0=-1.05, wz=-0.2,
                         m_nu=u.Quantity([0.001, 0.01, 0.015], u.eV))
    expected = 'w0wzCDM(H0=55 km / (Mpc s), Om0=0.4, Ode0=0.8, w0=-1.05, wz=-0.2 Tcmb0=2.725 K, Neff=3.04, m_nu=[ 0.001  0.01   0.015] eV)'
    assert "{0:s}".format(str(cosmo)) == expected

@pytest.mark.skipif('not HAS_SCIPY')
def test_flat_z1():
    """ Test a flat cosmology at z=1 against several other on-line
    calculators.
    """
    cosmo = core.FlatLambdaCDM(H0=70, Om0=0.27, Tcmb0=0.0)
    z = 1

    # Test values were taken from the following web cosmology
    # calculators on 27th Feb 2012:

    # Wright: http://www.astro.ucla.edu/~wright/CosmoCalc.html
    #         (http://adsabs.harvard.edu/abs/2006PASP..118.1711W)
    # Kempner: http://www.kempner.net/cosmic.php
    # iCosmos: http://www.icosmos.co.uk/index.html

    # The order of values below is Wright, Kempner, iCosmos'
    assert np.allclose(cosmo.comoving_distance(z).value,
                       [3364.5, 3364.8, 3364.7988], rtol=1e-4)
    assert np.allclose(cosmo.angular_diameter_distance(z).value,
                       [1682.3, 1682.4, 1682.3994], rtol=1e-4)
    assert np.allclose(cosmo.luminosity_distance(z).value,
                       [6729.2, 6729.6, 6729.5976], rtol=1e-4)
    assert np.allclose(cosmo.lookback_time(z).value,
                       [7.841, 7.84178, 7.843], rtol=1e-3)


# This class is to test whether the routines work correctly
# if one only overloads w(z)
class test_cos_sub(core.FLRW):

    def __init__(self):
        core.FLRW.__init__(self, 70.0, 0.27, 0.73, Tcmb0=0.0, name="test_cos")
        self._w0 = -0.9

    def w(self, z):
        return self._w0 * np.ones_like(z)


@pytest.mark.skipif('not HAS_SCIPY')
def test_de_subclass():
    # This is the comparison object
    z = [0.2, 0.4, 0.6, 0.9]
    cosmo = core.wCDM(H0=70, Om0=0.27, Ode0=0.73, w0=-0.9, Tcmb0=0.0)
    # Values taken from Ned Wrights advanced cosmo calcluator, Aug 17 2012
    assert np.allclose(cosmo.luminosity_distance(z).value,
                       [975.5, 2158.2, 3507.3, 5773.1], rtol=1e-3)
    # Now try the subclass that only gives w(z)
    cosmo = test_cos_sub()
    assert np.allclose(cosmo.luminosity_distance(z).value,
                       [975.5, 2158.2, 3507.3, 5773.1], rtol=1e-3)
    # Test efunc
    assert np.allclose(cosmo.efunc(1.0), 1.7489240754, rtol=1e-5)
    assert np.allclose(cosmo.efunc([0.5, 1.0]),
                       [1.31744953, 1.7489240754], rtol=1e-5)
    assert np.allclose(cosmo.inv_efunc([0.5, 1.0]),
                       [0.75904236, 0.57178011], rtol=1e-5)
    # Test de_density_scale
    assert np.allclose(cosmo.de_density_scale(1.0), 1.23114444, rtol=1e-4)
    assert np.allclose(cosmo.de_density_scale([0.5, 1.0]),
                       [1.12934694, 1.23114444], rtol=1e-4)


@pytest.mark.skipif('not HAS_SCIPY')
def test_varyde_lumdist_mathematica():
    """Tests a few varying dark energy EOS models against a mathematica
    computation"""

    # w0wa models
    z = np.array([0.2, 0.4, 0.9, 1.2])
    cosmo = core.w0waCDM(H0=70, Om0=0.2, Ode0=0.8, w0=-1.1, wa=0.2, Tcmb0=0.0)
    assert np.allclose(cosmo.luminosity_distance(z).value,
                       [1004.0, 2268.62, 6265.76, 9061.84], rtol=1e-4)
    assert np.allclose(cosmo.de_density_scale(0.0), 1.0, rtol=1e-5)
    assert np.allclose(cosmo.de_density_scale([0.0, 0.5, 1.5]),
                       [1.0, 0.9246310669529021, 0.9184087000251957])

    cosmo = core.w0waCDM(H0=70, Om0=0.3, Ode0=0.7, w0=-0.9, wa=0.0, Tcmb0=0.0)
    assert np.allclose(cosmo.luminosity_distance(z).value,
                       [971.667, 2141.67, 5685.96, 8107.41], rtol=1e-4)
    cosmo = core.w0waCDM(H0=70, Om0=0.3, Ode0=0.7, w0=-0.9, wa=-0.5, Tcmb0=0.0)
    assert np.allclose(cosmo.luminosity_distance(z).value,
                       [974.087, 2157.08, 5783.92, 8274.08], rtol=1e-4)

    # wpwa models
    cosmo = core.wpwaCDM(H0=70, Om0=0.2, Ode0=0.8, wp=-1.1, wa=0.2, zp=0.5,
                         Tcmb0=0.0)
    assert np.allclose(cosmo.luminosity_distance(z).value,
                       [1010.81, 2294.45, 6369.45, 9218.95], rtol=1e-4)
    cosmo = core.wpwaCDM(H0=70, Om0=0.2, Ode0=0.8, wp=-1.1, wa=0.2, zp=0.9,
                         Tcmb0=0.0)
    assert np.allclose(cosmo.luminosity_distance(z).value,
                       [1013.68, 2305.3, 6412.37, 9283.33], rtol=1e-4)


@pytest.mark.skipif('not HAS_SCIPY')
def test_omatter():
    # Test Om evolution
    tcos = core.FlatLambdaCDM(70.0, 0.3)
    assert np.allclose(tcos.Om0, 0.3)
    assert np.allclose(tcos.H0.value, 70.0)
    assert np.allclose(tcos.Om(0), 0.3)
    z = np.array([0.0, 0.5, 1.0, 2.0])
    assert np.allclose(tcos.Om(z), [0.3, 0.59112134, 0.77387435, 0.91974179],
                       rtol=1e-4)


@pytest.mark.skipif('not HAS_SCIPY')
def test_ocurv():
    # Test Ok evolution
    # Flat, boring case
    tcos = core.FlatLambdaCDM(70.0, 0.3)
    assert np.allclose(tcos.Ok0, 0.0)
    assert np.allclose(tcos.Ok(0), 0.0)
    z = np.array([0.0, 0.5, 1.0, 2.0])
    assert np.allclose(tcos.Ok(z), [0.0, 0.0, 0.0, 0.0],
                       rtol=1e-6)

    # Not flat
    tcos = core.LambdaCDM(70.0, 0.3, 0.5, Tcmb0=u.Quantity(0.0, u.K))
    assert np.allclose(tcos.Ok0, 0.2)
    assert np.allclose(tcos.Ok(0), 0.2)
    assert np.allclose(tcos.Ok(z), [0.2, 0.22929936, 0.21621622, 0.17307692],
                       rtol=1e-4)

    # Test the sum; note that Ogamma/Onu are 0
    assert np.allclose(tcos.Ok(z) + tcos.Om(z) + tcos.Ode(z),
                       [1.0, 1.0, 1.0, 1.0], rtol=1e-5)


@pytest.mark.skipif('not HAS_SCIPY')
def test_ode():
    # Test Ode evolution, turn off neutrinos, cmb
    tcos = core.FlatLambdaCDM(70.0, 0.3, Tcmb0=0)
    assert np.allclose(tcos.Ode0, 0.7)
    assert np.allclose(tcos.Ode(0), 0.7)
    z = np.array([0.0, 0.5, 1.0, 2.0])
    assert np.allclose(tcos.Ode(z), [0.7, 0.408759, 0.2258065, 0.07954545],
                       rtol=1e-5)


@pytest.mark.skipif('not HAS_SCIPY')
def test_ogamma():
    """Tests the effects of changing the temperature of the CMB"""

    # Tested against Ned Wright's advanced cosmology calculator,
    # Sep 7 2012.  The accuracy of our comparision is limited by
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
    cosmo = core.FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=0, Neff=3)
    assert np.allclose(cosmo.angular_diameter_distance(z).value,
                       [1651.9, 858.2, 26.855, 13.642], rtol=5e-4)
    cosmo = core.FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725, Neff=3)
    assert np.allclose(cosmo.angular_diameter_distance(z).value,
                       [1651.8, 857.9, 26.767, 13.582], rtol=5e-4)
    cosmo = core.FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=4.0, Neff=3)
    assert np.allclose(cosmo.angular_diameter_distance(z).value,
                       [1651.4, 856.6, 26.489, 13.405], rtol=5e-4)

    # Next compare with doing the integral numerically in Mathematica,
    # which allows more precision in the test.  It is at least as
    # good as 0.01%, possibly better
    cosmo = core.FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=0, Neff=3.04)
    assert np.allclose(cosmo.angular_diameter_distance(z).value,
                       [1651.91, 858.205, 26.8586, 13.6469], rtol=1e-5)
    cosmo = core.FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725, Neff=3.04)
    assert np.allclose(cosmo.angular_diameter_distance(z).value,
                       [1651.76, 857.817, 26.7688, 13.5841], rtol=1e-5)
    cosmo = core.FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=4.0, Neff=3.04)
    assert np.allclose(cosmo.angular_diameter_distance(z).value,
                       [1651.21, 856.411, 26.4845, 13.4028], rtol=1e-5)

    # Just to be really sure, we also do a version where the integral
    # is analytic, which is a Ode = 0 flat universe.  In this case
    # Integrate(1/E(x),{x,0,z}) = 2 ( sqrt((1+Or z)/(1+z)) - 1 )/(Or - 1)
    # Recall that c/H0 * Integrate(1/E) is FLRW.comoving_distance.
    Ogamma0h2 = 4 * 5.670373e-8 / 299792458.0 ** 3 * 2.725 ** 4 / 1.87837e-26
    Onu0h2 = Ogamma0h2 * 7.0 / 8.0 * (4.0 / 11.0) ** (4.0 / 3.0) * 3.04
    Or0 = (Ogamma0h2 + Onu0h2) / 0.7 ** 2
    Om0 = 1.0 - Or0
    hubdis = 299792.458 / 70.0
    cosmo = core.FlatLambdaCDM(H0=70, Om0=Om0, Tcmb0=2.725, Neff=3.04)
    targvals = 2.0 * hubdis * \
        (np.sqrt((1.0 + Or0 * z) / (1.0 + z)) - 1.0) / (Or0 - 1.0)
    assert np.allclose(cosmo.comoving_distance(z).value, targvals, rtol=1e-5)

    # Try Tcmb0 = 4
    Or0 *= (4.0 / 2.725) ** 4
    Om0 = 1.0 - Or0
    cosmo = core.FlatLambdaCDM(H0=70, Om0=Om0, Tcmb0=4.0, Neff=3.04)
    targvals = 2.0 * hubdis * \
        (np.sqrt((1.0 + Or0 * z) / (1.0 + z)) - 1.0) / (Or0 - 1.0)
    assert np.allclose(cosmo.comoving_distance(z).value, targvals, rtol=1e-5)


@pytest.mark.skipif('not HAS_SCIPY')
def test_tcmb():
    cosmo = core.FlatLambdaCDM(70.4, 0.272, Tcmb0=3.0)
    assert np.allclose(cosmo.Tcmb0.value, 3.0)
    assert np.allclose(cosmo.Tcmb(2).value, 9.0)
    z = [0.0, 1.0, 2.0, 3.0, 9.0]
    assert np.allclose(cosmo.Tcmb(z).value,
                       [3.0, 6.0, 9.0, 12.0, 30.0], rtol=1e-6)


@pytest.mark.skipif('not HAS_SCIPY')
def test_tnu():
    cosmo = core.FlatLambdaCDM(70.4, 0.272, Tcmb0=3.0)
    assert np.allclose(cosmo.Tnu0.value, 2.1412975665108247, rtol=1e-6)
    assert np.allclose(cosmo.Tnu(2).value, 6.423892699532474, rtol=1e-6)
    z = [0.0, 1.0, 2.0, 3.0]
    assert np.allclose(cosmo.Tnu(z), [2.14129757, 4.28259513,
                                      6.4238927, 8.56519027], rtol=1e-6)


def test_efunc_vs_invefunc():
    # Test that efunc and inv_efunc give the same values
    z0 = 0.5
    z = np.array([0.5, 1.0, 2.0, 5.0])

    # Below are the 'standard' included cosmologies
    # We do the non-standard case in test_efunc_vs_invefunc_flrw,
    # since it requires scipy
    cosmo = core.LambdaCDM(70, 0.3, 0.5)
    assert np.allclose(cosmo.efunc(z0), 1.0 / cosmo.inv_efunc(z0))
    assert np.allclose(cosmo.efunc(z), 1.0 / cosmo.inv_efunc(z))
    cosmo = core.LambdaCDM(70, 0.3, 0.5, m_nu=u.Quantity(0.01, u.eV))
    assert np.allclose(cosmo.efunc(z0), 1.0 / cosmo.inv_efunc(z0))
    assert np.allclose(cosmo.efunc(z), 1.0 / cosmo.inv_efunc(z))
    cosmo = core.FlatLambdaCDM(50.0, 0.27)
    assert np.allclose(cosmo.efunc(z0), 1.0 / cosmo.inv_efunc(z0))
    assert np.allclose(cosmo.efunc(z), 1.0 / cosmo.inv_efunc(z))
    cosmo = core.wCDM(60.0, 0.27, 0.6, w0=-0.8)
    assert np.allclose(cosmo.efunc(z0), 1.0 / cosmo.inv_efunc(z0))
    assert np.allclose(cosmo.efunc(z), 1.0 / cosmo.inv_efunc(z))
    cosmo = core.FlatwCDM(65.0, 0.27, w0=-0.6)
    assert np.allclose(cosmo.efunc(z0), 1.0 / cosmo.inv_efunc(z0))
    assert np.allclose(cosmo.efunc(z), 1.0 / cosmo.inv_efunc(z))
    cosmo = core.w0waCDM(60.0, 0.25, 0.4, w0=-0.6, wa=0.1)
    assert np.allclose(cosmo.efunc(z0), 1.0 / cosmo.inv_efunc(z0))
    assert np.allclose(cosmo.efunc(z), 1.0 / cosmo.inv_efunc(z))
    cosmo = core.Flatw0waCDM(55.0, 0.35, w0=-0.9, wa=-0.2)
    assert np.allclose(cosmo.efunc(z0), 1.0 / cosmo.inv_efunc(z0))
    assert np.allclose(cosmo.efunc(z), 1.0 / cosmo.inv_efunc(z))
    cosmo = core.wpwaCDM(50.0, 0.3, 0.3, wp=-0.9, wa=-0.2, zp=0.3)
    assert np.allclose(cosmo.efunc(z0), 1.0 / cosmo.inv_efunc(z0))
    assert np.allclose(cosmo.efunc(z), 1.0 / cosmo.inv_efunc(z))
    cosmo = core.w0wzCDM(55.0, 0.4, 0.8, w0=-1.05, wz=-0.2)
    assert np.allclose(cosmo.efunc(z0), 1.0 / cosmo.inv_efunc(z0))
    assert np.allclose(cosmo.efunc(z), 1.0 / cosmo.inv_efunc(z))


@pytest.mark.skipif('not HAS_SCIPY')
def test_efunc_vs_invefunc_flrw():
    z0 = 0.5
    z = np.array([0.5, 1.0, 2.0, 5.0])

    # FLRW is abstract, so requires test_cos_sub defined earlier
    # This requires scipy, unlike the built-ins
    cosmo = test_cos_sub()
    assert np.allclose(cosmo.efunc(z0), 1.0 / cosmo.inv_efunc(z0))
    assert np.allclose(cosmo.efunc(z), 1.0 / cosmo.inv_efunc(z))


@pytest.mark.skipif('not HAS_SCIPY')
def test_kpc_methods():
    cosmo = core.FlatLambdaCDM(70.4, 0.272, Tcmb0=0.0)
    assert np.allclose(cosmo.arcsec_per_kpc_comoving(3).value, 0.0317179)
    assert np.allclose(cosmo.arcsec_per_kpc_proper(3).value, 0.1268716668)
    assert np.allclose(cosmo.kpc_comoving_per_arcmin(3).value, 1891.6753126)
    assert np.allclose(cosmo.kpc_proper_per_arcmin(3).value, 472.918828)


@pytest.mark.skipif('not HAS_SCIPY')
def test_convenience():
    # these are all for WMAP7 with Tcmb = 0
    tcos = core.FlatLambdaCDM(70.4, 0.272, Tcmb0=0.0)
    core.set_current(tcos)

    # scalars
    assert np.allclose(funcs.arcsec_per_kpc_comoving(3).value, 0.0317179)
    assert funcs.arcsec_per_kpc_comoving(3).unit == u.arcsec / u.kpc
    assert np.allclose(funcs.arcsec_per_kpc_proper(3).value, 0.1268716668)
    assert funcs.arcsec_per_kpc_proper(3).unit == u.arcsec / u.kpc
    assert np.allclose(funcs.kpc_comoving_per_arcmin(3).value, 1891.6753126)
    assert funcs.kpc_comoving_per_arcmin(3).unit == u.kpc / u.arcmin
    assert np.allclose(funcs.kpc_proper_per_arcmin(3).value, 472.918828)
    assert funcs.kpc_proper_per_arcmin(3).unit == u.kpc / u.arcmin
    assert np.allclose(funcs.distmod(3).value, 47.075902)
    assert funcs.distmod(3).unit == u.mag
    assert np.allclose(funcs.H(3).value, 299.80813491298068)
    assert funcs.H(3).unit == u.km / (u.Mpc * u.s)
    assert np.allclose(funcs.scale_factor(3), 0.25)
    assert np.allclose(funcs.scale_factor([3, 4]), [0.25, 0.2])
    assert np.allclose(funcs.critical_density(3).value, 1.6884621680232328e-28)
    assert funcs.critical_density(3).unit == u.g / u.cm ** 3
    assert np.allclose(funcs.lookback_time(3).value, 11.555469926558361)
    assert funcs.lookback_time(3).unit == u.Gyr
    assert np.allclose(funcs.lookback_time([3, 4]).value,
                       [11.555469927, 12.17718555], rtol=1e-5)
    assert np.allclose(funcs.comoving_distance(3).value, 6503.100697385924)
    assert funcs.comoving_distance(3).unit == u.Mpc
    assert np.allclose(funcs.angular_diameter_distance(3).value,
                       1625.775174346481)
    assert funcs.angular_diameter_distance(3).unit == u.Mpc
    assert np.allclose(funcs.luminosity_distance(3).value, 26012.402789543696)
    assert funcs.luminosity_distance(3).unit == u.Mpc
    assert np.allclose(funcs.age(3).value, 2.20407604076062)
    assert funcs.age(3).unit == u.Gyr
    assert np.allclose(funcs.comoving_volume(3).value, 1151993546079.626)
    assert funcs.comoving_volume(3).unit == u.Mpc**3

    # arrays
    assert np.allclose(funcs.arcsec_per_kpc_comoving([0.1, 0.5]).value,
                       [0.4946986, 0.10876163])
    assert np.allclose(funcs.arcsec_per_kpc_proper([0.1, 0.5]).value,
                       [0.54416846354697479, 0.16314245192751084])
    assert np.allclose(funcs.kpc_comoving_per_arcmin([0.1, 0.5]).value,
                       [121.2859701, 551.66511804])
    assert np.allclose(funcs.kpc_proper_per_arcmin([0.1, 0.5]).value,
                       [110.25997282, 367.77674536])
    assert np.allclose(funcs.distmod([0.1, 0.5]).value,
                       [38.30738567, 42.27020333])


@pytest.mark.skipif('not HAS_SCIPY')
def test_comoving_volume():

    c_flat = core.LambdaCDM(H0=70, Om0=0.27, Ode0=0.73, Tcmb0=0.0)
    c_open = core.LambdaCDM(H0=70, Om0=0.27, Ode0=0.0, Tcmb0=0.0)
    c_closed = core.LambdaCDM(H0=70, Om0=2, Ode0=0.0, Tcmb0=0.0)

    # test against ned wright's calculator (cubic Gpc)
    redshifts = np.array([0.5, 1, 2, 3, 5, 9])
    wright_flat = np.array([29.123, 159.529, 630.427, 1178.531, 2181.485,
                            3654.802]) * 1e9  # convert to Mpc**3
    wright_open = np.array([20.501, 99.019, 380.278, 747.049, 1558.363,
                            3123.814]) * 1e9
    wright_closed = np.array([12.619, 44.708, 114.904, 173.709, 258.82,
                              358.992]) * 1e9
    # The wright calculator isn't very accurate, so we use a rather
    # modest precision
    assert np.allclose(c_flat.comoving_volume(redshifts).value, wright_flat,
                       rtol=1e-2)
    assert np.allclose(c_open.comoving_volume(redshifts).value, wright_open,
                       rtol=1e-2)
    assert np.allclose(c_closed.comoving_volume(redshifts).value, wright_closed,
                       rtol=1e-2)


@pytest.mark.skipif('not HAS_SCIPY')
def test_flat_open_closed_icosmo():
    """ Test against the tabulated values generated from icosmo.org
    with three example cosmologies (flat, open and closed).
    """

    cosmo_flat = """\
# from icosmo (icosmo.org)
# Om 0.3 w -1 h 0.7 Ol 0.7
# z     comoving_transvers_dist   angular_diameter_dist  luminosity_dist
       0.0000000       0.0000000       0.0000000         0.0000000
      0.16250000       669.77536       576.15085         778.61386
      0.32500000       1285.5964       970.26143         1703.4152
      0.50000000       1888.6254       1259.0836         2832.9381
      0.66250000       2395.5489       1440.9317         3982.6000
      0.82500000       2855.5732       1564.6976         5211.4210
       1.0000000       3303.8288       1651.9144         6607.6577
       1.1625000       3681.1867       1702.2829         7960.5663
       1.3250000       4025.5229       1731.4077         9359.3408
       1.5000000       4363.8558       1745.5423         10909.640
       1.6625000       4651.4830       1747.0359         12384.573
       1.8250000       4916.5970       1740.3883         13889.387
       2.0000000       5179.8621       1726.6207         15539.586
       2.1625000       5406.0204       1709.4136         17096.540
       2.3250000       5616.5075       1689.1752         18674.888
       2.5000000       5827.5418       1665.0120         20396.396
       2.6625000       6010.4886       1641.0890         22013.414
       2.8250000       6182.1688       1616.2533         23646.796
       3.0000000       6355.6855       1588.9214         25422.742
       3.1625000       6507.2491       1563.3031         27086.425
       3.3250000       6650.4520       1537.6768         28763.205
       3.5000000       6796.1499       1510.2555         30582.674
       3.6625000       6924.2096       1485.0852         32284.127
       3.8250000       7045.8876       1460.2876         33996.408
       4.0000000       7170.3664       1434.0733         35851.832
       4.1625000       7280.3423       1410.2358         37584.767
       4.3250000       7385.3277       1386.9160         39326.870
       4.5000000       7493.2222       1362.4040         41212.722
       4.6625000       7588.9589       1340.2135         42972.480
"""

    cosmo_open = """\
# from icosmo (icosmo.org)
# Om 0.3 w -1 h 0.7 Ol 0.1
# z     comoving_transvers_dist   angular_diameter_dist  luminosity_dist
       0.0000000       0.0000000       0.0000000       0.0000000
      0.16250000       643.08185       553.18868       747.58265
      0.32500000       1200.9858       906.40441       1591.3062
      0.50000000       1731.6262       1154.4175       2597.4393
      0.66250000       2174.3252       1307.8648       3614.8157
      0.82500000       2578.7616       1413.0201       4706.2399
       1.0000000       2979.3460       1489.6730       5958.6920
       1.1625000       3324.2002       1537.2024       7188.5829
       1.3250000       3646.8432       1568.5347       8478.9104
       1.5000000       3972.8407       1589.1363       9932.1017
       1.6625000       4258.1131       1599.2913       11337.226
       1.8250000       4528.5346       1603.0211       12793.110
       2.0000000       4804.9314       1601.6438       14414.794
       2.1625000       5049.2007       1596.5852       15968.097
       2.3250000       5282.6693       1588.7727       17564.875
       2.5000000       5523.0914       1578.0261       19330.820
       2.6625000       5736.9813       1566.4113       21011.694
       2.8250000       5942.5803       1553.6158       22730.370
       3.0000000       6155.4289       1538.8572       24621.716
       3.1625000       6345.6997       1524.4924       26413.975
       3.3250000       6529.3655       1509.6799       28239.506
       3.5000000       6720.2676       1493.3928       30241.204
       3.6625000       6891.5474       1478.0799       32131.840
       3.8250000       7057.4213       1462.6780       34052.058
       4.0000000       7230.3723       1446.0745       36151.862
       4.1625000       7385.9998       1430.7021       38130.224
       4.3250000       7537.1112       1415.4199       40135.117
       4.5000000       7695.0718       1399.1040       42322.895
       4.6625000       7837.5510       1384.1150       44380.133
"""

    cosmo_closed = """\
# from icosmo (icosmo.org)
# Om 2 w -1 h 0.7 Ol 0.1
# z     comoving_transvers_dist   angular_diameter_dist  luminosity_dist
       0.0000000       0.0000000       0.0000000       0.0000000
      0.16250000       601.80160       517.67879       699.59436
      0.32500000       1057.9502       798.45297       1401.7840
      0.50000000       1438.2161       958.81076       2157.3242
      0.66250000       1718.6778       1033.7912       2857.3019
      0.82500000       1948.2400       1067.5288       3555.5381
       1.0000000       2152.7954       1076.3977       4305.5908
       1.1625000       2312.3427       1069.2914       5000.4410
       1.3250000       2448.9755       1053.3228       5693.8681
       1.5000000       2575.6795       1030.2718       6439.1988
       1.6625000       2677.9671       1005.8092       7130.0873
       1.8250000       2768.1157       979.86398       7819.9270
       2.0000000       2853.9222       951.30739       8561.7665
       2.1625000       2924.8116       924.84161       9249.7167
       2.3250000       2988.5333       898.80701       9936.8732
       2.5000000       3050.3065       871.51614       10676.073
       2.6625000       3102.1909       847.01459       11361.774
       2.8250000       3149.5043       823.39982       12046.854
       3.0000000       3195.9966       798.99915       12783.986
       3.1625000       3235.5334       777.30533       13467.908
       3.3250000       3271.9832       756.52790       14151.327
       3.5000000       3308.1758       735.15017       14886.791
       3.6625000       3339.2521       716.19347       15569.263
       3.8250000       3368.1489       698.06195       16251.319
       4.0000000       3397.0803       679.41605       16985.401
       4.1625000       3422.1142       662.87926       17666.664
       4.3250000       3445.5542       647.05243       18347.576
       4.5000000       3469.1805       630.76008       19080.493
       4.6625000       3489.7534       616.29199       19760.729
"""

    redshifts, dm, da, dl = np.loadtxt(StringIO(cosmo_flat), unpack=1)
    cosmo = core.LambdaCDM(H0=70, Om0=0.3, Ode0=0.70, Tcmb0=0.0)
    assert np.allclose(cosmo.comoving_transverse_distance(redshifts).value, dm)
    assert np.allclose(cosmo.angular_diameter_distance(redshifts).value, da)
    assert np.allclose(cosmo.luminosity_distance(redshifts).value, dl)

    redshifts, dm, da, dl = np.loadtxt(StringIO(cosmo_open), unpack=1)
    cosmo = core.LambdaCDM(H0=70, Om0=0.3, Ode0=0.1, Tcmb0=0.0)
    assert np.allclose(cosmo.comoving_transverse_distance(redshifts).value, dm)
    assert np.allclose(cosmo.angular_diameter_distance(redshifts).value, da)
    assert np.allclose(cosmo.luminosity_distance(redshifts).value, dl)

    redshifts, dm, da, dl = np.loadtxt(StringIO(cosmo_closed), unpack=1)
    cosmo = core.LambdaCDM(H0=70, Om0=2, Ode0=0.1, Tcmb0=0.0)
    assert np.allclose(cosmo.comoving_transverse_distance(redshifts).value, dm)
    assert np.allclose(cosmo.angular_diameter_distance(redshifts).value, da)
    assert np.allclose(cosmo.luminosity_distance(redshifts).value, dl)


def test_current():
    core.set_current('WMAP7')
    cosmo = core.get_current()
    assert cosmo == core.WMAP7
    core.set_current('WMAP5')
    assert core.get_current() == core.WMAP5
    core.set_current('WMAP9')
    assert core.get_current() == core.WMAP9
    core.set_current('Planck13')
    assert core.get_current() == core.Planck13
    core.set_current(cosmo)
    assert core.get_current() == cosmo


def test_wz():
    cosmo = core.LambdaCDM(H0=70, Om0=0.3, Ode0=0.70)
    assert np.allclose(cosmo.w([0.1, 0.2, 0.5, 1.5, 2.5, 11.5]),
                       [-1., -1, -1, -1, -1, -1])
    cosmo = core.wCDM(H0=70, Om0=0.3, Ode0=0.70, w0=-0.5)
    assert np.allclose(cosmo.w([0.1, 0.2, 0.5, 1.5, 2.5, 11.5]),
                       [-0.5, -0.5, -0.5, -0.5, -0.5, -0.5])
    cosmo = core.w0wzCDM(H0=70, Om0=0.3, Ode0=0.70, w0=-1, wz=0.5)
    assert np.allclose(cosmo.w([0.0, 0.5, 1.0, 1.5, 2.3]),
                       [-1.0, -0.75, -0.5, -0.25, 0.15])
    cosmo = core.w0waCDM(H0=70, Om0=0.3, Ode0=0.70, w0=-1, wa=-0.5)
    assert np.allclose(cosmo.w([0.0, 0.5, 1.0, 1.5, 2.3]),
                       [-1, -1.16666667, -1.25, -1.3, -1.34848485])
    cosmo = core.wpwaCDM(H0=70, Om0=0.3, Ode0=0.70, wp=-0.9,
                         wa=0.2, zp=0.5)
    assert np.allclose(cosmo.w([0.1, 0.2, 0.5, 1.5, 2.5, 11.5]),
                       [-0.94848485, -0.93333333, -0.9, -0.84666667,
                        -0.82380952, -0.78266667])


@pytest.mark.skipif('not HAS_SCIPY')
def test_de_densityscale():
    cosmo = core.LambdaCDM(H0=70, Om0=0.3, Ode0=0.70)
    z = np.array([0.1, 0.2, 0.5, 1.5, 2.5])
    assert np.allclose(cosmo.de_density_scale(z),
                       [1.0, 1.0, 1.0, 1.0, 1.0])
    cosmo = core.wCDM(H0=70, Om0=0.3, Ode0=0.60, w0=-0.5)
    assert np.allclose(cosmo.de_density_scale(z),
                       [1.15369, 1.31453, 1.83712, 3.95285, 6.5479],
                       rtol=1e-4)
    cosmo = core.w0wzCDM(H0=70, Om0=0.3, Ode0=0.50, w0=-1, wz=0.5)
    assert np.allclose(cosmo.de_density_scale(z),
                       [0.746048, 0.5635595, 0.25712378, 0.026664129,
                        0.0035916468], rtol=1e-4)
    cosmo = core.w0waCDM(H0=70, Om0=0.3, Ode0=0.70, w0=-1, wa=-0.5)
    assert np.allclose(cosmo.de_density_scale(z),
                       [0.9934201, 0.9767912, 0.897450,
                        0.622236, 0.4458753], rtol=1e-4)
    cosmo = core.wpwaCDM(H0=70, Om0=0.3, Ode0=0.70, wp=-0.9,
                         wa=0.2, zp=0.5)
    assert np.allclose(cosmo.de_density_scale(z),
                       [1.012246048, 1.0280102, 1.087439,
                        1.324988, 1.565746], rtol=1e-4)


@pytest.mark.skipif('not HAS_SCIPY')
def test_age():
    # WMAP7 but with Omega_relativisitic = 0
    tcos = core.FlatLambdaCDM(70.4, 0.272, Tcmb0=0.0)
    assert np.allclose(tcos.age([1, 5]).value, [5.97113193, 1.20553129])


@pytest.mark.skipif('not HAS_SCIPY')
def test_distmod():
    # WMAP7 but with Omega_relativisitic = 0
    tcos = core.FlatLambdaCDM(70.4, 0.272, Tcmb0=0.0)
    core.set_current(tcos)
    assert np.allclose(tcos.distmod([1, 5]).value, [44.124857, 48.40167258])
    assert np.allclose(funcs.distmod([1, 5], cosmo=tcos).value,
                       [44.124857, 48.40167258])


@pytest.mark.skipif('not HAS_SCIPY')
def test_critical_density():
    # WMAP7 but with Omega_relativisitic = 0
    tcos = core.FlatLambdaCDM(70.4, 0.272, Tcmb0=0.0)
    assert np.allclose(tcos.critical_density([1, 5]).value,
                       [2.70362491e-29, 5.53758986e-28])


@pytest.mark.skipif('not HAS_SCIPY')
def test_angular_diameter_distance_z1z2():
    tcos = core.FlatLambdaCDM(70.4, 0.272, Tcmb0=0.0)
    assert np.allclose(tcos.angular_diameter_distance_z1z2(1, 2).value,
                       646.22968662822018)
    z1 = 0, 0, 1, 0.5, 1
    z2 = 2, 1, 2, 2.5, 1.1
    results = (1760.0628637762106,
               1670.7497657219858,
               646.22968662822018,
               1159.0970895962193,
               115.72768186186921)

    assert np.allclose(tcos.angular_diameter_distance_z1z2(z1, z2).value,
                       results)


@pytest.mark.skipif('not HAS_SCIPY')
def test_absorption_distance():
    tcos = core.FlatLambdaCDM(70.4, 0.272, Tcmb0=0.0)
    assert np.allclose(tcos.absorption_distance([1, 3]),
                       [1.72576635, 7.98685853])
    assert np.allclose(tcos.absorption_distance(3), 7.98685853)


@pytest.mark.skipif('not HAS_SCIPY')
def test_massivenu_basic():
    # Test no neutrinos case
    tcos = core.FlatLambdaCDM(70.4, 0.272, Neff=4.05, m_nu=u.Quantity(0, u.eV))
    assert np.allclose(tcos.Neff, 4.05)
    assert not tcos.has_massive_nu
    mnu = tcos.m_nu
    assert len(mnu) == 4
    assert mnu.unit == u.eV
    assert np.allclose(mnu.value, [0.0, 0.0, 0.0, 0.0])
    assert np.allclose(tcos.nu_relative_density(1.0), 0.22710731766 * 4.05,
                       rtol=1e-6)

    # Test basic setting, retrieval of values
    tcos = core.FlatLambdaCDM(70.4, 0.272,
                              m_nu=u.Quantity([0.0, 0.01, 0.02], u.eV))
    assert tcos.has_massive_nu
    mnu = tcos.m_nu
    assert len(mnu) == 3
    assert mnu.unit == u.eV
    assert np.allclose(mnu.value, [0.0, 0.01, 0.02])

    # All massive neutrinos case
    tcos = core.FlatLambdaCDM(70.4, 0.272, m_nu=u.Quantity(0.1, u.eV),
                              Neff=3.1)
    assert np.allclose(tcos.Neff, 3.1)
    assert tcos.has_massive_nu
    mnu = tcos.m_nu
    assert len(mnu) == 3
    assert mnu.unit == u.eV
    assert np.allclose(mnu.value, [0.1, 0.1, 0.1])


@pytest.mark.skipif('not HAS_SCIPY')
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
    tcos = core.FlatLambdaCDM(75.0, 0.25, Tcmb0=3.0, Neff=3,
                              m_nu=u.Quantity(100.0, u.eV))
    assert tcos.has_massive_nu
    assert tcos.Neff == 3
    nurel_exp = nuprefac * tcos.Neff * np.array([171969, 85984.5, 57323,
                                                 15633.5, 171.801])
    assert np.allclose(tcos.nu_relative_density(ztest), nurel_exp, rtol=5e-3)
    assert np.allclose(tcos.efunc([0.0, 1.0]), [1.0, 7.46144727668], rtol=5e-3)

    # Next, slightly less massive
    tcos = core.FlatLambdaCDM(75.0, 0.25, Tcmb0=3.0, Neff=3,
                              m_nu=u.Quantity(0.25, u.eV))
    nurel_exp = nuprefac * tcos.Neff * np.array([429.924, 214.964, 143.312,
                                                 39.1005, 1.11086])
    assert np.allclose(tcos.nu_relative_density(ztest), nurel_exp,
                       rtol=5e-3)

    # For this one also test Onu directly
    onu_exp = np.array([0.01890217, 0.05244681, 0.0638236,
                        0.06999286, 0.1344951])
    assert np.allclose(tcos.Onu(ztest), onu_exp, rtol=5e-3)

    # And fairly light
    tcos = core.FlatLambdaCDM(80.0, 0.30, Tcmb0=3.0, Neff=3,
                              m_nu=u.Quantity(0.01, u.eV))

    nurel_exp = nuprefac * tcos.Neff * np.array([17.2347, 8.67345, 5.84348,
                                                 1.90671, 1.00021])
    assert np.allclose(tcos.nu_relative_density(ztest), nurel_exp,
                       rtol=5e-3)
    onu_exp = np.array([0.00066599, 0.00172677, 0.0020732,
                        0.00268404, 0.0978313])
    assert np.allclose(tcos.Onu(ztest), onu_exp, rtol=5e-3)
    assert np.allclose(tcos.efunc([1.0, 2.0]), [1.76225893, 2.97022048],
                       rtol=1e-4)
    assert np.allclose(tcos.inv_efunc([1.0, 2.0]), [0.5674535, 0.33667534],
                       rtol=1e-4)

    # Now a mixture of neutrino masses, with non-integer Neff
    tcos = core.FlatLambdaCDM(80.0, 0.30, Tcmb0=3.0, Neff=3.04,
                              m_nu=u.Quantity([0.0, 0.01, 0.25], u.eV))
    nurel_exp = nuprefac * tcos.Neff * np.array([149.386233, 74.87915, 50.0518,
                                                 14.002403, 1.03702333])
    assert np.allclose(tcos.nu_relative_density(ztest), nurel_exp,
                       rtol=5e-3)
    onu_exp = np.array([0.00584959, 0.01493142, 0.01772291,
                        0.01963451, 0.10227728])
    assert np.allclose(tcos.Onu(ztest), onu_exp, rtol=5e-3)


def test_default_reset():
    # Check that the default is being reset after tests. This test should be
    # updated if the default cosmology is updated.
    assert core.get_current() == core.WMAP9
