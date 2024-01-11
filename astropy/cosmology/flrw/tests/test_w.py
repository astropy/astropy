# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.flrw.lambdacdm`."""

import numpy as np
import pytest

import astropy.units as u
from astropy.cosmology import FLRW, wCDM
from astropy.utils.compat.optional_deps import HAS_SCIPY

##############################################################################
# TYPES


class W1(FLRW):
    """
    This class is to test whether the routines work correctly if one only overloads w(z).
    """

    def __init__(self):
        super().__init__(70.0, 0.27, 0.73, Tcmb0=0.0, name="test_cos")
        self._w0 = -0.9

    def w(self, z):
        return self._w0 * np.ones_like(z)


class W1nu(FLRW):
    """Similar, but with neutrinos."""

    def __init__(self):
        super().__init__(
            70.0, 0.27, 0.73, Tcmb0=3.0, m_nu=0.1 * u.eV, name="test_cos_nu"
        )
        self._w0 = -0.8

    def w(self, z):
        return self._w0 * np.ones_like(z)


##############################################################################
# TESTS
##############################################################################


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_de_subclass():
    z = [0.2, 0.4, 0.6, 0.9]

    # This is the comparison object
    cosmo = wCDM(H0=70, Om0=0.27, Ode0=0.73, w0=-0.9, Tcmb0=0.0)
    # Values taken from Ned Wrights advanced cosmo calculator, Aug 17 2012
    assert u.allclose(
        cosmo.luminosity_distance(z), [975.5, 2158.2, 3507.3, 5773.1] * u.Mpc, rtol=1e-3
    )

    # Now try the subclass that only gives w(z)
    cosmo = W1()
    assert u.allclose(
        cosmo.luminosity_distance(z), [975.5, 2158.2, 3507.3, 5773.1] * u.Mpc, rtol=1e-3
    )
    # Test efunc
    assert u.allclose(cosmo.efunc(1.0), 1.7489240754, rtol=1e-5)
    assert u.allclose(cosmo.efunc([0.5, 1.0]), [1.31744953, 1.7489240754], rtol=1e-5)
    assert u.allclose(cosmo.inv_efunc([0.5, 1.0]), [0.75904236, 0.57178011], rtol=1e-5)
    # Test de_density_scale
    assert u.allclose(cosmo.de_density_scale(1.0), 1.23114444, rtol=1e-4)
    assert u.allclose(
        cosmo.de_density_scale([0.5, 1.0]), [1.12934694, 1.23114444], rtol=1e-4
    )


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_efunc_vs_invefunc_flrw():
    """Test that efunc and inv_efunc give inverse values"""
    z0 = 0.5
    z = np.array([0.5, 1.0, 2.0, 5.0])

    # FLRW is abstract, so requires W1 defined earlier
    # This requires scipy, unlike the built-ins, because it
    # calls de_density_scale, which has an integral in it
    cosmo = W1()
    assert u.allclose(cosmo.efunc(z0), 1.0 / cosmo.inv_efunc(z0))
    assert u.allclose(cosmo.efunc(z), 1.0 / cosmo.inv_efunc(z))
    # Add neutrinos
    cosmo = W1nu()
    assert u.allclose(cosmo.efunc(z0), 1.0 / cosmo.inv_efunc(z0))
    assert u.allclose(cosmo.efunc(z), 1.0 / cosmo.inv_efunc(z))
