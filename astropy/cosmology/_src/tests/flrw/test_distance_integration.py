# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Accuracy of shared comoving-distance integration for redshift catalogues."""

import numpy as np
import pytest

import astropy.units as u
from astropy.cosmology import WMAP9, LambdaCDM, Planck18, w0waCDM
from astropy.cosmology import units as cu
from astropy.utils.compat.optional_deps import HAS_SCIPY

pytestmark = pytest.mark.skipif(not HAS_SCIPY, reason="scipy required")


@pytest.mark.parametrize(
    "cosmo",
    [
        Planck18,
        WMAP9,
        LambdaCDM(70, 0.3, 0.6, Tcmb0=2.7),
        w0waCDM(70, 0.3, 0.6, w0=-0.9, wa=0.2, Tcmb0=2.7),
    ],
)
@pytest.mark.parametrize("origin", [0.0, 1.0])
def test_catalogue_distance_matches_independent_intervals(cosmo, origin):
    # Repeated, unsorted redshifts in multiple dimensions, on both sides of
    # the observer. Compare to independent quadrature, not an interpolation.
    redshifts = np.random.default_rng(5281).uniform(-0.5, 10, (16, 32))
    redshifts[0] = origin
    redshifts[1] = redshifts[2]
    expected = cosmo.hubble_distance * cosmo._integral_comoving_distance_z1z2_scalar(
        origin, redshifts
    )
    actual = cosmo.comoving_distance(origin * cu.redshift, redshifts * cu.redshift)
    assert actual.shape == redshifts.shape
    assert actual.unit == u.Mpc
    np.testing.assert_allclose(actual.value, expected.value, rtol=1e-9, atol=1e-8)
    np.testing.assert_array_equal(actual[0].value, 0)
    np.testing.assert_array_equal(cosmo.comoving_distance(redshifts, origin), -actual)


@pytest.mark.parametrize("origin", [1.0, 10000.0])
def test_catalogue_distance_near_observer(origin):
    # A cumulative integral from zero followed by subtraction loses these
    # distances. Integrating from the actual observer retains their precision.
    step = np.spacing(origin)
    redshifts = origin + np.arange(-256, 256) * step
    # Over a few hundred ulps the integrand is constant to 1e-12 relative;
    # this local expansion avoids QUADPACK's roundoff warnings on tiny intervals.
    expected = (
        Planck18.hubble_distance * Planck18.inv_efunc(origin) * (redshifts - origin)
    )
    actual = Planck18.comoving_distance(origin, redshifts)
    np.testing.assert_allclose(actual.value, expected.value, rtol=1e-11, atol=0)


def test_catalogue_pairwise_distance():
    # Independent observer/source pairs must keep direct interval integration.
    origins = np.linspace(0.5, 10, 512)
    redshifts = origins + 4 * np.spacing(origins)
    expected = (
        Planck18.hubble_distance
        * Planck18._integral_comoving_distance_z1z2_scalar(origins, redshifts)
    )
    np.testing.assert_array_equal(
        Planck18.comoving_distance(origins, redshifts), expected
    )


def test_catalogue_unsigned_redshifts():
    redshifts = np.arange(512, dtype=np.uint64)
    expected = (
        Planck18.hubble_distance
        * Planck18._integral_comoving_distance_z1z2_scalar(256, redshifts)
    )
    np.testing.assert_allclose(
        Planck18.comoving_distance(256, redshifts).value,
        expected.value,
        rtol=1e-9,
        atol=1e-8,
    )


def test_catalogue_integer_redshifts_coincident_as_doubles():
    # Distinct integer redshifts can coincide at the integrator's precision.
    # Deduplication must use that precision to keep t_eval strictly ordered.
    origin = 10**16
    redshifts = origin + np.arange(-256, 256, dtype=np.int64)
    expected = (
        Planck18.hubble_distance
        * Planck18.inv_efunc(float(origin))
        * (redshifts.astype(float) - float(origin))
    )
    np.testing.assert_allclose(
        Planck18.comoving_distance(origin, redshifts).value,
        expected.value,
        rtol=1e-11,
        atol=0,
    )


def test_catalogue_infinite_endpoint():
    redshifts = np.resize([0, 0.1, 1, np.inf], 512)
    expected = (
        Planck18.hubble_distance
        * Planck18._integral_comoving_distance_z1z2_scalar(0, redshifts)
    )
    np.testing.assert_array_equal(Planck18.comoving_distance(redshifts), expected)


def test_catalogue_scalar_override():
    class CustomCosmology(LambdaCDM):
        def _integral_comoving_distance_z1z2_scalar(self, z1, z2, /):
            return np.asarray(z2) - z1

    cosmo = CustomCosmology(70, 0.3, 0.6, Tcmb0=2.7)
    redshifts = np.linspace(0, 10, 512)
    np.testing.assert_array_equal(
        cosmo.comoving_distance(1, redshifts), cosmo.hubble_distance * (redshifts - 1)
    )
