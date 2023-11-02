"""Testing :mod:`astropy.cosmology.units`."""

##############################################################################
# IMPORTS

import pytest

import astropy.cosmology.units as cu
import astropy.units as u
from astropy.cosmology import Planck13, default_cosmology
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.compat.optional_deps import HAS_SCIPY

##############################################################################
# TESTS
##############################################################################


def test_littleh():
    """Test :func:`astropy.cosmology.units.with_H0`."""
    H0_70 = 70 * u.km / u.s / u.Mpc
    h70dist = 70 * u.Mpc / cu.littleh

    assert_quantity_allclose(h70dist.to(u.Mpc, cu.with_H0(H0_70)), 100 * u.Mpc)

    # make sure using the default cosmology works
    cosmodist = default_cosmology.get().H0.value * u.Mpc / cu.littleh
    assert_quantity_allclose(cosmodist.to(u.Mpc, cu.with_H0()), 100 * u.Mpc)

    # Now try a luminosity scaling
    h1lum = 0.49 * u.Lsun * cu.littleh**-2
    assert_quantity_allclose(h1lum.to(u.Lsun, cu.with_H0(H0_70)), 1 * u.Lsun)

    # And the trickiest one: magnitudes.  Using H0=10 here for the round numbers
    H0_10 = 10 * u.km / u.s / u.Mpc
    # assume the "true" magnitude M = 12.
    # Then M - 5*log_10(h)  = M + 5 = 17
    withlittlehmag = 17 * (u.mag - u.MagUnit(cu.littleh**2))
    assert_quantity_allclose(withlittlehmag.to(u.mag, cu.with_H0(H0_10)), 12 * u.mag)


@pytest.mark.skipif(not HAS_SCIPY, reason="Cosmology needs scipy")
def test_dimensionless_redshift():
    """Test :func:`astropy.cosmology.units.dimensionless_redshift`."""
    z = 3 * cu.redshift
    val = 3 * u.one

    # show units not equal
    assert z.unit == cu.redshift
    assert z.unit != u.one
    assert u.get_physical_type(z) == "redshift"

    # test equivalency enabled by default
    assert z == val

    # also test that it works for powers
    assert (3 * cu.redshift**3) == val

    # and in composite units
    assert (3 * u.km / cu.redshift**3) == 3 * u.km

    # test it also works as an equivalency
    with u.set_enabled_equivalencies([]):  # turn off default equivalencies
        assert z.to(u.one, equivalencies=cu.dimensionless_redshift()) == val

        with pytest.raises(ValueError):
            z.to(u.one)

    # if this fails, something is really wrong
    with u.add_enabled_equivalencies(cu.dimensionless_redshift()):
        assert z == val


@pytest.mark.skipif(not HAS_SCIPY, reason="Cosmology needs scipy")
def test_redshift_temperature():
    """Test :func:`astropy.cosmology.units.redshift_temperature`."""
    cosmo = Planck13.clone(Tcmb0=3 * u.K)
    default_cosmo = default_cosmology.get()
    z = 15 * cu.redshift
    Tcmb = cosmo.Tcmb(z)

    # 1) Default (without specifying the cosmology)
    with default_cosmology.set(cosmo):
        equivalency = cu.redshift_temperature()
        assert_quantity_allclose(z.to(u.K, equivalency), Tcmb)
        assert_quantity_allclose(Tcmb.to(cu.redshift, equivalency), z)

    # showing the answer changes if the cosmology changes
    # this test uses the default cosmology
    equivalency = cu.redshift_temperature()
    assert_quantity_allclose(z.to(u.K, equivalency), default_cosmo.Tcmb(z))
    assert default_cosmo.Tcmb(z) != Tcmb

    # 2) Specifying the cosmology
    equivalency = cu.redshift_temperature(cosmo)
    assert_quantity_allclose(z.to(u.K, equivalency), Tcmb)
    assert_quantity_allclose(Tcmb.to(cu.redshift, equivalency), z)

    # Test `atzkw`
    equivalency = cu.redshift_temperature(cosmo, ztol=1e-10)
    assert_quantity_allclose(Tcmb.to(cu.redshift, equivalency), z)


@pytest.mark.skipif(not HAS_SCIPY, reason="Cosmology needs scipy")
def test_redshift_hubble():
    """Test :func:`astropy.cosmology.units.redshift_hubble`."""
    unit = u.km / u.s / u.Mpc
    cosmo = Planck13.clone(H0=100 * unit)
    default_cosmo = default_cosmology.get()
    z = 15 * cu.redshift
    H = cosmo.H(z)
    h = H.to_value(u.km / u.s / u.Mpc) / 100 * cu.littleh

    # 1) Default (without specifying the cosmology)
    with default_cosmology.set(cosmo):
        equivalency = cu.redshift_hubble()
        # H
        assert_quantity_allclose(z.to(unit, equivalency), H)
        assert_quantity_allclose(H.to(cu.redshift, equivalency), z)
        # little-h
        assert_quantity_allclose(z.to(cu.littleh, equivalency), h)
        assert_quantity_allclose(h.to(cu.redshift, equivalency), z)

    # showing the answer changes if the cosmology changes
    # this test uses the default cosmology
    equivalency = cu.redshift_hubble()
    assert_quantity_allclose(z.to(unit, equivalency), default_cosmo.H(z))
    assert default_cosmo.H(z) != H

    # 2) Specifying the cosmology
    equivalency = cu.redshift_hubble(cosmo)
    # H
    assert_quantity_allclose(z.to(unit, equivalency), H)
    assert_quantity_allclose(H.to(cu.redshift, equivalency), z)
    # little-h
    assert_quantity_allclose(z.to(cu.littleh, equivalency), h)
    assert_quantity_allclose(h.to(cu.redshift, equivalency), z)

    # Test `atzkw`
    equivalency = cu.redshift_hubble(cosmo, ztol=1e-10)
    assert_quantity_allclose(H.to(cu.redshift, equivalency), z)  # H
    assert_quantity_allclose(h.to(cu.redshift, equivalency), z)  # little-h


@pytest.mark.skipif(not HAS_SCIPY, reason="Cosmology needs scipy")
@pytest.mark.parametrize(
    "kind",
    [cu.redshift_distance.__defaults__[-1], "comoving", "lookback", "luminosity"],
)
def test_redshift_distance(kind):
    """Test :func:`astropy.cosmology.units.redshift_distance`."""
    z = 15 * cu.redshift
    d = getattr(Planck13, kind + "_distance")(z)

    equivalency = cu.redshift_distance(cosmology=Planck13, kind=kind)

    # properties of Equivalency
    assert equivalency.name[0] == "redshift_distance"
    assert equivalency.kwargs[0]["cosmology"] == Planck13
    assert equivalency.kwargs[0]["distance"] == kind

    # roundtrip
    assert_quantity_allclose(z.to(u.Mpc, equivalency), d)
    assert_quantity_allclose(d.to(cu.redshift, equivalency), z)


def test_redshift_distance_wrong_kind():
    """Test :func:`astropy.cosmology.units.redshift_distance` wrong kind."""
    with pytest.raises(ValueError, match="`kind`"):
        cu.redshift_distance(kind=None)


@pytest.mark.skipif(not HAS_SCIPY, reason="Cosmology needs scipy")
class Test_with_redshift:
    """Test `astropy.cosmology.units.with_redshift`."""

    @pytest.fixture(scope="class")
    def cosmo(self):
        """Test cosmology."""
        return Planck13.clone(Tcmb0=3 * u.K)

    # ===========================================

    def test_cosmo_different(self, cosmo):
        """The default is different than the test cosmology."""
        default_cosmo = default_cosmology.get()
        assert default_cosmo != cosmo  # shows changing default

    def test_no_equivalency(self, cosmo):
        """Test the equivalency  ``with_redshift`` without any enabled."""
        equivalency = cu.with_redshift(distance=None, hubble=False, Tcmb=False)
        assert len(equivalency) == 0

    # -------------------------------------------

    def test_temperature_off(self, cosmo):
        """Test ``with_redshift`` with the temperature off."""
        z = 15 * cu.redshift
        err_msg = (
            r"^'redshift' \(redshift\) and 'K' \(temperature\) are not convertible$"
        )

        # 1) Default (without specifying the cosmology)
        with default_cosmology.set(cosmo):
            equivalency = cu.with_redshift(Tcmb=False)
            with pytest.raises(u.UnitConversionError, match=err_msg):
                z.to(u.K, equivalency)

        # 2) Specifying the cosmology
        equivalency = cu.with_redshift(cosmo, Tcmb=False)
        with pytest.raises(u.UnitConversionError, match=err_msg):
            z.to(u.K, equivalency)

    def test_temperature(self, cosmo):
        """Test temperature equivalency component."""
        default_cosmo = default_cosmology.get()
        z = 15 * cu.redshift
        Tcmb = cosmo.Tcmb(z)

        # 1) Default (without specifying the cosmology)
        with default_cosmology.set(cosmo):
            equivalency = cu.with_redshift(Tcmb=True)
            assert_quantity_allclose(z.to(u.K, equivalency), Tcmb)
            assert_quantity_allclose(Tcmb.to(cu.redshift, equivalency), z)

        # showing the answer changes if the cosmology changes
        # this test uses the default cosmology
        equivalency = cu.with_redshift(Tcmb=True)
        assert_quantity_allclose(z.to(u.K, equivalency), default_cosmo.Tcmb(z))
        assert default_cosmo.Tcmb(z) != Tcmb

        # 2) Specifying the cosmology
        equivalency = cu.with_redshift(cosmo, Tcmb=True)
        assert_quantity_allclose(z.to(u.K, equivalency), Tcmb)
        assert_quantity_allclose(Tcmb.to(cu.redshift, equivalency), z)

        # Test `atzkw`
        # this is really just a test that 'atzkw' doesn't fail
        equivalency = cu.with_redshift(cosmo, Tcmb=True, atzkw={"ztol": 1e-10})
        assert_quantity_allclose(Tcmb.to(cu.redshift, equivalency), z)

    # -------------------------------------------

    def test_hubble_off(self, cosmo):
        """Test ``with_redshift`` with Hubble off."""
        unit = u.km / u.s / u.Mpc
        z = 15 * cu.redshift
        err_msg = (
            r"^'redshift' \(redshift\) and 'km / \(Mpc s\)' \(frequency\) are not "
            "convertible$"
        )

        # 1) Default (without specifying the cosmology)
        with default_cosmology.set(cosmo):
            equivalency = cu.with_redshift(hubble=False)
            with pytest.raises(u.UnitConversionError, match=err_msg):
                z.to(unit, equivalency)

        # 2) Specifying the cosmology
        equivalency = cu.with_redshift(cosmo, hubble=False)
        with pytest.raises(u.UnitConversionError, match=err_msg):
            z.to(unit, equivalency)

    def test_hubble(self, cosmo):
        """Test Hubble equivalency component."""
        unit = u.km / u.s / u.Mpc
        default_cosmo = default_cosmology.get()
        z = 15 * cu.redshift
        H = cosmo.H(z)
        h = H.to_value(u.km / u.s / u.Mpc) / 100 * cu.littleh

        # 1) Default (without specifying the cosmology)
        with default_cosmology.set(cosmo):
            equivalency = cu.with_redshift(hubble=True)
            # H
            assert_quantity_allclose(z.to(unit, equivalency), H)
            assert_quantity_allclose(H.to(cu.redshift, equivalency), z)
            # little-h
            assert_quantity_allclose(z.to(cu.littleh, equivalency), h)
            assert_quantity_allclose(h.to(cu.redshift, equivalency), z)

        # showing the answer changes if the cosmology changes
        # this test uses the default cosmology
        equivalency = cu.with_redshift(hubble=True)
        assert_quantity_allclose(z.to(unit, equivalency), default_cosmo.H(z))
        assert default_cosmo.H(z) != H

        # 2) Specifying the cosmology
        equivalency = cu.with_redshift(cosmo, hubble=True)
        # H
        assert_quantity_allclose(z.to(unit, equivalency), H)
        assert_quantity_allclose(H.to(cu.redshift, equivalency), z)
        # little-h
        assert_quantity_allclose(z.to(cu.littleh, equivalency), h)
        assert_quantity_allclose(h.to(cu.redshift, equivalency), z)

        # Test `atzkw`
        # this is really just a test that 'atzkw' doesn't fail
        equivalency = cu.with_redshift(cosmo, hubble=True, atzkw={"ztol": 1e-10})
        assert_quantity_allclose(H.to(cu.redshift, equivalency), z)  # H
        assert_quantity_allclose(h.to(cu.redshift, equivalency), z)  # h

    # -------------------------------------------

    def test_distance_off(self, cosmo):
        """Test ``with_redshift`` with the distance off."""
        z = 15 * cu.redshift
        err_msg = r"^'redshift' \(redshift\) and 'Mpc' \(length\) are not convertible$"

        # 1) Default (without specifying the cosmology)
        with default_cosmology.set(cosmo):
            equivalency = cu.with_redshift(distance=None)
            with pytest.raises(u.UnitConversionError, match=err_msg):
                z.to(u.Mpc, equivalency)

        # 2) Specifying the cosmology
        equivalency = cu.with_redshift(cosmo, distance=None)
        with pytest.raises(u.UnitConversionError, match=err_msg):
            z.to(u.Mpc, equivalency)

    def test_distance_default(self):
        """Test distance equivalency default."""
        z = 15 * cu.redshift
        d = default_cosmology.get().comoving_distance(z)

        equivalency = cu.with_redshift()
        assert_quantity_allclose(z.to(u.Mpc, equivalency), d)
        assert_quantity_allclose(d.to(cu.redshift, equivalency), z)

    def test_distance_wrong_kind(self):
        """Test distance equivalency, but the wrong kind."""
        with pytest.raises(ValueError, match="`kind`"):
            cu.with_redshift(distance=ValueError)

    @pytest.mark.parametrize("kind", ["comoving", "lookback", "luminosity"])
    def test_distance(self, kind):
        """Test distance equivalency."""
        cosmo = Planck13
        z = 15 * cu.redshift
        dist = getattr(cosmo, kind + "_distance")(z)

        default_cosmo = default_cosmology.get()
        assert default_cosmo != cosmo  # shows changing default

        # 1) without specifying the cosmology
        with default_cosmology.set(cosmo):
            equivalency = cu.with_redshift(distance=kind)
            assert_quantity_allclose(z.to(u.Mpc, equivalency), dist)

        # showing the answer changes if the cosmology changes
        # this test uses the default cosmology
        equivalency = cu.with_redshift(distance=kind)
        assert_quantity_allclose(
            z.to(u.Mpc, equivalency), getattr(default_cosmo, kind + "_distance")(z)
        )
        assert not u.allclose(getattr(default_cosmo, kind + "_distance")(z), dist)

        # 2) Specifying the cosmology
        equivalency = cu.with_redshift(cosmo, distance=kind)
        assert_quantity_allclose(z.to(u.Mpc, equivalency), dist)
        assert_quantity_allclose(dist.to(cu.redshift, equivalency), z)

        # Test atzkw
        # this is really just a test that 'atzkw' doesn't fail
        equivalency = cu.with_redshift(cosmo, distance=kind, atzkw={"ztol": 1e-10})
        assert_quantity_allclose(dist.to(cu.redshift, equivalency), z)


def test_equivalency_context_manager():
    base_registry = u.get_current_unit_registry()

    # check starting with only the dimensionless_redshift equivalency.
    assert len(base_registry.equivalencies) == 1
    assert str(base_registry.equivalencies[0][0]) == "redshift"
