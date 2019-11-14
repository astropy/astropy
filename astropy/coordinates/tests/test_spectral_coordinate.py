import pytest

from astropy.tests.helper import assert_quantity_allclose
from astropy.coordinates import (SkyCoord, EarthLocation, ICRS, GCRS, Galactic,
                                 CartesianDifferential, SphericalRepresentation,
                                 RadialDifferential, get_body_barycentric_posvel)
from astropy import time
import numpy as np
import astropy.units as u

from ..spectra.spectral_coordinate import SpectralCoord

# Local Standard of Rest
LSRD = Galactic(u=0 * u.km, v=0 * u.km, w=0 * u.km,
                U=9 * u.km / u.s, V=12 * u.km / u.s, W=7 * u.km / u.s,
                representation_type='cartesian', differential_type='cartesian')

# Target located in direction of motion of LSRD
LSRD_DIR_STATIONARY = Galactic(u=9 * u.km, v=12 * u.km, w=7 * u.km,
                               representation_type='cartesian')


def assert_frame_allclose(frame1, frame2,
                          pos_rtol=1e-7, pos_atol=1 * u.m,
                          vel_rtol=1e-7, vel_atol=1 * u.mm / u.s):
    # TODO: write a function that checks that:
    # - the positions are equal to within some tolerance (the relative tolerance
    #   should be dimensionless, the absolute tolerance should be a distance)
    # - either both or nether frame has velocities
    # - if velocities are present, they are equal to some tolerance
    # Ideally this should accept both frames and SkyCoords
    assert frame1 is frame2


# GENERAL TESTS


class TestCreateSpectralCoord:

    def setup_class(self):
        self.wav = [100, 200, 300] * u.nm

    def test_simple(self):
        # Initialize a SpectralCoord without any observer or target
        coord = SpectralCoord(self.wav)
        assert coord.observer is None
        assert coord.target is None
        assert coord.rest is None
        assert coord.velocity_convention is None
        assert coord.radial_velocity is None
        assert coord.redshift is None

    @pytest.mark.parametrize('observer_skycoord', [False, True])
    def test_with_frame_observer(self, observer_skycoord):
        # Ensure we can initialize a SpectralCoord with an observer
        observer = LSRD
        if observer_skycoord:
            observer = SkyCoord(observer)
        coord = SpectralCoord(self.wav, observer=LSRD)
        assert coord.observer is LSRD
        assert coord.target is None
        assert coord.rest is None
        assert coord.velocity_convention is None
        assert coord.radial_velocity is None
        assert coord.redshift is None

    @pytest.mark.parametrize('observer_skycoord', [False, True])
    @pytest.mark.parametrize('target_skycoord', [False, True])
    def test_with_skycoord_observer_skycoord_target(self, observer_skycoord, target_skycoord):
        # Ensure we can initialize a SpectralCoord with an observer and a target
        observer = LSRD
        if observer_skycoord:
            observer = SkyCoord(observer)
        target = Galactic(u=9 * u.km, v=12 * u.km, w=7 * u.km,
                          representation_type='cartesian')
        if target_skycoord:
            target = SkyCoord(target)
        coord = SpectralCoord(self.wav, observer=observer, target=target)
        assert coord.observer is skycoord
        assert coord.target is target
        assert coord.rest is None
        assert coord.velocity_convention is None
        assert_quantity_allclose(coord.radial_velocity,
                                 (9 * 9 + 12 * 12 + 7 * 7) ** 0.5 * u.km / u.s)
        assert coord.redshift is None


def test_create_spectral_coord():
    keck = EarthLocation.of_site('keck')
    obstime = time.Time('2018-12-13 9:00')

    observer_gcrs = keck.get_gcrs(obstime)

    spectral_axis = np.linspace(500, 2500, 1000) * u.AA
    spec_coord = SpectralCoord(spectral_axis)

    assert isinstance(spec_coord, u.Quantity)
    assert len(spec_coord) == 1000

    spec_coord.rest = 1215 * u.AA


# SCIENCE USE CASE TESTS
