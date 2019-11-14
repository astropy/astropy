import pytest

from astropy.tests.helper import assert_quantity_allclose
from astropy.coordinates import (SkyCoord, EarthLocation, ICRS, GCRS, Galactic,
                                 CartesianDifferential, SphericalRepresentation,
                                 RadialDifferential, get_body_barycentric_posvel,
                                 FK5)
from astropy import time
import numpy as np
import astropy.units as u

from ..spectra.spectral_coordinate import SpectralCoord


def assert_frame_allclose(frame1, frame2,
                          pos_rtol=1e-7, pos_atol=1 * u.m,
                          vel_rtol=1e-7, vel_atol=1 * u.mm / u.s):
    # TODO: write a function that checks that:
    # - the positions are equal to within some tolerance (the relative tolerance
    #   should be dimensionless, the absolute tolerance should be a distance)
    # - either both or nether frame has velocities, or if one has no velocities
    #   the other one can have zero velocities
    # - if velocities are present, they are equal to some tolerance
    # Ideally this should accept both frames and SkyCoords
    pass


def test_create_spectral_coord_orig():

    # TODO: decide whether this test is still needed once the rest is implemented

    keck = EarthLocation.of_site('keck')
    obstime = time.Time('2018-12-13 9:00')

    observer_gcrs = keck.get_gcrs(obstime)

    spectral_axis = np.linspace(500, 2500, 1000) * u.AA
    spec_coord = SpectralCoord(spectral_axis)

    assert isinstance(spec_coord, u.Quantity)
    assert len(spec_coord) == 1000

    spec_coord.rest = 1215 * u.AA


# GENERAL TESTS

# We first run through a series of cases to test different ways of initializing
# the observer and target for SpectralCoord, including for example frames,
# SkyCoords, and making sure that SpectralCoord is not sensitive to the actual
# frame or representation class.

# Local Standard of Rest
LSRD = Galactic(u=0 * u.km, v=0 * u.km, w=0 * u.km,
                U=9 * u.km / u.s, V=12 * u.km / u.s, W=7 * u.km / u.s,
                representation_type='cartesian', differential_type='cartesian')

LSRD_EQUIV = [
              LSRD,
              SkyCoord(LSRD),  # as a SkyCoord
              LSRD.transform_to(ICRS),  # different frame
              LSRD.transform_to(ICRS).transform_to(Galactic)  # different representation
              ]


@pytest.fixture(params=[None] + LSRD_EQUIV)
def observer(request):
    return request.param


# Target located in direction of motion of LSRD with no velocities
LSRD_DIR_STATIONARY = Galactic(u=9 * u.km, v=12 * u.km, w=7 * u.km,
                               representation_type='cartesian')

LSRD_DIR_STATIONARY_EQUIV = [
                             LSRD_DIR_STATIONARY,
                             SkyCoord(LSRD_DIR_STATIONARY),  # as a SkyCoord
                             LSRD_DIR_STATIONARY.transform_to(FK5()),  # different frame
                             LSRD_DIR_STATIONARY.transform_to(ICRS()).transform_to(Galactic())  # different representation
                            ]


@pytest.fixture(params=[None] + LSRD_DIR_STATIONARY_EQUIV)
def target(request):
    return request.param


def test_create_spectral_coord_observer_target(observer, target):

    coord = SpectralCoord([100, 200, 300] * u.nm, observer=observer, target=target)

    if observer is None:
        assert coord.observer is None
    else:
        assert_frame_allclose(observer, coord.observer)

    if target is None:
        assert coord.target is None
    else:
        assert_frame_allclose(target, coord.target)

    assert coord.rest is None
    assert coord.velocity_convention is None

    if observer is None or target is None:
        assert coord.redshift is None
        assert coord.radial_velocity is None
    elif observer in LSRD_EQUIV and target in LSRD_DIR_STATIONARY_EQUIV:
        assert_quantity_allclose(coord.redshift, ((1 + 274 / 299792.458**2) / (1 - 274 / 299792.458**2)) ** 0.5 - 1)
        assert_quantity_allclose(coord.radial_velocity, 274 ** 0.5 * u.km / u.s)
    else:
        raise NotImplementedError()


# SCIENCE USE CASE TESTS
