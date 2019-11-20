import pytest

from astropy.tests.helper import assert_quantity_allclose, quantity_allclose
from astropy.coordinates import (SkyCoord, EarthLocation, ICRS, GCRS, Galactic,
                                 CartesianDifferential, SphericalRepresentation,
                                 RadialDifferential, get_body_barycentric_posvel,
                                 FK5, CartesianRepresentation)
from astropy import time
import numpy as np
import astropy.units as u
from astropy.constants import c

from ..spectra.spectral_coordinate import SpectralCoord


def assert_frame_allclose(frame1, frame2,
                          pos_rtol=1e-7, pos_atol=1 * u.m,
                          vel_rtol=1e-7, vel_atol=1 * u.mm / u.s):
    # checks that:
    # - the positions are equal to within some tolerance (the relative tolerance
    #   should be dimensionless, the absolute tolerance should be a distance).
    #   note that these are the tolerances *in 3d*
    # - either both or nether frame has velocities, or if one has no velocities
    #   the other one can have zero velocities
    # - if velocities are present, they are equal to some tolerance
    # Ideally this should accept both frames and SkyCoords
    if hasattr(frame1, 'frame'):  # SkyCoord-like
        frame1 = frame1.frame
    if hasattr(frame2, 'frame'):  # SkyCoord-like
        frame2 = frame2.frame

    assert (frame1.differentials and frame2.differentials or
            (not frame1.differentials and not frame2.differentials))

    frame2_in_1 = frame2.transform_to(frame1)

    assert_quantity_allclose(0, frame1.separation_3d(frame2_in_1), rtol=pos_rtol, atol=pos_atol)

    if frame1.differentials:
        d1 = frame1.data.represent_as(CartesianRepresentation, CartesianDifferential).differentials['s']
        d2 = frame2_in_1.data.represent_as(CartesianRepresentation, CartesianDifferential).differentials['s']

        assert_quantity_allclose(0, d1.norm(d2), rtol=pos_rtol, atol=pos_atol)



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

def test_spectral_coord_jupiter():
    """
    Checks that jupiter yields an RV consistent with the solar system
    """
    obstime = time.Time('2018-12-13 9:00')
    obs = EarthLocation.of_site('greenwich').get_gcrs(obstime)  # always built-in, no download required

    # jupiter = get_body('jupiter', obstime)  # not supported by astropy yet,
    # but this is the eventual goal, which should return:
    pos, vel = get_body_barycentric_posvel('jupiter', obstime)
    jupiter = SkyCoord(pos.with_differentials(CartesianDifferential(vel.xyz)), obstime=obstime)

    spc = SpectralCoord([100, 200, 300] * u.nm, observer=obs, target=jupiter)

    # "magic number" of 45 is actually  ~43 + a bit extra, because the
    # maximum possible speed should be the earth-jupiter relative velocity
    assert np.abs(spc.radial_velocity) < (45*u.km/u.s)


def test_spectral_coord_alphacen():
    """
    Checks that a nearby star yields a reasonable RV
    """
    obstime = time.Time('2018-12-13 9:00')
    obs = EarthLocation.of_site('greenwich').get_gcrs(obstime)  # always built-in, no download required

    # acen = SkyCoord.from_name('alpha cen')  # coordinates are from here,
    # but hard-coded below so no download required
    acen = SkyCoord(ra=219.90085*u.deg, dec=-60.83562*u.deg, frame='icrs',
                    distance=4.37*u.lightyear, radial_velocity=-18.*u.km/u.s)

    spc = SpectralCoord([100, 200, 300] * u.nm, observer=obs, target=acen)

    # "magic number" of 50 is actually  ~18 + 30 + a bit extra, because the
    # maximum possible speed should be the star + earth + earth's surface
    assert np.abs(spc.radial_velocity) < (50*u.km/u.s)


def test_spectral_coord_m31():
    """
    Checks that a nearby star yields a reasonable RV
    """
    obstime = time.Time('2018-12-13 9:00')
    obs = EarthLocation.of_site('greenwich').get_gcrs(obstime)  # always built-in, no download required

    # acen = SkyCoord.from_name('alpha cen')  # coordinates are from here, but
    # hard-coded below so no download required
    m31 = SkyCoord(ra=10.6847*u.deg, dec=41.269*u.deg,
                   distance=710*u.kpc, radial_velocity=-300*u.km/u.s)

    spc = SpectralCoord([100, 200, 300] * u.nm, observer=obs, target=m31)

    # "magic number"  is actually  ~300 + 30 + a bit extra km/s, because the
    # maximum possible speed should be M31 + earth + earth's surface.  Then
    # transform to an approximate redshift bound
    assert np.abs(spc.redshift) < 0.00112


def test_shift_to_rest_galaxy():
    """
    This tests storing a spectral coordinate with a specific redshift, and then
    doing basic rest-to-observed-and-back transformations
    """
    z = 5
    rest_line_wls = [5007, 6563]*u.angstrom

    observed_spc = SpectralCoord(rest_line_wls*(z+1), redshift=z)
    rest_spc = observed_spc.to_rest()
    # alternatively:
    # rest_spc = observed_spc.with_observer(observed_spec.target)
    # although then it would have to be clearly documented, or the `to_rest`
    # implemented in Spectrum1D?

    assert_quantity_allclose(rest_spc, rest_line_wls)
    assert_frame_allclose(rest_spc.observer, rest_spc.target)

    with pytest.raises(ValueError):
        # *any* observer shift should fail, since we didn't specify one at the
        # outset
        rest_spc.with_observer(ICRS(CartesianRepresentation([0,0,0]*u.au)))

    # note: it may be an acceptable fallback for the next part to onle work on
    # spectrum1D but not SpectralCoord - the thinking being that the shift to
    # rest might require dropping the "observer" information from spectralcoord
    # but maybe not? So including the test here.

    roundtrip_obs_spc = rest_spc.to_observed()
    assert_quantity_allclose(roundtrip_obs_spc, rest_line_wls*(z+1))


def test_shift_to_rest_star_withobserver():
    rv = -18*u.km/u.s
    rest_line_wls = [5007, 6563]*u.angstrom

    obstime = time.Time('2018-12-13 9:00')
    eloc = EarthLocation.of_site('greenwich')
    obs = eloc.get_gcrs(obstime)  # always built-in, no download required
    acen = SkyCoord(ra=219.90085*u.deg, dec=-60.83562*u.deg, frame='icrs',
                    distance=4.37*u.lightyear)
    # Note that above the rv is missing from the SkyCoord.
    # That's intended, as it will instead be set in the `SpectralCoord`.  But
    # the SpectralCoord machinery should yield something comparable to test_
    # spectral_coord_alphacen

    observed_spc = SpectralCoord(rest_line_wls*(rv/c + 1),
                                 observer=obs, target=acen, radial_velocity=rv)

    rest_spc = observed_spc.to_rest()
    assert_quantity_allclose(rest_spc, rest_line_wls)

    barycentric_spc = observed_spc.with_observer(ICRS(CartesianRepresentation([0,0,0]*u.au)))
    baryrest_spc = barycentric_spc.to_rest()
    assert not quantity_allclose(baryrest_spc, rest_line_wls)

    # now make sure the change the barycentric shift did is comparable to the
    # offset rv_correction produces
    barytarg = SkyCoord(barycentric_spc.target)
    vcorr = barytarg.radial_velocity_correction(kind='barycentric',
                                                obstime=obstime, location=eloc)

    drv = baryrest_spc.radial_velocity - observed_spc.radial_velocity

    # note this probably will not work on the first try, but it's ok if this is
    # "good enough", where good enough is estimated below.  But that could be
    # adjusted if we think that's too agressive of a precision target for what
    # the machinery can handle
    assert_quantity_allclose(vcorr, drv, atol=10*u.m/u.s)
