import astropy.units as u
import numpy as np
import pytest
from astropy import time
from astropy.constants import c
from astropy.coordinates import (SkyCoord, EarthLocation, ICRS, GCRS, Galactic,
                                 CartesianDifferential,
                                 get_body_barycentric_posvel,
                                 FK5, CartesianRepresentation)
from astropy.tests.helper import assert_quantity_allclose, quantity_allclose

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

    # assert (frame1.data.differentials and frame2.data.differentials or
    #         (not frame1.data.differentials and not frame2.data.differentials))
    assert frame1.is_equivalent_frame(frame2)

    frame2_in_1 = frame2.transform_to(frame1)

    assert_quantity_allclose(0 * u.m, frame1.separation_3d(frame2_in_1), rtol=pos_rtol, atol=pos_atol)

    if frame1.data.differentials:
        d1 = frame1.data.represent_as(CartesianRepresentation, CartesianDifferential).differentials['s']
        d2 = frame2_in_1.data.represent_as(CartesianRepresentation, CartesianDifferential).differentials['s']

        assert_quantity_allclose(d1.norm(d1), d1.norm(d2), rtol=vel_rtol, atol=vel_atol)


def get_greenwich_earthlocation():
    """
    A helper function to get an EarthLocation for greenwich (without trying to
    do a download)
    """
    site_registry = EarthLocation._get_site_registry(force_builtin=True)
    return site_registry.get('greenwich')

def test_create_spectral_coord_orig():

    # TODO: decide whether this test is still needed once the rest is
    #  implemented

    site = get_greenwich_earthlocation()
    obstime = time.Time('2018-12-13 9:00')

    observer_gcrs = site.get_gcrs(obstime)

    spectral_axis = np.linspace(500, 2500, 1000) * u.AA
    spec_coord = SpectralCoord(spectral_axis, observer=observer_gcrs)

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

    assert coord.doppler_rest is None
    assert coord.doppler_convention is None

    if observer is None or target is None:
        assert quantity_allclose(coord.redshift, 0)
        assert quantity_allclose(coord.radial_velocity, 0 * u.km/u.s)
    elif observer in LSRD_EQUIV and target in LSRD_DIR_STATIONARY_EQUIV:
        assert_quantity_allclose(coord.radial_velocity, -274 ** 0.5 * u.km / u.s)
        assert_quantity_allclose(coord.redshift, -274 ** 0.5 / 299792.458)
    else:
        raise NotImplementedError()

def test_create_from_spectral_coord(observer, target):
    """
    Checks that parameters are correctly copied to the new SpectralCoord object
    """
    spec_coord1 = SpectralCoord([100, 200, 300] * u.nm, observer=observer,
            target=target, radial_velocity=u.Quantity(1000, 'km/s'),
            doppler_convention = 'optical', doppler_rest = 6000*u.AA)
    spec_coord2 = SpectralCoord(spec_coord1)
    assert spec_coord1.observer == spec_coord2.observer
    assert spec_coord2.target == spec_coord2.target
    assert spec_coord2.radial_velocity == spec_coord2.radial_velocity
    assert spec_coord2.doppler_convention == spec_coord2.doppler_convention
    assert spec_coord2.doppler_rest == spec_coord2.doppler_rest

# SCIENCE USE CASE TESTS

def test_spectral_coord_jupiter():
    """
    Checks that jupiter yields an RV consistent with the solar system
    """
    obstime = time.Time('2018-12-13 9:00')
    obs = get_greenwich_earthlocation().get_gcrs(obstime)  # always built-in, no download required

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
    obs = get_greenwich_earthlocation().get_gcrs(obstime)  # always built-in, no download required

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
    obs = get_greenwich_earthlocation().get_gcrs(obstime)  # always built-in, no download required

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

    # No frames are explicitly defined, so to the user, the observer and
    #  target are not set.
    with pytest.raises(AttributeError):
        assert_frame_allclose(rest_spc.observer, rest_spc.target)

    with pytest.raises(ValueError):
        # *any* observer shift should fail, since we didn't specify one at the
        # outset
        rest_spc._change_observer_to(ICRS(CartesianRepresentation([0, 0, 0] * u.au)))

    # note: it may be an acceptable fallback for the next part to onle work on
    # spectrum1D but not SpectralCoord - the thinking being that the shift to
    # rest might require dropping the "observer" information from spectralcoord
    # but maybe not? So including the test here.

    roundtrip_obs_spc = rest_spc.to_observed()
    assert_quantity_allclose(roundtrip_obs_spc, rest_line_wls*(z+1))


def test_shift_to_rest_star_withobserver():
    rv = -8.3283011*u.km/u.s
    rest_line_wls = [5007, 6563]*u.angstrom

    obstime = time.Time('2018-12-13 9:00')
    eloc = get_greenwich_earthlocation()
    obs = eloc.get_gcrs(obstime)  # always built-in, no download required
    acen = SkyCoord(ra=219.90085*u.deg, dec=-60.83562*u.deg, frame='icrs',
                    distance=4.37*u.lightyear)
    # Note that above the rv is missing from the SkyCoord.
    # That's intended, as it will instead be set in the `SpectralCoord`.  But
    # the SpectralCoord machinery should yield something comparable to test_
    # spectral_coord_alphacen

    observed_spc = SpectralCoord(rest_line_wls*(rv/c + 1),
                                 observer=obs, target=acen)

    rest_spc = observed_spc.to_rest()
    assert_quantity_allclose(rest_spc, rest_line_wls)

    barycentric_spc = observed_spc._change_observer_to(ICRS(CartesianRepresentation([0, 0, 0] * u.au)))
    baryrest_spc = barycentric_spc.to_rest()
    assert quantity_allclose(baryrest_spc, rest_line_wls)

    # now make sure the change the barycentric shift did is comparable to the
    # offset rv_correction produces
    # barytarg = SkyCoord(barycentric_spc.target.frame)  # should be this but that doesn't work for unclear reasons
    barytarg = SkyCoord(barycentric_spc.target.data.without_differentials(),
                        frame=barycentric_spc.target.realize_frame(None))
    vcorr = barytarg.radial_velocity_correction(kind='barycentric',
                                                obstime=obstime, location=eloc)

    drv = baryrest_spc.radial_velocity - observed_spc.radial_velocity

    # note this probably will not work on the first try, but it's ok if this is
    # "good enough", where good enough is estimated below.  But that could be
    # adjusted if we think that's too agressive of a precision target for what
    # the machinery can handle
    # with pytest.raises(AssertionError):
    assert_quantity_allclose(vcorr, drv, atol=10*u.m/u.s)


def test_observer_init_rv_behavior():
    """
    Test basic initialization behavior or observer/target and redshift/rv
    """
    sc_init = SpectralCoord([4000, 5000]*u.angstrom,
                            observer=None, target=None,
                            radial_velocity=100*u.km/u.s)
    assert sc_init.observer is None
    assert sc_init.target is None
    assert_quantity_allclose(sc_init.radial_velocity, 100*u.km/u.s)

    sc_init.observer = GCRS(CartesianRepresentation([0*u.km, 0*u.km, 0*u.km]))
    assert sc_init.observer is not None
    assert_quantity_allclose(sc_init.radial_velocity, 100*u.km/u.s)

    sc_init.target = SkyCoord(CartesianRepresentation([0*u.km, 0*u.km, 0*u.km]),
                              frame='icrs')
    assert sc_init.target is not None
    assert_quantity_allclose(sc_init.radial_velocity, 0.20502225*u.km/u.s)

    with pytest.raises(ValueError):
        # cannot reset after setting once as above
        sc_init.observer = GCRS(CartesianRepresentation([0*u.km, 1*u.km, 0*u.km]))

    with pytest.raises(ValueError):
        # same with target
        sc_init.target = GCRS(CartesianRepresentation([0*u.km, 1*u.km, 0*u.km]))


def test_rv_redshift_initialization():
    sc_init = SpectralCoord([4000, 5000]*u.angstrom, redshift=1)
    assert_quantity_allclose(sc_init.redshift, 1*u.dimensionless_unscaled)
    assert_quantity_allclose(sc_init.radial_velocity, c)

    sc_init2 = SpectralCoord([4000, 5000]*u.angstrom, radial_velocity=c)
    assert_quantity_allclose(sc_init2.redshift, 1*u.dimensionless_unscaled)
    assert_quantity_allclose(sc_init2.radial_velocity, c)

    sc_init3 = SpectralCoord([4000, 5000]*u.angstrom, redshift=1*u.dimensionless_unscaled)
    assert sc_init.redshift == sc_init3.redshift

    with pytest.raises(ValueError):
        # can't set conflicting rv and redshift values
        SpectralCoord([4000, 5000]*u.angstrom,
                      radial_velocity=10*u.km/u.s,
                      redshift=2)


def test_with_rvredshift():
    sc_init = SpectralCoord([4000, 5000]*u.angstrom, redshift=1)

    sc_set_rv = sc_init.with_redshift(.5)
    assert_quantity_allclose(sc_set_rv.radial_velocity, c/2)

    sc_set_rv = sc_init.with_radial_velocity(c/2)
    assert_quantity_allclose(sc_set_rv.redshift, .5)

    gcrs_origin = GCRS(CartesianRepresentation([0*u.km, 0*u.km, 0*u.km]))
    sc_init2 = SpectralCoord([4000, 5000]*u.angstrom, redshift=1,
                             observer=gcrs_origin)
    sc_init2.with_redshift(.5)

    sc_init3 = SpectralCoord([4000, 5000]*u.angstrom, redshift=1,
                             target=gcrs_origin)
    sc_init3.with_redshift(.5)

    sc_init4 = SpectralCoord([4000, 5000]*u.angstrom, redshift=1,
                             observer=gcrs_origin, target=gcrs_origin)
    with pytest.raises(ValueError):
        # fails if both observer and target are set
        sc_init4.with_redshift(.5)


gcrs_origin = GCRS(CartesianRepresentation([0*u.km, 0*u.km, 0*u.km]))
gcrs_not_origin = GCRS(CartesianRepresentation([1*u.km, 0*u.km, 0*u.km]))
@pytest.mark.parametrize("sc_kwargs", [
                         dict(radial_velocity=0*u.km/u.s),
                         dict(observer=gcrs_origin, radial_velocity=0*u.km/u.s),
                         dict(target=gcrs_origin, radial_velocity=0*u.km/u.s),
                         dict(observer=gcrs_origin, target=gcrs_not_origin)])
def test_los_shift(sc_kwargs):
    wl = [4000, 5000]*u.angstrom
    sc_init = SpectralCoord(wl, **sc_kwargs)

    # these should always work in *all* cases because it's unambiguous that
    # a target shift should behave this way
    new_sc1 = sc_init.with_los_shift(.1)
    assert_quantity_allclose(new_sc1, wl*1.1)
    new_sc2 = sc_init.with_los_shift(.1*u.dimensionless_unscaled)  # interpret at redshift
    assert_quantity_allclose(new_sc1, new_sc2)

    new_sc3 = sc_init.with_los_shift(-100*u.km/u.s)
    assert_quantity_allclose(new_sc3, wl*(1 + (-100*u.km/u.s / c)))

    # now try the cases where observer is specified as well/instead
    if sc_init.observer is None or sc_init.target is None:
        with pytest.raises(ValueError):
            # both must be specified if you're going to mess with observer
            sc_init.with_los_shift(observer_shift=.1)

    if sc_init.observer is not None and sc_init.target is not None:
        # redshifting the observer should *blushift* the LOS velocity since
        # its the observer-to-target vector that matters
        new_sc4 = sc_init.with_los_shift(observer_shift=.1)
        assert_quantity_allclose(new_sc4, wl*.9)

        # an equal shift in both should produce no offset at all
        new_sc5 = sc_init.with_los_shift(target_shift=.1, observer_shift=.1)
        assert_quantity_allclose(new_sc5, wl)


def test_asteroid_velocity_frame_shifts():
    """
    This test mocks up the use case of observing a spectrum of an asteroid
    at different times and from different observer locations.
    """
    time1 = time.Time('2018-12-13 9:00')
    dt = 12*u.hour
    time2 = time1 + dt

    # make the silly but simplifying assumption that the astroid is moving along
    # the x-axis of GCRS, and makes a 10 earth-radius closest approach

    v_ast = [5, 0, 0]*u.km/u.s
    x1 = -v_ast[0]*dt / 2
    x2 = v_ast[0]*dt / 2
    z = 10*u.Rearth

    cdiff = CartesianDifferential(v_ast)

    asteroid_loc1 = GCRS(CartesianRepresentation(x1.to(u.km),
                                                 0*u.km,
                                                 z.to(u.km),
                                                 differentials=cdiff),
                                                 obstime=time1)
    asteroid_loc2 = GCRS(CartesianRepresentation(x2.to(u.km),
                                                 0*u.km,
                                                 z.to(u.km),
                                                 differentials=cdiff),
                                                 obstime=time2)

    # assume satellites that are essentially fixed in geostationary orbit on
    # opposite sides of the earth
    observer1 = GCRS(CartesianRepresentation([0*u.km, 35000*u.km, 0*u.km]),
                                             obstime=time1)
    observer2 = GCRS(CartesianRepresentation([0*u.km, -35000*u.km, 0*u.km]),
                                             obstime=time2)

    wls = np.linspace(4000, 7000, 100) * u.angstrom
    spec_coord1 = SpectralCoord(wls, observer=observer1, target=asteroid_loc1)

    assert spec_coord1.radial_velocity < 0*u.km/u.s
    assert spec_coord1.radial_velocity > -5*u.km/u.s

    spec_coord2 = SpectralCoord(wls, observer=observer2, target=asteroid_loc2)

    assert spec_coord2.radial_velocity > 0*u.km/u.s
    assert spec_coord2.radial_velocity < 5*u.km/u.s

    # now check the behavior of in_observer_velocity_frame: we shift each coord
    # into the velocity frame of its *own* target.  That would then be a
    # spectralcoord that would allow direct physical comparison of the two
    # diffferent spec_corrds.  There's no way to test that, without
    # actual data, though.

    # spec_coord2 is redshifted, so we test that it behaves the way "shifting
    # to rest frame" should - the as-observed spectral coordinate should become
    # the rest frame, so something that starts out red should become bluer
    target_sc2 = spec_coord2.in_observer_velocity_frame(spec_coord2.target)
    assert np.all(target_sc2 < spec_coord2)
    # rv/redshift should be 0 since the observer and target velocities should
    # be the same
    assert_quantity_allclose(target_sc2.radial_velocity, 0*u.km/u.s,
                             atol=1e-7 * u.km / u.s)

    # check that the same holds for spec_coord1, but be more specific: it
    # should follow the standard redshift formula (which in this case yields
    # a blueshift, although the formula is the same as 1+z)
    target_sc1 = spec_coord1.in_observer_velocity_frame(spec_coord1.target)
    assert_quantity_allclose(target_sc1, spec_coord1/(1+spec_coord1.redshift))

    # TODO: Figure out what is meant by the below use case
    # ensure the "target-rest" use gives the same answer
    # target_sc1_alt = spec_coord1.in_observer_velocity_frame('target-rest')
    # assert_quantity_allclose(target_sc1, target_sc1_alt)
