from contextlib import nullcontext

import astropy.units as u
import numpy as np
from numpy.testing import assert_allclose

import pytest
from astropy import time
from astropy.constants import c
from astropy.table import Table
from astropy.time import Time
from astropy.utils import iers
from astropy.coordinates import (SkyCoord, EarthLocation, ICRS, GCRS, Galactic,
                                 CartesianDifferential,
                                 get_body_barycentric_posvel,
                                 FK5, CartesianRepresentation,
                                 SpectralQuantity)
from astropy.tests.helper import assert_quantity_allclose, quantity_allclose
from astropy.utils.exceptions import AstropyUserWarning, AstropyWarning
from astropy.utils.data import get_pkg_data_filename

from astropy.coordinates.spectral_coordinate import SpectralCoord, _apply_relativistic_doppler_shift
from astropy.wcs.wcsapi.fitswcs import VELOCITY_FRAMES as FITSWCS_VELOCITY_FRAMES


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


# GENERAL TESTS

# We first run through a series of cases to test different ways of initializing
# the observer and target for SpectralCoord, including for example frames,
# SkyCoords, and making sure that SpectralCoord is not sensitive to the actual
# frame or representation class.

# Local Standard of Rest
LSRD = Galactic(u=0.1 * u.km, v=0.1 * u.km, w=0.1 * u.km,
                U=9 * u.km / u.s, V=12 * u.km / u.s, W=7 * u.km / u.s,
                representation_type='cartesian', differential_type='cartesian')

LSRD_EQUIV = [
              LSRD,
              SkyCoord(LSRD),  # as a SkyCoord
              LSRD.transform_to(ICRS()),  # different frame
              LSRD.transform_to(ICRS()).transform_to(Galactic())  # different representation
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

    with nullcontext() if target is None else pytest.warns(AstropyUserWarning, match='No velocity defined on frame'):
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
    elif (any(observer is lsrd for lsrd in LSRD_EQUIV)
          and any(target is lsrd for lsrd in LSRD_DIR_STATIONARY_EQUIV)):
        assert_quantity_allclose(coord.radial_velocity, -274 ** 0.5 * u.km / u.s, atol=1e-4 * u.km / u.s)
        assert_quantity_allclose(coord.redshift, -5.5213158163147646e-05, atol=1e-9)
    else:
        raise NotImplementedError()


def test_create_from_spectral_coord(observer, target):
    """
    Checks that parameters are correctly copied to the new SpectralCoord object
    """
    with nullcontext() if target is None else pytest.warns(AstropyUserWarning, match='No velocity defined on frame'):
        spec_coord1 = SpectralCoord([100, 200, 300] * u.nm, observer=observer,
                                    target=target, doppler_convention='optical',
                                    doppler_rest=6000*u.AA)
    spec_coord2 = SpectralCoord(spec_coord1)
    assert spec_coord1.observer == spec_coord2.observer
    assert spec_coord1.target == spec_coord2.target
    assert spec_coord1.radial_velocity == spec_coord2.radial_velocity
    assert spec_coord1.doppler_convention == spec_coord2.doppler_convention
    assert spec_coord1.doppler_rest == spec_coord2.doppler_rest


# INTERNAL FUNCTIONS TESTS


def test_apply_relativistic_doppler_shift():

    # Frequency
    sq1 = SpectralQuantity(1 * u.GHz)
    sq2 = _apply_relativistic_doppler_shift(sq1, 0.5 * c)
    assert_quantity_allclose(sq2, np.sqrt(1. / 3.) * u.GHz)

    # Wavelength
    sq3 = SpectralQuantity(500 * u.nm)
    sq4 = _apply_relativistic_doppler_shift(sq3, 0.5 * c)
    assert_quantity_allclose(sq4, np.sqrt(3) * 500 * u.nm)

    # Energy
    sq5 = SpectralQuantity(300 * u.eV)
    sq6 = _apply_relativistic_doppler_shift(sq5, 0.5 * c)
    assert_quantity_allclose(sq6, np.sqrt(1. / 3.) * 300 * u.eV)

    # Wavenumber
    sq7 = SpectralQuantity(0.01 / u.micron)
    sq8 = _apply_relativistic_doppler_shift(sq7, 0.5 * c)
    assert_quantity_allclose(sq8, np.sqrt(1. / 3.) * 0.01 / u.micron)

    # Velocity (doppler_convention='relativistic')
    sq9 = SpectralQuantity(200 * u.km / u.s, doppler_convention='relativistic', doppler_rest=1 * u.GHz)
    sq10 = _apply_relativistic_doppler_shift(sq9, 300 * u.km / u.s)
    assert_quantity_allclose(sq10, 499.999666 * u.km / u.s)
    assert sq10.doppler_convention == 'relativistic'

    # Velocity (doppler_convention='optical')
    sq11 = SpectralQuantity(200 * u.km / u.s, doppler_convention='radio', doppler_rest=1 * u.GHz)
    sq12 = _apply_relativistic_doppler_shift(sq11, 300 * u.km / u.s)
    assert_quantity_allclose(sq12, 499.650008 * u.km / u.s)
    assert sq12.doppler_convention == 'radio'

    # Velocity (doppler_convention='radio')
    sq13 = SpectralQuantity(200 * u.km / u.s, doppler_convention='optical', doppler_rest=1 * u.GHz)
    sq14 = _apply_relativistic_doppler_shift(sq13, 300 * u.km / u.s)
    assert_quantity_allclose(sq14, 500.350493 * u.km / u.s)
    assert sq14.doppler_convention == 'optical'

    # Velocity - check relativistic velocity addition
    sq13 = SpectralQuantity(0 * u.km / u.s, doppler_convention='relativistic', doppler_rest=1 * u.GHz)
    sq14 = _apply_relativistic_doppler_shift(sq13, 0.999 * c)
    assert_quantity_allclose(sq14, 0.999 * c)
    sq14 = _apply_relativistic_doppler_shift(sq14, 0.999 * c)
    assert_quantity_allclose(sq14, (0.999 * 2) / (1 + 0.999**2) * c)
    assert sq14.doppler_convention == 'relativistic'

    # Cases that should raise errors
    sq15 = SpectralQuantity(200 * u.km / u.s)
    with pytest.raises(ValueError, match='doppler_convention not set'):
        _apply_relativistic_doppler_shift(sq15, 300 * u.km / u.s)
    sq16 = SpectralQuantity(200 * u.km / u.s, doppler_rest=10 * u.GHz)
    with pytest.raises(ValueError, match='doppler_convention not set'):
        _apply_relativistic_doppler_shift(sq16, 300 * u.km / u.s)
    sq17 = SpectralQuantity(200 * u.km / u.s, doppler_convention='optical')
    with pytest.raises(ValueError, match='doppler_rest not set'):
        _apply_relativistic_doppler_shift(sq17, 300 * u.km / u.s)


# BASIC TESTS


def test_init_quantity():
    sc = SpectralCoord(10 * u.GHz)
    assert sc.value == 10.
    assert sc.unit is u.GHz
    assert sc.doppler_convention is None
    assert sc.doppler_rest is None
    assert sc.observer is None
    assert sc.target is None


def test_init_spectral_quantity():
    sc = SpectralCoord(SpectralQuantity(10 * u.GHz, doppler_convention='optical'))
    assert sc.value == 10.
    assert sc.unit is u.GHz
    assert sc.doppler_convention == 'optical'
    assert sc.doppler_rest is None
    assert sc.observer is None
    assert sc.target is None


def test_init_too_many_args():

    with pytest.raises(ValueError, match='Cannot specify radial velocity or redshift if both'):
        SpectralCoord(10 * u.GHz, observer=LSRD, target=SkyCoord(10, 20, unit='deg'),
                      radial_velocity=1 * u.km / u.s)

    with pytest.raises(ValueError, match='Cannot specify radial velocity or redshift if both'):
        SpectralCoord(10 * u.GHz, observer=LSRD, target=SkyCoord(10, 20, unit='deg'),
                      redshift=1)

    with pytest.raises(ValueError, match='Cannot set both a radial velocity and redshift'):
        SpectralCoord(10 * u.GHz, radial_velocity=1 * u.km / u.s, redshift=1)


def test_init_wrong_type():

    with pytest.raises(TypeError, match='observer must be a SkyCoord or coordinate frame instance'):
        SpectralCoord(10 * u.GHz, observer=3.4)

    with pytest.raises(TypeError, match='target must be a SkyCoord or coordinate frame instance'):
        SpectralCoord(10 * u.GHz, target=3.4)

    with pytest.raises(u.UnitsError, match="Argument 'radial_velocity' to function "
                                           "'__new__' must be in units convertible to 'km / s'"):
        SpectralCoord(10 * u.GHz, radial_velocity=1 * u.kg)

    with pytest.raises(TypeError, match="Argument 'radial_velocity' to function "
                                        "'__new__' has no 'unit' attribute. You should "
                                        "pass in an astropy Quantity instead."):
        SpectralCoord(10 * u.GHz, radial_velocity='banana')

    with pytest.raises(u.UnitsError, match='redshift should be dimensionless'):
        SpectralCoord(10 * u.GHz, redshift=1 * u.m)

    with pytest.raises(TypeError, match='Cannot parse "banana" as a Quantity. It does not start with a number.'):
        SpectralCoord(10 * u.GHz, redshift='banana')


def test_observer_init_rv_behavior():
    """
    Test basic initialization behavior or observer/target and redshift/rv
    """

    # Start off by specifying the radial velocity only
    sc_init = SpectralCoord([4000, 5000]*u.AA,
                            radial_velocity=100*u.km/u.s)
    assert sc_init.observer is None
    assert sc_init.target is None
    assert_quantity_allclose(sc_init.radial_velocity, 100*u.km/u.s)

    # Next, set the observer, and check that the radial velocity hasn't changed
    with pytest.warns(AstropyUserWarning, match='No velocity defined on frame'):
        sc_init.observer = ICRS(CartesianRepresentation([0*u.km, 0*u.km, 0*u.km]))
    assert sc_init.observer is not None
    assert_quantity_allclose(sc_init.radial_velocity, 100*u.km/u.s)

    # Setting the target should now cause the original radial velocity to be
    # dropped in favor of the automatically computed one
    sc_init.target = SkyCoord(CartesianRepresentation([1*u.km, 0*u.km, 0*u.km]),
                              frame='icrs', radial_velocity=30 * u.km / u.s)
    assert sc_init.target is not None
    assert_quantity_allclose(sc_init.radial_velocity, 30 * u.km / u.s)

    # The observer can only be set if originally None - now that it isn't
    # setting it again should fail
    with pytest.raises(ValueError, match='observer has already been set'):
        sc_init.observer = GCRS(CartesianRepresentation([0*u.km, 1*u.km, 0*u.km]))

    # And similarly, changing the target should not be possible
    with pytest.raises(ValueError, match='target has already been set'):
        sc_init.target = GCRS(CartesianRepresentation([0*u.km, 1*u.km, 0*u.km]))


def test_rv_redshift_initialization():

    # Check that setting the redshift sets the radial velocity appropriately,
    # and that the redshift can be recovered
    sc_init = SpectralCoord([4000, 5000]*u.AA, redshift=1)
    assert isinstance(sc_init.redshift, u.Quantity)
    assert_quantity_allclose(sc_init.redshift, 1*u.dimensionless_unscaled)
    assert_quantity_allclose(sc_init.radial_velocity, 0.6 * c)

    # Check that setting the same radial velocity produces the same redshift
    # and that the radial velocity can be recovered
    sc_init2 = SpectralCoord([4000, 5000]*u.AA, radial_velocity=0.6 * c)
    assert_quantity_allclose(sc_init2.redshift, 1*u.dimensionless_unscaled)
    assert_quantity_allclose(sc_init2.radial_velocity, 0.6 * c)

    # Check that specifying redshift as a quantity works
    sc_init3 = SpectralCoord([4000, 5000]*u.AA, redshift=1 * u.one)
    assert sc_init.redshift == sc_init3.redshift

    # Make sure that both redshift and radial velocity can't be specified at
    # the same time.
    with pytest.raises(ValueError, match='Cannot set both a radial velocity and redshift'):
        SpectralCoord([4000, 5000]*u.AA,
                      radial_velocity=10*u.km/u.s,
                      redshift=2)


def test_replicate():

    # The replicate method makes a new object with attributes updated, but doesn't
    # do any conversion

    sc_init = SpectralCoord([4000, 5000]*u.AA, redshift=2)

    sc_set_rv = sc_init.replicate(redshift=1)
    assert_quantity_allclose(sc_set_rv.radial_velocity, 0.6 * c)
    assert_quantity_allclose(sc_init, [4000, 5000] * u.AA)

    sc_set_rv = sc_init.replicate(radial_velocity=c / 2)
    assert_quantity_allclose(sc_set_rv.redshift, np.sqrt(3) - 1)
    assert_quantity_allclose(sc_init, [4000, 5000] * u.AA)

    gcrs_origin = GCRS(CartesianRepresentation([0*u.km, 0*u.km, 0*u.km]))
    with pytest.warns(AstropyUserWarning, match='No velocity defined on frame'):
        sc_init2 = SpectralCoord([4000, 5000]*u.AA, redshift=1,
                                 observer=gcrs_origin)
    with np.errstate(all='ignore'):
        sc_init2.replicate(redshift=.5)
    assert_quantity_allclose(sc_init2, [4000, 5000] * u.AA)

    with pytest.warns(AstropyUserWarning, match='No velocity defined on frame'):
        sc_init3 = SpectralCoord([4000, 5000]*u.AA, redshift=1,
                                 target=gcrs_origin)
    with np.errstate(all='ignore'):
        sc_init3.replicate(redshift=.5)
    assert_quantity_allclose(sc_init2, [4000, 5000] * u.AA)

    with pytest.warns(AstropyUserWarning, match='No velocity defined on frame'):
        sc_init4 = SpectralCoord([4000, 5000]*u.AA,
                                 observer=gcrs_origin, target=gcrs_origin)
    with pytest.raises(ValueError, match='Cannot specify radial velocity or redshift if both target and observer are specified'):
        sc_init4.replicate(redshift=.5)

    sc_init = SpectralCoord([4000, 5000]*u.AA, redshift=2)
    sc_init_copy = sc_init.replicate(copy=True)
    sc_init[0] = 6000 * u.AA
    assert_quantity_allclose(sc_init_copy, [4000, 5000] * u.AA)

    sc_init = SpectralCoord([4000, 5000]*u.AA, redshift=2)
    sc_init_ref = sc_init.replicate()
    sc_init[0] = 6000 * u.AA
    assert_quantity_allclose(sc_init_ref, [6000, 5000] * u.AA)


def test_with_observer_stationary_relative_to():

    # Simple tests of with_observer_stationary_relative_to to cover different
    # ways of calling it

    # The replicate method makes a new object with attributes updated, but doesn't
    # do any conversion

    sc1 = SpectralCoord([4000, 5000]*u.AA)
    with pytest.raises(ValueError, match='This method can only be used if both '
                                         'observer and target are defined on the '
                                         'SpectralCoord'):
        sc1.with_observer_stationary_relative_to('icrs')

    sc2 = SpectralCoord([4000, 5000] * u.AA,
                        observer=ICRS(0 * u.km, 0 * u.km, 0 * u.km,
                                      -1 * u.km / u.s, 0 * u.km / u.s, -1 * u.km / u.s,
                                      representation_type='cartesian',
                                      differential_type='cartesian'),
                        target=ICRS(0 * u.deg, 45 * u.deg, distance=1 * u.kpc, radial_velocity=2 * u.km / u.s))

    # Motion of observer is in opposite direction to target
    assert_quantity_allclose(sc2.radial_velocity, (2 + 2 ** 0.5) * u.km / u.s)

    # Change to observer that is stationary in ICRS
    sc3 = sc2.with_observer_stationary_relative_to('icrs')

    # Velocity difference is now pure radial velocity of target
    assert_quantity_allclose(sc3.radial_velocity, 2 * u.km / u.s)

    # Check setting the velocity in with_observer_stationary_relative_to
    sc4 = sc2.with_observer_stationary_relative_to('icrs', velocity=[-2**0.5, 0, -2**0.5] * u.km / u.s)

    # Observer once again moving away from target but faster
    assert_quantity_allclose(sc4.radial_velocity, 4 * u.km / u.s)

    # Check that we can also pass frame classes instead of names

    sc5 = sc2.with_observer_stationary_relative_to(ICRS, velocity=[-2**0.5, 0, -2**0.5] * u.km / u.s)
    assert_quantity_allclose(sc5.radial_velocity, 4 * u.km / u.s)

    # And make sure we can also pass instances of classes without data

    sc6 = sc2.with_observer_stationary_relative_to(ICRS(), velocity=[-2**0.5, 0, -2**0.5] * u.km / u.s)
    assert_quantity_allclose(sc6.radial_velocity, 4 * u.km / u.s)

    # And with data provided no velocities are present

    sc7 = sc2.with_observer_stationary_relative_to(ICRS(0 * u.km, 0 * u.km, 0 * u.km,
                                                        representation_type='cartesian'),
                                                   velocity=[-2**0.5, 0, -2**0.5] * u.km / u.s)
    assert_quantity_allclose(sc7.radial_velocity, 4 * u.km / u.s)

    # And also have the ability to pass frames with velocities already defined

    sc8 = sc2.with_observer_stationary_relative_to(ICRS(0 * u.km, 0 * u.km, 0 * u.km,
                                                        2**0.5 * u.km / u.s, 0 * u.km / u.s, 2**0.5 * u.km / u.s,
                                                        representation_type='cartesian',
                                                        differential_type='cartesian'))
    assert_quantity_allclose(sc8.radial_velocity, 0 * u.km / u.s, atol=1e-10 * u.km / u.s)

    # Make sure that things work properly if passing a SkyCoord

    sc9 = sc2.with_observer_stationary_relative_to(SkyCoord(ICRS(0 * u.km, 0 * u.km, 0 * u.km,
                                                                 representation_type='cartesian')),
                                                   velocity=[-2**0.5, 0, -2**0.5] * u.km / u.s)
    assert_quantity_allclose(sc9.radial_velocity, 4 * u.km / u.s)

    sc10 = sc2.with_observer_stationary_relative_to(SkyCoord(ICRS(0 * u.km, 0 * u.km, 0 * u.km,
                                                                  2**0.5 * u.km / u.s, 0 * u.km / u.s, 2**0.5 * u.km / u.s,
                                                                  representation_type='cartesian',
                                                                  differential_type='cartesian')))
    assert_quantity_allclose(sc10.radial_velocity, 0 * u.km / u.s, atol=1e-10 * u.km / u.s)

    # But we shouldn't be able to pass both a frame with velocities, and explicit velocities

    with pytest.raises(ValueError, match='frame already has differentials, cannot also specify velocity'):
        sc2.with_observer_stationary_relative_to(ICRS(0 * u.km, 0 * u.km, 0 * u.km,
                                                 2**0.5 * u.km / u.s, 0 * u.km / u.s, 2**0.5 * u.km / u.s,
                                                 representation_type='cartesian',
                                                 differential_type='cartesian'),
                                                 velocity=[-2**0.5, 0, -2**0.5] * u.km / u.s)

    # And velocities should have three elements

    with pytest.raises(ValueError, match='velocity should be a Quantity vector with 3 elements'):
        sc2.with_observer_stationary_relative_to(ICRS, velocity=[-2**0.5, 0, -2**0.5, -3] * u.km / u.s)

    # Make sure things don't change depending on what frame class is used for reference
    sc11 = sc2.with_observer_stationary_relative_to(SkyCoord(ICRS(0 * u.km, 0 * u.km, 0 * u.km,
                                                                  2**0.5 * u.km / u.s, 0 * u.km / u.s, 2**0.5 * u.km / u.s,
                                                                  representation_type='cartesian',
                                                                  differential_type='cartesian')).transform_to(Galactic))
    assert_quantity_allclose(sc11.radial_velocity, 0 * u.km / u.s, atol=1e-10 * u.km / u.s)

    # Check that it is possible to preserve the observer frame
    sc12 = sc2.with_observer_stationary_relative_to(LSRD)
    sc13 = sc2.with_observer_stationary_relative_to(LSRD, preserve_observer_frame=True)

    assert isinstance(sc12.observer, Galactic)
    assert isinstance(sc13.observer, ICRS)


def test_los_shift_radial_velocity():

    # Tests to make sure that with_radial_velocity_shift correctly calculates
    # the new radial velocity

    # First check case where observer and/or target aren't specified

    sc1 = SpectralCoord(500 * u.nm, radial_velocity=1 * u.km / u.s)

    sc2 = sc1.with_radial_velocity_shift(1 * u.km / u.s)
    assert_quantity_allclose(sc2.radial_velocity, 2 * u.km / u.s)

    sc3 = sc1.with_radial_velocity_shift(-3 * u.km / u.s)
    assert_quantity_allclose(sc3.radial_velocity, -2 * u.km / u.s)

    with pytest.warns(AstropyUserWarning, match='No velocity defined on frame'):
        sc4 = SpectralCoord(500 * u.nm, radial_velocity=1 * u.km / u.s, observer=gcrs_not_origin)

    sc5 = sc4.with_radial_velocity_shift(1 * u.km / u.s)
    assert_quantity_allclose(sc5.radial_velocity, 2 * u.km / u.s)

    sc6 = sc4.with_radial_velocity_shift(-3 * u.km / u.s)
    assert_quantity_allclose(sc6.radial_velocity, -2 * u.km / u.s)

    with pytest.warns(AstropyUserWarning, match='No velocity defined on frame'):
        sc7 = SpectralCoord(500 * u.nm, radial_velocity=1 * u.km / u.s, target=ICRS(10 * u.deg, 20 * u.deg))

    sc8 = sc7.with_radial_velocity_shift(1 * u.km / u.s)
    assert_quantity_allclose(sc8.radial_velocity, 2 * u.km / u.s)

    sc9 = sc7.with_radial_velocity_shift(-3 * u.km / u.s)
    assert_quantity_allclose(sc9.radial_velocity, -2 * u.km / u.s)

    # Check that things still work when both observer and target are specified

    with pytest.warns(AstropyUserWarning, match='No velocity defined on frame'):
        sc10 = SpectralCoord(500 * u.nm,
                             observer=ICRS(0 * u.deg, 0 * u.deg, distance=1 * u.m),
                             target=ICRS(10 * u.deg, 20 * u.deg,
                                         radial_velocity=1 * u.km / u.s,
                                         distance=10 * u.kpc))

    sc11 = sc10.with_radial_velocity_shift(1 * u.km / u.s)
    assert_quantity_allclose(sc11.radial_velocity, 2 * u.km / u.s)

    sc12 = sc10.with_radial_velocity_shift(-3 * u.km / u.s)
    assert_quantity_allclose(sc12.radial_velocity, -2 * u.km / u.s)

    # Check that things work if radial_velocity wasn't specified at all

    sc13 = SpectralCoord(500 * u.nm)
    sc14 = sc13.with_radial_velocity_shift(1 * u.km / u.s)
    assert_quantity_allclose(sc14.radial_velocity, 1 * u.km / u.s)

    sc15 = sc1.with_radial_velocity_shift()
    assert_quantity_allclose(sc15.radial_velocity, 1 * u.km / u.s)

    # Check that units are verified

    with pytest.raises(u.UnitsError, match="Argument must have unit physical "
                                           "type 'speed' for radial velocty or "
                                           "'dimensionless' for redshift."):
        sc1.with_radial_velocity_shift(target_shift=1 * u.kg)


@pytest.mark.xfail
def test_relativistic_radial_velocity():

    # Test for when both observer and target have relativistic velocities.
    # This is not yet supported, so the test is xfailed for now.

    sc = SpectralCoord(500 * u.nm,
                       observer=ICRS(0 * u.km, 0 * u.km, 0 * u.km,
                                     -0.5 * c, -0.5 * c, -0.5 * c,
                                     representation_type='cartesian',
                                     differential_type='cartesian'),
                       target=ICRS(1 * u.kpc, 1 * u.kpc, 1 * u.kpc,
                                   0.5 * c, 0.5 * c, 0.5 * c,
                                   representation_type='cartesian',
                                   differential_type='cartesian'))

    assert_quantity_allclose(sc.radial_velocity, 0.989743318610787 * u.km / u.s)


# SCIENCE USE CASE TESTS


def test_spectral_coord_jupiter():
    """
    Checks radial velocity between Earth and Jupiter
    """
    obstime = time.Time('2018-12-13 9:00')
    obs = get_greenwich_earthlocation().get_gcrs(obstime)

    pos, vel = get_body_barycentric_posvel('jupiter', obstime)
    jupiter = SkyCoord(pos.with_differentials(CartesianDifferential(vel.xyz)), obstime=obstime)

    spc = SpectralCoord([100, 200, 300] * u.nm, observer=obs, target=jupiter)

    # The velocity should be less than ~43 + a bit extra, which is the
    # maximum possible earth-jupiter relative velocity. We check the exact
    # value here (determined from SpectralCoord, so this serves as a test to
    # check that this value doesn't change - the value is not a ground truth)
    assert_quantity_allclose(spc.radial_velocity, -7.35219854 * u.km / u.s)


def test_spectral_coord_alphacen():
    """
    Checks radial velocity between Earth and Alpha Centauri
    """
    obstime = time.Time('2018-12-13 9:00')
    obs = get_greenwich_earthlocation().get_gcrs(obstime)

    # Coordinates were obtained from the following then hard-coded to avoid download
    # acen = SkyCoord.from_name('alpha cen')
    acen = SkyCoord(ra=219.90085*u.deg, dec=-60.83562*u.deg, frame='icrs',
                    distance=4.37*u.lightyear, radial_velocity=-18.*u.km/u.s)

    spc = SpectralCoord([100, 200, 300] * u.nm, observer=obs, target=acen)

    # The velocity should be less than ~18 + 30 + a bit extra, which is the
    # maximum possible relative velocity.  We check the exact value here
    # (determined from SpectralCoord, so this serves as a test to check that
    # this value doesn't change - the value is not a ground truth)
    assert_quantity_allclose(spc.radial_velocity, -26.328301 * u.km / u.s)


def test_spectral_coord_m31():
    """
    Checks radial velocity between Earth and M31
    """
    obstime = time.Time('2018-12-13 9:00')
    obs = get_greenwich_earthlocation().get_gcrs(obstime)

    # Coordinates were obtained from the following then hard-coded to avoid download
    # m31 = SkyCoord.from_name('M31')
    m31 = SkyCoord(ra=10.6847*u.deg, dec=41.269*u.deg,
                   distance=710*u.kpc, radial_velocity=-300*u.km/u.s)

    spc = SpectralCoord([100, 200, 300] * u.nm, observer=obs, target=m31)

    # The velocity should be less than ~300 + 30 + a bit extra in km/s, which
    # is the maximum possible relative velocity.  We check the exact values
    # here (determined from SpectralCoord, so this serves as a test to check
    # that this value doesn't change - the value is not a ground truth)
    assert_quantity_allclose(spc.radial_velocity, -279.755128 * u.km / u.s)
    assert_allclose(spc.redshift, -0.0009327276702120191)


def test_shift_to_rest_galaxy():
    """
    This tests storing a spectral coordinate with a specific redshift, and then
    doing basic rest-to-observed-and-back transformations
    """
    z = 5
    rest_line_wls = [5007, 6563]*u.AA

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


def test_shift_to_rest_star_withobserver():
    rv = -8.3283011*u.km/u.s
    rest_line_wls = [5007, 6563]*u.AA

    obstime = time.Time('2018-12-13 9:00')
    eloc = get_greenwich_earthlocation()
    obs = eloc.get_gcrs(obstime)
    acen = SkyCoord(ra=219.90085*u.deg, dec=-60.83562*u.deg, frame='icrs',
                    distance=4.37*u.lightyear)
    # Note that above the rv is missing from the SkyCoord.
    # That's intended, as it will instead be set in the `SpectralCoord`.  But
    # the SpectralCoord machinery should yield something comparable to test_
    # spectral_coord_alphacen

    with pytest.warns(AstropyUserWarning, match='No velocity defined on frame'):
        observed_spc = SpectralCoord(rest_line_wls*(rv/c + 1),
                                     observer=obs, target=acen)

    rest_spc = observed_spc.to_rest()
    assert_quantity_allclose(rest_spc, rest_line_wls)

    barycentric_spc = observed_spc.with_observer_stationary_relative_to('icrs')
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
    # adjusted if we think that's too aggressive of a precision target for what
    # the machinery can handle
    # with pytest.raises(AssertionError):
    assert_quantity_allclose(vcorr, drv, atol=10*u.m/u.s)


gcrs_origin = GCRS(CartesianRepresentation([0*u.km, 0*u.km, 0*u.km]))
gcrs_not_origin = GCRS(CartesianRepresentation([1*u.km, 0*u.km, 0*u.km]))


@pytest.mark.parametrize("sc_kwargs", [
                         dict(radial_velocity=0*u.km/u.s),
                         dict(observer=gcrs_origin, radial_velocity=0*u.km/u.s),
                         dict(target=gcrs_origin, radial_velocity=0*u.km/u.s),
                         dict(observer=gcrs_origin, target=gcrs_not_origin)])
def test_los_shift(sc_kwargs):
    wl = [4000, 5000]*u.AA
    with nullcontext() if 'observer' not in sc_kwargs and 'target' not in sc_kwargs else pytest.warns(AstropyUserWarning, match='No velocity defined on frame'):
        sc_init = SpectralCoord(wl, **sc_kwargs)

    # these should always work in *all* cases because it's unambiguous that
    # a target shift should behave this way
    new_sc1 = sc_init.with_radial_velocity_shift(.1)
    assert_quantity_allclose(new_sc1, wl*1.1)
    new_sc2 = sc_init.with_radial_velocity_shift(.1*u.dimensionless_unscaled)  # interpret at redshift
    assert_quantity_allclose(new_sc1, new_sc2)

    new_sc3 = sc_init.with_radial_velocity_shift(-100*u.km/u.s)
    assert_quantity_allclose(new_sc3, wl*(1 + (-100*u.km/u.s / c)))

    # now try the cases where observer is specified as well/instead
    if sc_init.observer is None or sc_init.target is None:
        with pytest.raises(ValueError):
            # both must be specified if you're going to mess with observer
            sc_init.with_radial_velocity_shift(observer_shift=.1)

    if sc_init.observer is not None and sc_init.target is not None:
        # redshifting the observer should *blushift* the LOS velocity since
        # its the observer-to-target vector that matters
        new_sc4 = sc_init.with_radial_velocity_shift(observer_shift=.1)
        assert_quantity_allclose(new_sc4, wl/1.1)

        # an equal shift in both should produce no offset at all
        new_sc5 = sc_init.with_radial_velocity_shift(target_shift=.1, observer_shift=.1)
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

    wls = np.linspace(4000, 7000, 100) * u.AA
    with pytest.warns(AstropyUserWarning, match='No velocity defined on frame'):
        spec_coord1 = SpectralCoord(wls, observer=observer1, target=asteroid_loc1)

    assert spec_coord1.radial_velocity < 0*u.km/u.s
    assert spec_coord1.radial_velocity > -5*u.km/u.s

    with pytest.warns(AstropyUserWarning, match='No velocity defined on frame'):
        spec_coord2 = SpectralCoord(wls, observer=observer2, target=asteroid_loc2)

    assert spec_coord2.radial_velocity > 0*u.km/u.s
    assert spec_coord2.radial_velocity < 5*u.km/u.s

    # now check the behavior of with_observer_stationary_relative_to: we shift each coord
    # into the velocity frame of its *own* target.  That would then be a
    # spectralcoord that would allow direct physical comparison of the two
    # different spec_corrds.  There's no way to test that, without
    # actual data, though.

    # spec_coord2 is redshifted, so we test that it behaves the way "shifting
    # to rest frame" should - the as-observed spectral coordinate should become
    # the rest frame, so something that starts out red should become bluer
    target_sc2 = spec_coord2.with_observer_stationary_relative_to(spec_coord2.target)
    assert np.all(target_sc2 < spec_coord2)
    # rv/redshift should be 0 since the observer and target velocities should
    # be the same
    assert_quantity_allclose(target_sc2.radial_velocity, 0*u.km/u.s,
                             atol=1e-7 * u.km / u.s)

    # check that the same holds for spec_coord1, but be more specific: it
    # should follow the standard redshift formula (which in this case yields
    # a blueshift, although the formula is the same as 1+z)
    target_sc1 = spec_coord1.with_observer_stationary_relative_to(spec_coord1.target)
    assert_quantity_allclose(target_sc1, spec_coord1/(1+spec_coord1.redshift))

    # TODO: Figure out what is meant by the below use case
    # ensure the "target-rest" use gives the same answer
    # target_sc1_alt = spec_coord1.with_observer_stationary_relative_to('target-rest')
    # assert_quantity_allclose(target_sc1, target_sc1_alt)


def test_spectral_coord_from_sky_coord_without_distance():
    # see https://github.com/astropy/specutils/issues/658 for issue context
    obs = SkyCoord(0 * u.m, 0 * u.m, 0 * u.m, representation_type='cartesian')
    with pytest.warns(AstropyUserWarning, match='No velocity defined on frame'):
        coord = SpectralCoord([1, 2, 3] * u.micron, observer=obs)
    # coord.target = SkyCoord.from_name('m31')  # <- original issue, but below is the same but requires no remote data access
    with pytest.warns(AstropyUserWarning, match='Distance on coordinate object is dimensionless'):
        coord.target = SkyCoord(ra=10.68470833*u.deg, dec=41.26875*u.deg)


EXPECTED_VELOCITY_FRAMES = {'geocent': 'gcrs',
                            'heliocent': 'hcrs',
                            'lsrk': 'lsrk',
                            'lsrd': 'lsrd',
                            'galactoc': FITSWCS_VELOCITY_FRAMES['GALACTOC'],
                            'localgrp': FITSWCS_VELOCITY_FRAMES['LOCALGRP']}


@pytest.mark.parametrize('specsys', list(EXPECTED_VELOCITY_FRAMES))
@pytest.mark.slow
def test_spectralcoord_accuracy(specsys):

    # This is a test to check the numerical results of transformations between
    # different velocity frames in SpectralCoord. This compares the velocity
    # shifts determined with SpectralCoord to those determined from the rv
    # package in Starlink.

    velocity_frame = EXPECTED_VELOCITY_FRAMES[specsys]

    reference_filename = get_pkg_data_filename('accuracy/data/rv.ecsv')
    reference_table = Table.read(reference_filename, format='ascii.ecsv')

    rest = 550 * u.nm

    with iers.conf.set_temp('auto_download', False):

        for row in reference_table:

            observer = EarthLocation.from_geodetic(-row['obslon'], row['obslat']).get_itrs(obstime=row['obstime'])

            with pytest.warns(AstropyUserWarning, match='No velocity defined on frame'):
                sc_topo = SpectralCoord(545 * u.nm, observer=observer, target=row['target'])

            # FIXME: A warning is emitted for dates after MJD=57754.0 even
            # though the leap second table should be valid until the end of
            # 2020.
            with nullcontext() if row['obstime'].mjd < 57754 else pytest.warns(AstropyWarning, match='Tried to get polar motions'):
                sc_final = sc_topo.with_observer_stationary_relative_to(velocity_frame)

            delta_vel = (sc_topo.to(u.km / u.s, doppler_convention='relativistic', doppler_rest=rest) -
                         sc_final.to(u.km / u.s, doppler_convention='relativistic', doppler_rest=rest))

            if specsys == 'galactoc':
                assert_allclose(delta_vel.to_value(u.km / u.s), row[specsys.lower()], atol=30)
            else:
                assert_allclose(delta_vel.to_value(u.km / u.s), row[specsys.lower()], atol=0.02, rtol=0.002)


# TODO: add test when target is not ICRS
# TODO: add test when SpectralCoord is in velocity to start with
