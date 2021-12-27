# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from copy import deepcopy

import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from astropy import units as u
from astropy.tests.helper import (assert_quantity_allclose as
                                  assert_allclose_quantity)
from astropy.utils import isiterable
from astropy.utils.exceptions import DuplicateRepresentationWarning
from astropy.coordinates.angles import Longitude, Latitude, Angle
from astropy.coordinates.distances import Distance
from astropy.coordinates.matrix_utilities import rotation_matrix
from astropy.coordinates.representation import (
    REPRESENTATION_CLASSES, DIFFERENTIAL_CLASSES, DUPLICATE_REPRESENTATIONS,
    BaseRepresentation, SphericalRepresentation, UnitSphericalRepresentation,
    SphericalCosLatDifferential, CartesianRepresentation, RadialRepresentation,
    RadialDifferential, CylindricalRepresentation,
    PhysicsSphericalRepresentation, CartesianDifferential,
    SphericalDifferential, CylindricalDifferential,
    PhysicsSphericalDifferential, UnitSphericalDifferential,
    UnitSphericalCosLatDifferential)


# create matrices for use in testing ``.transform()`` methods
matrices = {
    "rotation": rotation_matrix(-10, "z", u.deg),
    "general": np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
}


# Preserve the original REPRESENTATION_CLASSES dict so that importing
#   the test file doesn't add a persistent test subclass (LogDRepresentation)
def setup_function(func):
    func.REPRESENTATION_CLASSES_ORIG = deepcopy(REPRESENTATION_CLASSES)
    func.DUPLICATE_REPRESENTATIONS_ORIG = deepcopy(DUPLICATE_REPRESENTATIONS)


def teardown_function(func):
    REPRESENTATION_CLASSES.clear()
    REPRESENTATION_CLASSES.update(func.REPRESENTATION_CLASSES_ORIG)
    DUPLICATE_REPRESENTATIONS.clear()
    DUPLICATE_REPRESENTATIONS.update(func.DUPLICATE_REPRESENTATIONS_ORIG)


def components_equal(rep1, rep2):
    result = True
    if type(rep1) is not type(rep2):
        return False
    for component in rep1.components:
        result &= getattr(rep1, component) == getattr(rep2, component)
    return result


def components_allclose(rep1, rep2):
    result = True
    if type(rep1) is not type(rep2):
        return False
    for component in rep1.components:
        result &= u.allclose(getattr(rep1, component), getattr(rep2, component))
    return result


def representation_equal(rep1, rep2):
    result = True
    if type(rep1) is not type(rep2):
        return False
    if getattr(rep1, '_differentials', False):
        if rep1._differentials.keys() != rep2._differentials.keys():
            return False
        for key, diff1 in rep1._differentials.items():
            result &= components_equal(diff1, rep2._differentials[key])
    elif getattr(rep2, '_differentials', False):
        return False

    return result & components_equal(rep1, rep2)


def representation_equal_up_to_angular_type(rep1, rep2):
    result = True
    if type(rep1) is not type(rep2):
        return False
    if getattr(rep1, '_differentials', False):
        if rep1._differentials.keys() != rep2._differentials.keys():
            return False
        for key, diff1 in rep1._differentials.items():
            result &= components_allclose(diff1, rep2._differentials[key])
    elif getattr(rep2, '_differentials', False):
        return False

    return result & components_allclose(rep1, rep2)


class TestRadialRepresentation:

    def test_transform(self):
        """Test the ``transform`` method. Only multiplication matrices pass."""
        rep = RadialRepresentation(distance=10 * u.kpc)

        # a rotation matrix does not work
        matrix = rotation_matrix(10 * u.deg)
        with pytest.raises(ValueError, match="scaled identity matrix"):
            rep.transform(matrix)

        # only a scaled identity matrix
        matrix = 3 * np.identity(3)
        newrep = rep.transform(matrix)
        assert newrep.distance == 30 * u.kpc

        # let's also check with differentials
        dif = RadialDifferential(d_distance=-3 * u.km / u.s)
        rep = rep.with_differentials(dict(s=dif))

        newrep = rep.transform(matrix)
        assert newrep.distance == 30 * u.kpc
        assert newrep.differentials["s"].d_distance == -9 * u.km / u.s


class TestSphericalRepresentation:

    def test_name(self):
        assert SphericalRepresentation.get_name() == 'spherical'
        assert SphericalRepresentation.get_name() in REPRESENTATION_CLASSES

    def test_empty_init(self):
        with pytest.raises(TypeError) as exc:
            s = SphericalRepresentation()

    def test_init_quantity(self):

        s3 = SphericalRepresentation(lon=8 * u.hourangle, lat=5 * u.deg, distance=10 * u.kpc)
        assert s3.lon == 8. * u.hourangle
        assert s3.lat == 5. * u.deg
        assert s3.distance == 10 * u.kpc

        assert isinstance(s3.lon, Longitude)
        assert isinstance(s3.lat, Latitude)
        assert isinstance(s3.distance, Distance)

    def test_init_no_mutate_input(self):

        lon = -1 * u.hourangle
        s = SphericalRepresentation(lon=lon, lat=-1 * u.deg, distance=1 * u.kpc, copy=True)

        # The longitude component should be wrapped at 24 hours
        assert_allclose_quantity(s.lon, 23 * u.hourangle)

        # The input should not have been mutated by the constructor
        assert_allclose_quantity(lon, -1 * u.hourangle)

    def test_init_lonlat(self):

        s2 = SphericalRepresentation(Longitude(8, u.hour),
                                     Latitude(5, u.deg),
                                     Distance(10, u.kpc))

        assert s2.lon == 8. * u.hourangle
        assert s2.lat == 5. * u.deg
        assert s2.distance == 10. * u.kpc

        assert isinstance(s2.lon, Longitude)
        assert isinstance(s2.lat, Latitude)
        assert isinstance(s2.distance, Distance)

        # also test that wrap_angle is preserved
        s3 = SphericalRepresentation(Longitude(-90, u.degree,
                                               wrap_angle=180*u.degree),
                                     Latitude(-45, u.degree),
                                     Distance(1., u.Rsun))
        assert s3.lon == -90. * u.degree
        assert s3.lon.wrap_angle == 180 * u.degree

    def test_init_subclass(self):
        class Longitude180(Longitude):
            _default_wrap_angle = 180*u.degree

        s = SphericalRepresentation(Longitude180(-90, u.degree),
                                    Latitude(-45, u.degree),
                                    Distance(1., u.Rsun))
        assert isinstance(s.lon, Longitude180)
        assert s.lon == -90. * u.degree
        assert s.lon.wrap_angle == 180 * u.degree

    def test_init_array(self):

        s1 = SphericalRepresentation(lon=[8, 9] * u.hourangle,
                                     lat=[5, 6] * u.deg,
                                     distance=[1, 2] * u.kpc)

        assert_allclose(s1.lon.degree, [120, 135])
        assert_allclose(s1.lat.degree, [5, 6])
        assert_allclose(s1.distance.kpc, [1, 2])

        assert isinstance(s1.lon, Longitude)
        assert isinstance(s1.lat, Latitude)
        assert isinstance(s1.distance, Distance)

    def test_init_array_nocopy(self):

        lon = Longitude([8, 9] * u.hourangle)
        lat = Latitude([5, 6] * u.deg)
        distance = Distance([1, 2] * u.kpc)

        s1 = SphericalRepresentation(lon=lon, lat=lat, distance=distance, copy=False)

        lon[:] = [1, 2] * u.rad
        lat[:] = [3, 4] * u.arcmin
        distance[:] = [8, 9] * u.Mpc

        assert_allclose_quantity(lon, s1.lon)
        assert_allclose_quantity(lat, s1.lat)
        assert_allclose_quantity(distance, s1.distance)

    def test_init_float32_array(self):
        """Regression test against #2983"""
        lon = Longitude(np.float32([1., 2.]), u.degree)
        lat = Latitude(np.float32([3., 4.]), u.degree)
        s1 = UnitSphericalRepresentation(lon=lon, lat=lat, copy=False)
        assert s1.lon.dtype == np.float32
        assert s1.lat.dtype == np.float32
        assert s1._values['lon'].dtype == np.float32
        assert s1._values['lat'].dtype == np.float32

    def test_reprobj(self):

        s1 = SphericalRepresentation(lon=8 * u.hourangle, lat=5 * u.deg, distance=10 * u.kpc)

        s2 = SphericalRepresentation.from_representation(s1)

        assert_allclose_quantity(s2.lon, 8. * u.hourangle)
        assert_allclose_quantity(s2.lat, 5. * u.deg)
        assert_allclose_quantity(s2.distance, 10 * u.kpc)

        s3 = SphericalRepresentation(s1)

        assert representation_equal(s1, s3)

    def test_broadcasting(self):

        s1 = SphericalRepresentation(lon=[8, 9] * u.hourangle,
                                     lat=[5, 6] * u.deg,
                                     distance=10 * u.kpc)

        assert_allclose_quantity(s1.lon, [120, 135] * u.degree)
        assert_allclose_quantity(s1.lat, [5, 6] * u.degree)
        assert_allclose_quantity(s1.distance, [10, 10] * u.kpc)

    def test_broadcasting_mismatch(self):

        with pytest.raises(ValueError) as exc:
            s1 = SphericalRepresentation(lon=[8, 9, 10] * u.hourangle,
                                         lat=[5, 6] * u.deg,
                                         distance=[1, 2] * u.kpc)
        assert exc.value.args[0] == "Input parameters lon, lat, and distance cannot be broadcast"

    def test_broadcasting_and_nocopy(self):

        s1 = SphericalRepresentation(lon=[200] * u.deg,
                                     lat=[0] * u.deg,
                                     distance=[0] * u.kpc,
                                     copy=False)
        # With no copying, we should be able to modify the wrap angle of the longitude component
        s1.lon.wrap_angle = 180 * u.deg

        s2 = SphericalRepresentation(lon=[200] * u.deg,
                                     lat=0 * u.deg,
                                     distance=0 * u.kpc,
                                     copy=False)
        # We should be able to modify the wrap angle of the longitude component even if other
        # components need to be broadcasted
        s2.lon.wrap_angle = 180 * u.deg

    def test_readonly(self):

        s1 = SphericalRepresentation(lon=8 * u.hourangle,
                                     lat=5 * u.deg,
                                     distance=1. * u.kpc)

        with pytest.raises(AttributeError):
            s1.lon = 1. * u.deg

        with pytest.raises(AttributeError):
            s1.lat = 1. * u.deg

        with pytest.raises(AttributeError):
            s1.distance = 1. * u.kpc

    def test_getitem_len_iterable(self):

        s = SphericalRepresentation(lon=np.arange(10) * u.deg,
                                    lat=-np.arange(10) * u.deg,
                                    distance=1 * u.kpc)

        s_slc = s[2:8:2]

        assert_allclose_quantity(s_slc.lon, [2, 4, 6] * u.deg)
        assert_allclose_quantity(s_slc.lat, [-2, -4, -6] * u.deg)
        assert_allclose_quantity(s_slc.distance, [1, 1, 1] * u.kpc)

        assert len(s) == 10
        assert isiterable(s)

    def test_getitem_len_iterable_scalar(self):

        s = SphericalRepresentation(lon=1 * u.deg,
                                    lat=-2 * u.deg,
                                    distance=3 * u.kpc)

        with pytest.raises(TypeError):
            s_slc = s[0]
        with pytest.raises(TypeError):
            len(s)
        assert not isiterable(s)

    def test_setitem(self):
        s = SphericalRepresentation(lon=np.arange(5) * u.deg,
                                    lat=-np.arange(5) * u.deg,
                                    distance=1 * u.kpc)
        s[:2] = SphericalRepresentation(lon=10.*u.deg, lat=2.*u.deg,
                                        distance=5.*u.kpc)
        assert_allclose_quantity(s.lon, [10, 10, 2, 3, 4] * u.deg)
        assert_allclose_quantity(s.lat, [2, 2, -2, -3, -4] * u.deg)
        assert_allclose_quantity(s.distance, [5, 5, 1, 1, 1] * u.kpc)

    def test_negative_distance(self):
        """Only allowed if explicitly passed on."""
        with pytest.raises(ValueError, match='allow_negative'):
            SphericalRepresentation(10*u.deg, 20*u.deg, -10*u.m)

        s1 = SphericalRepresentation(10*u.deg, 20*u.deg,
                                     Distance(-10*u.m, allow_negative=True))

        assert s1.distance == -10.*u.m

    def test_nan_distance(self):
        """ This is a regression test: calling represent_as() and passing in the
        same class as the object shouldn't round-trip through cartesian.
        """

        sph = SphericalRepresentation(1*u.deg, 2*u.deg, np.nan*u.kpc)
        new_sph = sph.represent_as(SphericalRepresentation)
        assert_allclose_quantity(new_sph.lon, sph.lon)
        assert_allclose_quantity(new_sph.lat, sph.lat)

        dif = SphericalCosLatDifferential(1*u.mas/u.yr, 2*u.mas/u.yr,
                                          3*u.km/u.s)
        sph = sph.with_differentials(dif)
        new_sph = sph.represent_as(SphericalRepresentation)
        assert_allclose_quantity(new_sph.lon, sph.lon)
        assert_allclose_quantity(new_sph.lat, sph.lat)

    def test_raise_on_extra_arguments(self):
        with pytest.raises(TypeError, match='got multiple values'):
            SphericalRepresentation(1*u.deg, 2*u.deg, 1.*u.kpc, lat=10)

        with pytest.raises(TypeError, match='unexpected keyword.*parrot'):
            SphericalRepresentation(1*u.deg, 2*u.deg, 1.*u.kpc, parrot=10)

    def test_representation_shortcuts(self):
        """Test that shortcuts in ``represent_as`` don't fail."""
        difs = SphericalCosLatDifferential(4*u.mas/u.yr,5*u.mas/u.yr,6*u.km/u.s)
        sph = SphericalRepresentation(1*u.deg, 2*u.deg, 3*u.kpc,
                                      differentials={'s': difs})

        got = sph.represent_as(PhysicsSphericalRepresentation,
                               PhysicsSphericalDifferential)
        assert np.may_share_memory(sph.lon, got.phi)
        assert np.may_share_memory(sph.distance, got.r)
        expected = BaseRepresentation.represent_as(
            sph, PhysicsSphericalRepresentation, PhysicsSphericalDifferential)
        # equal up to angular type
        assert representation_equal_up_to_angular_type(got, expected)

        got = sph.represent_as(UnitSphericalRepresentation,
                               UnitSphericalDifferential)
        assert np.may_share_memory(sph.lon, got.lon)
        assert np.may_share_memory(sph.lat, got.lat)
        expected = BaseRepresentation.represent_as(
            sph, UnitSphericalRepresentation, UnitSphericalDifferential)
        assert representation_equal_up_to_angular_type(got, expected)

    def test_transform(self):
        """Test ``.transform()`` on rotation and general matrices."""
        # set up representation
        ds1 = SphericalDifferential(
            d_lon=[1, 2] * u.mas / u.yr, d_lat=[3, 4] * u.mas / u.yr,
            d_distance=[-5, 6] * u.km / u.s)
        s1 = SphericalRepresentation(lon=[1, 2] * u.deg, lat=[3, 4] * u.deg,
                                     distance=[5, 6] * u.kpc, differentials=ds1)

        # transform representation & get comparison (thru CartesianRep)
        s2 = s1.transform(matrices["rotation"])
        ds2 =  s2.differentials["s"]

        dexpected = SphericalDifferential.from_cartesian(
            ds1.to_cartesian(base=s1).transform(matrices["rotation"]), base=s2)

        assert_allclose_quantity(s2.lon, s1.lon + 10 * u.deg)
        assert_allclose_quantity(s2.lat, s1.lat)
        assert_allclose_quantity(s2.distance, s1.distance)
        # check differentials. they shouldn't have changed.
        assert_allclose_quantity(ds2.d_lon, ds1.d_lon)
        assert_allclose_quantity(ds2.d_lat, ds1.d_lat)
        assert_allclose_quantity(ds2.d_distance, ds1.d_distance)
        assert_allclose_quantity(ds2.d_lon, dexpected.d_lon)
        assert_allclose_quantity(ds2.d_lat, dexpected.d_lat)
        assert_allclose_quantity(ds2.d_distance, dexpected.d_distance)

        # now with a non rotation matrix
        # transform representation & get comparison (thru CartesianRep)
        s3 = s1.transform(matrices["general"])
        ds3 = s3.differentials["s"]

        expected = (s1.represent_as(CartesianRepresentation,
                                    CartesianDifferential)
                    .transform(matrices["general"])
                    .represent_as(SphericalRepresentation,
                                  SphericalDifferential))
        dexpected = expected.differentials["s"]

        assert_allclose_quantity(s3.lon, expected.lon)
        assert_allclose_quantity(s3.lat, expected.lat)
        assert_allclose_quantity(s3.distance, expected.distance)
        assert_allclose_quantity(ds3.d_lon, dexpected.d_lon)
        assert_allclose_quantity(ds3.d_lat, dexpected.d_lat)
        assert_allclose_quantity(ds3.d_distance, dexpected.d_distance)

    def test_transform_with_NaN(self):
        # all over again, but with a NaN in the distance

        ds1 = SphericalDifferential(
            d_lon=[1, 2] * u.mas / u.yr, d_lat=[3, 4] * u.mas / u.yr,
            d_distance=[-5, 6] * u.km / u.s)
        s1 = SphericalRepresentation(lon=[1, 2] * u.deg, lat=[3, 4] * u.deg,
                                     distance=[5, np.nan] * u.kpc,
                                     differentials=ds1)

        # transform representation & get comparison (thru CartesianRep)
        s2 = s1.transform(matrices["rotation"])
        ds2 =  s2.differentials["s"]

        dexpected = SphericalDifferential.from_cartesian(
            ds1.to_cartesian(base=s1).transform(matrices["rotation"]), base=s2)

        assert_allclose_quantity(s2.lon, s1.lon + 10 * u.deg)
        assert_allclose_quantity(s2.lat, s1.lat)
        assert_allclose_quantity(s2.distance, s1.distance)
        assert_allclose_quantity(ds2.d_lon, dexpected.d_lon)
        assert_allclose_quantity(ds2.d_lat, dexpected.d_lat)
        assert_allclose_quantity(ds2.d_distance, dexpected.d_distance)
        # the 2nd component is NaN since the 2nd distance is NaN
        # TODO! this will change when ``.transform`` skips Cartesian
        assert_array_equal(np.isnan(ds2.d_lon), (False, True))
        assert_array_equal(np.isnan(ds2.d_lat), (False, True))
        assert_array_equal(np.isnan(ds2.d_distance), (False, True))

        # now with a non rotation matrix
        s3 = s1.transform(matrices["general"])
        ds3 = s3.differentials["s"]

        thruC = (s1.represent_as(CartesianRepresentation,
                                CartesianDifferential)
                    .transform(matrices["general"])
                    .represent_as(SphericalRepresentation,
                              differential_class=SphericalDifferential))
        dthruC = thruC.differentials["s"]

        # s3 should not propagate Nan.
        assert_array_equal(np.isnan(s3.lon), (False, False))
        assert_array_equal(np.isnan(s3.lat), (False, False))
        assert_array_equal(np.isnan(s3.distance), (False, True))
        # ds3 does b/c currently aren't any shortcuts on the transform
        assert_array_equal(np.isnan(ds3.d_lon), (False, True))
        assert_array_equal(np.isnan(ds3.d_lat), (False, True))
        assert_array_equal(np.isnan(ds3.d_distance), (False, True))

        # through Cartesian should
        assert_array_equal(np.isnan(thruC.lon), (False, True))
        assert_array_equal(np.isnan(thruC.lat), (False, True))
        assert_array_equal(np.isnan(thruC.distance), (False, True))
        assert_array_equal(np.isnan(dthruC.d_lon), (False, True))
        assert_array_equal(np.isnan(dthruC.d_lat), (False, True))
        assert_array_equal(np.isnan(dthruC.d_distance), (False, True))
        # test that they are close on the first value
        assert_allclose_quantity(s3.lon[0], thruC.lon[0])
        assert_allclose_quantity(s3.lat[0], thruC.lat[0])
        assert_allclose_quantity(ds3.d_lon[0], dthruC.d_lon[0])
        assert_allclose_quantity(ds3.d_lat[0], dthruC.d_lat[0])


class TestUnitSphericalRepresentation:

    def test_name(self):
        assert UnitSphericalRepresentation.get_name() == 'unitspherical'
        assert UnitSphericalRepresentation.get_name() in REPRESENTATION_CLASSES

    def test_empty_init(self):
        with pytest.raises(TypeError) as exc:
            s = UnitSphericalRepresentation()

    def test_init_quantity(self):

        s3 = UnitSphericalRepresentation(lon=8 * u.hourangle, lat=5 * u.deg)
        assert s3.lon == 8. * u.hourangle
        assert s3.lat == 5. * u.deg

        assert isinstance(s3.lon, Longitude)
        assert isinstance(s3.lat, Latitude)

    def test_init_lonlat(self):

        s2 = UnitSphericalRepresentation(Longitude(8, u.hour),
                                         Latitude(5, u.deg))

        assert s2.lon == 8. * u.hourangle
        assert s2.lat == 5. * u.deg

        assert isinstance(s2.lon, Longitude)
        assert isinstance(s2.lat, Latitude)

    def test_init_array(self):

        s1 = UnitSphericalRepresentation(lon=[8, 9] * u.hourangle,
                                         lat=[5, 6] * u.deg)

        assert_allclose(s1.lon.degree, [120, 135])
        assert_allclose(s1.lat.degree, [5, 6])

        assert isinstance(s1.lon, Longitude)
        assert isinstance(s1.lat, Latitude)

    def test_init_array_nocopy(self):

        lon = Longitude([8, 9] * u.hourangle)
        lat = Latitude([5, 6] * u.deg)

        s1 = UnitSphericalRepresentation(lon=lon, lat=lat, copy=False)

        lon[:] = [1, 2] * u.rad
        lat[:] = [3, 4] * u.arcmin

        assert_allclose_quantity(lon, s1.lon)
        assert_allclose_quantity(lat, s1.lat)

    def test_reprobj(self):

        s1 = UnitSphericalRepresentation(lon=8 * u.hourangle, lat=5 * u.deg)

        s2 = UnitSphericalRepresentation.from_representation(s1)

        assert_allclose_quantity(s2.lon, 8. * u.hourangle)
        assert_allclose_quantity(s2.lat, 5. * u.deg)

        s3 = UnitSphericalRepresentation(s1)

        assert representation_equal(s3, s1)

    def test_broadcasting(self):

        s1 = UnitSphericalRepresentation(lon=[8, 9] * u.hourangle,
                                         lat=[5, 6] * u.deg)

        assert_allclose_quantity(s1.lon, [120, 135] * u.degree)
        assert_allclose_quantity(s1.lat, [5, 6] * u.degree)

    def test_broadcasting_mismatch(self):

        with pytest.raises(ValueError) as exc:
            s1 = UnitSphericalRepresentation(lon=[8, 9, 10] * u.hourangle,
                                             lat=[5, 6] * u.deg)
        assert exc.value.args[0] == "Input parameters lon and lat cannot be broadcast"

    def test_readonly(self):

        s1 = UnitSphericalRepresentation(lon=8 * u.hourangle,
                                         lat=5 * u.deg)

        with pytest.raises(AttributeError):
            s1.lon = 1. * u.deg

        with pytest.raises(AttributeError):
            s1.lat = 1. * u.deg

    def test_getitem(self):

        s = UnitSphericalRepresentation(lon=np.arange(10) * u.deg,
                                        lat=-np.arange(10) * u.deg)

        s_slc = s[2:8:2]

        assert_allclose_quantity(s_slc.lon, [2, 4, 6] * u.deg)
        assert_allclose_quantity(s_slc.lat, [-2, -4, -6] * u.deg)

    def test_getitem_scalar(self):

        s = UnitSphericalRepresentation(lon=1 * u.deg,
                                        lat=-2 * u.deg)

        with pytest.raises(TypeError):
            s_slc = s[0]

    def test_representation_shortcuts(self):
        """Test that shortcuts in ``represent_as`` don't fail."""
        # TODO! representation transformations with differentials cannot
        # (currently) be implemented due to a mismatch between the UnitSpherical
        # expected keys (e.g. "s") and that expected in the other class
        # (here "s / m"). For more info, see PR #11467
        # We leave the test code commented out for future use.
        # diffs = UnitSphericalCosLatDifferential(4*u.mas/u.yr, 5*u.mas/u.yr,
        #                                         6*u.km/u.s)
        sph = UnitSphericalRepresentation(1*u.deg, 2*u.deg)
                                          # , differentials={'s': diffs}
        got = sph.represent_as(PhysicsSphericalRepresentation)
                               # , PhysicsSphericalDifferential)
        assert np.may_share_memory(sph.lon, got.phi)
        expected = BaseRepresentation.represent_as(
            sph, PhysicsSphericalRepresentation)  # PhysicsSphericalDifferential
        assert representation_equal_up_to_angular_type(got, expected)

        got = sph.represent_as(SphericalRepresentation)
                               # , SphericalDifferential)
        assert np.may_share_memory(sph.lon, got.lon)
        assert np.may_share_memory(sph.lat, got.lat)
        expected = BaseRepresentation.represent_as(
            sph, SphericalRepresentation) # , SphericalDifferential)
        assert representation_equal_up_to_angular_type(got, expected)

    def test_transform(self):
        """Test ``.transform()`` on rotation and general matrices."""
        # set up representation
        ds1 = UnitSphericalDifferential(d_lon=[1, 2] * u.mas / u.yr,
                                        d_lat=[3, 4] * u.mas / u.yr,)
        s1 = UnitSphericalRepresentation(lon=[1, 2] * u.deg, lat=[3, 4] * u.deg,
                                         differentials=ds1)

        # transform representation & get comparison (thru CartesianRep)
        s2 = s1.transform(matrices["rotation"])
        ds2 =  s2.differentials["s"]

        dexpected = UnitSphericalDifferential.from_cartesian(
            ds1.to_cartesian(base=s1).transform(matrices["rotation"]), base=s2)

        assert_allclose_quantity(s2.lon, s1.lon + 10 * u.deg)
        assert_allclose_quantity(s2.lat, s1.lat)
        # compare differentials. they should be unchanged (ds1).
        assert_allclose_quantity(ds2.d_lon, ds1.d_lon)
        assert_allclose_quantity(ds2.d_lat, ds1.d_lat)
        assert_allclose_quantity(ds2.d_lon, dexpected.d_lon)
        assert_allclose_quantity(ds2.d_lat, dexpected.d_lat)
        assert not hasattr(ds2, "d_distance")

        # now with a non rotation matrix
        # note that the result will be a Spherical, not UnitSpherical
        s3 = s1.transform(matrices["general"])
        ds3 = s3.differentials["s"]

        expected = (s1.represent_as(CartesianRepresentation,
                                    CartesianDifferential)
                    .transform(matrices["general"])
                    .represent_as(SphericalRepresentation,
                                  differential_class=SphericalDifferential))
        dexpected = expected.differentials["s"]

        assert_allclose_quantity(s3.lon, expected.lon)
        assert_allclose_quantity(s3.lat, expected.lat)
        assert_allclose_quantity(s3.distance, expected.distance)
        assert_allclose_quantity(ds3.d_lon, dexpected.d_lon)
        assert_allclose_quantity(ds3.d_lat, dexpected.d_lat)
        assert_allclose_quantity(ds3.d_distance, dexpected.d_distance)


class TestPhysicsSphericalRepresentation:

    def test_name(self):
        assert PhysicsSphericalRepresentation.get_name() == 'physicsspherical'
        assert PhysicsSphericalRepresentation.get_name() in REPRESENTATION_CLASSES

    def test_empty_init(self):
        with pytest.raises(TypeError) as exc:
            s = PhysicsSphericalRepresentation()

    def test_init_quantity(self):

        s3 = PhysicsSphericalRepresentation(phi=8 * u.hourangle, theta=5 * u.deg, r=10 * u.kpc)
        assert s3.phi == 8. * u.hourangle
        assert s3.theta == 5. * u.deg
        assert s3.r == 10 * u.kpc

        assert isinstance(s3.phi, Angle)
        assert isinstance(s3.theta, Angle)
        assert isinstance(s3.r, Distance)

    def test_init_phitheta(self):

        s2 = PhysicsSphericalRepresentation(Angle(8, u.hour),
                                            Angle(5, u.deg),
                                            Distance(10, u.kpc))

        assert s2.phi == 8. * u.hourangle
        assert s2.theta == 5. * u.deg
        assert s2.r == 10. * u.kpc

        assert isinstance(s2.phi, Angle)
        assert isinstance(s2.theta, Angle)
        assert isinstance(s2.r, Distance)

    def test_init_array(self):

        s1 = PhysicsSphericalRepresentation(phi=[8, 9] * u.hourangle,
                                            theta=[5, 6] * u.deg,
                                            r=[1, 2] * u.kpc)

        assert_allclose(s1.phi.degree, [120, 135])
        assert_allclose(s1.theta.degree, [5, 6])
        assert_allclose(s1.r.kpc, [1, 2])

        assert isinstance(s1.phi, Angle)
        assert isinstance(s1.theta, Angle)
        assert isinstance(s1.r, Distance)

    def test_init_array_nocopy(self):

        phi = Angle([8, 9] * u.hourangle)
        theta = Angle([5, 6] * u.deg)
        r = Distance([1, 2] * u.kpc)

        s1 = PhysicsSphericalRepresentation(phi=phi, theta=theta, r=r, copy=False)

        phi[:] = [1, 2] * u.rad
        theta[:] = [3, 4] * u.arcmin
        r[:] = [8, 9] * u.Mpc

        assert_allclose_quantity(phi, s1.phi)
        assert_allclose_quantity(theta, s1.theta)
        assert_allclose_quantity(r, s1.r)

    def test_reprobj(self):

        s1 = PhysicsSphericalRepresentation(phi=8 * u.hourangle, theta=5 * u.deg, r=10 * u.kpc)

        s2 = PhysicsSphericalRepresentation.from_representation(s1)

        assert_allclose_quantity(s2.phi, 8. * u.hourangle)
        assert_allclose_quantity(s2.theta, 5. * u.deg)
        assert_allclose_quantity(s2.r, 10 * u.kpc)

        s3 = PhysicsSphericalRepresentation(s1)

        assert representation_equal(s3, s1)

    def test_broadcasting(self):

        s1 = PhysicsSphericalRepresentation(phi=[8, 9] * u.hourangle,
                                            theta=[5, 6] * u.deg,
                                            r=10 * u.kpc)

        assert_allclose_quantity(s1.phi, [120, 135] * u.degree)
        assert_allclose_quantity(s1.theta, [5, 6] * u.degree)
        assert_allclose_quantity(s1.r, [10, 10] * u.kpc)

    def test_broadcasting_mismatch(self):

        with pytest.raises(ValueError) as exc:
            s1 = PhysicsSphericalRepresentation(phi=[8, 9, 10] * u.hourangle,
                                                theta=[5, 6] * u.deg,
                                                r=[1, 2] * u.kpc)
        assert exc.value.args[0] == "Input parameters phi, theta, and r cannot be broadcast"

    def test_readonly(self):

        s1 = PhysicsSphericalRepresentation(phi=[8, 9] * u.hourangle,
                                            theta=[5, 6] * u.deg,
                                            r=[10, 20] * u.kpc)

        with pytest.raises(AttributeError):
            s1.phi = 1. * u.deg

        with pytest.raises(AttributeError):
            s1.theta = 1. * u.deg

        with pytest.raises(AttributeError):
            s1.r = 1. * u.kpc

    def test_getitem(self):

        s = PhysicsSphericalRepresentation(phi=np.arange(10) * u.deg,
                                           theta=np.arange(5, 15) * u.deg,
                                           r=1 * u.kpc)

        s_slc = s[2:8:2]

        assert_allclose_quantity(s_slc.phi, [2, 4, 6] * u.deg)
        assert_allclose_quantity(s_slc.theta, [7, 9, 11] * u.deg)
        assert_allclose_quantity(s_slc.r, [1, 1, 1] * u.kpc)

    def test_getitem_scalar(self):

        s = PhysicsSphericalRepresentation(phi=1 * u.deg,
                                           theta=2 * u.deg,
                                           r=3 * u.kpc)

        with pytest.raises(TypeError):
            s_slc = s[0]

    def test_representation_shortcuts(self):
        """Test that shortcuts in ``represent_as`` don't fail."""
        difs = PhysicsSphericalDifferential(4*u.mas/u.yr,5*u.mas/u.yr,6*u.km/u.s)
        sph = PhysicsSphericalRepresentation(1*u.deg, 2*u.deg, 3*u.kpc,
                                             differentials={'s': difs})

        got = sph.represent_as(SphericalRepresentation,
                               SphericalDifferential)
        assert np.may_share_memory(sph.phi, got.lon)
        assert np.may_share_memory(sph.r, got.distance)
        expected = BaseRepresentation.represent_as(
            sph, SphericalRepresentation, SphericalDifferential)
        assert representation_equal_up_to_angular_type(got, expected)

        got = sph.represent_as(UnitSphericalRepresentation,
                               UnitSphericalDifferential)
        assert np.may_share_memory(sph.phi, got.lon)
        expected = BaseRepresentation.represent_as(
            sph, UnitSphericalRepresentation, UnitSphericalDifferential)
        assert representation_equal_up_to_angular_type(got, expected)

    def test_initialize_with_nan(self):
        # Regression test for gh-11558: initialization used to fail.
        psr = PhysicsSphericalRepresentation([1., np.nan]*u.deg, [np.nan, 2.]*u.deg,
                                             [3., np.nan]*u.m)
        assert_array_equal(np.isnan(psr.phi), [False, True])
        assert_array_equal(np.isnan(psr.theta), [True, False])
        assert_array_equal(np.isnan(psr.r), [False, True])

    def test_transform(self):
        """Test ``.transform()`` on rotation and general transform matrices."""
        # set up representation
        ds1 = PhysicsSphericalDifferential(
            d_phi=[1, 2] * u.mas / u.yr, d_theta=[3, 4] * u.mas / u.yr,
            d_r=[-5, 6] * u.km / u.s)
        s1 = PhysicsSphericalRepresentation(
            phi=[1, 2] * u.deg, theta=[3, 4] * u.deg, r=[5, 6] * u.kpc,
            differentials=ds1)

        # transform representation & get comparison (thru CartesianRep)
        s2 = s1.transform(matrices["rotation"])
        ds2 = s2.differentials["s"]

        dexpected = PhysicsSphericalDifferential.from_cartesian(
            ds1.to_cartesian(base=s1).transform(matrices["rotation"]), base=s2)

        assert_allclose_quantity(s2.phi, s1.phi + 10 * u.deg)
        assert_allclose_quantity(s2.theta, s1.theta)
        assert_allclose_quantity(s2.r, s1.r)
        # compare differentials. should be unchanged (ds1).
        assert_allclose_quantity(ds2.d_phi, ds1.d_phi)
        assert_allclose_quantity(ds2.d_theta, ds1.d_theta)
        assert_allclose_quantity(ds2.d_r, ds1.d_r)
        assert_allclose_quantity(ds2.d_phi, dexpected.d_phi)
        assert_allclose_quantity(ds2.d_theta, dexpected.d_theta)
        assert_allclose_quantity(ds2.d_r, dexpected.d_r)

        # now with a non rotation matrix
        # transform representation & get comparison (thru CartesianRep)
        s3 = s1.transform(matrices["general"])
        ds3 = s3.differentials["s"]

        expected = (s1.represent_as(CartesianRepresentation,
                                    CartesianDifferential)
                    .transform(matrices["general"])
                    .represent_as(PhysicsSphericalRepresentation,
                                  PhysicsSphericalDifferential))
        dexpected = expected.differentials["s"]

        assert_allclose_quantity(s3.phi, expected.phi)
        assert_allclose_quantity(s3.theta, expected.theta)
        assert_allclose_quantity(s3.r, expected.r)
        assert_allclose_quantity(ds3.d_phi, dexpected.d_phi)
        assert_allclose_quantity(ds3.d_theta, dexpected.d_theta)
        assert_allclose_quantity(ds3.d_r, dexpected.d_r)

    def test_transform_with_NaN(self):
        # all over again, but with a NaN in the distance

        ds1 = PhysicsSphericalDifferential(
            d_phi=[1, 2] * u.mas / u.yr, d_theta=[3, 4] * u.mas / u.yr,
            d_r=[-5, 6] * u.km / u.s)
        s1 = PhysicsSphericalRepresentation(
            phi=[1, 2] * u.deg, theta=[3, 4] * u.deg, r=[5, np.nan] * u.kpc,
            differentials=ds1)

        # transform representation & get comparison (thru CartesianRep)
        s2 = s1.transform(matrices["rotation"])
        ds2 =  s2.differentials["s"]

        dexpected = PhysicsSphericalDifferential.from_cartesian(
            ds1.to_cartesian(base=s1).transform(matrices["rotation"]), base=s2)

        assert_allclose_quantity(s2.phi, s1.phi + 10 * u.deg)
        assert_allclose_quantity(s2.theta, s1.theta)
        assert_allclose_quantity(s2.r, s1.r)
        assert_allclose_quantity(ds2.d_phi, dexpected.d_phi)
        assert_allclose_quantity(ds2.d_theta, dexpected.d_theta)
        assert_allclose_quantity(ds2.d_r, dexpected.d_r)

        # now with a non rotation matrix
        s3 = s1.transform(matrices["general"])
        ds3 = s3.differentials["s"]

        thruC = (s1.represent_as(CartesianRepresentation,
                                 CartesianDifferential)
                    .transform(matrices["general"])
                    .represent_as(PhysicsSphericalRepresentation,
                                  PhysicsSphericalDifferential))
        dthruC = thruC.differentials["s"]

        # s3 should not propagate Nan.
        assert_array_equal(np.isnan(s3.phi), (False, False))
        assert_array_equal(np.isnan(s3.theta), (False, False))
        assert_array_equal(np.isnan(s3.r), (False, True))
        # ds3 does b/c currently aren't any shortcuts on the transform
        assert_array_equal(np.isnan(ds3.d_phi), (False, True))
        assert_array_equal(np.isnan(ds3.d_theta), (False, True))
        assert_array_equal(np.isnan(ds3.d_r), (False, True))

        # through Cartesian does
        assert_array_equal(np.isnan(thruC.phi), (False, True))
        assert_array_equal(np.isnan(thruC.theta), (False, True))
        assert_array_equal(np.isnan(thruC.r), (False, True))
        # so only test on the first value
        assert_allclose_quantity(s3.phi[0], thruC.phi[0])
        assert_allclose_quantity(s3.theta[0], thruC.theta[0])
        assert_allclose_quantity(ds3.d_phi[0], dthruC.d_phi[0])
        assert_allclose_quantity(ds3.d_theta[0], dthruC.d_theta[0])


class TestCartesianRepresentation:

    def test_name(self):
        assert CartesianRepresentation.get_name() == 'cartesian'
        assert CartesianRepresentation.get_name() in REPRESENTATION_CLASSES

    def test_empty_init(self):
        with pytest.raises(TypeError) as exc:
            s = CartesianRepresentation()

    def test_init_quantity(self):

        s1 = CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc)

        assert s1.x.unit is u.kpc
        assert s1.y.unit is u.kpc
        assert s1.z.unit is u.kpc

        assert_allclose(s1.x.value, 1)
        assert_allclose(s1.y.value, 2)
        assert_allclose(s1.z.value, 3)

    def test_init_singleunit(self):

        s1 = CartesianRepresentation(x=1, y=2, z=3, unit=u.kpc)

        assert s1.x.unit is u.kpc
        assert s1.y.unit is u.kpc
        assert s1.z.unit is u.kpc

        assert_allclose(s1.x.value, 1)
        assert_allclose(s1.y.value, 2)
        assert_allclose(s1.z.value, 3)

    def test_init_array(self):

        s1 = CartesianRepresentation(x=[1, 2, 3] * u.pc,
                                     y=[2, 3, 4] * u.Mpc,
                                     z=[3, 4, 5] * u.kpc)

        assert s1.x.unit is u.pc
        assert s1.y.unit is u.Mpc
        assert s1.z.unit is u.kpc

        assert_allclose(s1.x.value, [1, 2, 3])
        assert_allclose(s1.y.value, [2, 3, 4])
        assert_allclose(s1.z.value, [3, 4, 5])

    def test_init_one_array(self):

        s1 = CartesianRepresentation(x=[1, 2, 3] * u.pc)

        assert s1.x.unit is u.pc
        assert s1.y.unit is u.pc
        assert s1.z.unit is u.pc

        assert_allclose(s1.x.value, 1)
        assert_allclose(s1.y.value, 2)
        assert_allclose(s1.z.value, 3)

        r = np.arange(27.).reshape(3, 3, 3) * u.kpc
        s2 = CartesianRepresentation(r, xyz_axis=0)
        assert s2.shape == (3, 3)
        assert s2.x.unit == u.kpc
        assert np.all(s2.x == r[0])
        assert np.all(s2.xyz == r)
        assert np.all(s2.get_xyz(xyz_axis=0) == r)
        s3 = CartesianRepresentation(r, xyz_axis=1)
        assert s3.shape == (3, 3)
        assert np.all(s3.x == r[:, 0])
        assert np.all(s3.y == r[:, 1])
        assert np.all(s3.z == r[:, 2])
        assert np.all(s3.get_xyz(xyz_axis=1) == r)
        s4 = CartesianRepresentation(r, xyz_axis=2)
        assert s4.shape == (3, 3)
        assert np.all(s4.x == r[:, :, 0])
        assert np.all(s4.get_xyz(xyz_axis=2) == r)
        s5 = CartesianRepresentation(r, unit=u.pc)
        assert s5.x.unit == u.pc
        assert np.all(s5.xyz == r)
        s6 = CartesianRepresentation(r.value, unit=u.pc, xyz_axis=2)
        assert s6.x.unit == u.pc
        assert np.all(s6.get_xyz(xyz_axis=2).value == r.value)

    def test_init_one_array_size_fail(self):
        with pytest.raises(ValueError) as exc:
            CartesianRepresentation(x=[1, 2, 3, 4] * u.pc)
        assert exc.value.args[0].startswith("too many values to unpack")

    def test_init_xyz_but_more_than_one_array_fail(self):
        with pytest.raises(ValueError) as exc:
            CartesianRepresentation(x=[1, 2, 3] * u.pc, y=[2, 3, 4] * u.pc,
                                    z=[3, 4, 5] * u.pc, xyz_axis=0)
        assert 'xyz_axis should only be set' in str(exc.value)

    def test_init_one_array_yz_fail(self):
        with pytest.raises(ValueError) as exc:
            CartesianRepresentation(x=[1, 2, 3, 4] * u.pc, y=[1, 2] * u.pc)
        assert exc.value.args[0] == ("x, y, and z are required to instantiate "
                                     "CartesianRepresentation")

    def test_init_array_nocopy(self):

        x = [8, 9, 10] * u.pc
        y = [5, 6, 7] * u.Mpc
        z = [2, 3, 4] * u.kpc

        s1 = CartesianRepresentation(x=x, y=y, z=z, copy=False)

        x[:] = [1, 2, 3] * u.kpc
        y[:] = [9, 9, 8] * u.kpc
        z[:] = [1, 2, 1] * u.kpc

        assert_allclose_quantity(x, s1.x)
        assert_allclose_quantity(y, s1.y)
        assert_allclose_quantity(z, s1.z)

    def test_xyz_is_view_if_possible(self):
        xyz = np.arange(1., 10.).reshape(3, 3)
        s1 = CartesianRepresentation(xyz, unit=u.kpc, copy=False)
        s1_xyz = s1.xyz
        assert s1_xyz.value[0, 0] == 1.
        xyz[0, 0] = 0.
        assert s1.x[0] == 0.
        assert s1_xyz.value[0, 0] == 0.
        # Not possible: we don't check that tuples are from the same array
        xyz = np.arange(1., 10.).reshape(3, 3)
        s2 = CartesianRepresentation(*xyz, unit=u.kpc, copy=False)
        s2_xyz = s2.xyz
        assert s2_xyz.value[0, 0] == 1.
        xyz[0, 0] = 0.
        assert s2.x[0] == 0.
        assert s2_xyz.value[0, 0] == 1.

    def test_reprobj(self):

        s1 = CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc)

        s2 = CartesianRepresentation.from_representation(s1)

        assert s2.x == 1 * u.kpc
        assert s2.y == 2 * u.kpc
        assert s2.z == 3 * u.kpc

        s3 = CartesianRepresentation(s1)

        assert representation_equal(s3, s1)

    def test_broadcasting(self):

        s1 = CartesianRepresentation(x=[1, 2] * u.kpc, y=[3, 4] * u.kpc, z=5 * u.kpc)

        assert s1.x.unit == u.kpc
        assert s1.y.unit == u.kpc
        assert s1.z.unit == u.kpc

        assert_allclose(s1.x.value, [1, 2])
        assert_allclose(s1.y.value, [3, 4])
        assert_allclose(s1.z.value, [5, 5])

    def test_broadcasting_mismatch(self):

        with pytest.raises(ValueError) as exc:
            s1 = CartesianRepresentation(x=[1, 2] * u.kpc, y=[3, 4] * u.kpc, z=[5, 6, 7] * u.kpc)
        assert exc.value.args[0] == "Input parameters x, y, and z cannot be broadcast"

    def test_readonly(self):

        s1 = CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc)

        with pytest.raises(AttributeError):
            s1.x = 1. * u.kpc

        with pytest.raises(AttributeError):
            s1.y = 1. * u.kpc

        with pytest.raises(AttributeError):
            s1.z = 1. * u.kpc

    def test_xyz(self):

        s1 = CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc)

        assert isinstance(s1.xyz, u.Quantity)
        assert s1.xyz.unit is u.kpc

        assert_allclose(s1.xyz.value, [1, 2, 3])

    def test_unit_mismatch(self):

        q_len = u.Quantity([1], u.km)
        q_nonlen = u.Quantity([1], u.kg)

        with pytest.raises(u.UnitsError) as exc:
            s1 = CartesianRepresentation(x=q_nonlen, y=q_len, z=q_len)
        assert exc.value.args[0] == "x, y, and z should have matching physical types"

        with pytest.raises(u.UnitsError) as exc:
            s1 = CartesianRepresentation(x=q_len, y=q_nonlen, z=q_len)
        assert exc.value.args[0] == "x, y, and z should have matching physical types"

        with pytest.raises(u.UnitsError) as exc:
            s1 = CartesianRepresentation(x=q_len, y=q_len, z=q_nonlen)
        assert exc.value.args[0] == "x, y, and z should have matching physical types"

    def test_unit_non_length(self):

        s1 = CartesianRepresentation(x=1 * u.kg, y=2 * u.kg, z=3 * u.kg)

        s2 = CartesianRepresentation(x=1 * u.km / u.s, y=2 * u.km / u.s, z=3 * u.km / u.s)

        banana = u.def_unit('banana')
        s3 = CartesianRepresentation(x=1 * banana, y=2 * banana, z=3 * banana)

    def test_getitem(self):

        s = CartesianRepresentation(x=np.arange(10) * u.m,
                                    y=-np.arange(10) * u.m,
                                    z=3 * u.km)

        s_slc = s[2:8:2]

        assert_allclose_quantity(s_slc.x, [2, 4, 6] * u.m)
        assert_allclose_quantity(s_slc.y, [-2, -4, -6] * u.m)
        assert_allclose_quantity(s_slc.z, [3, 3, 3] * u.km)

    def test_getitem_scalar(self):

        s = CartesianRepresentation(x=1 * u.m,
                                    y=-2 * u.m,
                                    z=3 * u.km)

        with pytest.raises(TypeError):
            s_slc = s[0]

    def test_transform(self):

        ds1 = CartesianDifferential(d_x=[1, 2] * u.km / u.s,
                                    d_y=[3, 4] * u.km / u.s,
                                    d_z=[5, 6] * u.km / u.s)
        s1 = CartesianRepresentation(x=[1, 2] * u.kpc, y=[3, 4] * u.kpc,
                                     z=[5, 6] * u.kpc, differentials=ds1)

        # transform representation & get comparison (thru CartesianRep)
        s2 = s1.transform(matrices["general"])
        ds2 =  s2.differentials["s"]

        dexpected = CartesianDifferential.from_cartesian(
            ds1.to_cartesian(base=s1).transform(matrices["general"]), base=s2)

        assert_allclose_quantity(ds2.d_x, dexpected.d_x)
        assert_allclose_quantity(ds2.d_y, dexpected.d_y)
        assert_allclose_quantity(ds2.d_z, dexpected.d_z)

        # also explicitly calculate, since we can
        assert_allclose(s2.x.value, [1 * 1 + 2 * 3 + 3 * 5, 1 * 2 + 2 * 4 + 3 * 6])
        assert_allclose(s2.y.value, [4 * 1 + 5 * 3 + 6 * 5, 4 * 2 + 5 * 4 + 6 * 6])
        assert_allclose(s2.z.value, [7 * 1 + 8 * 3 + 9 * 5, 7 * 2 + 8 * 4 + 9 * 6])
        assert_allclose(ds2.d_x.value, [1 * 1 + 2 * 3 + 3 * 5, 1 * 2 + 2 * 4 + 3 * 6])
        assert_allclose(ds2.d_y.value, [4 * 1 + 5 * 3 + 6 * 5, 4 * 2 + 5 * 4 + 6 * 6])
        assert_allclose(ds2.d_z.value, [7 * 1 + 8 * 3 + 9 * 5, 7 * 2 + 8 * 4 + 9 * 6])

        assert s2.x.unit is u.kpc
        assert s2.y.unit is u.kpc
        assert s2.z.unit is u.kpc
        assert ds2.d_x.unit == u.km / u.s
        assert ds2.d_y.unit == u.km / u.s
        assert ds2.d_z.unit == u.km / u.s


class TestCylindricalRepresentation:

    def test_name(self):
        assert CylindricalRepresentation.get_name() == 'cylindrical'
        assert CylindricalRepresentation.get_name() in REPRESENTATION_CLASSES

    def test_empty_init(self):
        with pytest.raises(TypeError) as exc:
            s = CylindricalRepresentation()

    def test_init_quantity(self):

        s1 = CylindricalRepresentation(rho=1 * u.kpc, phi=2 * u.deg, z=3 * u.kpc)

        assert s1.rho.unit is u.kpc
        assert s1.phi.unit is u.deg
        assert s1.z.unit is u.kpc

        assert_allclose(s1.rho.value, 1)
        assert_allclose(s1.phi.value, 2)
        assert_allclose(s1.z.value, 3)

    def test_init_array(self):

        s1 = CylindricalRepresentation(rho=[1, 2, 3] * u.pc,
                                       phi=[2, 3, 4] * u.deg,
                                       z=[3, 4, 5] * u.kpc)

        assert s1.rho.unit is u.pc
        assert s1.phi.unit is u.deg
        assert s1.z.unit is u.kpc

        assert_allclose(s1.rho.value, [1, 2, 3])
        assert_allclose(s1.phi.value, [2, 3, 4])
        assert_allclose(s1.z.value, [3, 4, 5])

    def test_init_array_nocopy(self):

        rho = [8, 9, 10] * u.pc
        phi = [5, 6, 7] * u.deg
        z = [2, 3, 4] * u.kpc

        s1 = CylindricalRepresentation(rho=rho, phi=phi, z=z, copy=False)

        rho[:] = [9, 2, 3] * u.kpc
        phi[:] = [1, 2, 3] * u.arcmin
        z[:] = [-2, 3, 8] * u.kpc

        assert_allclose_quantity(rho, s1.rho)
        assert_allclose_quantity(phi, s1.phi)
        assert_allclose_quantity(z, s1.z)

    def test_reprobj(self):

        s1 = CylindricalRepresentation(rho=1 * u.kpc, phi=2 * u.deg, z=3 * u.kpc)

        s2 = CylindricalRepresentation.from_representation(s1)

        assert s2.rho == 1 * u.kpc
        assert s2.phi == 2 * u.deg
        assert s2.z == 3 * u.kpc

        s3 = CylindricalRepresentation(s1)

        assert representation_equal(s3, s1)

    def test_broadcasting(self):

        s1 = CylindricalRepresentation(rho=[1, 2] * u.kpc, phi=[3, 4] * u.deg, z=5 * u.kpc)

        assert s1.rho.unit == u.kpc
        assert s1.phi.unit == u.deg
        assert s1.z.unit == u.kpc

        assert_allclose(s1.rho.value, [1, 2])
        assert_allclose(s1.phi.value, [3, 4])
        assert_allclose(s1.z.value, [5, 5])

    def test_broadcasting_mismatch(self):

        with pytest.raises(ValueError) as exc:
            s1 = CylindricalRepresentation(rho=[1, 2] * u.kpc, phi=[3, 4] * u.deg, z=[5, 6, 7] * u.kpc)
        assert exc.value.args[0] == "Input parameters rho, phi, and z cannot be broadcast"

    def test_readonly(self):

        s1 = CylindricalRepresentation(rho=1 * u.kpc,
                                       phi=20 * u.deg,
                                       z=3 * u.kpc)

        with pytest.raises(AttributeError):
            s1.rho = 1. * u.kpc

        with pytest.raises(AttributeError):
            s1.phi = 20 * u.deg

        with pytest.raises(AttributeError):
            s1.z = 1. * u.kpc

    def unit_mismatch(self):

        q_len = u.Quantity([1], u.kpc)
        q_nonlen = u.Quantity([1], u.kg)

        with pytest.raises(u.UnitsError) as exc:
            s1 = CylindricalRepresentation(rho=q_nonlen, phi=10 * u.deg, z=q_len)
        assert exc.value.args[0] == "rho and z should have matching physical types"

        with pytest.raises(u.UnitsError) as exc:
            s1 = CylindricalRepresentation(rho=q_len, phi=10 * u.deg, z=q_nonlen)
        assert exc.value.args[0] == "rho and z should have matching physical types"

    def test_getitem(self):

        s = CylindricalRepresentation(rho=np.arange(10) * u.pc,
                                      phi=-np.arange(10) * u.deg,
                                      z=1 * u.kpc)

        s_slc = s[2:8:2]

        assert_allclose_quantity(s_slc.rho, [2, 4, 6] * u.pc)
        assert_allclose_quantity(s_slc.phi, [-2, -4, -6] * u.deg)
        assert_allclose_quantity(s_slc.z, [1, 1, 1] * u.kpc)

    def test_getitem_scalar(self):

        s = CylindricalRepresentation(rho=1 * u.pc,
                                      phi=-2 * u.deg,
                                      z=3 * u.kpc)

        with pytest.raises(TypeError):
            s_slc = s[0]

    def test_transform(self):

        s1 = CylindricalRepresentation(phi=[1, 2] * u.deg, z=[3, 4] * u.pc,
                                       rho=[5, 6] * u.kpc)

        s2 = s1.transform(matrices["rotation"])

        assert_allclose_quantity(s2.phi, s1.phi + 10 * u.deg)
        assert_allclose_quantity(s2.z, s1.z)
        assert_allclose_quantity(s2.rho, s1.rho)

        assert s2.phi.unit is u.rad
        assert s2.z.unit is u.kpc
        assert s2.rho.unit is u.kpc

        # now with a non rotation matrix
        s3 = s1.transform(matrices["general"])
        expected = (s1.to_cartesian().transform(matrices["general"])
                    ).represent_as(CylindricalRepresentation)

        assert_allclose_quantity(s3.phi, expected.phi)
        assert_allclose_quantity(s3.z, expected.z)
        assert_allclose_quantity(s3.rho, expected.rho)


class TestUnitSphericalCosLatDifferential:

    @pytest.mark.parametrize("matrix", list(matrices.values()))
    def test_transform(self, matrix):
        """Test ``.transform()`` on rotation and general matrices."""
        # set up representation
        ds1 = UnitSphericalCosLatDifferential(d_lon_coslat=[1, 2] * u.mas / u.yr,
                                              d_lat=[3, 4] * u.mas / u.yr,)
        s1 = UnitSphericalRepresentation(lon=[1, 2] * u.deg, lat=[3, 4] * u.deg)

        # transform representation & get comparison (thru CartesianRep)
        s2 = s1.transform(matrix)
        ds2 = ds1.transform(matrix, s1, s2)

        dexpected = UnitSphericalCosLatDifferential.from_cartesian(
            ds1.to_cartesian(base=s1).transform(matrix), base=s2)

        assert_allclose_quantity(ds2.d_lon_coslat, dexpected.d_lon_coslat)
        assert_allclose_quantity(ds2.d_lat, dexpected.d_lat)


def test_cartesian_spherical_roundtrip():

    s1 = CartesianRepresentation(x=[1, 2000.] * u.kpc,
                                 y=[3000., 4.] * u.pc,
                                 z=[5., 6000.] * u.pc)

    s2 = SphericalRepresentation.from_representation(s1)

    s3 = CartesianRepresentation.from_representation(s2)

    s4 = SphericalRepresentation.from_representation(s3)

    assert_allclose_quantity(s1.x, s3.x)
    assert_allclose_quantity(s1.y, s3.y)
    assert_allclose_quantity(s1.z, s3.z)

    assert_allclose_quantity(s2.lon, s4.lon)
    assert_allclose_quantity(s2.lat, s4.lat)
    assert_allclose_quantity(s2.distance, s4.distance)


def test_cartesian_setting_with_other():

    s1 = CartesianRepresentation(x=[1, 2000.] * u.kpc,
                                 y=[3000., 4.] * u.pc,
                                 z=[5., 6000.] * u.pc)
    s1[0] = SphericalRepresentation(0.*u.deg, 0.*u.deg, 1*u.kpc)
    assert_allclose_quantity(s1.x, [1., 2000.] * u.kpc)
    assert_allclose_quantity(s1.y, [0., 4.] * u.pc)
    assert_allclose_quantity(s1.z, [0., 6000.] * u.pc)

    with pytest.raises(ValueError, match='loss of information'):
        s1[1] = UnitSphericalRepresentation(0.*u.deg, 10.*u.deg)


def test_cartesian_physics_spherical_roundtrip():

    s1 = CartesianRepresentation(x=[1, 2000.] * u.kpc,
                                 y=[3000., 4.] * u.pc,
                                 z=[5., 6000.] * u.pc)

    s2 = PhysicsSphericalRepresentation.from_representation(s1)

    s3 = CartesianRepresentation.from_representation(s2)

    s4 = PhysicsSphericalRepresentation.from_representation(s3)

    assert_allclose_quantity(s1.x, s3.x)
    assert_allclose_quantity(s1.y, s3.y)
    assert_allclose_quantity(s1.z, s3.z)

    assert_allclose_quantity(s2.phi, s4.phi)
    assert_allclose_quantity(s2.theta, s4.theta)
    assert_allclose_quantity(s2.r, s4.r)


def test_spherical_physics_spherical_roundtrip():

    s1 = SphericalRepresentation(lon=3 * u.deg, lat=4 * u.deg, distance=3 * u.kpc)

    s2 = PhysicsSphericalRepresentation.from_representation(s1)

    s3 = SphericalRepresentation.from_representation(s2)

    s4 = PhysicsSphericalRepresentation.from_representation(s3)

    assert_allclose_quantity(s1.lon, s3.lon)
    assert_allclose_quantity(s1.lat, s3.lat)
    assert_allclose_quantity(s1.distance, s3.distance)

    assert_allclose_quantity(s2.phi, s4.phi)
    assert_allclose_quantity(s2.theta, s4.theta)
    assert_allclose_quantity(s2.r, s4.r)

    assert_allclose_quantity(s1.lon, s4.phi)
    assert_allclose_quantity(s1.lat, 90. * u.deg - s4.theta)
    assert_allclose_quantity(s1.distance, s4.r)


def test_cartesian_cylindrical_roundtrip():

    s1 = CartesianRepresentation(x=np.array([1., 2000.]) * u.kpc,
                                 y=np.array([3000., 4.]) * u.pc,
                                 z=np.array([5., 600.]) * u.cm)

    s2 = CylindricalRepresentation.from_representation(s1)

    s3 = CartesianRepresentation.from_representation(s2)

    s4 = CylindricalRepresentation.from_representation(s3)

    assert_allclose_quantity(s1.x, s3.x)
    assert_allclose_quantity(s1.y, s3.y)
    assert_allclose_quantity(s1.z, s3.z)

    assert_allclose_quantity(s2.rho, s4.rho)
    assert_allclose_quantity(s2.phi, s4.phi)
    assert_allclose_quantity(s2.z, s4.z)


def test_unit_spherical_roundtrip():

    s1 = UnitSphericalRepresentation(lon=[10., 30.] * u.deg,
                                     lat=[5., 6.] * u.arcmin)

    s2 = CartesianRepresentation.from_representation(s1)

    s3 = SphericalRepresentation.from_representation(s2)

    s4 = UnitSphericalRepresentation.from_representation(s3)

    assert_allclose_quantity(s1.lon, s4.lon)
    assert_allclose_quantity(s1.lat, s4.lat)


def test_no_unnecessary_copies():

    s1 = UnitSphericalRepresentation(lon=[10., 30.] * u.deg,
                                     lat=[5., 6.] * u.arcmin)
    s2 = s1.represent_as(UnitSphericalRepresentation)
    assert s2 is s1
    assert np.may_share_memory(s1.lon, s2.lon)
    assert np.may_share_memory(s1.lat, s2.lat)
    s3 = s1.represent_as(SphericalRepresentation)
    assert np.may_share_memory(s1.lon, s3.lon)
    assert np.may_share_memory(s1.lat, s3.lat)
    s4 = s1.represent_as(CartesianRepresentation)
    s5 = s4.represent_as(CylindricalRepresentation)
    assert np.may_share_memory(s5.z, s4.z)


def test_representation_repr():
    r1 = SphericalRepresentation(lon=1 * u.deg, lat=2.5 * u.deg, distance=1 * u.kpc)
    assert repr(r1) == ('<SphericalRepresentation (lon, lat, distance) in (deg, deg, kpc)\n'
                        '    (1., 2.5, 1.)>')

    r2 = CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc)
    assert repr(r2) == ('<CartesianRepresentation (x, y, z) in kpc\n'
                        '    (1., 2., 3.)>')

    r3 = CartesianRepresentation(x=[1, 2, 3] * u.kpc, y=4 * u.kpc, z=[9, 10, 11] * u.kpc)
    assert repr(r3) == ('<CartesianRepresentation (x, y, z) in kpc\n'
                        '    [(1., 4.,  9.), (2., 4., 10.), (3., 4., 11.)]>')


def test_representation_repr_multi_d():
    """Regression test for #5889."""
    cr = CartesianRepresentation(np.arange(27).reshape(3, 3, 3), unit='m')
    assert repr(cr) == (
        '<CartesianRepresentation (x, y, z) in m\n'
        '    [[(0.,  9., 18.), (1., 10., 19.), (2., 11., 20.)],\n'
        '     [(3., 12., 21.), (4., 13., 22.), (5., 14., 23.)],\n'
        '     [(6., 15., 24.), (7., 16., 25.), (8., 17., 26.)]]>')
    # This was broken before.
    assert repr(cr.T) == (
        '<CartesianRepresentation (x, y, z) in m\n'
        '    [[(0.,  9., 18.), (3., 12., 21.), (6., 15., 24.)],\n'
        '     [(1., 10., 19.), (4., 13., 22.), (7., 16., 25.)],\n'
        '     [(2., 11., 20.), (5., 14., 23.), (8., 17., 26.)]]>')


def test_representation_str():
    r1 = SphericalRepresentation(lon=1 * u.deg, lat=2.5 * u.deg, distance=1 * u.kpc)
    assert str(r1) == '(1., 2.5, 1.) (deg, deg, kpc)'

    r2 = CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc)
    assert str(r2) == '(1., 2., 3.) kpc'

    r3 = CartesianRepresentation(x=[1, 2, 3] * u.kpc, y=4 * u.kpc, z=[9, 10, 11] * u.kpc)
    assert str(r3) == '[(1., 4.,  9.), (2., 4., 10.), (3., 4., 11.)] kpc'


def test_representation_str_multi_d():
    """Regression test for #5889."""
    cr = CartesianRepresentation(np.arange(27).reshape(3, 3, 3), unit='m')
    assert str(cr) == (
        '[[(0.,  9., 18.), (1., 10., 19.), (2., 11., 20.)],\n'
        ' [(3., 12., 21.), (4., 13., 22.), (5., 14., 23.)],\n'
        ' [(6., 15., 24.), (7., 16., 25.), (8., 17., 26.)]] m')
    # This was broken before.
    assert str(cr.T) == (
        '[[(0.,  9., 18.), (3., 12., 21.), (6., 15., 24.)],\n'
        ' [(1., 10., 19.), (4., 13., 22.), (7., 16., 25.)],\n'
        ' [(2., 11., 20.), (5., 14., 23.), (8., 17., 26.)]] m')


def test_subclass_representation():
    from astropy.coordinates.builtin_frames import ICRS

    class Longitude180(Longitude):
        def __new__(cls, angle, unit=None, wrap_angle=180 * u.deg, **kwargs):
            self = super().__new__(cls, angle, unit=unit, wrap_angle=wrap_angle,
                                   **kwargs)
            return self

    class SphericalWrap180Representation(SphericalRepresentation):
        attr_classes = {'lon': Longitude180,
                        'lat': Latitude,
                        'distance': u.Quantity}

    class ICRSWrap180(ICRS):
        frame_specific_representation_info = ICRS._frame_specific_representation_info.copy()
        frame_specific_representation_info[SphericalWrap180Representation] = \
            frame_specific_representation_info[SphericalRepresentation]
        default_representation = SphericalWrap180Representation

    c = ICRSWrap180(ra=-1 * u.deg, dec=-2 * u.deg, distance=1 * u.m)
    assert c.ra.value == -1
    assert c.ra.unit is u.deg
    assert c.dec.value == -2
    assert c.dec.unit is u.deg


def test_minimal_subclass():
    # Basically to check what we document works;
    # see doc/coordinates/representations.rst
    class LogDRepresentation(BaseRepresentation):
        attr_classes = {'lon': Longitude,
                        'lat': Latitude,
                        'logd': u.Dex}

        def to_cartesian(self):
            d = self.logd.physical
            x = d * np.cos(self.lat) * np.cos(self.lon)
            y = d * np.cos(self.lat) * np.sin(self.lon)
            z = d * np.sin(self.lat)
            return CartesianRepresentation(x=x, y=y, z=z, copy=False)

        @classmethod
        def from_cartesian(cls, cart):
            s = np.hypot(cart.x, cart.y)
            r = np.hypot(s, cart.z)
            lon = np.arctan2(cart.y, cart.x)
            lat = np.arctan2(cart.z, s)
            return cls(lon=lon, lat=lat, logd=u.Dex(r), copy=False)

    ld1 = LogDRepresentation(90.*u.deg, 0.*u.deg, 1.*u.dex(u.kpc))
    ld2 = LogDRepresentation(lon=90.*u.deg, lat=0.*u.deg, logd=1.*u.dex(u.kpc))
    assert np.all(ld1.lon == ld2.lon)
    assert np.all(ld1.lat == ld2.lat)
    assert np.all(ld1.logd == ld2.logd)
    c = ld1.to_cartesian()
    assert_allclose_quantity(c.xyz, [0., 10., 0.] * u.kpc, atol=1.*u.npc)
    ld3 = LogDRepresentation.from_cartesian(c)
    assert np.all(ld3.lon == ld2.lon)
    assert np.all(ld3.lat == ld2.lat)
    assert np.all(ld3.logd == ld2.logd)
    s = ld1.represent_as(SphericalRepresentation)
    assert_allclose_quantity(s.lon, ld1.lon)
    assert_allclose_quantity(s.distance, 10.*u.kpc)
    assert_allclose_quantity(s.lat, ld1.lat)

    with pytest.raises(TypeError):
        LogDRepresentation(0.*u.deg, 1.*u.deg)
    with pytest.raises(TypeError):
        LogDRepresentation(0.*u.deg, 1.*u.deg, 1.*u.dex(u.kpc), lon=1.*u.deg)
    with pytest.raises(TypeError):
        LogDRepresentation(0.*u.deg, 1.*u.deg, 1.*u.dex(u.kpc), True, False)
    with pytest.raises(TypeError):
        LogDRepresentation(0.*u.deg, 1.*u.deg, 1.*u.dex(u.kpc), foo='bar')

    # if we define it a second time, even the qualnames are the same,
    # so we raise
    with pytest.raises(ValueError):
        class LogDRepresentation(BaseRepresentation):
            attr_classes = {'lon': Longitude,
                            'lat': Latitude,
                            'logr': u.Dex}


def test_duplicate_warning():
    from astropy.coordinates.representation import DUPLICATE_REPRESENTATIONS
    from astropy.coordinates.representation import REPRESENTATION_CLASSES

    with pytest.warns(DuplicateRepresentationWarning):
        class UnitSphericalRepresentation(BaseRepresentation):
            attr_classes = {'lon': Longitude,
                            'lat': Latitude}

    assert 'unitspherical' in DUPLICATE_REPRESENTATIONS
    assert 'unitspherical' not in REPRESENTATION_CLASSES
    assert 'astropy.coordinates.representation.UnitSphericalRepresentation' in REPRESENTATION_CLASSES
    assert __name__ + '.test_duplicate_warning.<locals>.UnitSphericalRepresentation' in REPRESENTATION_CLASSES


class TestCartesianRepresentationWithDifferential:

    def test_init_differential(self):

        diff = CartesianDifferential(d_x=1 * u.km/u.s,
                                     d_y=2 * u.km/u.s,
                                     d_z=3 * u.km/u.s)

        # Check that a single differential gets turned into a 1-item dict.
        s1 = CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc,
                                     differentials=diff)

        assert s1.x.unit is u.kpc
        assert s1.y.unit is u.kpc
        assert s1.z.unit is u.kpc
        assert len(s1.differentials) == 1
        assert s1.differentials['s'] is diff

        # can also pass in an explicit dictionary
        s1 = CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc,
                                     differentials={'s': diff})
        assert len(s1.differentials) == 1
        assert s1.differentials['s'] is diff

        # using the wrong key will cause it to fail
        with pytest.raises(ValueError):
            s1 = CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc,
                                         differentials={'1 / s2': diff})

        # make sure other kwargs are handled properly
        s1 = CartesianRepresentation(x=1, y=2, z=3,
                                     differentials=diff, copy=False, unit=u.kpc)
        assert len(s1.differentials) == 1
        assert s1.differentials['s'] is diff

        with pytest.raises(TypeError):  # invalid type passed to differentials
            CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc,
                                    differentials='garmonbozia')

        # And that one can add it to another representation.
        s1 = CartesianRepresentation(
            CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc),
            differentials=diff)
        assert len(s1.differentials) == 1
        assert s1.differentials['s'] is diff

        # make sure differentials can't accept differentials
        with pytest.raises(TypeError):
            CartesianDifferential(d_x=1 * u.km/u.s, d_y=2 * u.km/u.s,
                                  d_z=3 * u.km/u.s, differentials=diff)

    def test_init_differential_compatible(self):
        # TODO: more extensive checking of this

        # should fail - representation and differential not compatible
        diff = SphericalDifferential(d_lon=1 * u.mas/u.yr,
                                     d_lat=2 * u.mas/u.yr,
                                     d_distance=3 * u.km/u.s)
        with pytest.raises(TypeError):
            CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc,
                                    differentials=diff)

        # should succeed - representation and differential are compatible
        diff = SphericalCosLatDifferential(d_lon_coslat=1 * u.mas/u.yr,
                                           d_lat=2 * u.mas/u.yr,
                                           d_distance=3 * u.km/u.s)

        r1 = SphericalRepresentation(lon=15*u.deg, lat=21*u.deg,
                                     distance=1*u.pc,
                                     differentials=diff)

    def test_init_differential_multiple_equivalent_keys(self):
        d1 = CartesianDifferential(*[1, 2, 3] * u.km/u.s)
        d2 = CartesianDifferential(*[4, 5, 6] * u.km/u.s)

        # verify that the check against expected_unit validates against passing
        # in two different but equivalent keys
        with pytest.raises(ValueError):
            r1 = CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc,
                                         differentials={'s': d1, 'yr': d2})

    def test_init_array_broadcasting(self):

        arr1 = np.arange(8).reshape(4, 2) * u.km/u.s
        diff = CartesianDifferential(d_x=arr1, d_y=arr1, d_z=arr1)

        # shapes aren't compatible
        arr2 = np.arange(27).reshape(3, 9) * u.kpc
        with pytest.raises(ValueError):
            rep = CartesianRepresentation(x=arr2, y=arr2, z=arr2,
                                          differentials=diff)

        arr2 = np.arange(8).reshape(4, 2) * u.kpc
        rep = CartesianRepresentation(x=arr2, y=arr2, z=arr2,
                                      differentials=diff)

        assert rep.x.unit is u.kpc
        assert rep.y.unit is u.kpc
        assert rep.z.unit is u.kpc
        assert len(rep.differentials) == 1
        assert rep.differentials['s'] is diff

        assert rep.xyz.shape == rep.differentials['s'].d_xyz.shape

    def test_reprobj(self):

        # should succeed - representation and differential are compatible
        diff = SphericalCosLatDifferential(d_lon_coslat=1 * u.mas/u.yr,
                                           d_lat=2 * u.mas/u.yr,
                                           d_distance=3 * u.km/u.s)

        r1 = SphericalRepresentation(lon=15*u.deg, lat=21*u.deg,
                                     distance=1*u.pc,
                                     differentials=diff)

        r2 = CartesianRepresentation.from_representation(r1)
        assert r2.get_name() == 'cartesian'
        assert not r2.differentials

        r3 = SphericalRepresentation(r1)
        assert r3.differentials
        assert representation_equal(r3, r1)

    def test_readonly(self):

        s1 = CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc)

        with pytest.raises(AttributeError):  # attribute is not settable
            s1.differentials = 'thing'

    def test_represent_as(self):

        diff = CartesianDifferential(d_x=1 * u.km/u.s,
                                     d_y=2 * u.km/u.s,
                                     d_z=3 * u.km/u.s)
        rep1 = CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc,
                                       differentials=diff)

        # Only change the representation, drop the differential
        new_rep = rep1.represent_as(SphericalRepresentation)
        assert new_rep.get_name() == 'spherical'
        assert not new_rep.differentials  # dropped

        # Pass in separate classes for representation, differential
        new_rep = rep1.represent_as(SphericalRepresentation,
                                    SphericalCosLatDifferential)
        assert new_rep.get_name() == 'spherical'
        assert new_rep.differentials['s'].get_name() == 'sphericalcoslat'

        # Pass in a dictionary for the differential classes
        new_rep = rep1.represent_as(SphericalRepresentation,
                                    {'s': SphericalCosLatDifferential})
        assert new_rep.get_name() == 'spherical'
        assert new_rep.differentials['s'].get_name() == 'sphericalcoslat'

        # make sure represent_as() passes through the differentials
        for name in REPRESENTATION_CLASSES:
            if name == 'radial':
                # TODO: Converting a CartesianDifferential to a
                #       RadialDifferential fails, even on `main`
                continue
            elif name.endswith("geodetic"):
                # TODO: Geodetic representations do not have differentials yet
                continue
            new_rep = rep1.represent_as(REPRESENTATION_CLASSES[name],
                                        DIFFERENTIAL_CLASSES[name])
            assert new_rep.get_name() == name
            assert len(new_rep.differentials) == 1
            assert new_rep.differentials['s'].get_name() == name

        with pytest.raises(ValueError) as excinfo:
            rep1.represent_as('name')
        assert 'use frame object' in str(excinfo.value)

    @pytest.mark.parametrize('sph_diff,usph_diff', [
        (SphericalDifferential, UnitSphericalDifferential),
        (SphericalCosLatDifferential, UnitSphericalCosLatDifferential)])
    def test_represent_as_unit_spherical_with_diff(self, sph_diff, usph_diff):
        """Test that differential angles are correctly reduced."""
        diff = CartesianDifferential(d_x=1 * u.km/u.s,
                                     d_y=2 * u.km/u.s,
                                     d_z=3 * u.km/u.s)
        rep = CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc,
                                      differentials=diff)
        sph = rep.represent_as(SphericalRepresentation, sph_diff)
        usph = rep.represent_as(UnitSphericalRepresentation, usph_diff)
        assert components_equal(usph, sph.represent_as(UnitSphericalRepresentation))
        assert components_equal(usph.differentials['s'],
                                sph.differentials['s'].represent_as(usph_diff))
        # Just to be sure components_equal and the represent_as work as advertised,
        # a sanity check: d_lat is always defined and should be the same.
        assert_array_equal(sph.differentials['s'].d_lat,
                           usph.differentials['s'].d_lat)

    def test_getitem(self):

        d = CartesianDifferential(d_x=np.arange(10) * u.m/u.s,
                                  d_y=-np.arange(10) * u.m/u.s,
                                  d_z=1. * u.m/u.s)
        s = CartesianRepresentation(x=np.arange(10) * u.m,
                                    y=-np.arange(10) * u.m,
                                    z=3 * u.km,
                                    differentials=d)

        s_slc = s[2:8:2]
        s_dif = s_slc.differentials['s']

        assert_allclose_quantity(s_slc.x, [2, 4, 6] * u.m)
        assert_allclose_quantity(s_slc.y, [-2, -4, -6] * u.m)
        assert_allclose_quantity(s_slc.z, [3, 3, 3] * u.km)

        assert_allclose_quantity(s_dif.d_x, [2, 4, 6] * u.m/u.s)
        assert_allclose_quantity(s_dif.d_y, [-2, -4, -6] * u.m/u.s)
        assert_allclose_quantity(s_dif.d_z, [1, 1, 1] * u.m/u.s)

    def test_setitem(self):
        d = CartesianDifferential(d_x=np.arange(5) * u.m/u.s,
                                  d_y=-np.arange(5) * u.m/u.s,
                                  d_z=1. * u.m/u.s)
        s = CartesianRepresentation(x=np.arange(5) * u.m,
                                    y=-np.arange(5) * u.m,
                                    z=3 * u.km,
                                    differentials=d)
        s[:2] = s[2]
        assert_array_equal(s.x, [2, 2, 2, 3, 4] * u.m)
        assert_array_equal(s.y, [-2, -2, -2, -3, -4] * u.m)
        assert_array_equal(s.z, [3, 3, 3, 3, 3] * u.km)
        assert_array_equal(s.differentials['s'].d_x,
                           [2, 2, 2, 3, 4] * u.m/u.s)
        assert_array_equal(s.differentials['s'].d_y,
                           [-2, -2, -2, -3, -4] * u.m/u.s)
        assert_array_equal(s.differentials['s'].d_z,
                           [1, 1, 1, 1, 1] * u.m/u.s)

        s2 = s.represent_as(SphericalRepresentation,
                            SphericalDifferential)

        s[0] = s2[3]
        assert_allclose_quantity(s.x, [3, 2, 2, 3, 4] * u.m)
        assert_allclose_quantity(s.y, [-3, -2, -2, -3, -4] * u.m)
        assert_allclose_quantity(s.z, [3, 3, 3, 3, 3] * u.km)
        assert_allclose_quantity(s.differentials['s'].d_x,
                                 [3, 2, 2, 3, 4] * u.m/u.s)
        assert_allclose_quantity(s.differentials['s'].d_y,
                                 [-3, -2, -2, -3, -4] * u.m/u.s)
        assert_allclose_quantity(s.differentials['s'].d_z,
                                 [1, 1, 1, 1, 1] * u.m/u.s)

        s3 = CartesianRepresentation(s.xyz, differentials={
            's': d,
            's2': CartesianDifferential(np.ones((3, 5))*u.m/u.s**2)})
        with pytest.raises(ValueError, match='same differentials'):
            s[0] = s3[2]

        s4 = SphericalRepresentation(0.*u.deg, 0.*u.deg, 1.*u.kpc,
                                     differentials=RadialDifferential(
                                         10*u.km/u.s))
        with pytest.raises(ValueError, match='loss of information'):
            s[0] = s4

    def test_transform(self):
        d1 = CartesianDifferential(d_x=[1, 2] * u.km/u.s,
                                   d_y=[3, 4] * u.km/u.s,
                                   d_z=[5, 6] * u.km/u.s)
        r1 = CartesianRepresentation(x=[1, 2] * u.kpc,
                                     y=[3, 4] * u.kpc,
                                     z=[5, 6] * u.kpc,
                                     differentials=d1)

        r2 = r1.transform(matrices["general"])
        d2 = r2.differentials['s']
        assert_allclose_quantity(d2.d_x, [22., 28]*u.km/u.s)
        assert_allclose_quantity(d2.d_y, [49, 64]*u.km/u.s)
        assert_allclose_quantity(d2.d_z, [76, 100.]*u.km/u.s)

    def test_with_differentials(self):
        # make sure with_differential correctly creates a new copy with the same
        # differential
        cr = CartesianRepresentation([1, 2, 3]*u.kpc)
        diff = CartesianDifferential([.1, .2, .3]*u.km/u.s)
        cr2 = cr.with_differentials(diff)
        assert cr.differentials != cr2.differentials
        assert cr2.differentials['s'] is diff

        # make sure it works even if a differential is present already
        diff2 = CartesianDifferential([.1, .2, .3]*u.m/u.s)
        cr3 = CartesianRepresentation([1, 2, 3]*u.kpc, differentials=diff)
        cr4 = cr3.with_differentials(diff2)
        assert cr4.differentials['s'] != cr3.differentials['s']
        assert cr4.differentials['s'] == diff2

        # also ensure a *scalar* differential will works
        cr5 = cr.with_differentials(diff)
        assert len(cr5.differentials) == 1
        assert cr5.differentials['s'] == diff

        # make sure we don't update the original representation's dict
        d1 = CartesianDifferential(*np.random.random((3, 5)), unit=u.km/u.s)
        d2 = CartesianDifferential(*np.random.random((3, 5)), unit=u.km/u.s**2)
        r1 = CartesianRepresentation(*np.random.random((3, 5)), unit=u.pc,
                                     differentials=d1)

        r2 = r1.with_differentials(d2)
        assert r1.differentials['s'] is r2.differentials['s']
        assert 's2' not in r1.differentials
        assert 's2' in r2.differentials


def test_repr_with_differentials():
    diff = CartesianDifferential([.1, .2, .3]*u.km/u.s)
    cr = CartesianRepresentation([1, 2, 3]*u.kpc, differentials=diff)
    assert "has differentials w.r.t.: 's'" in repr(cr)


def test_to_cartesian():
    """
    Test that to_cartesian drops the differential.
    """
    sd = SphericalDifferential(d_lat=1*u.deg, d_lon=2*u.deg, d_distance=10*u.m)
    sr = SphericalRepresentation(lat=1*u.deg, lon=2*u.deg, distance=10*u.m,
                                 differentials=sd)

    cart = sr.to_cartesian()
    assert cart.get_name() == 'cartesian'
    assert not cart.differentials


@pytest.fixture
def unitphysics():
    """
    This fixture is used
    """
    had_unit = False
    if hasattr(PhysicsSphericalRepresentation, '_unit_representation'):
        orig = PhysicsSphericalRepresentation._unit_representation
        had_unit = True

    class UnitPhysicsSphericalRepresentation(BaseRepresentation):
        attr_classes = {'phi': Angle,
                        'theta': Angle}

        def __init__(self, *args, copy=True, **kwargs):
            super().__init__(*args, copy=copy, **kwargs)

            # Wrap/validate phi/theta
            if copy:
                self._phi = self._phi.wrap_at(360 * u.deg)
            else:
                # necessary because the above version of `wrap_at` has to be a copy
                self._phi.wrap_at(360 * u.deg, inplace=True)

            if np.any(self._theta < 0.*u.deg) or np.any(self._theta > 180.*u.deg):
                raise ValueError('Inclination angle(s) must be within '
                                 '0 deg <= angle <= 180 deg, '
                                 'got {}'.format(self._theta.to(u.degree)))

        @property
        def phi(self):
            return self._phi

        @property
        def theta(self):
            return self._theta

        def unit_vectors(self):
            sinphi, cosphi = np.sin(self.phi), np.cos(self.phi)
            sintheta, costheta = np.sin(self.theta), np.cos(self.theta)
            return {
                'phi': CartesianRepresentation(-sinphi, cosphi, 0., copy=False),
                'theta': CartesianRepresentation(costheta*cosphi,
                                                 costheta*sinphi,
                                                 -sintheta, copy=False)}

        def scale_factors(self):
            sintheta = np.sin(self.theta)
            l = np.broadcast_to(1.*u.one, self.shape, subok=True)
            return {'phi', sintheta,
                    'theta', l}

        def to_cartesian(self):
            x = np.sin(self.theta) * np.cos(self.phi)
            y = np.sin(self.theta) * np.sin(self.phi)
            z = np.cos(self.theta)

            return CartesianRepresentation(x=x, y=y, z=z, copy=False)

        @classmethod
        def from_cartesian(cls, cart):
            """
            Converts 3D rectangular cartesian coordinates to spherical polar
            coordinates.
            """
            s = np.hypot(cart.x, cart.y)

            phi = np.arctan2(cart.y, cart.x)
            theta = np.arctan2(s, cart.z)

            return cls(phi=phi, theta=theta, copy=False)

        def norm(self):
            return u.Quantity(np.ones(self.shape), u.dimensionless_unscaled,
                              copy=False)

    PhysicsSphericalRepresentation._unit_representation = UnitPhysicsSphericalRepresentation
    yield UnitPhysicsSphericalRepresentation

    if had_unit:
        PhysicsSphericalRepresentation._unit_representation = orig
    else:
        del PhysicsSphericalRepresentation._unit_representation

    # remove from the module-level representations, if present
    REPRESENTATION_CLASSES.pop(UnitPhysicsSphericalRepresentation.get_name(), None)


def test_unitphysics(unitphysics):
    obj = unitphysics(phi=0*u.deg, theta=10*u.deg)
    objkw = unitphysics(phi=0*u.deg, theta=10*u.deg)
    assert objkw.phi == obj.phi
    assert objkw.theta == obj.theta

    asphys = obj.represent_as(PhysicsSphericalRepresentation)
    assert asphys.phi == obj.phi
    assert_allclose(asphys.theta, obj.theta)
    assert_allclose_quantity(asphys.r, 1*u.dimensionless_unscaled)

    assph = obj.represent_as(SphericalRepresentation)
    assert assph.lon == obj.phi
    assert assph.lat == 80*u.deg
    assert_allclose_quantity(assph.distance, 1*u.dimensionless_unscaled)

    with pytest.raises(TypeError, match='got multiple values'):
        unitphysics(1*u.deg, 2*u.deg, theta=10)

    with pytest.raises(TypeError, match='unexpected keyword.*parrot'):
        unitphysics(1*u.deg, 2*u.deg, parrot=10)


def test_distance_warning(recwarn):
    SphericalRepresentation(1*u.deg, 2*u.deg, 1*u.kpc)
    with pytest.raises(ValueError) as excinfo:
        SphericalRepresentation(1*u.deg, 2*u.deg, -1*u.kpc)
    assert 'Distance must be >= 0' in str(excinfo.value)
    # second check is because the "originating" ValueError says the above,
    # while the representation one includes the below
    assert 'you must explicitly pass' in str(excinfo.value)


def test_dtype_preservation_in_indexing():
    # Regression test for issue #8614 (fixed in #8876)
    xyz = np.array([[1, 0, 0], [0.9, 0.1, 0]], dtype='f4')
    cr = CartesianRepresentation(xyz, xyz_axis=-1, unit="km")
    assert cr.xyz.dtype == xyz.dtype
    cr0 = cr[0]
    # This used to fail.
    assert cr0.xyz.dtype == xyz.dtype


class TestInfo:
    def setup_class(cls):
        cls.rep = SphericalRepresentation([0, 1]*u.deg, [2, 3]*u.deg,
                                          10*u.pc)
        cls.diff = SphericalDifferential([10, 20]*u.mas/u.yr,
                                         [30, 40]*u.mas/u.yr,
                                         [50, 60]*u.km/u.s)
        cls.rep_w_diff = SphericalRepresentation(cls.rep,
                                                 differentials=cls.diff)

    def test_info_unit(self):
        assert self.rep.info.unit == 'deg, deg, pc'
        assert self.diff.info.unit == 'mas / yr, mas / yr, km / s'
        assert self.rep_w_diff.info.unit == 'deg, deg, pc'

    @pytest.mark.parametrize('item', ['rep', 'diff', 'rep_w_diff'])
    def test_roundtrip(self, item):
        rep_or_diff = getattr(self, item)
        as_dict = rep_or_diff.info._represent_as_dict()
        new = rep_or_diff.__class__.info._construct_from_dict(as_dict)
        assert np.all(representation_equal(new, rep_or_diff))


@pytest.mark.parametrize('cls',
                         [SphericalDifferential,
                          SphericalCosLatDifferential,
                          CylindricalDifferential,
                          PhysicsSphericalDifferential,
                          UnitSphericalDifferential,
                          UnitSphericalCosLatDifferential])
def test_differential_norm_noncartesian(cls):
    # The norm of a non-Cartesian differential without specifying `base` should error
    rep = cls(0, 0, 0)
    with pytest.raises(ValueError, match=r"`base` must be provided .* " + cls.__name__):
        rep.norm()


def test_differential_norm_radial():
    # Unlike most non-Cartesian differentials, the norm of a radial differential does not require `base`
    rep = RadialDifferential(1*u.km/u.s)
    assert_allclose_quantity(rep.norm(), 1*u.km/u.s)
