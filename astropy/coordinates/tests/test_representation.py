# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from copy import deepcopy
from collections import OrderedDict

import pytest
import numpy as np
from numpy.testing import assert_allclose

from ... import units as u
from ...tests.helper import (assert_quantity_allclose as
                             assert_allclose_quantity)
from ...utils import isiterable
from ...utils.compat import NUMPY_LT_1_14
from ..angles import Longitude, Latitude, Angle
from ..distances import Distance
from ..representation import (REPRESENTATION_CLASSES,
                              DIFFERENTIAL_CLASSES,
                              BaseRepresentation,
                              SphericalRepresentation,
                              UnitSphericalRepresentation,
                              SphericalCosLatDifferential,
                              CartesianRepresentation,
                              CylindricalRepresentation,
                              PhysicsSphericalRepresentation,
                              CartesianDifferential,
                              SphericalDifferential,
                              _combine_xyz)


# Preserve the original REPRESENTATION_CLASSES dict so that importing
#   the test file doesn't add a persistent test subclass (LogDRepresentation)
def setup_function(func):
    func.REPRESENTATION_CLASSES_ORIG = deepcopy(REPRESENTATION_CLASSES)


def teardown_function(func):
    REPRESENTATION_CLASSES.clear()
    REPRESENTATION_CLASSES.update(func.REPRESENTATION_CLASSES_ORIG)


class TestSphericalRepresentation(object):

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


class TestUnitSphericalRepresentation(object):

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


class TestPhysicsSphericalRepresentation(object):

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


class TestCartesianRepresentation(object):

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
        assert 'xyz_axis should only be set' in str(exc)

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

    def test_reprobj(self):

        s1 = CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc)

        s2 = CartesianRepresentation.from_representation(s1)

        assert s2.x == 1 * u.kpc
        assert s2.y == 2 * u.kpc
        assert s2.z == 3 * u.kpc

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

        s1 = CartesianRepresentation(x=[1, 2] * u.kpc, y=[3, 4] * u.kpc, z=[5, 6] * u.kpc)

        matrix = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

        s2 = s1.transform(matrix)

        assert_allclose(s2.x.value, [1 * 1 + 2 * 3 + 3 * 5, 1 * 2 + 2 * 4 + 3 * 6])
        assert_allclose(s2.y.value, [4 * 1 + 5 * 3 + 6 * 5, 4 * 2 + 5 * 4 + 6 * 6])
        assert_allclose(s2.z.value, [7 * 1 + 8 * 3 + 9 * 5, 7 * 2 + 8 * 4 + 9 * 6])

        assert s2.x.unit is u.kpc
        assert s2.y.unit is u.kpc
        assert s2.z.unit is u.kpc


class TestCylindricalRepresentation(object):

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
                        '    ({})>').format(' 1.,  2.5,  1.' if NUMPY_LT_1_14
                                            else '1., 2.5, 1.')

    r2 = CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc)
    assert repr(r2) == ('<CartesianRepresentation (x, y, z) in kpc\n'
                        '    ({})>').format(' 1.,  2.,  3.' if NUMPY_LT_1_14
                                            else '1., 2., 3.')

    r3 = CartesianRepresentation(x=[1, 2, 3] * u.kpc, y=4 * u.kpc, z=[9, 10, 11] * u.kpc)
    if NUMPY_LT_1_14:
        assert repr(r3) == ('<CartesianRepresentation (x, y, z) in kpc\n'
                            '    [( 1.,  4.,   9.), ( 2.,  4.,  10.), ( 3.,  4.,  11.)]>')
    else:
        assert repr(r3) == ('<CartesianRepresentation (x, y, z) in kpc\n'
                            '    [(1., 4.,  9.), (2., 4., 10.), (3., 4., 11.)]>')


def test_representation_repr_multi_d():
    """Regression test for #5889."""
    cr = CartesianRepresentation(np.arange(27).reshape(3, 3, 3), unit='m')
    if NUMPY_LT_1_14:
        assert repr(cr) == (
            '<CartesianRepresentation (x, y, z) in m\n'
            '    [[( 0.,   9.,  18.), ( 1.,  10.,  19.), ( 2.,  11.,  20.)],\n'
            '     [( 3.,  12.,  21.), ( 4.,  13.,  22.), ( 5.,  14.,  23.)],\n'
            '     [( 6.,  15.,  24.), ( 7.,  16.,  25.), ( 8.,  17.,  26.)]]>')
    else:
        assert repr(cr) == (
            '<CartesianRepresentation (x, y, z) in m\n'
            '    [[(0.,  9., 18.), (1., 10., 19.), (2., 11., 20.)],\n'
            '     [(3., 12., 21.), (4., 13., 22.), (5., 14., 23.)],\n'
            '     [(6., 15., 24.), (7., 16., 25.), (8., 17., 26.)]]>')
    # This was broken before.
    if NUMPY_LT_1_14:
        assert repr(cr.T) == (
            '<CartesianRepresentation (x, y, z) in m\n'
            '    [[( 0.,   9.,  18.), ( 3.,  12.,  21.), ( 6.,  15.,  24.)],\n'
            '     [( 1.,  10.,  19.), ( 4.,  13.,  22.), ( 7.,  16.,  25.)],\n'
            '     [( 2.,  11.,  20.), ( 5.,  14.,  23.), ( 8.,  17.,  26.)]]>')
    else:
        assert repr(cr.T) == (
            '<CartesianRepresentation (x, y, z) in m\n'
            '    [[(0.,  9., 18.), (3., 12., 21.), (6., 15., 24.)],\n'
            '     [(1., 10., 19.), (4., 13., 22.), (7., 16., 25.)],\n'
            '     [(2., 11., 20.), (5., 14., 23.), (8., 17., 26.)]]>')


def test_representation_str():
    r1 = SphericalRepresentation(lon=1 * u.deg, lat=2.5 * u.deg, distance=1 * u.kpc)
    assert str(r1) == ('( 1.,  2.5,  1.) (deg, deg, kpc)' if NUMPY_LT_1_14 else
                       '(1., 2.5, 1.) (deg, deg, kpc)')
    r2 = CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc)
    assert str(r2) == ('( 1.,  2.,  3.) kpc' if NUMPY_LT_1_14 else
                       '(1., 2., 3.) kpc')
    r3 = CartesianRepresentation(x=[1, 2, 3] * u.kpc, y=4 * u.kpc, z=[9, 10, 11] * u.kpc)
    assert str(r3) == ('[( 1.,  4.,   9.), ( 2.,  4.,  10.), ( 3.,  4.,  11.)] kpc'
                       if NUMPY_LT_1_14 else
                       '[(1., 4.,  9.), (2., 4., 10.), (3., 4., 11.)] kpc')


def test_representation_str_multi_d():
    """Regression test for #5889."""
    cr = CartesianRepresentation(np.arange(27).reshape(3, 3, 3), unit='m')
    if NUMPY_LT_1_14:
        assert str(cr) == (
            '[[( 0.,   9.,  18.), ( 1.,  10.,  19.), ( 2.,  11.,  20.)],\n'
            ' [( 3.,  12.,  21.), ( 4.,  13.,  22.), ( 5.,  14.,  23.)],\n'
            ' [( 6.,  15.,  24.), ( 7.,  16.,  25.), ( 8.,  17.,  26.)]] m')
    else:
        assert str(cr) == (
            '[[(0.,  9., 18.), (1., 10., 19.), (2., 11., 20.)],\n'
            ' [(3., 12., 21.), (4., 13., 22.), (5., 14., 23.)],\n'
            ' [(6., 15., 24.), (7., 16., 25.), (8., 17., 26.)]] m')
    # This was broken before.
    if NUMPY_LT_1_14:
        assert str(cr.T) == (
            '[[( 0.,   9.,  18.), ( 3.,  12.,  21.), ( 6.,  15.,  24.)],\n'
            ' [( 1.,  10.,  19.), ( 4.,  13.,  22.), ( 7.,  16.,  25.)],\n'
            ' [( 2.,  11.,  20.), ( 5.,  14.,  23.), ( 8.,  17.,  26.)]] m')
    else:
        assert str(cr.T) == (
            '[[(0.,  9., 18.), (3., 12., 21.), (6., 15., 24.)],\n'
            ' [(1., 10., 19.), (4., 13., 22.), (7., 16., 25.)],\n'
            ' [(2., 11., 20.), (5., 14., 23.), (8., 17., 26.)]] m')


def test_subclass_representation():
    from ..builtin_frames import ICRS

    class Longitude180(Longitude):
        def __new__(cls, angle, unit=None, wrap_angle=180 * u.deg, **kwargs):
            self = super(Longitude180, cls).__new__(cls, angle, unit=unit,
                                                    wrap_angle=wrap_angle, **kwargs)
            return self

    class SphericalWrap180Representation(SphericalRepresentation):
        attr_classes = OrderedDict([('lon', Longitude180),
                                    ('lat', Latitude),
                                    ('distance', u.Quantity)])
        recommended_units = {'lon': u.deg, 'lat': u.deg}

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
        attr_classes = OrderedDict([('lon', Longitude),
                                    ('lat', Latitude),
                                    ('logd', u.Dex)])

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

    with pytest.raises(ValueError):
        # check we cannot redefine an existing class.
        class LogDRepresentation(BaseRepresentation):
            attr_classes = OrderedDict([('lon', Longitude),
                                        ('lat', Latitude),
                                        ('logr', u.Dex)])


def test_combine_xyz():

    x, y, z = np.arange(27).reshape(3, 9) * u.kpc
    xyz = _combine_xyz(x, y, z, xyz_axis=0)
    assert xyz.shape == (3, 9)
    assert np.all(xyz[0] == x)
    assert np.all(xyz[1] == y)
    assert np.all(xyz[2] == z)

    x, y, z = np.arange(27).reshape(3, 3, 3) * u.kpc
    xyz = _combine_xyz(x, y, z, xyz_axis=0)
    assert xyz.ndim == 3
    assert np.all(xyz[0] == x)
    assert np.all(xyz[1] == y)
    assert np.all(xyz[2] == z)

    xyz = _combine_xyz(x, y, z, xyz_axis=1)
    assert xyz.ndim == 3
    assert np.all(xyz[:, 0] == x)
    assert np.all(xyz[:, 1] == y)
    assert np.all(xyz[:, 2] == z)

    xyz = _combine_xyz(x, y, z, xyz_axis=-1)
    assert xyz.ndim == 3
    assert np.all(xyz[..., 0] == x)
    assert np.all(xyz[..., 1] == y)
    assert np.all(xyz[..., 2] == z)


class TestCartesianRepresentationWithDifferential(object):

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
                #       RadialDifferential fails, even on `master`
                continue
            new_rep = rep1.represent_as(REPRESENTATION_CLASSES[name],
                                        DIFFERENTIAL_CLASSES[name])
            assert new_rep.get_name() == name
            assert len(new_rep.differentials) == 1
            assert new_rep.differentials['s'].get_name() == name

        with pytest.raises(ValueError) as excinfo:
            rep1.represent_as('name')
        assert 'use frame object' in str(excinfo.value)

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

    def test_transform(self):
        d1 = CartesianDifferential(d_x=[1, 2] * u.km/u.s,
                                   d_y=[3, 4] * u.km/u.s,
                                   d_z=[5, 6] * u.km/u.s)
        r1 = CartesianRepresentation(x=[1, 2] * u.kpc,
                                     y=[3, 4] * u.kpc,
                                     z=[5, 6] * u.kpc,
                                     differentials=d1)

        matrix = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

        r2 = r1.transform(matrix)
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
