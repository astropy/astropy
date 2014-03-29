import numpy as np
from numpy.testing import assert_allclose

from ... import units as u
from ...tests.helper import pytest
from ..angles import Longitude, Latitude
from ..distances import Distance
from ..representation import (SphericalRepresentation,
                              CartesianRepresentation)


def assert_allclose_quantity(q1, q2):
    assert_allclose(q1.value, q2.to(q1.unit).value)

class TestSphericalRepresentation(object):

    def test_empty_init(self):
        with pytest.raises(ValueError) as exc:
            s = SphericalRepresentation()
        assert exc.value.args[0] == "lon and lat are required to instantiate SphericalRepresentation"

    def test_init_quantity(self):

        s1 = SphericalRepresentation(lon=8 * u.hour, lat=5 * u.deg)
        assert s1.lon == 8. * u.hourangle
        assert s1.lat == 5. * u.deg
        assert s1.distance is None

        assert isinstance(s1.lon, Longitude)
        assert isinstance(s1.lat, Latitude)

        s2 = SphericalRepresentation(lon=8 * u.hourangle, lat=5 * u.deg)
        assert s2.lon == 8. * u.hourangle
        assert s2.lat == 5. * u.deg
        assert s2.distance is None

        assert isinstance(s2.lon, Longitude)
        assert isinstance(s2.lat, Latitude)

        s3 = SphericalRepresentation(lon=8 * u.hourangle, lat=5 * u.deg, distance=10 * u.kpc)
        assert s3.lon == 8. * u.hourangle
        assert s3.lat == 5. * u.deg
        assert s3.distance == 10 * u.kpc

        assert isinstance(s3.lon, Longitude)
        assert isinstance(s3.lat, Latitude)
        assert isinstance(s3.distance, Distance)

    def test_init_lonlat(self):

        s1 = SphericalRepresentation(Longitude(8, u.hour),
                                     Latitude(5, u.deg))

        assert s1.lon == 8. * u.hourangle
        assert s1.lat == 5. * u.deg
        assert s1.distance is None

        assert isinstance(s1.lon, Longitude)
        assert isinstance(s1.lat, Latitude)

        s2 = SphericalRepresentation(Longitude(8, u.hour),
                                     Latitude(5, u.deg),
                                     Distance(10, u.kpc))

        assert s2.lon == 8. * u.hourangle
        assert s2.lat == 5. * u.deg
        assert s2.distance == 10. * u.kpc

        assert isinstance(s2.lon, Longitude)
        assert isinstance(s2.lat, Latitude)
        assert isinstance(s2.distance, Distance)

    def test_init_array(self):

        s1 = SphericalRepresentation(lon=[8, 9] * u.hourangle,
                                     lat=[5, 6] * u.deg)

        assert_allclose(s1.lon.degree, [120, 135])
        assert_allclose(s1.lat.degree, [5, 6])

        assert isinstance(s1.lon, Longitude)
        assert isinstance(s1.lat, Latitude)

    def test_init_array_nocopy(self):

        lon = [8, 9] * u.hourangle
        lat = [5, 6] * u.deg

        s1 = SphericalRepresentation(lon=lon, lat=lat, copy=False)

        # Not sure if we can expect this to work
        assert s1.lon is lon
        assert s1.lat is lat

    def test_init_str(self):

        s1 = SphericalRepresentation(lon='2h6m3.3s', lat='0.1rad')
        assert_allclose(s1.lon.degree, 31.513749999999995)
        assert_allclose(s1.lat.degree, 5.729577951308233)

    def test_reprobj(self):

        s1 = SphericalRepresentation(lon=8 * u.hourangle, lat=5 * u.deg, distance=10 * u.kpc)

        s2 = SphericalRepresentation(representation=s1)

        assert_allclose_quantity(s2.lon, 8. * u.hourangle)
        assert_allclose_quantity(s2.lat, 5. * u.deg)
        assert_allclose_quantity(s2.distance, 10 * u.kpc)

    def test_reprobj_invalid(self):

        s1 = SphericalRepresentation(lon=8 * u.hourangle, lat=5 * u.deg, distance=10 * u.kpc)

        with pytest.raises(ValueError) as exc:
            s2 = SphericalRepresentation(lon=8 * u.hourangle, lat=5 * u.deg, representation=s1)
        assert exc.value.args[0] == "If representation is passed, no other arguments can be passed"

    def test_broadcasting(self):

        s1 = SphericalRepresentation(lon=[8, 9]*u.hourangle,
                                     lat=[5, 6]*u.deg,
                                     distance=10*u.kpc)

        assert_allclose_quantity(s1.lon, [120, 135] * u.degree)
        assert_allclose_quantity(s1.lat, [5, 6] * u.degree)
        assert_allclose_quantity(s1.distance, [10, 10] * u.kpc)

    def test_broadcasting_mismatch(self):

        with pytest.raises(ValueError) as exc:
            s1 = SphericalRepresentation(lon=[8, 9, 10]*u.hourangle,
                                         lat=[5, 6]*u.deg)
        assert exc.value.args[0] == "Input parameters lon, lat, and distance cannot be broadcast"

    def test_mixed_units(self):

        # It's also possible to pass in scalar quantity lists with mixed
        # units. These are converted to array quantities following the same
        # rule as `Quantity`: all elements are converted to match the first
        # element's units.

        s1 = SphericalRepresentation(lon=[8*u.hourangle, 135*u.deg],
                                     lat=[5*u.deg, (6*np.pi/180)*u.rad])
        assert s1.lat.unit == u.deg
        assert s1.lon.unit == u.hourangle
        assert_allclose(s1.lon.value, [8,9])
        assert_allclose(s1.lat.value, [5,6])

    def test_readonly(self):

        s1 = SphericalRepresentation(lon=[8*u.hourangle, 135*u.deg],
                                     lat=[5*u.deg, (6*np.pi/180)*u.rad])

        with pytest.raises(AttributeError):
            s1.lon = 1. * u.deg

        with pytest.raises(AttributeError):
            s1.lat = 1. * u.deg

        with pytest.raises(AttributeError):
            s1.distance = 1. * u.kpc


class TestCartesianRepresentation(object):

    def test_empty_init(self):
        with pytest.raises(ValueError) as exc:
            s = CartesianRepresentation()
        assert exc.value.args[0] == "x, y, and z are required to instantiate CartesianRepresentation"

    def test_init_quantity(self):

        s1 = CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc)

        assert s1.x.unit is u.kpc
        assert s1.y.unit is u.kpc
        assert s1.z.unit is u.kpc

        assert_allclose(s1.x.value,1)
        assert_allclose(s1.y.value,2)
        assert_allclose(s1.z.value,3)

    def test_init_singleunit(self):

        s1 = CartesianRepresentation(x=1,y=2,z=3,unit=u.kpc)

        assert s1.x.unit is u.kpc
        assert s1.y.unit is u.kpc
        assert s1.z.unit is u.kpc

        assert_allclose(s1.x.value,1)
        assert_allclose(s1.y.value,2)
        assert_allclose(s1.z.value,3)

    def test_init_override_unit(self):

        s1 = CartesianRepresentation(x=1 * u.pc,y=2 * u.Mpc,z=3 * u.kpc,unit=u.kpc)

        assert s1.x.unit is u.kpc
        assert s1.y.unit is u.kpc
        assert s1.z.unit is u.kpc

        assert_allclose(s1.x.value,0.001)
        assert_allclose(s1.y.value,2000)
        assert_allclose(s1.z.value,3)

    def test_init_array(self):

        s1 = CartesianRepresentation(x=[1,2,3] * u.pc,
                                     y=[2,3,4]* u.Mpc,
                                     z=[3,4,5] * u.kpc)

        assert s1.x.unit is u.pc
        assert s1.y.unit is u.Mpc
        assert s1.z.unit is u.kpc

        assert_allclose(s1.x.value,[1,2,3])
        assert_allclose(s1.y.value,[2,3,4])
        assert_allclose(s1.z.value,[3,4,5])

    def test_init_array_nocopy(self):

        x = [8, 9, 10] * u.pc
        y = [5, 6, 7] * u.Mpc
        z = [2, 3, 4] * u.kpc

        s1 = CartesianRepresentation(x=x, y=y, z=z, copy=False)

        # Not sure if we can expect this to work
        assert s1.x is x
        assert s1.y is y
        assert s1.z is z

    def test_reprobj(self):

        s1 = CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc)

        s2 = CartesianRepresentation(representation=s1)

        assert s2.x == 1 * u.kpc
        assert s2.y == 2 * u.kpc
        assert s2.z == 3 * u.kpc

    def test_reprobj_invalid(self):

        s1 = CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc)

        with pytest.raises(ValueError) as exc:
            s1 = CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc, representation=s1)
        assert exc.value.args[0] == "If representation is passed, no other arguments can be passed"

    def test_broadcasting(self):

        s1 = CartesianRepresentation(x=[1,2] * u.kpc, y=[3,4] * u.kpc, z=5 * u.kpc)

        assert s1.x.unit == u.kpc
        assert s1.y.unit == u.kpc
        assert s1.z.unit == u.kpc

        assert_allclose(s1.x.value, [1, 2])
        assert_allclose(s1.y.value, [3, 4])
        assert_allclose(s1.z.value, [5, 5])

    def test_broadcasting_mismatch(self):

        with pytest.raises(ValueError) as exc:
            s1 = CartesianRepresentation(x=[1,2] * u.kpc, y=[3,4] * u.kpc, z=[5,6,7] * u.kpc)
        assert exc.value.args[0] == "Input parameters x, y, and z cannot be broadcast"

    def test_mixed_units(self):

        # It's also possible to pass in scalar quantity lists with mixed
        # units. These are converted to array quantities following the same
        # rule as `Quantity`: all elements are converted to match the first
        # element's units.

        s1 = CartesianRepresentation(x=[1 * u.kpc,2 * u.Mpc],
                                     y=[3 * u.kpc, 4 * u.pc] ,
                                     z=[5. * u.cm, 6 * u.m])

        assert s1.x.unit == u.kpc
        assert s1.y.unit == u.kpc
        assert s1.z.unit == u.cm
        assert_allclose(s1.x.value, [1, 2000])
        assert_allclose(s1.y.value, [3, 0.004])
        assert_allclose(s1.z.value, [5, 600])

    def test_readonly(self):

        s1 = CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc)

        with pytest.raises(AttributeError):
            s1.x = 1. * u.kpc

        with pytest.raises(AttributeError):
            s1.y = 1. * u.kpc

        with pytest.raises(AttributeError):
            s1.z = 1. * u.kpc


def test_cartesian_spherical_roundtrip():

    s1 = CartesianRepresentation(x=[1 * u.kpc,2 * u.Mpc],
                                 y=[3 * u.kpc, 4 * u.pc] ,
                                 z=[5. * u.cm, 6 * u.m])

    s2 = SphericalRepresentation(representation=s1)

    s3 = CartesianRepresentation(representation=s2)

    s4 = SphericalRepresentation(representation=s3)

    assert_allclose_quantity(s1.x, s3.x)
    assert_allclose_quantity(s1.y, s3.y)
    assert_allclose_quantity(s1.z, s3.z)

    assert_allclose_quantity(s2.lon, s4.lon)
    assert_allclose_quantity(s2.lat, s4.lat)
    assert_allclose_quantity(s2.distance, s4.distance)
