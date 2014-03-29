from numpy.testing import assert_allclose

from ... import units as u
from ...tests.helper import pytest
from ..angles import Longitude, Latitude
from ..distances import Distance
from ..representation import (SphericalRepresentation,
                              CartesianRepresentation)


class TestSphericalRepresentation(object):

    def test_empty_init(self):
        with pytest.raises(ValueError) as exc:
            s = SphericalRepresentation()
        assert exc.value.args[0] == "lon and lat are required to instantiate a spherical representation"

    def test_init_quantity(self):

        s1 = SphericalRepresentation(lon=8 * u.hour, lat=5 * u.deg)
        assert s1.lon == 8. * u.hourangle
        assert s1.lat == 5. * u.deg
        assert s1.distance is None

        s2 = SphericalRepresentation(lon=8 * u.hourangle, lat=5 * u.deg)
        assert s2.lon == 8. * u.hourangle
        assert s2.lat == 5. * u.deg
        assert s2.distance is None

        s3 = SphericalRepresentation(lon=8 * u.hourangle, lat=5 * u.deg, distance=10 * u.kpc)
        assert s3.lon == 8. * u.hourangle
        assert s3.lat == 5. * u.deg
        assert s3.distance == 10 * u.kpc

    def test_init_lonlat(self):

        s1 = SphericalRepresentation(Longitude(8, u.hour),
                                     Latitude(5, u.deg))

        assert s1.lon == 8. * u.hourangle
        assert s1.lat == 5. * u.deg
        assert s1.distance is None

        s2 = SphericalRepresentation(Longitude(8, u.hour),
                                     Latitude(5, u.deg),
                                     Distance(10, u.kpc))

        assert s1.lon == 8. * u.hourangle
        assert s1.lat == 5. * u.deg
        assert s1.distance == 10. * u.kpc

    def test_init_array(self):

        s1 = SphericalRepresentation(lon=[8, 9] * u.hourangle,
                                     lat=[5, 6] * u.deg)

        assert_allclose(s1.lon.degree, [120, 135])
        assert_allclose(s1.lat.degree, [5, 6])

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
