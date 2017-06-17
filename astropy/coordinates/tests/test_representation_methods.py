# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

import pytest
import numpy as np

from ... import units as u
from .. import (SphericalRepresentation, Longitude, Latitude,
                SphericalDifferential)
from ...utils.compat.numpy import broadcast_to as np_broadcast_to


class TestManipulation():
    """Manipulation of Representation shapes.

    Checking that attributes are manipulated correctly.

    Even more exhaustive tests are done in time.tests.test_methods
    """

    def setup(self):
        lon = Longitude(np.arange(0, 24, 4), u.hourangle)
        lat = Latitude(np.arange(-90, 91, 30), u.deg)

        # With same-sized arrays
        self.s0 = SphericalRepresentation(
            lon[:, np.newaxis] * np.ones(lat.shape),
            lat * np.ones(lon.shape)[:, np.newaxis],
            np.ones(lon.shape + lat.shape) * u.kpc)

        self.diff = SphericalDifferential(
            d_lon=np.ones(self.s0.shape)*u.mas/u.yr,
            d_lat=np.ones(self.s0.shape)*u.mas/u.yr,
            d_distance=np.ones(self.s0.shape)*u.km/u.s)
        self.s0 = self.s0.with_differentials(self.diff)

        # With unequal arrays -> these will be broadcasted.
        self.s1 = SphericalRepresentation(lon[:, np.newaxis], lat, 1. * u.kpc,
                                          differentials=self.diff)

        # For completeness on some tests, also a cartesian one
        self.c0 = self.s0.to_cartesian()

    def test_ravel(self):
        s0_ravel = self.s0.ravel()
        assert type(s0_ravel) is type(self.s0)
        assert s0_ravel.shape == (self.s0.size,)
        assert np.all(s0_ravel.lon == self.s0.lon.ravel())
        assert np.may_share_memory(s0_ravel.lon, self.s0.lon)
        assert np.may_share_memory(s0_ravel.lat, self.s0.lat)
        assert np.may_share_memory(s0_ravel.distance, self.s0.distance)
        assert s0_ravel.differentials['s'].shape == (self.s0.size,)

        # Since s1 was broadcast, ravel needs to make a copy.
        s1_ravel = self.s1.ravel()
        assert type(s1_ravel) is type(self.s1)
        assert s1_ravel.shape == (self.s1.size,)
        assert s1_ravel.differentials['s'].shape == (self.s1.size,)
        assert np.all(s1_ravel.lon == self.s1.lon.ravel())
        assert not np.may_share_memory(s1_ravel.lat, self.s1.lat)

    def test_copy(self):
        s0_copy = self.s0.copy()
        s0_copy_diff = s0_copy.differentials['s']
        assert s0_copy.shape == self.s0.shape
        assert np.all(s0_copy.lon == self.s0.lon)
        assert np.all(s0_copy.lat == self.s0.lat)

        # Check copy was made of internal data.
        assert not np.may_share_memory(s0_copy.distance, self.s0.distance)
        assert not np.may_share_memory(s0_copy_diff.d_lon, self.diff.d_lon)

    def test_flatten(self):
        s0_flatten = self.s0.flatten()
        s0_diff = s0_flatten.differentials['s']
        assert s0_flatten.shape == (self.s0.size,)
        assert s0_diff.shape == (self.s0.size,)
        assert np.all(s0_flatten.lon == self.s0.lon.flatten())
        assert np.all(s0_diff.d_lon == self.diff.d_lon.flatten())

        # Flatten always copies.
        assert not np.may_share_memory(s0_flatten.distance, self.s0.distance)
        assert not np.may_share_memory(s0_diff.d_lon, self.diff.d_lon)

        s1_flatten = self.s1.flatten()
        assert s1_flatten.shape == (self.s1.size,)
        assert np.all(s1_flatten.lon == self.s1.lon.flatten())
        assert not np.may_share_memory(s1_flatten.lat, self.s1.lat)

    def test_transpose(self):
        s0_transpose = self.s0.transpose()
        s0_diff = s0_transpose.differentials['s']
        assert s0_transpose.shape == (7, 6)
        assert s0_diff.shape == s0_transpose.shape
        assert np.all(s0_transpose.lon == self.s0.lon.transpose())
        assert np.all(s0_diff.d_lon == self.diff.d_lon.transpose())
        assert np.may_share_memory(s0_transpose.distance, self.s0.distance)
        assert np.may_share_memory(s0_diff.d_lon, self.diff.d_lon)

        s1_transpose = self.s1.transpose()
        s1_diff = s1_transpose.differentials['s']
        assert s1_transpose.shape == (7, 6)
        assert s1_diff.shape == s1_transpose.shape
        assert np.all(s1_transpose.lat == self.s1.lat.transpose())
        assert np.all(s1_diff.d_lon == self.diff.d_lon.transpose())
        assert np.may_share_memory(s1_transpose.lat, self.s1.lat)
        assert np.may_share_memory(s1_diff.d_lon, self.diff.d_lon)

        # Only one check on T, since it just calls transpose anyway.
        # Doing it on the CartesianRepresentation just for variety's sake.
        c0_T = self.c0.T
        assert c0_T.shape == (7, 6)
        assert np.all(c0_T.x == self.c0.x.T)
        assert np.may_share_memory(c0_T.y, self.c0.y)

    def test_diagonal(self):
        s0_diagonal = self.s0.diagonal()
        s0_diff = s0_diagonal.differentials['s']
        assert s0_diagonal.shape == (6,)
        assert s0_diff.shape == s0_diagonal.shape
        assert np.all(s0_diagonal.lat == self.s0.lat.diagonal())
        assert np.all(s0_diff.d_lon == self.diff.d_lon.diagonal())
        assert np.may_share_memory(s0_diagonal.lat, self.s0.lat)
        assert np.may_share_memory(s0_diff.d_lon, self.diff.d_lon)

    def test_swapaxes(self):
        s1_swapaxes = self.s1.swapaxes(0, 1)
        s1_diff = s1_swapaxes.differentials['s']
        assert s1_swapaxes.shape == (7, 6)
        assert s1_diff.shape == s1_swapaxes.shape
        assert np.all(s1_swapaxes.lat == self.s1.lat.swapaxes(0, 1))
        assert np.all(s1_diff.d_lon == self.diff.d_lon.swapaxes(0, 1))
        assert np.may_share_memory(s1_swapaxes.lat, self.s1.lat)
        assert np.may_share_memory(s1_diff.d_lon, self.diff.d_lon)

    def test_reshape(self):
        s0_reshape = self.s0.reshape(2, 3, 7)
        s0_diff = s0_reshape.differentials['s']
        assert s0_reshape.shape == (2, 3, 7)
        assert s0_diff.shape == s0_reshape.shape
        assert np.all(s0_reshape.lon == self.s0.lon.reshape(2, 3, 7))
        assert np.all(s0_reshape.lat == self.s0.lat.reshape(2, 3, 7))
        assert np.all(s0_reshape.distance == self.s0.distance.reshape(2, 3, 7))
        assert np.may_share_memory(s0_reshape.lon, self.s0.lon)
        assert np.may_share_memory(s0_reshape.lat, self.s0.lat)
        assert np.may_share_memory(s0_reshape.distance, self.s0.distance)

        s1_reshape = self.s1.reshape(3, 2, 7)
        s1_diff = s1_reshape.differentials['s']
        assert s1_reshape.shape == (3, 2, 7)
        assert s1_diff.shape == s1_reshape.shape
        assert np.all(s1_reshape.lat == self.s1.lat.reshape(3, 2, 7))
        assert np.all(s1_diff.d_lon == self.diff.d_lon.reshape(3, 2, 7))
        assert np.may_share_memory(s1_reshape.lat, self.s1.lat)
        assert np.may_share_memory(s1_diff.d_lon, self.diff.d_lon)

        # For reshape(3, 14), copying is necessary for lon, lat, but not for d
        s1_reshape2 = self.s1.reshape(3, 14)
        assert s1_reshape2.shape == (3, 14)
        assert np.all(s1_reshape2.lon == self.s1.lon.reshape(3, 14))
        assert not np.may_share_memory(s1_reshape2.lon, self.s1.lon)
        assert s1_reshape2.distance.shape == (3, 14)
        assert np.may_share_memory(s1_reshape2.distance, self.s1.distance)

    def test_shape_setting(self):
        # Shape-setting should be on the object itself, since copying removes
        # zero-strides due to broadcasting.  We reset the objects at the end.
        self.s0.shape = (2, 3, 7)
        assert self.s0.shape == (2, 3, 7)
        assert self.s0.lon.shape == (2, 3, 7)
        assert self.s0.lat.shape == (2, 3, 7)
        assert self.s0.distance.shape == (2, 3, 7)
        assert self.diff.shape == (2, 3, 7)
        assert self.diff.d_lon.shape == (2, 3, 7)
        assert self.diff.d_lat.shape == (2, 3, 7)
        assert self.diff.d_distance.shape == (2, 3, 7)

        # this works with the broadcasting.
        self.s1.shape = (2, 3, 7)
        assert self.s1.shape == (2, 3, 7)
        assert self.s1.lon.shape == (2, 3, 7)
        assert self.s1.lat.shape == (2, 3, 7)
        assert self.s1.distance.shape == (2, 3, 7)
        assert self.s1.distance.strides == (0, 0, 0)

        # but this one does not.
        oldshape = self.s1.shape
        with pytest.raises(AttributeError):
            self.s1.shape = (42,)
        assert self.s1.shape == oldshape
        assert self.s1.lon.shape == oldshape
        assert self.s1.lat.shape == oldshape
        assert self.s1.distance.shape == oldshape

        # Finally, a more complicated one that checks that things get reset
        # properly if it is not the first component that fails.
        s2 = SphericalRepresentation(self.s1.lon.copy(), self.s1.lat,
                                     self.s1.distance, copy=False)
        assert 0 not in s2.lon.strides
        assert 0 in s2.lat.strides
        with pytest.raises(AttributeError):
            s2.shape = (42,)
        assert s2.shape == oldshape
        assert s2.lon.shape == oldshape
        assert s2.lat.shape == oldshape
        assert s2.distance.shape == oldshape
        assert 0 not in s2.lon.strides
        assert 0 in s2.lat.strides
        self.setup()

    def test_squeeze(self):
        s0_squeeze = self.s0.reshape(3, 1, 2, 1, 7).squeeze()
        s0_diff = s0_squeeze.differentials['s']
        assert s0_squeeze.shape == (3, 2, 7)
        assert s0_diff.shape == s0_squeeze.shape
        assert np.all(s0_squeeze.lat == self.s0.lat.reshape(3, 2, 7))
        assert np.all(s0_diff.d_lon == self.diff.d_lon.reshape(3, 2, 7))
        assert np.may_share_memory(s0_squeeze.lat, self.s0.lat)

    def test_add_dimension(self):
        s0_adddim = self.s0[:, np.newaxis, :]
        s0_diff = s0_adddim.differentials['s']
        assert s0_adddim.shape == (6, 1, 7)
        assert s0_diff.shape == s0_adddim.shape
        assert np.all(s0_adddim.lon == self.s0.lon[:, np.newaxis, :])
        assert np.all(s0_diff.d_lon == self.diff.d_lon[:, np.newaxis, :])
        assert np.may_share_memory(s0_adddim.lat, self.s0.lat)

    def test_take(self):
        s0_take = self.s0.take((5, 2))
        s0_diff = s0_take.differentials['s']
        assert s0_take.shape == (2,)
        assert s0_diff.shape == s0_take.shape
        assert np.all(s0_take.lon == self.s0.lon.take((5, 2)))
        assert np.all(s0_diff.d_lon == self.diff.d_lon.take((5, 2)))

    def test_broadcast_to(self):
        s0_broadcast = self.s0._apply(np_broadcast_to, (3, 6, 7), subok=True)
        s0_diff = s0_broadcast.differentials['s']
        assert type(s0_broadcast) is type(self.s0)
        assert s0_broadcast.shape == (3, 6, 7)
        assert s0_diff.shape == s0_broadcast.shape
        assert np.all(s0_broadcast.lon == self.s0.lon)
        assert np.all(s0_broadcast.lat == self.s0.lat)
        assert np.all(s0_broadcast.distance == self.s0.distance)
        assert np.may_share_memory(s0_broadcast.lon, self.s0.lon)
        assert np.may_share_memory(s0_broadcast.lat, self.s0.lat)
        assert np.may_share_memory(s0_broadcast.distance, self.s0.distance)

        s1_broadcast = self.s1._apply(np_broadcast_to, shape=(3, 6, 7),
                                      subok=True)
        s1_diff = s1_broadcast.differentials['s']
        assert s1_broadcast.shape == (3, 6, 7)
        assert s1_diff.shape == s1_broadcast.shape
        assert np.all(s1_broadcast.lat == self.s1.lat)
        assert np.all(s1_broadcast.lon == self.s1.lon)
        assert np.all(s1_broadcast.distance == self.s1.distance)
        assert s1_broadcast.distance.shape == (3, 6, 7)
        assert np.may_share_memory(s1_broadcast.lat, self.s1.lat)
        assert np.may_share_memory(s1_broadcast.lon, self.s1.lon)
        assert np.may_share_memory(s1_broadcast.distance, self.s1.distance)

        # A final test that "may_share_memory" equals "does_share_memory"
        # Do this on a copy, to keep self.s0 unchanged.
        sc = self.s0.copy()
        assert not np.may_share_memory(sc.lon, self.s0.lon)
        assert not np.may_share_memory(sc.lat, self.s0.lat)
        sc_broadcast = sc._apply(np_broadcast_to, (3, 6, 7), subok=True)
        assert np.may_share_memory(sc_broadcast.lon, sc.lon)
        # Can only write to copy, not to broadcast version.
        sc.lon[0, 0] = 22. * u.hourangle
        assert np.all(sc_broadcast.lon[:, 0, 0] == 22. * u.hourangle)
