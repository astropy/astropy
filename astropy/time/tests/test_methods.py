# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy
import itertools
import warnings

import numpy as np
import pytest

import astropy.units as u
from astropy.time import Time
from astropy.time.utils import day_frac
from astropy.units.quantity_helper.function_helpers import ARRAY_FUNCTION_ENABLED
from astropy.utils import iers

needs_array_function = pytest.mark.xfail(
    not ARRAY_FUNCTION_ENABLED, reason="Needs __array_function__ support"
)


def assert_time_all_equal(t1, t2):
    """Checks equality of shape and content."""
    assert t1.shape == t2.shape
    assert np.all(t1 == t2)


class ShapeSetup:
    def setup_class(cls):
        mjd = np.arange(50000, 50010)
        frac = np.arange(0.0, 0.999, 0.2)
        frac_masked = np.ma.array(frac)
        frac_masked[1] = np.ma.masked

        cls.t0 = {
            "not_masked": Time(mjd[:, np.newaxis] + frac, format="mjd", scale="utc"),
            "masked": Time(mjd[:, np.newaxis] + frac_masked, format="mjd", scale="utc"),
        }
        cls.t1 = {
            "not_masked": Time(
                mjd[:, np.newaxis] + frac,
                format="mjd",
                scale="utc",
                location=("45d", "50d"),
            ),
            "masked": Time(
                mjd[:, np.newaxis] + frac_masked,
                format="mjd",
                scale="utc",
                location=("45d", "50d"),
            ),
        }
        cls.t2 = {
            "not_masked": Time(
                mjd[:, np.newaxis] + frac,
                format="mjd",
                scale="utc",
                location=(np.arange(len(frac)), np.arange(len(frac))),
            ),
            "masked": Time(
                mjd[:, np.newaxis] + frac_masked,
                format="mjd",
                scale="utc",
                location=(np.arange(len(frac_masked)), np.arange(len(frac_masked))),
            ),
        }

    def create_data(self, use_mask):
        self.t0 = self.__class__.t0[use_mask]
        self.t1 = self.__class__.t1[use_mask]
        self.t2 = self.__class__.t2[use_mask]


@pytest.mark.parametrize("use_mask", ("masked", "not_masked"))
class TestManipulation(ShapeSetup):
    """Manipulation of Time objects, ensuring attributes are done correctly."""

    def test_ravel(self, use_mask):
        self.create_data(use_mask)

        t0_ravel = self.t0.ravel()
        assert t0_ravel.shape == (self.t0.size,)
        assert np.all(t0_ravel.jd1 == self.t0.jd1.ravel())
        assert np.may_share_memory(t0_ravel.jd1, self.t0.jd1)
        assert t0_ravel.location is None
        t1_ravel = self.t1.ravel()
        assert t1_ravel.shape == (self.t1.size,)
        assert np.all(t1_ravel.jd1 == self.t1.jd1.ravel())
        assert np.may_share_memory(t1_ravel.jd1, self.t1.jd1)
        assert t1_ravel.location is self.t1.location
        t2_ravel = self.t2.ravel()
        assert t2_ravel.shape == (self.t2.size,)
        assert np.all(t2_ravel.jd1 == self.t2.jd1.ravel())
        assert np.may_share_memory(t2_ravel.jd1, self.t2.jd1)
        assert t2_ravel.location.shape == t2_ravel.shape
        # Broadcasting and ravelling cannot be done without a copy.
        assert not np.may_share_memory(t2_ravel.location, self.t2.location)

    def test_flatten(self, use_mask):
        self.create_data(use_mask)

        t0_flatten = self.t0.flatten()
        assert t0_flatten.shape == (self.t0.size,)
        assert t0_flatten.location is None
        # Flatten always makes a copy.
        assert not np.may_share_memory(t0_flatten.jd1, self.t0.jd1)
        t1_flatten = self.t1.flatten()
        assert t1_flatten.shape == (self.t1.size,)
        assert not np.may_share_memory(t1_flatten.jd1, self.t1.jd1)
        assert t1_flatten.location is not self.t1.location
        assert t1_flatten.location == self.t1.location
        t2_flatten = self.t2.flatten()
        assert t2_flatten.shape == (self.t2.size,)
        assert not np.may_share_memory(t2_flatten.jd1, self.t2.jd1)
        assert t2_flatten.location.shape == t2_flatten.shape
        assert not np.may_share_memory(t2_flatten.location, self.t2.location)

    def test_transpose(self, use_mask):
        self.create_data(use_mask)

        t0_transpose = self.t0.transpose()
        assert t0_transpose.shape == (5, 10)
        assert np.all(t0_transpose.jd1 == self.t0.jd1.transpose())
        assert np.may_share_memory(t0_transpose.jd1, self.t0.jd1)
        assert t0_transpose.location is None
        t1_transpose = self.t1.transpose()
        assert t1_transpose.shape == (5, 10)
        assert np.all(t1_transpose.jd1 == self.t1.jd1.transpose())
        assert np.may_share_memory(t1_transpose.jd1, self.t1.jd1)
        assert t1_transpose.location is self.t1.location
        t2_transpose = self.t2.transpose()
        assert t2_transpose.shape == (5, 10)
        assert np.all(t2_transpose.jd1 == self.t2.jd1.transpose())
        assert np.may_share_memory(t2_transpose.jd1, self.t2.jd1)
        assert t2_transpose.location.shape == t2_transpose.shape
        assert np.may_share_memory(t2_transpose.location, self.t2.location)
        # Only one check on T, since it just calls transpose anyway.
        t2_T = self.t2.T
        assert t2_T.shape == (5, 10)
        assert np.all(t2_T.jd1 == self.t2.jd1.T)
        assert np.may_share_memory(t2_T.jd1, self.t2.jd1)
        assert t2_T.location.shape == t2_T.location.shape
        assert np.may_share_memory(t2_T.location, self.t2.location)

    def test_diagonal(self, use_mask):
        self.create_data(use_mask)

        t0_diagonal = self.t0.diagonal()
        assert t0_diagonal.shape == (5,)
        assert np.all(t0_diagonal.jd1 == self.t0.jd1.diagonal())
        assert t0_diagonal.location is None
        assert np.may_share_memory(t0_diagonal.jd1, self.t0.jd1)
        t1_diagonal = self.t1.diagonal()
        assert t1_diagonal.shape == (5,)
        assert np.all(t1_diagonal.jd1 == self.t1.jd1.diagonal())
        assert t1_diagonal.location is self.t1.location
        assert np.may_share_memory(t1_diagonal.jd1, self.t1.jd1)
        t2_diagonal = self.t2.diagonal()
        assert t2_diagonal.shape == (5,)
        assert np.all(t2_diagonal.jd1 == self.t2.jd1.diagonal())
        assert t2_diagonal.location.shape == t2_diagonal.shape
        assert np.may_share_memory(t2_diagonal.jd1, self.t2.jd1)
        assert np.may_share_memory(t2_diagonal.location, self.t2.location)

    def test_swapaxes(self, use_mask):
        self.create_data(use_mask)

        t0_swapaxes = self.t0.swapaxes(0, 1)
        assert t0_swapaxes.shape == (5, 10)
        assert np.all(t0_swapaxes.jd1 == self.t0.jd1.swapaxes(0, 1))
        assert np.may_share_memory(t0_swapaxes.jd1, self.t0.jd1)
        assert t0_swapaxes.location is None
        t1_swapaxes = self.t1.swapaxes(0, 1)
        assert t1_swapaxes.shape == (5, 10)
        assert np.all(t1_swapaxes.jd1 == self.t1.jd1.swapaxes(0, 1))
        assert np.may_share_memory(t1_swapaxes.jd1, self.t1.jd1)
        assert t1_swapaxes.location is self.t1.location
        t2_swapaxes = self.t2.swapaxes(0, 1)
        assert t2_swapaxes.shape == (5, 10)
        assert np.all(t2_swapaxes.jd1 == self.t2.jd1.swapaxes(0, 1))
        assert np.may_share_memory(t2_swapaxes.jd1, self.t2.jd1)
        assert t2_swapaxes.location.shape == t2_swapaxes.shape
        assert np.may_share_memory(t2_swapaxes.location, self.t2.location)

    def test_reshape(self, use_mask):
        self.create_data(use_mask)

        t0_reshape = self.t0.reshape(5, 2, 5)
        assert t0_reshape.shape == (5, 2, 5)
        assert np.all(t0_reshape.jd1 == self.t0._time.jd1.reshape(5, 2, 5))
        assert np.all(t0_reshape.jd2 == self.t0._time.jd2.reshape(5, 2, 5))
        assert np.may_share_memory(t0_reshape.jd1, self.t0.jd1)
        assert np.may_share_memory(t0_reshape.jd2, self.t0.jd2)
        assert t0_reshape.location is None
        t1_reshape = self.t1.reshape(2, 5, 5)
        assert t1_reshape.shape == (2, 5, 5)
        assert np.all(t1_reshape.jd1 == self.t1.jd1.reshape(2, 5, 5))
        assert np.may_share_memory(t1_reshape.jd1, self.t1.jd1)
        assert t1_reshape.location is self.t1.location
        # For reshape(5, 2, 5), the location array can remain the same.
        t2_reshape = self.t2.reshape(5, 2, 5)
        assert t2_reshape.shape == (5, 2, 5)
        assert np.all(t2_reshape.jd1 == self.t2.jd1.reshape(5, 2, 5))
        assert np.may_share_memory(t2_reshape.jd1, self.t2.jd1)
        assert t2_reshape.location.shape == t2_reshape.shape
        assert np.may_share_memory(t2_reshape.location, self.t2.location)
        # But for reshape(5, 5, 2), location has to be broadcast and copied.
        t2_reshape2 = self.t2.reshape(5, 5, 2)
        assert t2_reshape2.shape == (5, 5, 2)
        assert np.all(t2_reshape2.jd1 == self.t2.jd1.reshape(5, 5, 2))
        assert np.may_share_memory(t2_reshape2.jd1, self.t2.jd1)
        assert t2_reshape2.location.shape == t2_reshape2.shape
        assert not np.may_share_memory(t2_reshape2.location, self.t2.location)
        t2_reshape_t = self.t2.reshape(10, 5).T
        assert t2_reshape_t.shape == (5, 10)
        assert np.may_share_memory(t2_reshape_t.jd1, self.t2.jd1)
        assert t2_reshape_t.location.shape == t2_reshape_t.shape
        assert np.may_share_memory(t2_reshape_t.location, self.t2.location)
        # Finally, reshape in a way that cannot be a view.
        t2_reshape_t_reshape = t2_reshape_t.reshape(10, 5)
        assert t2_reshape_t_reshape.shape == (10, 5)
        assert not np.may_share_memory(t2_reshape_t_reshape.jd1, self.t2.jd1)
        assert t2_reshape_t_reshape.location.shape == t2_reshape_t_reshape.shape
        assert not np.may_share_memory(
            t2_reshape_t_reshape.location, t2_reshape_t.location
        )

    def test_squeeze(self, use_mask):
        self.create_data(use_mask)

        t0_squeeze = self.t0.reshape(5, 1, 2, 1, 5).squeeze()
        assert t0_squeeze.shape == (5, 2, 5)
        assert np.all(t0_squeeze.jd1 == self.t0.jd1.reshape(5, 2, 5))
        assert np.may_share_memory(t0_squeeze.jd1, self.t0.jd1)
        assert t0_squeeze.location is None
        t1_squeeze = self.t1.reshape(1, 5, 1, 2, 5).squeeze()
        assert t1_squeeze.shape == (5, 2, 5)
        assert np.all(t1_squeeze.jd1 == self.t1.jd1.reshape(5, 2, 5))
        assert np.may_share_memory(t1_squeeze.jd1, self.t1.jd1)
        assert t1_squeeze.location is self.t1.location
        t2_squeeze = self.t2.reshape(1, 1, 5, 2, 5, 1, 1).squeeze()
        assert t2_squeeze.shape == (5, 2, 5)
        assert np.all(t2_squeeze.jd1 == self.t2.jd1.reshape(5, 2, 5))
        assert np.may_share_memory(t2_squeeze.jd1, self.t2.jd1)
        assert t2_squeeze.location.shape == t2_squeeze.shape
        assert np.may_share_memory(t2_squeeze.location, self.t2.location)

    def test_add_dimension(self, use_mask):
        self.create_data(use_mask)

        t0_adddim = self.t0[:, np.newaxis, :]
        assert t0_adddim.shape == (10, 1, 5)
        assert np.all(t0_adddim.jd1 == self.t0.jd1[:, np.newaxis, :])
        assert np.may_share_memory(t0_adddim.jd1, self.t0.jd1)
        assert t0_adddim.location is None
        t1_adddim = self.t1[:, :, np.newaxis]
        assert t1_adddim.shape == (10, 5, 1)
        assert np.all(t1_adddim.jd1 == self.t1.jd1[:, :, np.newaxis])
        assert np.may_share_memory(t1_adddim.jd1, self.t1.jd1)
        assert t1_adddim.location is self.t1.location
        t2_adddim = self.t2[:, :, np.newaxis]
        assert t2_adddim.shape == (10, 5, 1)
        assert np.all(t2_adddim.jd1 == self.t2.jd1[:, :, np.newaxis])
        assert np.may_share_memory(t2_adddim.jd1, self.t2.jd1)
        assert t2_adddim.location.shape == t2_adddim.shape
        assert np.may_share_memory(t2_adddim.location, self.t2.location)

    def test_take(self, use_mask):
        self.create_data(use_mask)

        t0_take = self.t0.take((5, 2))
        assert t0_take.shape == (2,)
        assert np.all(t0_take.jd1 == self.t0._time.jd1.take((5, 2)))
        assert t0_take.location is None
        t1_take = self.t1.take((2, 4), axis=1)
        assert t1_take.shape == (10, 2)
        assert np.all(t1_take.jd1 == self.t1.jd1.take((2, 4), axis=1))
        assert t1_take.location is self.t1.location
        t2_take = self.t2.take((1, 3, 7), axis=0)
        assert t2_take.shape == (3, 5)
        assert np.all(t2_take.jd1 == self.t2.jd1.take((1, 3, 7), axis=0))
        assert t2_take.location.shape == t2_take.shape
        t2_take2 = self.t2.take((5, 15))
        assert t2_take2.shape == (2,)
        assert np.all(t2_take2.jd1 == self.t2.jd1.take((5, 15)))
        assert t2_take2.location.shape == t2_take2.shape

    def test_broadcast_via_apply(self, use_mask):
        """Test using a callable method."""
        self.create_data(use_mask)

        t0_broadcast = self.t0._apply(np.broadcast_to, shape=(3, 10, 5))
        assert t0_broadcast.shape == (3, 10, 5)
        assert np.all(t0_broadcast.jd1 == self.t0.jd1)
        assert np.may_share_memory(t0_broadcast.jd1, self.t0.jd1)
        assert t0_broadcast.location is None
        t1_broadcast = self.t1._apply(np.broadcast_to, shape=(3, 10, 5))
        assert t1_broadcast.shape == (3, 10, 5)
        assert np.all(t1_broadcast.jd1 == self.t1.jd1)
        assert np.may_share_memory(t1_broadcast.jd1, self.t1.jd1)
        assert t1_broadcast.location is self.t1.location
        t2_broadcast = self.t2._apply(np.broadcast_to, shape=(3, 10, 5))
        assert t2_broadcast.shape == (3, 10, 5)
        assert np.all(t2_broadcast.jd1 == self.t2.jd1)
        assert np.may_share_memory(t2_broadcast.jd1, self.t2.jd1)
        assert t2_broadcast.location.shape == t2_broadcast.shape
        assert np.may_share_memory(t2_broadcast.location, self.t2.location)


@pytest.mark.parametrize("use_mask", ("masked", "not_masked"))
class TestSetShape(ShapeSetup):
    def test_shape_setting(self, use_mask):
        # Shape-setting should be on the object itself, since copying removes
        # zero-strides due to broadcasting.  Hence, this should be the only
        # test in this class.
        self.create_data(use_mask)

        t0_reshape = self.t0.copy()
        mjd = t0_reshape.mjd  # Creates a cache of the mjd attribute
        t0_reshape.shape = (5, 2, 5)
        assert t0_reshape.shape == (5, 2, 5)
        assert mjd.shape != t0_reshape.mjd.shape  # Cache got cleared
        assert np.all(t0_reshape.jd1 == self.t0._time.jd1.reshape(5, 2, 5))
        assert np.all(t0_reshape.jd2 == self.t0._time.jd2.reshape(5, 2, 5))
        assert t0_reshape.location is None
        # But if the shape doesn't work, one should get an error.
        t0_reshape_t = t0_reshape.T
        with pytest.raises(ValueError):
            t0_reshape_t.shape = (12,)  # Wrong number of elements.
        with pytest.raises(AttributeError):
            t0_reshape_t.shape = (10, 5)  # Cannot be done without copy.
        # check no shape was changed.
        assert t0_reshape_t.shape == t0_reshape.T.shape
        assert t0_reshape_t.jd1.shape == t0_reshape.T.shape
        assert t0_reshape_t.jd2.shape == t0_reshape.T.shape
        t1_reshape = self.t1.copy()
        t1_reshape.shape = (2, 5, 5)
        assert t1_reshape.shape == (2, 5, 5)
        assert np.all(t1_reshape.jd1 == self.t1.jd1.reshape(2, 5, 5))
        # location is a single element, so its shape should not change.
        assert t1_reshape.location.shape == ()
        # For reshape(5, 2, 5), the location array can remain the same.
        # Note that we need to work directly on self.t2 here, since any
        # copy would cause location to have the full shape.
        self.t2.shape = (5, 2, 5)
        assert self.t2.shape == (5, 2, 5)
        assert self.t2.jd1.shape == (5, 2, 5)
        assert self.t2.jd2.shape == (5, 2, 5)
        assert self.t2.location.shape == (5, 2, 5)
        assert self.t2.location.strides == (0, 0, 24)
        # But for reshape(50), location would need to be copied, so this
        # should fail.
        oldshape = self.t2.shape
        with pytest.raises(AttributeError):
            self.t2.shape = (50,)
        # check no shape was changed.
        assert self.t2.jd1.shape == oldshape
        assert self.t2.jd2.shape == oldshape
        assert self.t2.location.shape == oldshape


@pytest.mark.parametrize("use_mask", ("masked", "not_masked"))
class TestShapeFunctions(ShapeSetup):
    @needs_array_function
    def test_broadcast(self, use_mask):
        """Test as supported numpy function."""
        self.create_data(use_mask)

        t0_broadcast = np.broadcast_to(self.t0, shape=(3, 10, 5))
        assert t0_broadcast.shape == (3, 10, 5)
        assert np.all(t0_broadcast.jd1 == self.t0.jd1)
        assert np.may_share_memory(t0_broadcast.jd1, self.t0.jd1)
        assert t0_broadcast.location is None
        t1_broadcast = np.broadcast_to(self.t1, shape=(3, 10, 5))
        assert t1_broadcast.shape == (3, 10, 5)
        assert np.all(t1_broadcast.jd1 == self.t1.jd1)
        assert np.may_share_memory(t1_broadcast.jd1, self.t1.jd1)
        assert t1_broadcast.location is self.t1.location
        t2_broadcast = np.broadcast_to(self.t2, shape=(3, 10, 5))
        assert t2_broadcast.shape == (3, 10, 5)
        assert np.all(t2_broadcast.jd1 == self.t2.jd1)
        assert np.may_share_memory(t2_broadcast.jd1, self.t2.jd1)
        assert t2_broadcast.location.shape == t2_broadcast.shape
        assert np.may_share_memory(t2_broadcast.location, self.t2.location)

    @needs_array_function
    def test_atleast_1d(self, use_mask):
        self.create_data(use_mask)

        t00 = self.t0.ravel()[0]
        assert t00.ndim == 0
        t00_1d = np.atleast_1d(t00)
        assert t00_1d.ndim == 1
        assert_time_all_equal(t00[np.newaxis], t00_1d)
        # Actual jd1 will not share memory, as cast to scalar.
        assert np.may_share_memory(t00_1d._time.jd1, t00._time.jd1)

    @needs_array_function
    def test_atleast_2d(self, use_mask):
        self.create_data(use_mask)

        t0r = self.t0.ravel()
        assert t0r.ndim == 1
        t0r_2d = np.atleast_2d(t0r)
        assert t0r_2d.ndim == 2
        assert_time_all_equal(t0r[np.newaxis], t0r_2d)
        assert np.may_share_memory(t0r_2d.jd1, t0r.jd1)

    @needs_array_function
    def test_atleast_3d(self, use_mask):
        self.create_data(use_mask)

        assert self.t0.ndim == 2
        t0_3d, t1_3d = np.atleast_3d(self.t0, self.t1)
        assert t0_3d.ndim == t1_3d.ndim == 3
        assert_time_all_equal(self.t0[:, :, np.newaxis], t0_3d)
        assert_time_all_equal(self.t1[:, :, np.newaxis], t1_3d)
        assert np.may_share_memory(t0_3d.jd2, self.t0.jd2)

    def test_move_axis(self, use_mask):
        # Goes via transpose so works without __array_function__ as well.
        self.create_data(use_mask)

        t0_10 = np.moveaxis(self.t0, 0, 1)
        assert t0_10.shape == (self.t0.shape[1], self.t0.shape[0])
        assert_time_all_equal(self.t0.T, t0_10)
        assert np.may_share_memory(t0_10.jd1, self.t0.jd1)

    def test_roll_axis(self, use_mask):
        # Goes via transpose so works without __array_function__ as well.
        self.create_data(use_mask)

        t0_10 = np.rollaxis(self.t0, 1)
        assert t0_10.shape == (self.t0.shape[1], self.t0.shape[0])
        assert_time_all_equal(self.t0.T, t0_10)
        assert np.may_share_memory(t0_10.jd1, self.t0.jd1)

    @needs_array_function
    def test_fliplr(self, use_mask):
        self.create_data(use_mask)

        t0_lr = np.fliplr(self.t0)
        assert_time_all_equal(self.t0[:, ::-1], t0_lr)
        assert np.may_share_memory(t0_lr.jd2, self.t0.jd2)

    @needs_array_function
    def test_rot90(self, use_mask):
        self.create_data(use_mask)

        t0_270 = np.rot90(self.t0, 3)
        assert_time_all_equal(self.t0.T[:, ::-1], t0_270)
        assert np.may_share_memory(t0_270.jd2, self.t0.jd2)

    @needs_array_function
    def test_roll(self, use_mask):
        self.create_data(use_mask)

        t0r = np.roll(self.t0, 1, axis=0)
        assert_time_all_equal(t0r[1:], self.t0[:-1])
        assert_time_all_equal(t0r[0], self.t0[-1])

    @needs_array_function
    def test_delete(self, use_mask):
        self.create_data(use_mask)

        t0d = np.delete(self.t0, [2, 3], axis=0)
        assert_time_all_equal(t0d[:2], self.t0[:2])
        assert_time_all_equal(t0d[2:], self.t0[4:])


@pytest.mark.parametrize("use_mask", ("masked", "not_masked"))
class TestArithmetic:
    """Arithmetic on Time objects, using both doubles."""

    kwargs = ({}, {"axis": None}, {"axis": 0}, {"axis": 1}, {"axis": 2})
    functions = ("min", "max", "sort")

    def setup_class(cls):
        mjd = np.arange(50000, 50100, 10).reshape(2, 5, 1)
        frac = np.array([0.1, 0.1 + 1.0e-15, 0.1 - 1.0e-15, 0.9 + 2.0e-16, 0.9])
        frac_masked = np.ma.array(frac)
        frac_masked[1] = np.ma.masked

        cls.t0 = {
            "not_masked": Time(mjd, frac, format="mjd", scale="utc"),
            "masked": Time(mjd, frac_masked, format="mjd", scale="utc"),
        }

        # Define arrays with same ordinal properties
        frac = np.array([1, 2, 0, 4, 3])
        frac_masked = np.ma.array(frac)
        frac_masked[1] = np.ma.masked

        cls.t1 = {
            "not_masked": Time(mjd + frac, format="mjd", scale="utc"),
            "masked": Time(mjd + frac_masked, format="mjd", scale="utc"),
        }
        cls.jd = {"not_masked": mjd + frac, "masked": mjd + frac_masked}

        cls.t2 = {
            "not_masked": Time(
                mjd + frac,
                format="mjd",
                scale="utc",
                location=(np.arange(len(frac)), np.arange(len(frac))),
            ),
            "masked": Time(
                mjd + frac_masked,
                format="mjd",
                scale="utc",
                location=(np.arange(len(frac_masked)), np.arange(len(frac_masked))),
            ),
        }

    def create_data(self, use_mask):
        self.t0 = self.__class__.t0[use_mask]
        self.t1 = self.__class__.t1[use_mask]
        self.t2 = self.__class__.t2[use_mask]
        self.jd = self.__class__.jd[use_mask]

    @pytest.mark.parametrize("kw, func", itertools.product(kwargs, functions))
    def test_argfuncs(self, kw, func, use_mask):
        """
        Test that ``np.argfunc(jd, **kw)`` is the same as ``t0.argfunc(**kw)``
        where ``jd`` is a similarly shaped array with the same ordinal properties
        but all integer values.  Also test the same for t1 which has the same
        integral values as jd.
        """
        self.create_data(use_mask)

        t0v = getattr(self.t0, "arg" + func)(**kw)
        t1v = getattr(self.t1, "arg" + func)(**kw)
        jdv = getattr(np, "arg" + func)(self.jd, **kw)

        if self.t0.masked and kw == {"axis": None} and func == "sort":
            t0v = np.ma.array(t0v, mask=self.t0.mask.reshape(t0v.shape)[t0v])
            t1v = np.ma.array(t1v, mask=self.t1.mask.reshape(t1v.shape)[t1v])
            jdv = np.ma.array(jdv, mask=self.jd.mask.reshape(jdv.shape)[jdv])

        assert np.all(t0v == jdv)
        assert np.all(t1v == jdv)
        assert t0v.shape == jdv.shape
        assert t1v.shape == jdv.shape

    @pytest.mark.parametrize("kw, func", itertools.product(kwargs, functions))
    def test_funcs(self, kw, func, use_mask):
        """
        Test that ``np.func(jd, **kw)`` is the same as ``t1.func(**kw)`` where
        ``jd`` is a similarly shaped array and the same integral values.
        """
        self.create_data(use_mask)

        t1v = getattr(self.t1, func)(**kw)
        jdv = getattr(np, func)(self.jd, **kw)
        assert np.all(t1v.value == jdv)
        assert t1v.shape == jdv.shape

    def test_argmin(self, use_mask):
        self.create_data(use_mask)

        assert self.t0.argmin() == 2
        assert np.all(self.t0.argmin(axis=0) == 0)
        assert np.all(self.t0.argmin(axis=1) == 0)
        assert np.all(self.t0.argmin(axis=2) == 2)

    def test_argmax(self, use_mask):
        self.create_data(use_mask)

        assert self.t0.argmax() == self.t0.size - 2
        if use_mask == "masked":
            # The 0 is where all entries are masked in that axis
            assert np.all(self.t0.argmax(axis=0) == [1, 0, 1, 1, 1])
            assert np.all(self.t0.argmax(axis=1) == [4, 0, 4, 4, 4])
        else:
            assert np.all(self.t0.argmax(axis=0) == 1)
            assert np.all(self.t0.argmax(axis=1) == 4)
        assert np.all(self.t0.argmax(axis=2) == 3)

    def test_argsort(self, use_mask):
        self.create_data(use_mask)

        order = [2, 0, 4, 3, 1] if use_mask == "masked" else [2, 0, 1, 4, 3]
        assert np.all(self.t0.argsort() == np.array(order))
        assert np.all(self.t0.argsort(axis=0) == np.arange(2).reshape(2, 1, 1))
        assert np.all(self.t0.argsort(axis=1) == np.arange(5).reshape(5, 1))
        assert np.all(self.t0.argsort(axis=2) == np.array(order))
        ravel = np.arange(50).reshape(-1, 5)[:, order].ravel()
        if use_mask == "masked":
            t0v = self.t0.argsort(axis=None)
            # Manually remove elements in ravel that correspond to masked
            # entries in self.t0.  This removes the 10 entries that are masked
            # which show up at the end of the list.
            mask = self.t0.mask.ravel()[ravel]
            ravel = ravel[~mask]
            assert np.all(t0v[:-10] == ravel)
        else:
            assert np.all(self.t0.argsort(axis=None) == ravel)

    @pytest.mark.parametrize("scale", Time.SCALES)
    def test_argsort_warning(self, use_mask, scale):
        self.create_data(use_mask)

        if scale == "utc":
            pytest.xfail()
        with warnings.catch_warnings(record=True) as wlist:
            Time([1, 2, 3], format="jd", scale=scale).argsort()
        assert len(wlist) == 0

    def test_min(self, use_mask):
        self.create_data(use_mask)

        assert self.t0.min() == self.t0[0, 0, 2]
        assert np.all(self.t0.min(0) == self.t0[0])
        assert np.all(self.t0.min(1) == self.t0[:, 0])
        assert np.all(self.t0.min(2) == self.t0[:, :, 2])
        assert self.t0.min(0).shape == (5, 5)
        assert self.t0.min(0, keepdims=True).shape == (1, 5, 5)
        assert self.t0.min(1).shape == (2, 5)
        assert self.t0.min(1, keepdims=True).shape == (2, 1, 5)
        assert self.t0.min(2).shape == (2, 5)
        assert self.t0.min(2, keepdims=True).shape == (2, 5, 1)

    def test_max(self, use_mask):
        self.create_data(use_mask)

        assert self.t0.max() == self.t0[-1, -1, -2]
        assert np.all(self.t0.max(0) == self.t0[1])
        assert np.all(self.t0.max(1) == self.t0[:, 4])
        assert np.all(self.t0.max(2) == self.t0[:, :, 3])
        assert self.t0.max(0).shape == (5, 5)
        assert self.t0.max(0, keepdims=True).shape == (1, 5, 5)

    def test_ptp(self, use_mask):
        self.create_data(use_mask)

        assert self.t0.ptp() == self.t0.max() - self.t0.min()
        assert np.all(self.t0.ptp(0) == self.t0.max(0) - self.t0.min(0))
        assert self.t0.ptp(0).shape == (5, 5)
        assert self.t0.ptp(0, keepdims=True).shape == (1, 5, 5)

    def test_sort(self, use_mask):
        self.create_data(use_mask)

        order = [2, 0, 4, 3, 1] if use_mask == "masked" else [2, 0, 1, 4, 3]
        assert np.all(self.t0.sort() == self.t0[:, :, order])
        assert np.all(self.t0.sort(0) == self.t0)
        assert np.all(self.t0.sort(1) == self.t0)
        assert np.all(self.t0.sort(2) == self.t0[:, :, order])
        if use_mask == "not_masked":
            assert np.all(self.t0.sort(None) == self.t0[:, :, order].ravel())
            # Bit superfluous, but good to check.
            assert np.all(self.t0.sort(-1)[:, :, 0] == self.t0.min(-1))
            assert np.all(self.t0.sort(-1)[:, :, -1] == self.t0.max(-1))

    @pytest.mark.parametrize("axis", [None, 0, 1, 2, (0, 1)])
    @pytest.mark.parametrize(
        "where", [True, np.array([True, False, True, True, False])[..., np.newaxis]]
    )
    @pytest.mark.parametrize("keepdims", [False, True])
    def test_mean(self, use_mask, axis, where, keepdims):
        self.create_data(use_mask)

        kwargs = dict(axis=axis, where=where, keepdims=keepdims)

        def is_consistent(time):
            where_expected = where & ~time.mask
            where_expected = np.broadcast_to(where_expected, time.shape)

            kw = kwargs.copy()
            kw["where"] = where_expected

            divisor = where_expected.sum(axis=axis, keepdims=keepdims)

            if np.any(divisor == 0):
                with pytest.raises(ValueError):
                    time.mean(**kwargs)

            else:
                time_mean = time.mean(**kwargs)
                time_expected = Time(
                    *day_frac(
                        val1=np.ma.getdata(time.tai.jd1).sum(**kw),
                        val2=np.ma.getdata(time.tai.jd2).sum(**kw),
                        divisor=divisor,
                    ),
                    format="jd",
                    scale="tai",
                )
                time_expected._set_scale(time.scale)
                assert np.all(time_mean == time_expected)

        is_consistent(self.t0)
        is_consistent(self.t1)

        axes_location_not_constant = [None, 2]
        if axis in axes_location_not_constant:
            with pytest.raises(ValueError):
                self.t2.mean(**kwargs)
        else:
            is_consistent(self.t2)

    def test_mean_precision(self, use_mask):
        scale = "tai"
        epsilon = 1 * u.ns

        t0 = Time("2021-07-27T00:00:00", scale=scale)
        t1 = Time("2022-07-27T00:00:00", scale=scale)
        t2 = Time("2023-07-27T00:00:00", scale=scale)

        t = Time([t0, t2 + epsilon])

        if use_mask == "masked":
            t[0] = np.ma.masked
            assert t.mean() == (t2 + epsilon)

        else:
            assert t.mean() == (t1 + epsilon / 2)

    def test_mean_dtype(self, use_mask):
        self.create_data(use_mask)
        with pytest.raises(ValueError):
            self.t0.mean(dtype=int)

    def test_mean_out(self, use_mask):
        self.create_data(use_mask)
        with pytest.raises(ValueError):
            self.t0.mean(out=Time(np.zeros_like(self.t0.jd1), format="jd"))

    def test_mean_leap_second(self, use_mask):
        # Check that leap second is dealt with correctly: for UTC, across a leap
        # second boundary, one cannot just average jd, but has to go through TAI.
        if use_mask == "not_masked":
            t = Time(["2012-06-30 23:59:60.000", "2012-07-01 00:00:01.000"])
            mean_expected = t[0] + (t[1] - t[0]) / 2
            mean_expected_explicit = Time("2012-07-01 00:00:00")
            mean_test = t.mean()
            assert mean_expected == mean_expected_explicit
            assert mean_expected == mean_test
            assert mean_test != Time(
                *day_frac(t.jd1.sum(), t.jd2.sum(), divisor=2), format="jd"
            )


def test_regression():
    # For #5225, where a time with a single-element delta_ut1_utc could not
    # be copied, flattened, or ravelled. (For copy, it is in test_basic.)
    with iers.conf.set_temp("auto_download", False):
        t = Time(49580.0, scale="tai", format="mjd")
        t_ut1 = t.ut1
        t_ut1_copy = copy.deepcopy(t_ut1)
        assert type(t_ut1_copy.delta_ut1_utc) is np.ndarray
        t_ut1_flatten = t_ut1.flatten()
        assert type(t_ut1_flatten.delta_ut1_utc) is np.ndarray
        t_ut1_ravel = t_ut1.ravel()
        assert type(t_ut1_ravel.delta_ut1_utc) is np.ndarray
        assert t_ut1_copy.delta_ut1_utc == t_ut1.delta_ut1_utc
