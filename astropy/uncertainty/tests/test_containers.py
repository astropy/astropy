# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test that Distribution works with classes other than ndarray and Quantity."""

import numpy as np
import pytest
from numpy.testing import assert_array_equal

import astropy.units as u
from astropy.coordinates import Angle, Latitude, Longitude
from astropy.uncertainty import Distribution


class TestAngles:
    @classmethod
    def setup_class(cls):
        cls.a = np.arange(27.0).reshape(3, 9)
        cls.d = Distribution(cls.a)
        cls.q = cls.a << u.deg
        cls.dq = Distribution(cls.q)

    @pytest.mark.parametrize("angle_cls", [Angle, Longitude, Latitude])
    def test_as_input_for_angle(self, angle_cls):
        da = angle_cls(self.dq)
        assert isinstance(da, angle_cls)
        assert isinstance(da, Distribution)
        assert_array_equal(da.distribution, angle_cls(self.q))

    @pytest.mark.parametrize("angle_cls", [Angle, Longitude, Latitude])
    def test_using_angle_as_input(self, angle_cls):
        a = angle_cls(self.q)
        da = Distribution(a)
        assert isinstance(da, angle_cls)
        assert isinstance(da, Distribution)

    # Parametrize the unit to check the various branches in Latitude._validate_angles
    @pytest.mark.parametrize("dtype", ["f8", "f4"])
    @pytest.mark.parametrize(
        "value", [90 * u.deg, np.pi / 2 * u.radian, 90 * 60 * u.arcmin]
    )
    def test_at_limit_for_latitude(self, value, dtype):
        q = u.Quantity(value, dtype=dtype).reshape(1)
        qd = Distribution(q)
        ld = Latitude(qd)
        assert_array_equal(ld.distribution, Latitude(q))

    # Parametrize the unit in case Longitude._wrap_at becomes unit-dependent.
    @pytest.mark.parametrize("dtype", ["f8", "f4"])
    @pytest.mark.parametrize(
        "value", [360 * u.deg, 2 * np.pi * u.radian, 360 * 60 * u.arcmin]
    )
    def test_at_wrap_angle_for_longitude(self, value, dtype):
        q = u.Quantity(value, dtype=dtype).reshape(1)
        qd = Distribution(q)
        ld = Longitude(qd)
        assert_array_equal(ld.distribution, Longitude(q))
        assert np.all(ld.distribution == 0)

    @pytest.mark.parametrize("angle_cls", [Longitude, Latitude])
    def test_operation_gives_correct_subclass(self, angle_cls):
        # Lon and Lat always fall back to Angle
        da = angle_cls(self.dq)
        da2 = da + da
        assert isinstance(da, Angle)
        assert isinstance(da, Distribution)

    @pytest.mark.parametrize("angle_cls", [Longitude, Latitude])
    def test_pdfstd_gives_correct_subclass(self, angle_cls):
        # Lon and Lat always fall back to Angle
        da = angle_cls(self.dq)
        std = da.pdf_std()
        assert isinstance(std, Angle)
        assert_array_equal(std, Angle(self.q.std(-1)))
