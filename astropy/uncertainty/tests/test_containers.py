# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test that Distribution works with classes other than ndarray and Quantity."""

import numpy as np
import pytest
from numpy.testing import assert_array_equal

import astropy.units as u
from astropy.coordinates import Angle, Latitude, Longitude
from astropy.uncertainty import Distribution
from astropy.utils import NumpyRNGContext


@pytest.mark.parametrize("angle_cls", [Angle, Longitude, Latitude])
class TestAngles:
    @classmethod
    def setup_class(cls):
        with NumpyRNGContext(12345):
            cls.a = np.random.normal(10, [[0.2], [1.5], [4], [1]], size=(4, 10))
        cls.d = Distribution(cls.a)
        cls.q = cls.a << u.deg
        cls.dq = Distribution(cls.q)

    def test_as_input_for_angle(self, angle_cls):
        da = angle_cls(self.dq)
        assert isinstance(da, angle_cls)
        assert isinstance(da, Distribution)
        assert_array_equal(da.distribution, angle_cls(self.q))

    def test_using_angle_as_input(self, angle_cls):
        a = angle_cls(self.q)
        da = Distribution(a)
        assert isinstance(da, angle_cls)
        assert isinstance(da, Distribution)

    def test_operation_gives_correct_subclass(self, angle_cls):
        # Lon and Lat always fall back to Angle
        da = angle_cls(self.dq)
        da2 = da + da
        assert isinstance(da, Angle)
        assert isinstance(da, Distribution)

    def test_pdfmean_gives_correct_subclass(self, angle_cls):
        # Lon and Lat always fall back to Angle
        da = angle_cls(self.dq)
        mean = da.pdf_mean()
        assert isinstance(mean, Angle)
        assert_array_equal(mean, Angle(self.q.mean(-1)))
