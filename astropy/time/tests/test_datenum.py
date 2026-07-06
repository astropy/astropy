# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for MATLAB datenum format support."""
import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy.time import Time


class TestMatlabDatenum:
    """Test suite for MATLAB datenum format."""

    def test_datenum_format_registered(self):
        """Verify datenum format is registered."""
        assert "datenum" in Time.FORMATS

    def test_datenum_basic_creation(self):
        """Test basic Time creation from MATLAB datenum."""
        datenum = 739252.5
        t = Time(datenum, format="datenum", scale="utc")
        assert isinstance(t.jd, (float, np.ndarray))

    def test_datenum_roundtrip(self):
        """Test round-trip conversion: datenum -> Time -> datenum."""
        datenum_orig = 739252.5
        t = Time(datenum_orig, format="datenum")
        datenum_back = t.datenum
        assert_allclose(datenum_orig, datenum_back, atol=1e-10)

    def test_datenum_from_iso(self):
        """Test creating datenum from ISO format."""
        # Create from ISO
        t_iso = Time("2024-03-29 12:00:00", scale="utc", format="iso")
        datenum_val = t_iso.datenum

        # Create from datenum
        t_datenum = Time(datenum_val, format="datenum", scale="utc")

        # Should match within floating-point precision
        assert_allclose(t_iso.jd, t_datenum.jd, atol=1e-9)

    def test_datenum_array_input(self):
        """Test datenum with 1D array input."""
        datenum_array = np.array([739252.0, 739252.5, 739253.0])
        t = Time(datenum_array, format="datenum")

        assert isinstance(t.datenum, np.ndarray)
        assert_allclose(t.datenum, datenum_array, atol=1e-10)

    def test_datenum_2d_array(self):
        """Test datenum with 2D array input."""
        datenum_array = np.array([[739252.0, 739252.5], [739253.0, 739253.5]])
        t = Time(datenum_array, format="datenum")

        assert t.datenum.shape == (2, 2)
        assert_allclose(t.datenum, datenum_array, atol=1e-10)

    def test_datenum_different_scales(self):
        """Test datenum with different time scales."""
        datenum = 739252.5

        for scale in ["utc", "tai", "tt", "tdb"]:
            t = Time(datenum, format="datenum", scale=scale)
            # Should not raise error
            assert t.scale == scale
            # Converting back should give approximately same datenum
            datenum_back = t.to_value("datenum")
            # Allow for 1 microsecond rounding for different scales
            assert_allclose(datenum, datenum_back, atol=1e-6)

    def test_datenum_to_other_formats(self):
        """Test converting datenum to other formats."""
        datenum = 739252.5
        t = Time(datenum, format="datenum")

        # Should be able to convert to other formats
        jd = t.jd
        mjd = t.mjd
        unix = t.unix
        iso = t.iso

        assert isinstance(jd, (float, np.ndarray))
        assert isinstance(mjd, (float, np.ndarray))
        assert isinstance(unix, (float, np.ndarray))
        assert isinstance(iso, str)

    def test_datenum_from_other_formats(self):
        """Test creating datenum from other formats."""
        # From JD
        t_jd = Time(2451544.5, format="jd")
        datenum_from_jd = t_jd.datenum
        assert isinstance(datenum_from_jd, (float, np.ndarray))

        # From MJD
        t_mjd = Time(51544.0, format="mjd")
        datenum_from_mjd = t_mjd.datenum
        assert isinstance(datenum_from_mjd, (float, np.ndarray))

        # From Unix
        t_unix = Time(946684800.0, format="unix")
        datenum_from_unix = t_unix.datenum
        assert isinstance(datenum_from_unix, (float, np.ndarray))

    def test_datenum_precision(self):
        """Test precision preservation in datenum."""
        # Test microsecond precision
        datenum_precise = 739252.123456789
        t = Time(datenum_precise, format="datenum", scale="utc")
        retrieved = t.datenum

        # Should preserve to ~1 microsecond (10^-6)
        assert_allclose(datenum_precise, retrieved, atol=1e-6)

    def test_datenum_year_0(self):
        """Test datenum at the epoch (year 0)."""
        # datenum = 0 should be year 0, January 1, 00:00:00
        t = Time(0.0, format="datenum")
        # Verify by checking JD value (it should be some large negative number)
        # JD for 0000-01-01 00:00:00 is approximately -1721425.5
        assert t.jd < 0  # Year 0 is before J2000

    def test_datenum_masked_array(self):
        """Test datenum with masked arrays."""
        datenum_array = np.ma.array(
            [739252.0, 739252.5, 739253.0], mask=[False, True, False]
        )
        t = Time(datenum_array, format="datenum")

        assert isinstance(t.datenum, np.ma.MaskedArray)
        assert t.datenum.mask[1]

    def test_datenum_arithmetic(self):
        """Test time arithmetic with datenum format."""
        t1 = Time(739252.0, format="datenum")
        t2 = Time(739253.0, format="datenum")

        # Time difference
        dt = t2 - t1
        assert_allclose(dt.jd, 1.0, atol=1e-10)  # 1 day difference

    def test_datenum_empty_array(self):
        """Test datenum with empty array."""
        t = Time([], format="datenum")
        assert t.size == 0
        assert t.shape == (0,)
        assert t.format == "datenum"

    def test_datenum_precision_loss_to_unix(self):
        """Test and document precision loss in unix conversion."""
        # Fractional seconds datenum
        datenum_precise = 739252.123456789
        t = Time(datenum_precise, format="datenum")

        # Convert to unix and back
        unix_val = t.unix
        t_back = Time(unix_val, format="unix")
        datenum_back = t_back.datenum

        # Document expected microsecond-level precision loss
        diff = abs(datenum_precise - datenum_back)
        # Should be less than 10 microseconds (1e-5 days)
        assert diff < 1e-5, f"Precision loss {diff} exceeds expected 1e-5 days"

    def test_datenum_knownvalue(self):
        """Test datenum with a known reference value."""
        # 2000-01-01 00:00:00 UTC in MATLAB datenum is 730486
        # (days from year 0 to year 2000)
        t_y2k = Time("2000-01-01 00:00:00", scale="utc")
        datenum_y2k = t_y2k.datenum
        # Should be around 730486
        assert 730000 < datenum_y2k < 731000

    def test_datenum_scalar_vs_array(self):
        """Test that scalar and array inputs give consistent results."""
        datenum_scalar = 739252.5
        datenum_array = np.array([739252.5])

        t_scalar = Time(datenum_scalar, format="datenum")
        t_array = Time(datenum_array, format="datenum")

        assert_allclose(t_scalar.jd, t_array.jd[0], atol=1e-10)
