# Licensed under a 3-clause BSD style license - see LICENSE.rst

import functools

import numpy as np
import numpy.testing as npt
import pytest

from astropy import units as u
from astropy.coordinates import EarthLocation
from astropy.table import QTable, Table, vstack
from astropy.time import Time, TimeDelta, conf
from astropy.utils import iers
from astropy.utils.compat.optional_deps import HAS_H5PY
from astropy.utils.masked import Masked

allclose_sec = functools.partial(
    np.allclose, rtol=2.0**-52, atol=2.0**-52 * 24 * 3600
)  # 20 ps atol
is_masked = np.ma.is_masked


def test_simple():
    t = Time([1, 2, 3], format="cxcsec")
    assert t.masked is False
    assert np.all(t.mask == [False, False, False])

    # Before masking, format output does not have a mask
    # (it is an ndarray like always)
    assert not hasattr(t.value, "mask")
    assert not hasattr(t.unix, "mask")

    t[2] = np.ma.masked
    assert t.masked is True
    assert np.all(t.mask == [False, False, True])
    assert allclose_sec(t.value[:2], [1, 2])
    assert is_masked(t.value[2])
    assert is_masked(t[2].value)

    # After masking format output has a mask.
    assert hasattr(t.value, "mask")
    assert hasattr(t.unix, "mask")

    # Can also unmask.
    t[2] = np.ma.nomask
    assert np.all(t.mask == [False, False, False])
    # But since the internal data are still masked, the instance stays masked too,
    # as does any output.
    assert t.masked
    assert hasattr(t.value, "mask")
    assert hasattr(t.unix, "mask")
    # Combo just for completeness
    t[1:] = np.ma.masked
    t[1] = np.ma.nomask
    assert np.all(t.mask == [False, False, True])
    assert t.masked


def test_unmasked():
    t = Time([1, 2, 3], format="cxcsec")
    t[2] = np.ma.masked
    assert t.masked
    t_unmasked = t.unmasked
    assert not t_unmasked.masked
    assert not hasattr(t_unmasked.value, "mask")
    assert not hasattr(t_unmasked.unix, "mask")
    assert (t_unmasked.value == t.value.unmasked).all()
    # Check that data is shared.
    assert np.may_share_memory(t_unmasked._time.jd1, t._time.jd1)


@pytest.mark.parametrize("fill_value", [4, Time(5, format="cxcsec")])
def test_filled(fill_value):
    t = Time([1, 2, 3], format="cxcsec")
    t[2] = np.ma.masked
    # Fill with something suitable.
    t_filled = t.filled(fill_value)
    assert not t_filled.masked
    assert not hasattr(t_filled.value, "mask")
    assert not hasattr(t_filled.unix, "mask")
    expected = t.value.filled(Time(fill_value, format="cxcsec").value)
    assert (t_filled.value == expected).all()
    # Check that data is not shared.
    assert not np.may_share_memory(t_filled._time.jd1, t._time.jd1)


def test_scalar_init():
    t = Time("2000:001")
    assert t.masked is False
    assert t.mask == np.array(False)


def test_mask_not_writeable():
    t = Time("2000:001")
    with pytest.raises(
        AttributeError,
        match=rf"property 'mask' of '{t.__class__.__name__}' object has no setter",
    ):
        t.mask = True

    t = Time(["2000:001"])
    with pytest.raises(ValueError, match="assignment destination is read-only"):
        t.mask[0] = True

    # But we can set it to masked directly.
    t[0] = np.ma.masked
    assert np.all(t.mask == [True])
    # After this, the mask should again not be writeable.
    with pytest.raises(ValueError, match="assignment destination is read-only"):
        t.mask[0] = False

    # Should also not be writeable if we initialize with masked elements.
    t = Time(Masked(["2000:001", "2000:002"], mask=[False, True]))
    with pytest.raises(ValueError, match="assignment destination is read-only"):
        t.mask[0] = True
    # But again we can set it to masked directly.
    t[0] = np.ma.masked
    assert np.all(t.mask == [True, True])
    # Check that the mask remains shared.
    assert np.may_share_memory(t._time.jd1.mask, t._time.jd2.mask)


def test_setitem_masked_value():
    # Regression test for gh-20173: assigning a masked Time into an unmasked
    # one silently dropped the mask, revealing the underlying value again.
    t = Time(["2000:001", "2000:002", "2000:003"])
    assert not t.masked

    value = Time(["2001:001", "2001:002"])
    value[1] = np.ma.masked

    t[:2] = value
    assert t.masked
    assert np.all(t.mask == [False, True, False])
    assert t.unmasked[0] == Time("2001:001")
    assert t.unmasked[2] == Time("2000:003")
    # jd1 and jd2 should share the mask, like they do elsewhere.
    assert np.may_share_memory(t._time.jd1.mask, t._time.jd2.mask)

    # Scalar assignment of a masked element should mask the target too.
    t2 = Time(["2000:001", "2000:002"])
    t2[0] = value[1]
    assert np.all(t2.mask == [True, False])

    # An unmasked value must not clear an existing mask elsewhere, and
    # assigning over a masked element unmasks it again.
    t[0] = np.ma.masked
    t[1:] = Time(["2002:001", "2002:002"])
    assert np.all(t.mask == [True, False, False])


def test_vstack_masked():
    # Regression test for gh-20173: vstack dropped the mask of Time columns.
    t1 = QTable()
    t1["time"] = Time(["2024-01-01T01:00:00", "2024-01-01T02:00:00"])
    t1["time"][1] = np.ma.masked

    t2 = QTable()
    t2["time"] = Time(["2024-01-02T03:00:00", "2024-01-02T04:00:00"])
    t2["time"][0] = np.ma.masked

    combined = vstack([t1, t2])
    assert combined["time"].masked
    assert np.all(combined["time"].mask == [False, True, True, False])
    assert np.all(combined["time"].unmasked[:2] == t1["time"].unmasked)
    assert np.all(combined["time"].unmasked[2:] == t2["time"].unmasked)


def test_str():
    t = Time(["2000:001", "2000:002"])
    t[1] = np.ma.masked
    assert str(t) == "['2000:001:00:00:00.000'                     ———]"
    assert (
        repr(t)
        == "<Time object: scale='utc' format='yday' value=['2000:001:00:00:00.000'                     ———]>"
    )

    expected = [
        "MaskedNDArray(['2000-01-01 00:00:00.000',                       ———],",
        "              dtype='<U23')",
    ]

    # Note that we need to take care to allow for big-endian platforms,
    # for which the dtype will be >U23 instead of <U23, which we do with
    # the call to replace().
    assert repr(t.iso).replace(">U23", "<U23").splitlines() == expected

    # Assign value to unmask, though the instance stays masked.
    t[1] = "2000:111"
    assert str(t) == "['2000:001:00:00:00.000' '2000:111:00:00:00.000']"
    assert t.masked


def test_transform():
    with iers.conf.set_temp("auto_download", False):
        t = Time(["2000:001", "2000:002"])
        t[1] = np.ma.masked

        # Change scale (this tests the ERFA machinery with masking as well)
        t_ut1 = t.ut1
        assert is_masked(t_ut1.value[1])
        assert not is_masked(t_ut1.value[0])
        assert np.all(t_ut1.mask == [False, True])
        # Check the mask is a copy, so we won't back-propagate changes.
        assert not np.may_share_memory(t_ut1.mask, t.mask)

        # Change format
        t_unix = t.unix
        assert is_masked(t_unix[1])
        assert not is_masked(t_unix[0])
        assert np.all(t_unix.mask == [False, True])
        # Check the mask is a copy.
        assert not np.may_share_memory(t_unix.mask, t.mask)


def test_masked_input():
    v0 = np.ma.MaskedArray([[1, 2], [3, 4]])  # No masked elements
    v1 = np.ma.MaskedArray([[1, 2], [3, 4]], mask=[[True, False], [False, False]])
    v2 = np.ma.MaskedArray([[10, 20], [30, 40]], mask=[[False, False], [False, True]])

    # Init from various combinations of masked arrays
    t = Time(v0, format="cxcsec")
    assert np.ma.allclose(t.value, v0)
    assert np.all(t.mask == [[False, False], [False, False]])
    assert t.masked is False

    t = Time(v1, format="cxcsec")
    assert np.ma.allclose(t.value, v1)
    assert np.all(t.mask == v1.mask)
    value = t.value
    assert np.all(value.mask == v1.mask)
    assert t.masked is True
    # Check that masked are not shared with input or output.
    assert not np.may_share_memory(t.mask, v1.mask)
    assert not np.may_share_memory(value.mask, t.mask)
    # But they should be shared in the private _jd1, _jd2.
    assert np.may_share_memory(t._time.jd1.mask, t._time.jd2.mask)

    t = Time(v1, v2, format="cxcsec")
    assert np.ma.allclose(t.value, v1 + v2)
    assert np.all(t.mask == (v1 + v2).mask)
    assert t.masked is True

    t = Time(v0, v1, format="cxcsec")
    assert np.ma.allclose(t.value, v0 + v1)
    assert np.all(t.mask == (v0 + v1).mask)
    assert t.masked is True

    t = Time(0, v2, format="cxcsec")
    assert np.ma.allclose(t.value, v2)
    assert np.all(t.mask == v2.mask)
    assert t.masked is True

    # Init from a string masked array
    t_iso = t.iso
    t2 = Time(t_iso)
    assert np.all(t2.value == t_iso)
    assert np.all(t2.mask == v2.mask)
    assert t2.masked is True


@pytest.mark.parametrize("masked_cls", [np.ma.MaskedArray, Masked])
@pytest.mark.parametrize("val", [0, np.nan, [0], [np.nan]])
def test_all_masked_input(masked_cls, val):
    """Fix for #9612"""
    # Test with jd=0 and jd=np.nan. Both triggered an exception prior to #9624
    # due to astropy.utils.exceptions.ErfaError.
    val = masked_cls(val, mask=True)
    t = Time(val, format="jd")
    if val.ndim:
        assert str(t.iso).endswith("———]")
    else:
        assert str(t.iso).endswith("———")


@pytest.mark.parametrize("mask_cls", [np.ma.MaskedArray, Masked])
@pytest.mark.parametrize("val", ["", [".."], b"", [b".."]])
def test_all_masked_input_str(mask_cls, val):
    dates = mask_cls(val, mask=True)
    t = Time(dates, format="iso", in_subfmt="date")
    assert np.all(t.mask)
    assert np.all(t.unmasked == "2000-01-01")


def test_some_masked_input_str_no_subfmt():
    """Test for cls.fill_value() being longer than other input strings."""
    dates = Masked(["", "2023:001"], mask=True)
    t = Time(dates, format="yday")
    assert np.all(t.mask)
    assert np.all(t.unmasked == "2000:001:12:00:00.000")


def test_masked_input_and_unmasked_array_location():
    dates = Masked(["", "2023:001"], mask=[True, False])
    locations = EarthLocation([10, 20] * u.deg, [30, 40] * u.deg)
    t = Time(dates, format="yday", location=locations)
    assert t.masked
    assert not isinstance(t.location, Masked)
    assert np.all(t.unmasked.location == locations)


def test_serialize_fits_masked(tmp_path):
    tm = Time([1, 2, 3], format="cxcsec")
    tm[1] = np.ma.masked

    fn = tmp_path / "tempfile.fits"
    t = Table([tm])
    t.write(fn)

    t2 = Table.read(fn, astropy_native=True)

    # Time FITS handling does not current round-trip format in FITS
    t2["col0"].format = tm.format

    assert t2["col0"].masked
    assert np.all(t2["col0"].mask == [False, True, False])
    assert np.all(t2["col0"].value == t["col0"].value)


@pytest.mark.skipif(not HAS_H5PY, reason="Needs h5py")
def test_serialize_hdf5_masked(tmp_path):
    tm = Time([1, 2, 3], format="cxcsec")
    tm[1] = np.ma.masked

    fn = tmp_path / "tempfile.hdf5"
    t = Table([tm])
    t.write(fn, path="root", serialize_meta=True)
    t2 = Table.read(fn)

    assert t2["col0"].masked
    assert np.all(t2["col0"].mask == [False, True, False])
    assert np.all(t2["col0"].value == t["col0"].value)


# Ignore warning in MIPS https://github.com/astropy/astropy/issues/9750
@pytest.mark.filterwarnings("ignore:invalid value encountered")
@pytest.mark.parametrize("serialize_method", ["jd1_jd2", "formatted_value"])
def test_serialize_ecsv_masked(serialize_method, tmp_path):
    tm = Time([1, 2, 3], format="cxcsec")
    tm[1] = np.ma.masked

    tm.info.serialize_method["ecsv"] = serialize_method

    fn = tmp_path / "tempfile.ecsv"
    t = Table([tm])
    t.write(fn)
    t2 = Table.read(fn)

    assert t2["col0"].masked
    assert np.all(t2["col0"].mask == [False, True, False])
    # Serializing formatted_value loses some precision.
    atol = 0.1 * u.us if serialize_method == "formatted_value" else 1 * u.ps
    assert np.all(abs(t2["col0"] - t["col0"]) <= atol)


@pytest.mark.parametrize("format_", Time.FORMATS)
@pytest.mark.parametrize("masked_cls", [np.ma.MaskedArray, Masked])
@pytest.mark.parametrize("masked_array_type", ["numpy", "astropy"])
def test_all_formats(format_, masked_cls, masked_array_type):
    mjd = np.array([55000.25, 55000.375, 55001.125])
    mask = np.array([True, False, False])
    mjdm = masked_cls(mjd, mask=mask)
    assert np.may_share_memory(mjdm.mask, mask)
    t = Time(mjd, format="mjd")
    tm = Time(mjdm, format="mjd")
    assert tm.masked and np.all(tm.mask == mask)
    assert not np.may_share_memory(tm.mask, mask)

    out_cls = np.ma.MaskedArray if masked_array_type == "numpy" else Masked
    with conf.set_temp("masked_array_type", masked_array_type):
        # Get values in the given format, check that these have the appropriate class
        # and are correct (ignoring masked ones, which get adjusted on Time
        # initialization, in core._check_for_masked_and_fill).
        t_format = getattr(t, format_)
        tm_format = getattr(tm, format_)
        assert isinstance(tm_format, out_cls)
        expected = t_format
        assert np.all(tm_format == expected)

        # Check masked scalar.
        tm0_format = getattr(tm[0], format_)
        assert isinstance(tm0_format, out_cls)
        comparison = tm0_format == tm_format[0]
        assert comparison.mask
        if out_cls is Masked:
            assert comparison.unmasked
        elif tm0_format.dtype.names:
            assert comparison.data
        else:
            assert comparison is np.ma.masked

        # Check unmasked scalar.
        tm1_format = getattr(tm[1], format_)
        assert tm1_format == tm_format[1]

    # Verify that configuration gets reset to "astropy".
    tm_format2 = getattr(tm, format_)
    assert isinstance(tm_format2, Masked)

    # Verify that we can also initialize with the format and that this gives
    # the right result and mask too.
    t2 = Time(t_format, format=format_)
    tm2 = Time(tm_format, format=format_)
    assert tm2.masked and np.all(tm2.mask == mask)
    assert np.all(tm2 == t2)


def test_datetime64_with_nat():
    dt64 = np.array(["nat", "2001-02-03", "2001-02-04"], dtype="datetime64[ns]")
    mdt64 = Masked(dt64, mask=[False, True, False])
    t = Time(mdt64)
    assert t.masked
    assert np.all(t.mask == [True, True, False])


def test_insert_masked():
    """Time.insert must preserve the mask of the original object (gh-20230)."""
    t = Time(["2001:001", "2001:002", "2001:003"], out_subfmt="date")
    t[0] = np.ma.masked
    t[2] = np.ma.masked

    out = t.insert(1, "1999:001")
    assert out.masked
    assert np.all(out.mask == [True, False, False, True])
    assert out.value[1] == "1999:001"
    assert out.value[2] == "2001:002"
    npt.assert_array_equal(
        out.unmasked.value, ["2001:001", "1999:001", "2001:002", "2001:003"]
    )
    # jd1 and jd2 must share the mask, as they do everywhere else.
    assert np.may_share_memory(out._time.jd1.mask, out._time.jd2.mask)

    # Insert more than one value.
    out = t.insert(1, ["1999:001", "1999:002"])
    npt.assert_array_equal(out.mask, [True, False, False, False, True])

    # An unmasked Time stays unmasked.
    t_unmasked = Time(["2001:001", "2001:002"])
    assert not t_unmasked.insert(1, "1999:001").masked


def test_insert_masked_values():
    """Inserting masked values into an unmasked Time keeps their mask (gh-20230)."""
    t = Time([1.0, 2.0, 3.0], format="cxcsec")
    values = Time([10.0, 20.0], format="cxcsec")
    values[1] = np.ma.masked

    out = t.insert(1, values)
    assert out.masked
    npt.assert_array_equal(out.mask, [False, False, True, False, False])
    npt.assert_array_equal(out.unmasked.value, [1.0, 10.0, 20.0, 2.0, 3.0])

    # A bare Masked array as the inserted values works too.
    out = t.insert(1, Masked(np.array([10.0, 20.0]), mask=[False, True]))
    npt.assert_array_equal(out.mask, [False, False, True, False, False])

    # A bare numpy masked array scalar as the inserted value works too.
    out = t.insert(1, np.ma.array(10.0, mask=True))
    npt.assert_array_equal(out.mask, [False, True, False, False])

    # A masked scalar Time.
    value = Time(10.0, format="cxcsec")
    value[()] = np.ma.masked
    npt.assert_array_equal(t.insert(1, value).mask, [False, True, False, False])

    # Both sides masked.
    t = t.copy()
    t[2] = np.ma.masked
    out = t.insert(1, values)
    npt.assert_array_equal(out.mask, [False, False, True, False, True])


def test_insert_masked_2d():
    t = Time(np.array(["2001:001", "2001:002", "2001:003", "2001:004"]).reshape(2, 2))
    t[0, 1] = np.ma.masked
    out = t.insert(1, "2010:001")
    assert out.shape == (3, 2)
    npt.assert_array_equal(out.mask, [[False, True], [False, False], [False, False]])


def test_insert_masked_timedelta():
    dt = TimeDelta([1.0, 2.0, 3.0], format="jd")
    dt[1] = np.ma.masked
    out = dt.insert(0, TimeDelta(9.0, format="jd"))
    assert out.masked
    npt.assert_array_equal(out.mask, [False, False, True, False])
    npt.assert_array_equal(out.value[[0, 1, 3]], [9.0, 1.0, 3.0])


def test_insert_masked_table_add_row():
    """Table.add_row goes through Time.insert (gh-20230)."""
    t = Time(["2001:001", "2001:002", "2001:003"])
    t[0] = np.ma.masked
    t[2] = np.ma.masked
    tbl = Table([t], names=["time"])
    tbl.add_row([Time("2010:001")])
    assert tbl["time"].masked
    npt.assert_array_equal(tbl["time"].mask, [True, False, True, False])
