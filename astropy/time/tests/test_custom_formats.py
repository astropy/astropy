from itertools import count

import pytest

import numpy as np

import astropy.units as u
from astropy._erfa import DJM0
from astropy.time import Time, TimeFormat
from astropy.time.utils import day_frac


class SpecificException(ValueError):
    pass


@pytest.fixture
def custom_format_name():
    for i in count():
        if not i:
            custom = f"custom_format_name"
        else:
            custom = f"custom_format_name_{i}"
        if custom not in Time.FORMATS:
            break
    yield custom
    Time.FORMATS.pop(custom, None)


def test_custom_time_format_set_jds_exception(custom_format_name):
    class Custom(TimeFormat):
        name = custom_format_name

        def set_jds(self, val, val2):
            raise SpecificException

    try:
        Time(7.0, format=custom_format_name)
    except ValueError as e:
        assert hasattr(e, "__cause__") and isinstance(e.__cause__, SpecificException)


def test_custom_time_format_val_type_exception(custom_format_name):
    class Custom(TimeFormat):
        name = custom_format_name

        def _check_val_type(self, val, val2):
            raise SpecificException

    try:
        Time(7.0, format=custom_format_name)
    except ValueError as e:
        assert hasattr(e, "__cause__") and isinstance(e.__cause__, SpecificException)


def test_custom_time_format_value_exception(custom_format_name):
    class Custom(TimeFormat):
        name = custom_format_name

        def set_jds(self, val, val2):
            self.jd1, self.jd2 = val, val2

        @property
        def value(self):
            raise SpecificException

    t = Time.now()
    with pytest.raises(SpecificException):
        getattr(t, custom_format_name)


def test_custom_time_format_fine(custom_format_name):
    class Custom(TimeFormat):
        name = custom_format_name

        def set_jds(self, val, val2):
            self.jd1, self.jd2 = val, val2

        @property
        def value(self):
            return self.jd1 + self.jd2

    t = Time.now()
    getattr(t, custom_format_name)
    t2 = Time(7, 9, format=custom_format_name)
    getattr(t2, custom_format_name)


def test_custom_time_format_forgot_property(custom_format_name):
    class Custom(TimeFormat):
        name = custom_format_name

        def set_jds(self, val, val2):
            self.jd1, self.jd2 = val, val2

        def value(self):
            return self.jd1, self.jd2

    t = Time.now()
    with pytest.raises(AttributeError):
        getattr(t, custom_format_name)

    t.format = custom_format_name
    with pytest.raises(AttributeError):
        t.value

    with pytest.raises(AttributeError):
        Time(7, 9, format=custom_format_name).value


def test_custom_time_format_problematic_name():
    assert "sort" not in Time.FORMATS, "problematic name in default FORMATS!"
    assert hasattr(Time, "sort")

    try:

        class Custom(TimeFormat):
            name = "sort"

            def set_jds(self, val, val2):
                self.jd1, self.jd2 = val, val2

            @property
            def value(self):
                return self.jd1, self.jd2

        t = Time.now()
        assert t.sort() == t, "bogus time format clobbers everyone's Time objects"

        t.format = "sort"
        if not isinstance(t.value, tuple):
            pytest.xfail("No good way to detect that `sort` is invalid")

        assert Time(7, 9, format="sort").value == (7, 9)

    finally:
        Time.FORMATS.pop("sort", None)


def test_mjd_longdouble_preserves_precision(custom_format_name):
    class CustomMJD(TimeFormat):
        name = custom_format_name

        def _check_val_type(self, val, val2):
            val = np.longdouble(val)
            if val2 is not None:
                raise ValueError("Only one value permitted")
            return val, 0

        def set_jds(self, val, val2):
            mjd1 = np.float64(np.floor(val))
            mjd2 = np.float64(val - mjd1)
            self.jd1, self.jd2 = day_frac(mjd1 + DJM0, mjd2)

        @property
        def value(self):
            mjd1, mjd2 = day_frac(self.jd1 - DJM0, self.jd2)
            return np.longdouble(mjd1) + np.longdouble(mjd2)

    m = 58000.0
    t = Time(m, format=custom_format_name)
    t2 = Time(m + 2 * m * np.finfo(np.longdouble).eps, format=custom_format_name)
    assert t != t2
    assert isinstance(getattr(t, custom_format_name), np.longdouble)
    assert getattr(t, custom_format_name) != getattr(t2, custom_format_name)


def test_mjd_unit_validation():
    with pytest.raises(u.UnitConversionError):
        Time(58000 * u.m, format="mjd")


def test_mjd_unit_conversion():
    assert Time(58000 * u.day, format="mjd") == Time(58000 * u.day, format="mjd")
    assert Time(58000 * u.day, format="mjd") != Time(58000 * u.s, format="mjd")
    assert Time(58000 * u.day, format="mjd") == Time(58000 * 86400 * u.s, format="mjd")


@pytest.mark.parametrize("f", ["mjd", "unix", "cxcsec"])
def test_existing_types_refuse_longdoubles(f):
    t = np.longdouble(getattr(Time(58000, format="mjd"), f))
    t2 = t + np.finfo(np.longdouble).eps * 2 * t
    try:
        tm = Time(np.longdouble(t), format=f)
    except ValueError:
        # Time processing makes it always ValueError not TypeError
        return
    else:
        # accepts long doubles, better preserve accuracy!
        assert Time(np.longdouble(t2), format=f) != tm
