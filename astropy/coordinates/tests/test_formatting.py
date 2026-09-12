"""
Tests the Angle string formatting capabilities.  SkyCoord formatting is in
test_sky_coord
"""

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates.angles import formats
from astropy.utils.compat.numpycompat import NUMPY_LT_2_1


def test_to_string_precision():
    # There are already some tests in test_api.py, but this is a regression
    # test for the bug in issue #1319 which caused incorrect formatting of the
    # seconds for precision=0

    angle = Angle(-1.23456789, unit=u.degree)

    assert angle.to_string(precision=3) == "-1d14m04.444s"
    assert angle.to_string(precision=1) == "-1d14m04.4s"
    assert angle.to_string(precision=0) == "-1d14m04s"

    angle2 = Angle(-1.23456789, unit=u.hourangle)

    assert angle2.to_string(precision=3, unit=u.hour) == "-1h14m04.444s"
    assert angle2.to_string(precision=1, unit=u.hour) == "-1h14m04.4s"
    assert angle2.to_string(precision=0, unit=u.hour) == "-1h14m04s"

    # Regression test for #7141
    angle3 = Angle(-0.5, unit=u.degree)
    assert angle3.to_string(precision=0, fields=3) == "-0d30m00s"
    assert angle3.to_string(precision=0, fields=2) == "-0d30m"
    assert angle3.to_string(precision=0, fields=1) == "-1d"


def test_to_string_decimal():
    # There are already some tests in test_api.py, but this is a regression
    # test for the bug in issue #1323 which caused decimal formatting to not
    # work

    angle1 = Angle(2.0, unit=u.degree)

    assert angle1.to_string(decimal=True, precision=3) == "2.000"
    assert angle1.to_string(decimal=True, precision=1) == "2.0"
    assert angle1.to_string(decimal=True, precision=0) == "2"

    angle2 = Angle(3.0, unit=u.hourangle)

    assert angle2.to_string(decimal=True, precision=3) == "3.000"
    assert angle2.to_string(decimal=True, precision=1) == "3.0"
    assert angle2.to_string(decimal=True, precision=0) == "3"

    angle3 = Angle(4.0, unit=u.radian)

    assert angle3.to_string(decimal=True, precision=3) == "4.000"
    assert angle3.to_string(decimal=True, precision=1) == "4.0"
    assert angle3.to_string(decimal=True, precision=0) == "4"


@pytest.mark.parametrize("sep", [":", ":.", "dms", "hms"])
@pytest.mark.parametrize(
    "angle",
    [
        Angle(np.pi / 12, "rad"),
        Angle(15, "deg"),
        Angle(15, "hourangle"),
    ],
)
def test_angle_to_string_decimal_with_sep_error(angle, sep):
    # see https://github.com/astropy/astropy/pull/16085#discussion_r1501177163
    with pytest.raises(
        ValueError,
        match=rf"With decimal=True, separator cannot be used \(got sep='{sep}'\)",
    ):
        angle.to_string(sep=sep, decimal=True)


def test_to_string_formats():
    a = Angle(1.113355, unit=u.deg)
    latex_str = r"$1^\circ06{}^\prime48.078{}^{\prime\prime}$"
    assert a.to_string(format="latex") == latex_str
    assert a.to_string(format="latex_inline") == latex_str
    assert a.to_string(format="unicode") == "1°06′48.078″"

    a = Angle(1.113355, unit=u.hour)
    latex_str = r"$1^{\mathrm{h}}06^{\mathrm{m}}48.078^{\mathrm{s}}$"
    assert a.to_string(format="latex") == latex_str
    assert a.to_string(format="latex_inline") == latex_str
    assert a.to_string(format="unicode") == "1ʰ06ᵐ48.078ˢ"

    a = Angle(1.113355, unit=u.radian)
    assert a.to_string(format="latex") == r"$1.11336\;\mathrm{rad}$"
    assert a.to_string(format="latex_inline") == r"$1.11336\;\mathrm{rad}$"
    assert a.to_string(format="unicode") == "1.11336 rad"


def test_to_string_decimal_formats():
    angle1 = Angle(2.0, unit=u.degree)

    assert angle1.to_string(decimal=True, format="generic") == "2 deg"
    assert angle1.to_string(decimal=True, format="latex") == "$2\\mathrm{{}^{\\circ}}$"
    assert angle1.to_string(decimal=True, format="unicode") == "2°"

    angle2 = Angle(3.0, unit=u.hourangle)
    assert angle2.to_string(decimal=True, format="generic") == "3 hourangle"
    assert angle2.to_string(decimal=True, format="latex") == "$3\\mathrm{{}^{h}}$"
    assert angle2.to_string(decimal=True, format="unicode") == "3ʰ"

    angle3 = Angle(4.0, unit=u.radian)

    assert angle3.to_string(decimal=True, format="generic") == "4 rad"
    assert angle3.to_string(decimal=True, format="latex") == "$4\\;\\mathrm{rad}$"
    assert angle3.to_string(decimal=True, format="unicode") == "4 rad"

    with pytest.raises(ValueError, match="Unknown format"):
        angle3.to_string(decimal=True, format="myformat")


def test_to_string_fields():
    a = Angle(1.113355, unit=u.deg)
    assert a.to_string(fields=1) == r"1d"
    assert a.to_string(fields=2) == r"1d07m"
    assert a.to_string(fields=3) == r"1d06m48.078s"


def test_to_string_padding():
    a = Angle(0.5653, unit=u.deg)
    assert a.to_string(unit="deg", sep=":", pad=True) == r"00:33:55.08"

    # Test to make sure negative angles are padded correctly
    a = Angle(-0.5653, unit=u.deg)
    assert a.to_string(unit="deg", sep=":", pad=True) == r"-00:33:55.08"


def test_sexagesimal_rounding_up():
    a = Angle(359.999999999999, unit=u.deg)

    assert a.to_string(precision=None) == "360d00m00s"
    assert a.to_string(precision=4) == "360d00m00.0000s"
    assert a.to_string(precision=5) == "360d00m00.00000s"
    assert a.to_string(precision=6) == "360d00m00.000000s"
    assert a.to_string(precision=7) == "360d00m00.0000000s"
    assert a.to_string(precision=8) == "360d00m00.00000000s"
    assert a.to_string(precision=9) == "359d59m59.999999996s"

    a = Angle(3.999999, unit=u.deg)
    assert a.to_string(fields=2, precision=None) == "4d00m"
    assert a.to_string(fields=2, precision=1) == "4d00m"
    assert a.to_string(fields=2, precision=5) == "4d00m"
    assert a.to_string(fields=1, precision=1) == "4d"
    assert a.to_string(fields=1, precision=5) == "4d"


@pytest.mark.parametrize(
    "angle, precision, expected",
    [
        # Just below half a unit in the last digit shown: no carry.
        ("1d2m59.49999s", 0, "1d02m59s"),
        ("1d2m59.94s", 1, "1d02m59.9s"),
        ("1d2m59.994s", 2, "1d02m59.99s"),
        # At or above it: the seconds carry into the minutes, and on into
        # the degrees.
        ("1d2m59.5s", 0, "1d03m00s"),
        ("1d2m59.96s", 1, "1d03m00.0s"),
        ("1d59m59.5s", 0, "2d00m00s"),
        ("-1d2m59.49999s", 0, "-1d02m59s"),
    ],
)
def test_sexagesimal_seconds_round_at_half_a_unit(angle, precision, expected):
    assert Angle(angle).to_string(precision=precision) == expected


@pytest.mark.parametrize(
    "angle, expected",
    [
        # Below half a degree of the next one: no carry.  The seconds must not
        # be rounded into the minutes first, or these would round twice.
        ("285d29m41.76s", "285d"),
        ("285d29m59.9s", "285d"),
        ("-285d29m41.76s", "-285d"),
        # Half a degree or more: carry.
        ("285d30m", "286d"),
        ("285d30m00.1s", "286d"),
    ],
)
def test_sexagesimal_degrees_round_once(angle, expected):
    assert Angle(angle).to_string(fields=1) == expected


def test_to_string_scalar():
    a = Angle(1.113355, unit=u.deg)
    assert isinstance(a.to_string(), str)


def test_to_string_radian_with_precision():
    """
    Regression test for a bug that caused ``to_string`` to crash for angles in
    radians when specifying the precision.
    """

    # Check that specifying the precision works
    a = Angle(3.0, unit=u.rad)
    assert a.to_string(precision=3, sep="fromunit") == "3.000 rad"


def test_sexagesimal_round_down():
    a1 = Angle(1, u.deg).to(u.hourangle)
    a2 = Angle(2, u.deg)
    assert a1.to_string() == "0h04m00s"
    assert a2.to_string() == "2d00m00s"


def test_to_string_fields_colon():
    a = Angle(1.113355, unit=u.deg)
    assert a.to_string(fields=2, sep=":") == "1:07"
    assert a.to_string(fields=3, sep=":") == "1:06:48.078"
    assert a.to_string(fields=1, sep=":") == "1"


# A spread of angles that exercises the sign, the carrying between fields, the
# rounding and the non-finite handling of the vectorized formatter.  Kept well
# above the size threshold so that ``to_string`` takes the fast array path.
TO_STRING_ANGLES = np.array(
    [
        0.0,
        -0.0,
        5.25,
        -5.25,
        12.3456789,
        -12.3456789,
        59.9999999,
        359.9999999,
        90.0,
        -90.0,
        123.4567,
        0.5 / 3600,
        1e-9,
        720.0,
        np.nan,
        np.inf,
        -np.inf,
    ]
)

# Representative keyword combinations covering every option that changes the
# output: separators, padding, number of fields, precision, sign and format.
TO_STRING_KWARGS = [
    {},
    {"sep": ":"},
    {"sep": ("-", ":")},
    {"sep": ""},
    {"pad": True},
    {"fields": 1},
    {"fields": 2},
    {"precision": 0},
    {"precision": 8},
    {"alwayssign": True},
    {"format": "latex"},
    {"format": "latex_inline"},
    {"format": "unicode"},
    {"precision": 2, "pad": True, "alwayssign": True, "sep": ":"},
]


@pytest.mark.parametrize("unit", [u.degree, u.hourangle])
@pytest.mark.parametrize("kwargs", TO_STRING_KWARGS)
def test_to_string_vectorized_matches_per_element(monkeypatch, unit, kwargs):
    """The fast array formatting of ``to_string`` matches the per-element path.

    The same call is run twice: once normally, and once with the vectorized
    helper disabled so it falls back to looping over the elements.  The two
    results must be identical.
    """
    angle = Angle(TO_STRING_ANGLES, unit=unit)
    fast = np.asarray(angle.to_string(**kwargs))

    monkeypatch.setattr(
        formats, "_decimal_to_sexagesimal_string_array", lambda *args, **kw: None
    )
    per_element = np.asarray(angle.to_string(**kwargs))

    assert_array_equal(fast, per_element)


@pytest.mark.skipif(NUMPY_LT_2_1, reason="vectorized formatter requires NumPy >= 2.1")
@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"sep": ("d", "m", "s")},
        {"sep": ("-", ":"), "pad": True},
        {"fields": 1},
        {"fields": 2, "precision": 0},
        {"precision": 5},
        {"sep": (), "precision": 2},
    ],
)
def test_decimal_to_sexagesimal_string_array(kwargs):
    """The array formatter matches ``_decimal_to_sexagesimal_string`` per value.

    ``Angle.to_string`` handles NaN itself, so it is excluded here. The broad
    option coverage lives in ``test_to_string_vectorized_matches_per_element``;
    this just checks the helper directly for a few representative cases.
    """
    values = np.array(
        [0.0, -0.0, 5.25, -5.25, 59.9999999, 12.3456789, 123.456, np.inf, -np.inf]
    )
    got = formats._decimal_to_sexagesimal_string_array(values, **kwargs)
    ref = np.array(
        [formats._decimal_to_sexagesimal_string(v, **kwargs) for v in values]
    )
    assert_array_equal(got, ref)


@pytest.mark.skipif(NUMPY_LT_2_1, reason="vectorized formatter requires NumPy >= 2.1")
def test_to_string_array_falls_back_for_huge_degrees():
    """Degrees too large for the int64 buffer defer to the per-element path."""
    values = np.array([1e19, 2e19, 3e19])
    assert formats._decimal_to_sexagesimal_string_array(values) is None
    # ``to_string`` still returns the correct result via the fallback.
    angle = Angle(values, unit=u.degree)
    expected = np.array(
        [formats._decimal_to_sexagesimal_string(v, sep=("d", "m", "s")) for v in values]
    )
    assert_array_equal(np.asarray(angle.to_string(sep="dms")), expected)


def test_to_string_array_falls_back_on_old_numpy(monkeypatch):
    """The vectorized helper defers to the per-element path on NumPy < 2.1."""
    monkeypatch.setattr(formats, "NUMPY_LT_2_1", True)
    assert formats._decimal_to_sexagesimal_string_array(np.arange(10.0)) is None
