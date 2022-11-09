from pytest import raises

from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from astropy.wcs import WCS
from astropy.wcs.wcsapi.utils import deserialize_class, wcs_info_str


def test_construct():
    result = deserialize_class(("astropy.units.Quantity", (10,), {"unit": "deg"}))
    assert_quantity_allclose(result, 10 * u.deg)


def test_noconstruct():
    result = deserialize_class(
        ("astropy.units.Quantity", (), {"unit": "deg"}), construct=False
    )
    assert result == (u.Quantity, (), {"unit": "deg"})


def test_invalid():
    with raises(ValueError) as exc:
        deserialize_class(("astropy.units.Quantity", (), {"unit": "deg"}, ()))
    assert exc.value.args[0] == "Expected a tuple of three values"


DEFAULT_1D_STR = """
WCS Transformation

This transformation has 1 pixel and 1 world dimensions

Array shape (Numpy order): None

Pixel Dim  Axis Name  Data size  Bounds
        0  None            None  None

World Dim  Axis Name  Physical Type  Units
        0  None       None           unknown

Correlation between pixel and world axes:

           Pixel Dim
World Dim    0
        0  yes
"""


def test_wcs_info_str():
    # The tests in test_sliced_low_level_wcs.py exercise wcs_info_str
    # extensively. This test is to ensure that the function exists and the
    # API of the function works as expected.

    wcs_empty = WCS(naxis=1)

    assert wcs_info_str(wcs_empty).strip() == DEFAULT_1D_STR.strip()
