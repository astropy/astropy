from dataclasses import dataclass
from typing import Any

import pytest

from astropy.wcs import _wcs


class MockUnit:
    def __init__(self, unit_str: str) -> None:
        self._unit_str = unit_str

    def to_string(self, format: str = "fits") -> str:
        return self._unit_str


class MockInvalidUnit:
    pass


@dataclass(kw_only=True, slots=True, frozen=True)
class UnitProxyMatrix:
    incoming_data: tuple[Any, ...]
    expected_outputs: tuple[str, ...] | None
    error_match: str | None


@pytest.mark.parametrize(
    "test_case",
    [
        pytest.param(
            UnitProxyMatrix(
                incoming_data=(MockUnit("m/s"), MockUnit("deg")),
                expected_outputs=("m/s", "deg"),
                error_match=None,
            ),
            id="valid_duck_typed_units",
        ),
        pytest.param(
            UnitProxyMatrix(
                incoming_data=("m/s", "deg"),
                expected_outputs=("m/s", "deg"),
                error_match=None,
            ),
            id="valid_raw_strings_now_supported",
        ),
    ],
)
def test_unit_list_proxy_valid(test_case: UnitProxyMatrix) -> None:
    wcs_struct = _wcs.Wcsprm()
    wcs_struct.cunit = test_case.incoming_data
    proxy = wcs_struct.cunit

    assert type(proxy).__name__ == "UnitListProxy"

    for actual, expected in zip(proxy, test_case.expected_outputs, strict=True):
        assert hasattr(actual, "to_string")
        assert actual == expected


@pytest.mark.parametrize(
    "test_case",
    [
        pytest.param(
            UnitProxyMatrix(
                incoming_data=(MockInvalidUnit(), MockInvalidUnit()),
                expected_outputs=None,
                error_match="cannot be converted to a unit",
            ),
            id="missing_to_string_method",
        ),
    ],
)
def test_unit_list_proxy_exceptions(test_case: UnitProxyMatrix) -> None:
    wcs_struct = _wcs.Wcsprm()

    with pytest.raises(TypeError) as exc_info:
        wcs_struct.cunit = test_case.incoming_data

    assert test_case.error_match in str(exc_info.value).lower()


def test_unit_list_proxy_setitem() -> None:
    wcs_struct = _wcs.Wcsprm()
    wcs_struct.cunit = (MockUnit("m/s"), MockUnit("m/s"))
    proxy = wcs_struct.cunit

    proxy[0] = MockUnit("Hz")
    assert proxy[0] == "Hz"

    with pytest.raises(TypeError) as exc_info:
        proxy[0] = MockInvalidUnit()

    assert "cannot be converted to a unit" in str(exc_info.value).lower()
