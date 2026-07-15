from dataclasses import dataclass
from typing import Any

import pytest

# Bypass high-level astropy.wcs.WCS wrappers.
# Import the compiled C-extension directly to isolate the C-slots.
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
    ],
)
def test_unit_list_proxy_valid(test_case: UnitProxyMatrix) -> None:
    # Initialize the raw C-struct engine to allocate the memory buffer
    wcs_struct = _wcs.Wcsprm(header="NAXIS = 2")
    proxy = wcs_struct.cunit

    # Trigger the C-level setitem/set_unit_list boundary
    proxy[:] = test_case.incoming_data

    # Strictly check exact return types to ensure subclasses are stripped
    assert type(proxy) is _wcs.UnitListProxy

    for actual, expected in zip(proxy, test_case.expected_outputs, strict=True):
        assert type(actual) is str
        assert actual == expected


@pytest.mark.parametrize(
    "test_case",
    [
        pytest.param(
            UnitProxyMatrix(
                incoming_data=(MockInvalidUnit(), MockInvalidUnit()),
                expected_outputs=None,
                error_match="object must be a unit-like object",
            ),
            id="missing_to_string_method",
        ),
        pytest.param(
            UnitProxyMatrix(
                incoming_data=("m/s", "deg"),
                expected_outputs=None,
                error_match="object must be a unit-like object",
            ),
            id="raw_strings_rejected",
        ),
    ],
)
def test_unit_list_proxy_exceptions(test_case: UnitProxyMatrix) -> None:
    wcs_struct = _wcs.Wcsprm(header="NAXIS = 2")
    proxy = wcs_struct.cunit

    with pytest.raises(TypeError) as exc_info:
        proxy[:] = test_case.incoming_data

    # Lowercase and substring match to avoid brittle C-memory/address matching
    assert test_case.error_match in str(exc_info.value).lower()


def test_unit_list_proxy_setitem() -> None:
    wcs_struct = _wcs.Wcsprm(header="NAXIS = 2")
    proxy = wcs_struct.cunit

    proxy[:] = (MockUnit("m/s"), MockUnit("m/s"))

    proxy[0] = MockUnit("Hz")
    assert proxy[0] == "Hz"

    with pytest.raises(TypeError) as exc_info:
        proxy[0] = "deg"

    assert "object must be a unit-like object" in str(exc_info.value).lower()
