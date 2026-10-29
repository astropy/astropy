# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Utilities to parse CODATA plain-text tables."""

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class CODATAEntry:
    """A single entry from a CODATA table."""

    name: str
    value: float
    uncertainty: float
    unit: str


def _parse_number_field(field: str):
    """Parse a CODATA numeric field."""
    normalized = field.replace(" ", "").replace("...", "")
    # if the field is empty or just contains ellipses, raise an error since a numeric value is expected
    if not normalized:
        raise ValueError("Empty numeric field.")
    return float(normalized)


def _parse_uncertainty_field(field: str):
    """Parse CODATA uncertainty field."""
    if "(exact)" in field:
        return 0.0

    normalized = field.replace(" ", "")
    # if the field is empty or just contains ellipses, treat it as zero uncertainty but not exact
    if not normalized:
        return 0.0

    return float(normalized)


def parse_codata_text(
    text: str,
    *,
    name_width: int = 60,
    value_width: int = 25,
    uncertainty_width: int = 25,
) -> dict[str, CODATAEntry]:
    """Parse CODATA fixed-width text into a dict keyed by constant name."""
    entries: dict[str, CODATAEntry] = {}
    in_table = False

    value_start = name_width
    uncertainty_start = value_start + value_width
    unit_start = uncertainty_start + uncertainty_width

    for raw_line in text.splitlines():
        line = raw_line.rstrip("\n")

        if not in_table:
            if set(line.strip()) == {"-"}:
                in_table = True
            continue

        if not line.strip():
            continue

        name = line[:name_width].rstrip()
        if not name:
            continue

        value_field = line[value_start:uncertainty_start]
        uncertainty_field = line[uncertainty_start:unit_start]
        unit = line[unit_start:].strip()

        value = _parse_number_field(value_field)
        uncertainty = _parse_uncertainty_field(uncertainty_field)

        entries[name] = CODATAEntry(
            name=name,
            value=value,
            uncertainty=uncertainty,
            unit=unit,
        )

    if not entries:
        raise ValueError("No CODATA entries were parsed from input text.")

    return entries


def parse_codata_file(path: str | Path) -> dict[str, CODATAEntry]:
    """Parse a CODATA table from a text file path."""
    text = Path(path).read_text(encoding="utf-8")
    return parse_codata_text(text)
