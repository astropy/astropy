# Licensed under a 3-clause BSD style license - see LICENSE.rst
#
import pytest

from astropy.io.ascii.mrt import MrtHeader

test_dat = [
    "names e d s i",
    "HD81809 1E-7 22.25608 +2 67",
    "HD103095 -31.6e5 +27.2500 -9E34 -30",
]


def test_find_map_note():
    good_twolines = (
        "Note (1)",
        ["", "    M = some description of M", "        with a continuation lines"],
    )
    text_with_equal = (
        "Note (1)",
        ["", "    M = 10.0 gives a magnitude and has", "    a continuation lines"],
    )
    # TODO: How should this case be treated? Only one key like this unlikely,
    # but a single-line note on a new line is also unlikely
    # I guess it boils down to which one is least likely
    single_compatible_line = ("Note (1)", ["", "    M = some flag"])

    hdr = MrtHeader()
    assert isinstance(hdr._find_map_note(*good_twolines), dict)
    assert (
        hdr._find_map_note(*good_twolines)["M"]
        == "some description of M with a continuation lines"
    )
    assert hdr._find_map_note(*text_with_equal) is None
    assert isinstance(hdr._find_map_note(*single_compatible_line), dict)


def test_extract_metadata_lines():
    pattern_not_enough_groups = "A = B"
    good_pattern = "([A-Z][a-z]+):(.*)$"
    lines_with_gaps = [
        "Title: A title",
        "Authors: Someone",
        "Table: A table",
        "" "Notes (1): A note",
    ]

    hdr = MrtHeader()
    with pytest.raises(ValueError):
        hdr._extract_metadata_lines([], pattern_not_enough_groups)
    hdr._extract_metadata_lines(lines_with_gaps, good_pattern, allow_gaps=True)
    hdr._extract_metadata_lines(lines_with_gaps[:3], good_pattern, allow_gaps=False)
    with pytest.raises(ValueError):
        hdr._extract_metadata_lines(lines_with_gaps, good_pattern, allow_gaps=False)
