import io

import pytest

from astropy.io.fits._utils import parse_header
from astropy.utils.data import get_pkg_data_filename


@pytest.mark.parametrize(
    "filename, expected_keywords",
    [
        ("data/blank.fits", ["SIMPLE", "BITPIX", "NAXIS"]),
        ("data/history_header.fits", ["SIMPLE", "BITPIX", "NAXIS"]),
        ("data/test0.fits", ["SIMPLE", "BITPIX", "NAXIS", "EXTEND"]),
        ("data/test1.fits", ["SIMPLE", "BITPIX", "NAXIS", "EXTEND"]),
        ("data/btable.fits", ["SIMPLE", "BITPIX", "NAXIS", "EXTEND"]),
        ("data/table.fits", ["SIMPLE", "BITPIX", "NAXIS", "EXTEND"]),
        ("data/ascii.fits", ["SIMPLE", "BITPIX", "NAXIS", "EXTEND"]),
        ("data/tb.fits", ["SIMPLE", "BITPIX", "NAXIS", "EXTEND", "NEXTEND"]),
        (
            "data/group.fits",
            ["SIMPLE", "BITPIX", "NAXIS", "GROUPS", "PCOUNT", "GCOUNT"],
        ),
        (
            "data/random_groups.fits",
            ["SIMPLE", "BITPIX", "NAXIS", "GROUPS", "PCOUNT", "GCOUNT"],
        ),
        (
            "data/invalid/group_invalid.fits",
            ["SIMPLE", "BITPIX", "NAXIS", "EXTEND", "GROUPS", "PCOUNT", "GCOUNT"],
        ),
        (
            "data/checksum.fits",
            [
                "SIMPLE",
                "BITPIX",
                "NAXIS",
                "NAXIS1",
                "NAXIS2",
                "EXTEND",
                "OBJECT",
                "TELESCOP",
                "EQUINOX",
                "CHECKSUM",
                "DATASUM",
            ],
        ),
    ],
)
def test_parse_header_expected(filename, expected_keywords):
    """Test that expected keywords are parsed from FITS files"""
    filepath = get_pkg_data_filename(filename)
    with open(filepath, "rb") as f:
        header_str, cards = parse_header(f)

    assert len(header_str) % 2880 == 0
    assert "END" in header_str
    for kw in expected_keywords:
        assert kw in cards, f"Expected keyword {kw} not found in {filename}"


@pytest.mark.parametrize(
    "filename,unexpected_keywords",
    [
        ("data/checksum.fits", ["COMMENT"]),
        ("data/history_header.fits", ["HISTORY"]),
    ],
)
def test_parse_header_unexpected(filename, unexpected_keywords):
    """Test that certain keywords are skipped (HISTORY, COMMENT, END)"""
    filepath = get_pkg_data_filename(filename)
    with open(filepath, "rb") as f:
        header_str, cards = parse_header(f)

    for kw in unexpected_keywords:
        assert kw not in cards, f"Unexpected keyword {kw} found in {filename}"


def test_parse_header_incomplete():
    """Test that block < 2880 with no END card raises exception"""
    with pytest.raises(Exception):
        parse_header(io.BytesIO(b"x" * 100))


def test_parse_header_no_end():
    """Test that exactly one block without END card raises exception"""
    with pytest.raises(Exception):
        parse_header(io.BytesIO(b" " * 2880))


def test_parse_header_short_keyword():
    """Test that keywords with separator within first 8 chars are parsed"""
    cards_list = [
        "A= 5                                                          ",
        "SIMPLE  =                    T                                      ",
        "BITPIX  =                  -32                                      ",
        "NAXIS   =                    0                                      ",
        "EXTEND  =                    T                                      ",
        "END" + " " * 77,
    ]
    header = "".join(card.ljust(80) for card in cards_list)
    header = header + " " * (2880 - len(header))
    header_bytes = header.encode("ascii")
    header_str, cards_result = parse_header(io.BytesIO(header_bytes))

    assert "A" in cards_result
    assert "SIMPLE" in cards_result
    assert "BITPIX" in cards_result
