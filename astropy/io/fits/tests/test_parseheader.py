import io
from collections import OrderedDict

import pytest

from astropy.io.fits._utils import parse_header
from astropy.utils.data import get_pkg_data_fileobj

required_keywords = ["SIMPLE", "BITPIX", "NAXIS"]


@pytest.mark.parametrize(
    "filename, expected_keywords",
    [
        ("data/blank.fits", required_keywords),
        ("data/history_header.fits", required_keywords),
        ("data/test0.fits", required_keywords + ["EXTEND"]),
        ("data/test1.fits", required_keywords + ["EXTEND"]),
        ("data/btable.fits", required_keywords + ["EXTEND"]),
        ("data/table.fits", required_keywords + ["EXTEND"]),
        ("data/ascii.fits", required_keywords + ["EXTEND"]),
        ("data/tb.fits", required_keywords + ["EXTEND"]),
        (
            "data/group.fits",
            required_keywords + ["GROUPS", "PCOUNT", "GCOUNT"],
        ),
        (
            "data/random_groups.fits",
            required_keywords + ["GROUPS", "PCOUNT", "GCOUNT"],
        ),
        (
            "data/invalid/group_invalid.fits",
            required_keywords + ["EXTEND", "GROUPS", "PCOUNT", "GCOUNT"],
        ),
        (
            "data/checksum.fits",
            required_keywords
            + [
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
def test_parse_header_keywords(filename, expected_keywords):
    """Test that parse_header correctly parses all the standard 8 character keywords and doesnt parse COMMENT, HISTORY, CONTINUE, HIERARCH or any keyword with length > 8"""
    unexpected_keywords = ["COMMENT", "HISTORY", "CONTINUE", "HIERARCH"]

    with get_pkg_data_fileobj(
        filename, package="astropy.io.fits.tests", encoding="binary"
    ) as f:
        header_str, cards = parse_header(f)

    # Blocks are always padded to a multiple of 2880, therefore the header length must be a multiple of 2880 as well but not 0
    assert len(header_str) % 2880 == 0 and len(header_str) > 0
    assert "END" in header_str

    for keyword in cards:
        assert len(keyword) <= 8
        assert keyword not in unexpected_keywords

    for keyword in expected_keywords:
        assert keyword in cards


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
    cards_dict = OrderedDict(
        [
            ("A", "A= 5".ljust(80)),
            ("SIMPLE", "SIMPLE  =                    T".ljust(80)),
            ("BITPIX", "BITPIX  =                  -32".ljust(80)),
            ("NAXIS", "NAXIS   =                    0".ljust(80)),
            ("EXTEND", "EXTEND  =                    T".ljust(80)),
        ]
    )
    end_card = "END" + " " * 77
    header = "".join(card.ljust(80) for card in cards_dict.values()) + end_card.ljust(
        80
    )
    header = header + " " * (2880 - len(header))
    header_bytes = header.encode("ascii")
    header_str, cards_result = parse_header(io.BytesIO(header_bytes))

    assert cards_result == cards_dict
    assert header == header_str
