# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests using the VOTable 1.6 unicode example tables by Mark Taylor.

The data files are copied unchanged from the top level of the VOTable
standard repository (https://github.com/ivoa-std/VOTable, commit 2836c14).
All three hold the same table in TABLEDATA, BINARY and BINARY2, except that
BINARY can only write a null array as an array of nulls.
"""

import io

import numpy as np
import pytest

from astropy.io.votable import parse
from astropy.table import Table
from astropy.utils.data import get_pkg_data_filename

FILES = [
    "data/unicode-example.vot",
    "data/unicode-example-binary.vot",
    "data/unicode-example-binary2.vot",
]

# Needs both the UTF-8 char handling from VOTable 1.6 and the char arrays
# from #19974. Remove once this branch is rebased on it.
needs_char_arrays = pytest.mark.xfail(
    reason="char arrays (arraysize like 10x2) need astropy#19974", strict=True
)

# The first three rows. The fourth row is all nulls.
EXPECTED = {
    "chr1": ["A", "B", "G"],
    "txt": ["Alpha", "Βητα", "Γ"],
    "word10": ["Alpha.....", "Βητα..", "Γꙮ.\U0001f600"],
    "words": [
        ["A123456789", "a123456789"],
        ["Βητα..", "beta......"],
        # In TABLEDATA the first string is padded with a blank, in BINARY
        # with a NUL. We treat both as padding, not part of the value, so
        # trailing blanks are stripped from each string of a TABLEDATA char
        # array. The spec doesn't say this explicitly yet.
        ["Γꙮ\U0001f600", "ABCDEFGHIJ"],
    ],
    "uchr1": ["A", "B", "Γ"],
    "utxt": ["Alpha", "Βητα", "Γαμμα"],
    "uword10": ["Alpha", "Βητα.", "Γꙮ..."],
    "uwords": [["A1234", "a1234"], ["Βητα.", "beta."], ["Γꙮ...", "ABCDE"]],
    "ialph": [1, 2, 3],
}

EXPECTED_PARAMS = {
    "soup1": ["oxtail"],
    "soup2": ["Борщ"],
    "philosophers": ["孔子", "Σωκρατης", "Wittgenstein"],
}


def _values(value):
    """A cell or PARAM value as a flat list of Python values."""
    return np.ma.asarray(value).ravel().tolist()


def _is_null(value, mask):
    # A null array may be masked, or (in BINARY) an array of empty strings.
    if np.all(mask):
        return True
    return all(v in ("", None) for v in _values(value))


def check_table(array, check_nulls=True):
    for name, expected in EXPECTED.items():
        for row, exp in enumerate(expected):
            assert _values(array[name][row]) == _values(exp), (name, row)
        if check_nulls:
            assert _is_null(array[name][3], array.mask[name][3]), name


def check_params(table):
    params = {p.name: _values(p.value) for p in table.params}
    for name, expected in EXPECTED_PARAMS.items():
        assert params[name] == expected, name


@needs_char_arrays
@pytest.mark.parametrize("filename", FILES)
def test_read(filename):
    votable = parse(get_pkg_data_filename(filename))
    assert votable.version == "1.6"
    assert votable.infos[0].value == "Ωμεγα ꙮ \U0001f600"

    table = votable.get_first_table()
    check_params(table)
    check_table(table.array)


@needs_char_arrays
@pytest.mark.parametrize("tabledata_format", ["tabledata", "binary", "binary2"])
@pytest.mark.parametrize("filename", FILES)
def test_rewrite(filename, tabledata_format):
    votable = parse(get_pkg_data_filename(filename))
    votable.get_first_table().format = tabledata_format

    out = io.BytesIO()
    votable.to_xml(out)
    out.seek(0)

    table = parse(out).get_first_table()
    check_params(table)
    check_table(table.array)


@needs_char_arrays
@pytest.mark.parametrize("tabledata_format", ["tabledata", "binary", "binary2"])
@pytest.mark.parametrize("filename", FILES)
def test_table_roundtrip(filename, tabledata_format):
    table = Table.read(get_pkg_data_filename(filename), format="votable")

    out = io.BytesIO()
    table.write(out, format="votable", tabledata_format=tabledata_format, version="1.6")
    out.seek(0)

    # unicodeChar columns are written as char in 1.6, so check the content
    # rather than the FIELD attributes.
    votable = parse(out)
    assert votable.version == "1.6"
    check_table(votable.get_first_table().array)
