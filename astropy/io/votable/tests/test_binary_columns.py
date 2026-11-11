# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests of reading a whole BINARY or BINARY2 stream a column at a time, in
`~astropy.io.votable.tree.TableElement._parse_binary_columns`.
"""

import base64
import io
import struct

import numpy as np
import pytest

from astropy.io.votable import converters, parse_single_table
from astropy.io.votable.exceptions import W46
from astropy.io.votable.tree import Field, Resource, TableElement, VOTableFile


def every_datatype_votable(tabledata_format):
    """
    Return a VOTable document, serialized with *tabledata_format*, holding a
    column of each datatype that the column reader handles.
    """
    votable = VOTableFile()
    resource = Resource()
    votable.resources.append(resource)
    table = TableElement(votable)
    resource.tables.append(table)
    table.fields.extend(
        [
            Field(votable, name="double", datatype="double"),
            Field(votable, name="float", datatype="float"),
            Field(votable, name="long", datatype="long"),
            Field(votable, name="int", datatype="int"),
            Field(votable, name="short", datatype="short"),
            Field(votable, name="unsignedByte", datatype="unsignedByte"),
            Field(votable, name="boolean", datatype="boolean"),
            Field(votable, name="bit", datatype="bit"),
            Field(votable, name="char", datatype="char", arraysize="*"),
            Field(votable, name="char4", datatype="char", arraysize="4"),
            Field(votable, name="char8max", datatype="char", arraysize="8*"),
            Field(votable, name="unicode", datatype="unicodeChar", arraysize="*"),
            Field(votable, name="unicode3", datatype="unicodeChar", arraysize="3"),
        ]
    )
    table.create_arrays(nrows=3)
    table.array[0] = (
        1.5,
        2.5,
        2**40,
        -7,
        300,
        255,
        True,
        True,
        "a",
        "abcd",
        "",
        "é",
        "€",
    )
    table.array[1] = (
        np.nan,
        -1.0,
        -1,
        0,
        -300,
        0,
        False,
        False,
        "",
        "ab",
        "12345678",
        "",
        "a€b",
    )
    table.array[2] = (
        -0.0,
        3.0,
        5,
        2**31 - 1,
        1,
        7,
        False,
        True,
        "xyz",
        "",
        "q",
        "abc",
        "",
    )
    # BINARY can only mask floats, booleans and strings
    table.array.mask[1] = (True, True) + (False,) * 4 + (True, False) + (True,) * 5
    table.format = tabledata_format

    document = io.BytesIO()
    votable.to_xml(document)
    document.seek(0)
    return document


def binary_votable(tag, fields, stream):
    """
    Return a VOTable document with the given FIELD elements, and *stream* as
    its data in the *tag* element, BINARY or BINARY2.
    """
    return io.BytesIO(
        b'<?xml version="1.0" encoding="utf-8"?>\n'
        b'<VOTABLE version="1.4" xmlns="http://www.ivoa.net/xml/VOTable/v1.3">\n'
        b"<RESOURCE><TABLE>" + fields.encode() + b"<DATA><" + tag.encode() + b">"
        b'<STREAM encoding="base64">' + base64.b64encode(stream) + b"</STREAM>"
        b"</" + tag.encode() + b"></DATA></TABLE></RESOURCE></VOTABLE>"
    )


def record_column_reader(monkeypatch):
    """
    Make `TableElement._parse_binary_columns` record what it returns, in the
    list returned here.
    """
    results = []
    column_reader = TableElement._parse_binary_columns

    def recording_column_reader(self, *args):
        results.append(column_reader(self, *args))
        return results[-1]

    monkeypatch.setattr(TableElement, "_parse_binary_columns", recording_column_reader)
    return results


def read_by_columns(document, monkeypatch, **kwargs):
    """
    Return the array of ``parse_single_table(document, **kwargs)``, and fail
    if the column reader leaves the table to the value-by-value reader.
    """
    results = record_column_reader(monkeypatch)
    array = parse_single_table(document, **kwargs).array
    assert results[0] is not None, "read value by value"
    return array


@pytest.mark.parametrize("tabledata_format", ["binary", "binary2"])
def test_columns_match_value_by_value(tabledata_format, monkeypatch):
    """
    parse_single_table on a table with a column of each datatype the column
    reader handles gives the same values and mask as reading value by value.
    """
    by_columns = parse_single_table(every_datatype_votable(tabledata_format)).array
    monkeypatch.setattr(TableElement, "_parse_binary_columns", lambda *args: None)
    by_values = parse_single_table(every_datatype_votable(tabledata_format)).array

    assert by_columns.dtype == by_values.dtype
    assert by_columns.tolist() == by_values.tolist()


def test_column_reader_handles_every_datatype(monkeypatch):
    """
    parse_single_table on a BINARY2 table with a column of each datatype the
    column reader handles reads it with the column reader.
    """
    results = record_column_reader(monkeypatch)

    parse_single_table(every_datatype_votable("binary2"))

    assert results[0] is not None


def test_bounded_string_over_its_limit_warns():
    """
    parse_single_table on a BINARY2 string longer than its arraysize limit
    warns W46, as reading value by value does.
    """
    document = binary_votable(
        "BINARY2",
        '<FIELD name="c" datatype="char" arraysize="4*"/>',
        b"\x00" + struct.pack(">I", 6) + b"abcdef",
    )

    with pytest.warns(W46):
        parse_single_table(document)


def test_partial_last_row_is_dropped():
    """
    parse_single_table on a BINARY2 stream that ends partway through its
    third row gives the two complete rows.
    """
    row = b"\x00" + struct.pack(">d", 1.5) + struct.pack(">I", 2) + b"ab"
    document = binary_votable(
        "BINARY2",
        '<FIELD name="d" datatype="double"/>'
        '<FIELD name="c" datatype="char" arraysize="*"/>',
        row + row + row[:-1],
    )

    array = parse_single_table(document).array

    assert array.tolist() == [(1.5, "ab"), (1.5, "ab")]


def test_bad_string_in_a_column_not_read_raises():
    """
    parse_single_table with columns=["d"] on a BINARY2 table whose other
    column holds a non-ASCII char string raises UnicodeDecodeError.

    Reading value by value decodes every column, so the column reader checks
    the strings it does not return too.
    """
    document = binary_votable(
        "BINARY2",
        '<FIELD name="c" datatype="char" arraysize="*"/>'
        '<FIELD name="d" datatype="double"/>',
        b"\x00" + struct.pack(">I", 2) + b"a\xe9" + struct.pack(">d", 2.0),
    )

    with pytest.raises(UnicodeDecodeError):
        parse_single_table(document, columns=["d"])


def test_empty_stream_gives_no_rows():
    """
    parse_single_table on a BINARY2 table with an empty stream gives an array
    with no rows.
    """
    document = binary_votable(
        "BINARY2",
        '<FIELD name="d" datatype="double"/>'
        '<FIELD name="c" datatype="char" arraysize="*"/>',
        b"",
    )

    array = parse_single_table(document).array

    assert len(array) == 0


def test_partial_last_row_with_a_bad_string_raises():
    """
    parse_single_table on a BINARY2 stream whose partial last row holds a
    unicodeChar string that does not decode raises UnicodeDecodeError, as
    reading value by value does.
    """
    row = b"\x00" + struct.pack(">I", 1) + "a".encode("utf_16_be")
    document = binary_votable(
        "BINARY2",
        '<FIELD name="u" datatype="unicodeChar" arraysize="*"/>'
        '<FIELD name="d" datatype="double"/>',
        row + struct.pack(">d", 1.0) + b"\x00" + struct.pack(">I", 1) + b"\xd8\x00",
    )

    with pytest.raises(UnicodeDecodeError):
        parse_single_table(document)


@pytest.mark.parametrize("datatype", ["char", "unicodeChar"])
def test_zero_length_string_is_read_by_columns(datatype, monkeypatch):
    """
    parse_single_table on a BINARY2 table with a string field of arraysize 0
    gives an empty string for it, with the column reader.
    """
    document = binary_votable(
        "BINARY2",
        f'<FIELD name="s" datatype="{datatype}" arraysize="0"/>'
        '<FIELD name="d" datatype="double"/>',
        b"\x00" + struct.pack(">d", 1.0),
    )

    array = read_by_columns(document, monkeypatch)

    assert array.tolist() == [("", 1.0)]


def test_rows_of_no_bytes_give_no_rows():
    """
    parse_single_table on a BINARY table whose only field is a char of
    arraysize 0, so that its rows take up no bytes, gives no rows.

    Reading value by value never stops on such rows.
    """
    document = binary_votable(
        "BINARY", '<FIELD name="c" datatype="char" arraysize="0"/>', b""
    )

    array = parse_single_table(document).array

    assert len(array) == 0


def test_unexpected_error_reads_value_by_value(monkeypatch):
    """
    parse_single_table on a BINARY2 table, when converting a column raises
    ValueError, gives the values read value by value.
    """

    def fail(self, stored):
        raise ValueError()

    monkeypatch.setattr(converters.Double, "_binparse_column", fail)
    document = binary_votable(
        "BINARY2",
        '<FIELD name="d" datatype="double"/>',
        b"\x00" + struct.pack(">d", 1.5),
    )

    array = parse_single_table(document).array

    assert array.tolist() == [(1.5,)]
