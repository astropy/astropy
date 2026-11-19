# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Arrays of fixed-length strings: ``char`` and ``unicodeChar`` fields whose
arraysize has more than one dimension, such as "10x4" or "10x*".
"""

import base64
import io

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from astropy.io.votable import parse, tree
from astropy.io.votable.converters import (
    get_converter,
    table_column_to_votable_datatype,
)
from astropy.io.votable.exceptions import E01, E02, W46
from astropy.table import Column, Table

DATATYPES = ["char", "unicodeChar"]
FORMATS = ["tabledata", "binary", "binary2"]


def _votable(datatype, arraysize, cells):
    """A VOTable file with one column *c*, and one row per item of *cells*."""
    rows = "".join(f"<TR><TD>{cell}</TD></TR>" for cell in cells)
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<VOTABLE version="1.4" xmlns="http://www.ivoa.net/xml/VOTable/v1.3">
 <RESOURCE>
  <TABLE>
   <FIELD ID="c" name="c" datatype="{datatype}" arraysize="{arraysize}"/>
   <DATA><TABLEDATA>{rows}</TABLEDATA></DATA>
  </TABLE>
 </RESOURCE>
</VOTABLE>""".encode()


def _binary_votable(datatype, arraysize, stream):
    """A VOTable file with one column *c*, holding the BINARY bytes *stream*."""
    encoded = base64.b64encode(stream).decode("ascii")
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<VOTABLE version="1.4" xmlns="http://www.ivoa.net/xml/VOTable/v1.3">
 <RESOURCE>
  <TABLE>
   <FIELD ID="c" name="c" datatype="{datatype}" arraysize="{arraysize}"/>
   <DATA><BINARY><STREAM encoding="base64">{encoded}</STREAM></BINARY></DATA>
  </TABLE>
 </RESOURCE>
</VOTABLE>""".encode()


def _column(data, **kwargs):
    """The column *c* of the VOTable file *data*, as a masked array."""
    return parse(io.BytesIO(data), **kwargs).get_first_table().array["c"]


def _roundtrip(data, tabledata_format):
    """Read the VOTable file *data*, write it in another format, read that."""
    votable = parse(io.BytesIO(data))
    votable.get_first_table().format = tabledata_format
    written = io.BytesIO()
    votable.to_xml(written)
    return written.getvalue()


@pytest.mark.parametrize("datatype", DATATYPES)
def test_read_fixed_array(datatype):
    """
    parse() of an arraysize="3x2" cell "ABCDEF" gives the strings "ABC" and
    "DEF".  This is the FIELD of issue 17098.
    """
    column = _column(_votable(datatype, "3x2", ["ABCDEF"]))

    assert column.dtype == "U3"
    assert column.shape == (1, 2)
    assert column.tolist() == [["ABC", "DEF"]]


@pytest.mark.parametrize("datatype", DATATYPES)
def test_read_variable_array(datatype):
    """
    parse() of arraysize="3x*" cells "ABCDEF" and "GHI" gives 2 strings for
    the first row and 1 for the second.
    """
    column = _column(_votable(datatype, "3x*", ["ABCDEF", "GHI"]))

    assert column.dtype == object
    assert column[0].tolist() == ["ABC", "DEF"]
    assert column[1].tolist() == ["GHI"]


@pytest.mark.parametrize("datatype", DATATYPES)
def test_read_fixed_2d_array(datatype):
    """
    parse() of an arraysize="2x3x2" cell "aabbccddeeff" gives 2 rows of 3
    strings.
    """
    column = _column(_votable(datatype, "2x3x2", ["aabbccddeeff"]))

    assert column.tolist() == [[["aa", "bb", "cc"], ["dd", "ee", "ff"]]]


@pytest.mark.parametrize("datatype", DATATYPES)
def test_read_variable_2d_array(datatype):
    """
    parse() of arraysize="2x3x*" cells "aabbccddeeff" and "gghhii" gives 2
    rows of 3 strings, then 1 row of 3 strings.
    """
    column = _column(_votable(datatype, "2x3x*", ["aabbccddeeff", "gghhii"]))

    assert column[0].tolist() == [["aa", "bb", "cc"], ["dd", "ee", "ff"]]
    assert column[1].tolist() == [["gg", "hh", "ii"]]


def test_read_param():
    """
    parse() of a PARAM with arraysize="5x3" and value "     ab   cd" gives the
    strings "", "ab" and "cd".  This is the PARAM of issue 17098, which
    Topcat writes for FITS COMMENT cards.
    """
    data = b"""<?xml version="1.0" encoding="UTF-8"?>
<VOTABLE version="1.4" xmlns="http://www.ivoa.net/xml/VOTable/v1.3">
 <RESOURCE>
  <TABLE>
   <PARAM ID="p" name="p" datatype="char" arraysize="5x3" value="     ab   cd"/>
   <FIELD ID="c" name="c" datatype="int"/>
   <DATA><TABLEDATA><TR><TD>0</TD></TR></TABLEDATA></DATA>
  </TABLE>
 </RESOURCE>
</VOTABLE>"""
    param = parse(io.BytesIO(data)).get_first_table().params[0]

    assert param.value.tolist() == ["", "ab", "cd"]


def test_read_strips_padding():
    """
    parse() of an arraysize="5x3" cell "ab   cdefgh" gives "ab", "cdefg" and
    "h": the blanks that pad a string to 5 characters are not part of it.
    """
    column = _column(_votable("char", "5x3", ["ab   cdefgh    "]))

    assert column.tolist() == [["ab", "cdefg", "h"]]


def test_read_too_few_strings():
    """
    parse() of an arraysize="3x2" cell "abc" gives "abc" and "".  The XML
    parser strips the blanks of an empty last string, so this is not an
    error.
    """
    column = _column(_votable("char", "3x2", ["abc"]))

    assert column.tolist() == [["abc", ""]]


def test_read_empty_cell():
    """parse() of an empty arraysize="3x2" cell gives 2 empty strings."""
    column = _column(_votable("char", "3x2", [""]))

    assert column.tolist() == [["", ""]]


def test_read_empty_variable_cell():
    """parse() of an empty arraysize="3x*" cell gives an array of 0 strings."""
    column = _column(_votable("char", "3x*", [""]))

    assert column[0].tolist() == []


def test_read_too_many_strings_warns():
    """
    parse() with verify="warn" of an arraysize="3x2" cell with 3 strings warns
    E02 and keeps the first 2 strings.
    """
    with pytest.warns(E02, match="Expected multiple of 2, got 3"):
        column = _column(_votable("char", "3x2", ["abcdefghi"]), verify="warn")

    assert column.tolist() == [["abc", "def"]]


def test_read_too_many_strings_raises():
    """
    parse() with verify="exception" of an arraysize="3x2" cell with 3 strings
    raises E02.
    """
    with pytest.raises(E02, match="Expected multiple of 2, got 3"):
        _column(_votable("char", "3x2", ["abcdefghi"]), verify="exception")


@pytest.mark.parametrize("arraysize", ["3x", "0x2", "3x0"])
def test_invalid_arraysize(arraysize):
    """parse() of a char FIELD with a malformed arraysize raises E01."""
    with pytest.raises(E01, match="Invalid size specifier"):
        _column(_votable("char", arraysize, []), verify="ignore")


def test_read_binary_char():
    """
    parse() of a BINARY stream for arraysize="3x2" reads 3 bytes per string
    and stops each string at its first NUL: b"ab\\0cde" gives "ab" and "cde".
    """
    column = _column(_binary_votable("char", "3x2", b"ab\0cde"))

    assert column.tolist() == [["ab", "cde"]]


def test_read_binary_unicodechar():
    """
    parse() of a BINARY stream for unicodeChar arraysize="2x2" reads 2 UTF-16
    code units per string, including U+0100 whose low byte is 0, and an empty
    string.
    """
    stream = "a\u0100".encode("utf-16-be") + b"\0\0\0\0"
    column = _column(_binary_votable("unicodeChar", "2x2", stream))

    assert column.tolist() == [["a\u0100", ""]]


def test_read_binary_variable_array():
    """
    parse() of a BINARY stream for arraysize="3x*" reads a count of
    characters, not of strings: a count of 6 is followed by 2 strings, and a
    count of 0 by none.
    """
    stream = b"\0\0\0\x06abcdef" + b"\0\0\0\0"
    column = _column(_binary_votable("char", "3x*", stream))

    assert column[0].tolist() == ["abc", "def"]
    assert column[1].tolist() == []


def test_read_binary_variable_2d_array():
    """
    parse() of a BINARY stream for arraysize="2x3x*" with a count of 12
    characters gives 2 rows of 3 strings.
    """
    stream = b"\0\0\0\x0caabbccddeeff"
    column = _column(_binary_votable("char", "2x3x*", stream))

    assert column[0].tolist() == [["aa", "bb", "cc"], ["dd", "ee", "ff"]]


@pytest.mark.parametrize("tabledata_format", FORMATS)
@pytest.mark.parametrize("datatype", DATATYPES)
def test_roundtrip_fixed_array(datatype, tabledata_format):
    """
    A fixed array with a short string and an empty string is the same after
    it is written in each format and read again.
    """
    data = _votable(datatype, "5x3", ["ab   cdefg", "12345"])

    column = _column(_roundtrip(data, tabledata_format))

    assert column.tolist() == [["ab", "cdefg", ""], ["12345", "", ""]]


@pytest.mark.parametrize("tabledata_format", FORMATS)
@pytest.mark.parametrize("datatype", DATATYPES)
def test_roundtrip_variable_array(datatype, tabledata_format):
    """
    A variable array with a short string and an empty cell is the same after
    it is written in each format and read again.
    """
    data = _votable(datatype, "5x*", ["ab   cdefg", ""])

    column = _column(_roundtrip(data, tabledata_format))

    assert column[0].tolist() == ["ab", "cdefg"]
    assert column[1].tolist() == []


@pytest.mark.parametrize("tabledata_format", FORMATS)
@pytest.mark.parametrize("datatype", DATATYPES)
def test_roundtrip_variable_2d_array(datatype, tabledata_format):
    """
    A variable array of string arrays is the same after it is written in each
    format and read again.
    """
    data = _votable(datatype, "2x3x*", ["aabbccddeeff", "gghhii"])

    column = _column(_roundtrip(data, tabledata_format))

    assert column[0].tolist() == [["aa", "bb", "cc"], ["dd", "ee", "ff"]]
    assert column[1].tolist() == [["gg", "hh", "ii"]]


def test_write_pads_strings():
    """
    to_xml() pads each string of an arraysize="3x2" cell to 3 characters, and
    escapes after padding: "a&" and "cde" are written as "a&amp; cde".
    """
    data = _votable("char", "3x2", ["a&amp; cde"])

    written = _roundtrip(data, "tabledata")

    assert b"<TD>a&amp; cde</TD>" in written


def test_write_masked_string():
    """
    to_xml() writes a masked string of an arraysize="3x3" cell as blanks, so
    the string after it stays in place.
    """
    votable = parse(io.BytesIO(_votable("char", "3x3", ["abcdefghi"])))
    votable.get_first_table().array.mask["c"][0][1] = True
    written = io.BytesIO()
    votable.to_xml(written)

    assert b"<TD>abc   ghi</TD>" in written.getvalue()
    assert _column(written.getvalue()).tolist() == [["abc", "", "ghi"]]


def test_write_masked_string_in_variable_array():
    """
    to_xml() writes a masked string of an arraysize="3x*" cell as blanks, as
    it does for a fixed-size array.
    """
    votable = parse(io.BytesIO(_votable("char", "3x*", ["abcdefghi"])))
    votable.get_first_table().array["c"][0][1] = np.ma.masked
    written = io.BytesIO()
    votable.to_xml(written)

    assert b"<TD>abc   ghi</TD>" in written.getvalue()
    assert _column(written.getvalue())[0].tolist() == ["abc", "", "ghi"]


def test_write_too_long_string_warns():
    """
    to_xml() of a 4-character string in an arraysize="3x2" FIELD warns W46
    and writes the first 3 characters.
    """
    votable = parse(io.BytesIO(_votable("char", "3x2", ["abcdef"])))
    table = votable.get_first_table()
    table.array = np.ma.array([(["abcd", "ef"],)], dtype=[("c", "U4", (2,))])
    written = io.BytesIO()

    with pytest.warns(W46, match="char value is too long for specified length of 3"):
        votable.to_xml(written)

    assert b"<TD>abcef </TD>" in written.getvalue()


@pytest.mark.parametrize("tabledata_format", FORMATS)
def test_table_write_string_array_column(tabledata_format):
    """
    Table.write() of a column of 2 strings per row writes a unicodeChar FIELD
    with arraysize="3x2", which reads back as the same strings.
    """
    table = Table({"c": np.array([["abc", "de"], ["f", "ghi"]])})
    written = io.BytesIO()
    table.write(written, format="votable", tabledata_format=tabledata_format)

    votable = parse(io.BytesIO(written.getvalue()))
    field = votable.get_first_table().fields[0]

    assert field.datatype == "unicodeChar"
    assert field.arraysize == "3x2"
    assert_array_equal(votable.get_first_table().array["c"], table["c"])


def test_table_write_2d_string_array_column():
    """
    Table.write() of a column with 2 rows of 3 strings per row writes
    arraysize="2x3x2", which reads back with the same shape.
    """
    table = Table({"c": np.array([[["aa", "bb", "cc"], ["dd", "ee", "ff"]]])})
    written = io.BytesIO()
    table.write(written, format="votable")

    votable = parse(io.BytesIO(written.getvalue()))

    assert votable.get_first_table().fields[0].arraysize == "2x3x2"
    assert_array_equal(votable.get_first_table().array["c"], table["c"])


@pytest.mark.parametrize("arraysize", ["3x2", "3x*", "3x2x*"])
def test_table_read_write_keeps_arraysize(arraysize):
    """
    Table.read() then Table.write() of a char string array FIELD writes a
    char FIELD with the same arraysize.
    """
    table = Table.read(io.BytesIO(_votable("char", arraysize, ["abcdef"])))
    written = io.BytesIO()
    table.write(written, format="votable")

    field = parse(io.BytesIO(written.getvalue())).get_first_table().fields[0]

    assert field.datatype == "char"
    assert field.arraysize == arraysize


def _v16(data):
    return data.replace(b'version="1.4"', b'version="1.6"')


def test_utf8_read_splits_by_bytes():
    # From 1.6 each string takes 3 bytes: "Λ" is 2 bytes plus one blank
    column = _column(_v16(_votable("char", "3x2", ["Λ ab"])))

    assert column.tolist() == [["Λ", "ab"]]


@pytest.mark.parametrize("tabledata_format", FORMATS)
def test_utf8_roundtrip(tabledata_format):
    t = Table({"c": np.array([["Λ", "ab"], ["ΛΛ", "c"]], dtype="U2")})
    written = io.BytesIO()
    t.write(written, format="votable", tabledata_format=tabledata_format, version="1.6")
    table = parse(io.BytesIO(written.getvalue())).get_first_table()

    assert table.fields[0].arraysize == "4x2"
    assert table.fields[0].width == 2
    assert table.array["c"].tolist() == [["Λ", "ab"], ["ΛΛ", "c"]]


@pytest.mark.parametrize("tabledata_format", FORMATS)
def test_utf8_roundtrip_variable_array(tabledata_format):
    # "ΛCDM" is 4 characters but fills all 5 bytes
    data = _v16(_votable("char", "5x*", ["ΛCDMab", ""]))

    column = _column(_roundtrip(data, tabledata_format))

    assert column[0].tolist() == ["ΛCDM", "ab"]
    assert column[1].tolist() == []


@pytest.mark.parametrize("tabledata_format", FORMATS)
def test_unicodechar_variable_array_written_as_utf8_char(tabledata_format):
    data = _votable("unicodeChar", "3x*", ["ΛΛΛab"])
    t = Table.read(io.BytesIO(data), format="votable")
    written = io.BytesIO()
    t.write(written, format="votable", tabledata_format=tabledata_format, version="1.6")
    table = parse(io.BytesIO(written.getvalue())).get_first_table()

    # 3 characters of unicodeChar become 6 bytes of char
    assert table.fields[0].datatype == "char"
    assert table.fields[0].arraysize == "6x*"
    assert table.array["c"][0].tolist() == ["ΛΛΛ", "ab"]


def test_utf8_output_truncates_by_bytes():
    config = {"version_1_6_or_later": True}
    field = tree.Field(None, name="c", datatype="char", arraysize="3x2", config=config)
    converter = get_converter(field, config=config)

    # "ΛΛ" is 4 bytes
    with pytest.warns(W46):
        assert converter.output(np.array(["ΛΛ", "a"]), False) == "Λ a  "


def test_utf8_datatype_of_string_array_columns():
    col = Column(
        np.array([["ab", "cd"]], dtype="U3"),
        meta={"_votable_string_dtype": "char", "_votable_arraysize": "3x2x*"},
    )
    assert table_column_to_votable_datatype(col, version="1.6") == {
        "datatype": "char",
        "arraysize": "3x2",
        "width": 3,
    }

    col = Column(
        [np.array(["ab", "cd"]), np.array(["ΛΛ", "e"])],
        dtype=object,
        meta={"_votable_string_dtype": "unicodeChar", "_votable_arraysize": "3x*"},
    )
    assert table_column_to_votable_datatype(col, version="1.6") == {
        "datatype": "char",
        "arraysize": "4x*",
    }


def test_widened_string_array_not_truncated():
    col = Column(
        np.array([["abcdefgh", "cd"]], dtype="U8"),
        meta={"_votable_string_dtype": "char", "_votable_arraysize": "3x2"},
    )
    assert table_column_to_votable_datatype(col) == {
        "datatype": "char",
        "arraysize": "8x2",
    }
