# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module tests some methods related to ``tdat`` format
reader/writer.
Requires `pyyaml <https://pyyaml.org/>`_ to be installed.
"""

import copy
import warnings
from collections import OrderedDict
from io import StringIO

import numpy as np
import pytest

from astropy.io import ascii
from astropy.io.ascii.tdat import TdatFormatError, TdatFormatWarning
from astropy.table import Table
from astropy.table.table_helpers import simple_table
from astropy.units import allclose as quantity_allclose

test_dat = [
    "<HEADER>",
    "#",
    "# TABLE: heasarc_simple",
    "# TOTAL ROWS: 7",
    "#",
    "table_name = heasarc_simple",
    'table_description = "Test table"',
    "table_security = public",
    "#",
    "# Table Parameters",
    "#",
    "field[record_number] = int4  [meta.id] (key) // Unique Identifier for Entry",
    "field[id] = int4  [meta.id] (index) // Source ID Number",
    "field[name] = char12  [meta.id;meta.main] (index) // String Name",
    "field[ra] = float8:.4f_degree [pos.eq.ra] (index) // Right Ascension",
    "field[dec] = float8:.4f_degree [pos.eq.dec] (index) // Declination",
    "field[empty] = float8:.4f // Empty // Comment",
    "#",
    "parameter_defaults = name ra dec",
    "#",
    "# Virtual Parameters",
    "#",
    "frequency_regime = Gamma-ray",
    "observatory_name = GAMMA-RAY BURSTS",
    "row_type = GRB",
    "table_author = Example et al.",
    "table_priority = 3.01",
    "table_type = Observation",
    "unique_key = record_number",
    "#",
    "# Data Format Specification",
    "#",
    "line[1] = record_number id name ra dec empty",
    "#",
    "<DATA>",
    "1|10|aaa|1.0|1.0||",
    "2|20|b|2.0|||",
    "3|30|c||3.0||",
    "4|20|||||",
    "5||||||",
    "|60|f|6.0|6.0||",
    "7| 70| g | 7.0 |7.0||",
    "<END>",
]
test_table = Table.read(test_dat, format="ascii.tdat")
# Corresponds to simple_table()
SIMPLE_LINES = [
    "<HEADER>",
    "table_name = astropy_table",
    "#",
    "# Table Parameters",
    "#",
    "field[a] = int4",
    "field[b] = float8",
    "field[c] = char1",
    "#",
    "# Data Format Specification",
    "#",
    "line[1] = a b c",
    "#",
    "<DATA>",
    "1|1.0|c|",
    "2|2.0|d|",
    "3|3.0|e|",
    "<END>",
]


def test_write_simple():
    """
    Write a simple table with common types.  This shows the compact version
    of serialization with one line per column.
    Table -> TDAT
    """
    t = simple_table()
    t.meta["table_name"] = "astropy_table"
    out = StringIO()
    t.write(out, format="ascii.tdat")

    # Explicit check that output matches expected
    assert out.getvalue().splitlines() == SIMPLE_LINES


def test_catch_format():
    """
    Ensure that a table handles header keywords appropriately.
    Table -> TDAT
    """
    t = simple_table()
    out = StringIO()

    # Missing table_name
    with pytest.warns(TdatFormatWarning, match="'table_name' must be specified"):
        t.write(out, format="ascii.tdat")

    # table_name too long (>20 characters)
    with pytest.warns(TdatFormatWarning, match="'table_name' is too long"):
        t.meta["table_name"] = "origin_12345678901234567890"
        t.write(out, format="ascii.tdat")

        # check if the name is truncated
        nt = Table.read(out.getvalue(), format="ascii.tdat")
        assert nt.meta["keywords"]["table_name"] == t.meta["table_name"][:20]

    t.meta["table_name"] = "astropy_table"

    # table_description too long (>80 characters)
    with pytest.warns(TdatFormatWarning, match="'table_description' is too long"):
        t.meta["table_description"] = (
            """\
        This is a description that exceeds the character limit allowed by the tdat
        format and it should be truncated before being written. A warning should
        pop up to inform the user of this behavior.
        """.replace("\n", " ")
            .replace("\t", "")
            .strip()
        )
        t.write(out, format="ascii.tdat")

        # check if the description is truncated
        nt = Table.read(out.getvalue(), format="ascii.tdat")
        assert (
            nt.meta["keywords"]["table_description"]
            == t.meta["table_description"][:80].strip()
        )


def test_read_tdat():
    """Ensure the Table is as it should be
    TDAT -> Table
    """
    assert test_table.meta["keywords"] == OrderedDict(
        [
            ("table_name", "heasarc_simple"),
            ("table_description", "Test table"),
            ("table_security", "public"),
            ("parameter_defaults", "name ra dec"),
            ("frequency_regime", "Gamma-ray"),
            ("observatory_name", "GAMMA-RAY BURSTS"),
            ("row_type", "GRB"),
            ("table_author", "Example et al."),
            ("table_priority", "3.01"),
            ("table_type", "Observation"),
            ("unique_key", "record_number"),
        ]
    )

    # Column checks
    descriptions = [
        "Unique Identifier for Entry",
        "Source ID Number",
        "String Name",
        "Right Ascension",
        "Declination",
        "Empty",
    ]
    dtypes = [int, int, "<U3", float, float, float]
    units = [None, None, None, "deg", "deg", None]
    meta = [
        {"ucd": "meta.id", "index": "key"},
        {"ucd": "meta.id", "index": "index"},
        {"ucd": "meta.id;meta.main", "index": "index"},
        {"ucd": "pos.eq.ra", "index": "index"},
        {"ucd": "pos.eq.dec", "index": "index"},
        {"comment": "Comment"},
    ]
    for i, col in enumerate(test_table.itercols()):
        assert col.description == descriptions[i]
        assert col.dtype == dtypes[i]
        assert col.unit == units[i]
        assert col.meta == meta[i]

    # data check
    ## Missing is masked
    assert isinstance(test_table["ra"][3], np.ma.core.MaskedConstant)
    ## Check data matches simple csv format
    test_data = [
        "record_number, id, name, ra, dec, empty",
        "1, 10, aaa, 1.0, 1.0, ",
        "2, 20,   b, 2.0,    , ",
        "3, 30,   c,    , 3.0, ",
        "4, 20,    ,    ,    , ",
        "5,   ,    ,    ,    , ",
        " , 60,   f, 6.0, 6.0, ",
        "7, 70,   g, 7.0, 7.0, ",
    ]
    table_data = Table.read(test_data, format="csv")
    assert all((table_data == test_table).data)


def test_full_table_content():
    """Check the table content matches expectation.
    TDAT -> Table
    """
    assert test_table.meta == {
        "comments": ["TABLE: heasarc_simple", "TOTAL ROWS: 7"],
        "keywords": OrderedDict(
            [
                ("table_name", "heasarc_simple"),
                ("table_description", "Test table"),
                ("table_security", "public"),
                ("parameter_defaults", "name ra dec"),
                ("frequency_regime", "Gamma-ray"),
                ("observatory_name", "GAMMA-RAY BURSTS"),
                ("row_type", "GRB"),
                ("table_author", "Example et al."),
                ("table_priority", "3.01"),
                ("table_type", "Observation"),
                ("unique_key", "record_number"),
            ]
        ),
    }
    assert len(test_table) == 7
    dtypes = [int, int, "<U3", float, float, float]
    descriptions = [
        "Unique Identifier for Entry",
        "Source ID Number",
        "String Name",
        "Right Ascension",
        "Declination",
        "Empty",
    ]
    ucds = ["meta.id", "meta.id", "meta.id;meta.main", "pos.eq.ra", "pos.eq.dec", None]
    indices = ["key", "index", "index", "index", "index", None]
    for i, col in enumerate(test_table.columns):
        assert test_table[col].dtype == dtypes[i]
        assert test_table[col].description == descriptions[i]
        assert test_table[col].meta.get("ucd", None) == ucds[i]
        assert test_table[col].meta.get("index", None) == indices[i]
    assert test_table["ra"].unit == "deg"
    assert test_table["dec"].unit == "deg"
    assert test_table["ra"].format == ".4f"
    assert test_table["dec"].format == ".4f"
    assert test_table["empty"].format == ".4f"


def test_write_full():
    """
    Write a full-featured table with common types and explicitly check output
    TDAT -> Table -> TDAT (formatted)

    Differences between `lines` and `test_dat`:
    - Empty comment lines are dropped (except when demarcating a section header
    like "Table Parameters").
    - The data type for string field "name" is downsized to char3 from char12
    reflecting the actual maximum size in the column.
    - Extraneous spaces in the data are stripped
    These differences reflect flexibility in reading in from a tdat file, and
    writing out in a standardized way.
    """
    t = test_table

    lines = [
        "<HEADER>",
        "# TABLE: heasarc_simple",
        "# TOTAL ROWS: 7",
        "table_name = heasarc_simple",
        "table_description = Test table",
        "table_security = public",
        "#",
        "# Table Parameters",
        "#",
        "field[record_number] = int4 [meta.id] (key) // Unique Identifier for Entry",
        "field[id] = int4 [meta.id] (index) // Source ID Number",
        "field[name] = char3 [meta.id;meta.main] (index) // String Name",
        "field[ra] = float8:.4f_deg [pos.eq.ra] (index) // Right Ascension",
        "field[dec] = float8:.4f_deg [pos.eq.dec] (index) // Declination",
        "field[empty] = float8:.4f // Empty // Comment",
        "#",
        "parameter_defaults = name ra dec",
        "#",
        "# Virtual Parameters",
        "#",
        "frequency_regime = Gamma-ray",
        "observatory_name = GAMMA-RAY BURSTS",
        "row_type = GRB",
        "table_author = Example et al.",
        "table_priority = 3.01",
        "table_type = Observation",
        "unique_key = record_number",
        "#",
        "# Data Format Specification",
        "#",
        "line[1] = record_number id name ra dec empty",
        "#",
        "<DATA>",
        "1|10|aaa|1.0000|1.0000||",
        "2|20|b|2.0000|||",
        "3|30|c||3.0000||",
        "4|20|||||",
        "5||||||",
        "|60|f|6.0000|6.0000||",
        "7|70|g|7.0000|7.0000||",
        "<END>",
    ]

    out = StringIO()
    t.write(out, format="ascii.tdat")
    assert out.getvalue().splitlines() == lines


def test_write_read_roundtrip():
    """
    Write a full-featured table with all types and see that it round-trips on
    readback.
    TDAT -> Table -> TDAT (formatted) -> Table
    Tables should have the same content.
    """
    t = test_table
    out = StringIO()
    t.write(out, format="ascii.tdat")

    t2s = [
        Table.read(out.getvalue(), format="ascii.tdat"),
        ascii.read(out.getvalue(), format="tdat"),
    ]
    for t2 in t2s:
        assert t.meta == t2.meta
        for name in t.colnames:
            assert t[name].attrs_equal(t2[name])
            assert np.all(t[name] == t2[name]) or np.all(t[name].mask == t2[name].mask)


def test_write_read_roundtrip_empty_table(tmp_path):
    """see https://github.com/astropy/astropy/issues/13191"""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=TdatFormatWarning)
        sfile = tmp_path / "x.tdat"
        Table().write(sfile)
        t = Table.read(sfile)
        assert len(t) == 0
        assert len(t.colnames) == 0


def test_keyword_quotes():
    lines = copy.copy(SIMPLE_LINES)
    # double quotes
    lines[1] = 'table_name = "astropy_table"'
    t = Table.read(lines, format="ascii.tdat")
    assert t.meta["keywords"]["table_name"] == "astropy_table"
    # single quotes
    lines[1] = "table_name = 'astropy_table'"
    t = Table.read(lines, format="ascii.tdat")
    assert t.meta["keywords"]["table_name"] == "astropy_table"
    # back quotes
    lines[1] = "table_name = `astropy_table`"
    t = Table.read(lines, format="ascii.tdat")
    assert t.meta["keywords"]["table_name"] == "astropy_table"
    # combination and multiple, nested properly
    lines[1] = "table_name = \"'`astropy_table`'\""
    t = Table.read(lines, format="ascii.tdat")
    assert t.meta["keywords"]["table_name"] == "astropy_table"

    # mismatched
    lines[1] = "table_name = \"astropy_table'"
    with pytest.raises(TdatFormatError, match="Mismatched"):
        t = Table.read(lines, format="ascii.tdat")
    # combination, nested improperly
    lines[1] = "table_name = \"'astropy_table\"'"
    with pytest.raises(TdatFormatError, match="Mismatched"):
        t = Table.read(lines, format="ascii.tdat")


def test_tablenames():
    """_summary_"""
    lines = copy.copy(SIMPLE_LINES)
    t = Table.read(lines, format="ascii.tdat")
    out = StringIO()

    with pytest.warns(TdatFormatWarning, match="'table_name' is too long"):
        t.meta["table_name"] = "origin_12345678901234567890"
        t.write(out, format="ascii.tdat")


def test_bad_delimiter():
    """
    Passing a delimiter other than | (pipe) gives an exception
    """
    out = StringIO()

    with pytest.warns(
        TdatFormatWarning, match="Delimiters other than the pipe character"
    ) as w:
        test_table.write(out, format="ascii.tdat", delimiter=",")


def test_bad_header_start():
    """
    Bad header without initial <HEADER>
    """
    lines = copy.copy(SIMPLE_LINES)
    lines[0] = "<DATA>"
    with pytest.raises(TdatFormatError) as err:
        Table.read("\n".join(lines), format="ascii.tdat")
        assert "<HEADER> not found in file." in str(err.value)


def test_bad_data_heading():
    """
    No <DATA> heading
    """
    lines = copy.copy(SIMPLE_LINES)
    lines[13] = "<DAYA>"
    with pytest.raises(ascii.tdat.TdatFormatError) as err:
        Table.read("\n".join(lines), format="ascii.tdat")
    assert "<DATA> not found in file." in str(err.value)


def test_bad_end_heading():
    """
    No <END> heading
    """
    lines = copy.copy(SIMPLE_LINES)
    lines[-1] = "<That's all folks>"
    Table.read("\n".join(lines), format="ascii.tdat")


def test_unrecognized_dtype():
    """Not all dtypes are supported by tdat files"""
    lines = copy.copy(SIMPLE_LINES)
    lines[5] = "field[a] = complex"
    with pytest.raises(TdatFormatError) as err:
        Table.read("\n".join(lines), format="ascii.tdat")
    assert "Unrecognized data type" in str(err.value)


def test_mismatch_line_field():
    """Not all dtypes are supported by tdat files"""
    lines = copy.copy(SIMPLE_LINES)
    lines[11] = "line[1] = a b c d"
    with pytest.raises(TdatFormatError) as err:
        Table.read("\n".join(lines), format="ascii.tdat")
    assert 'The columns "field" descriptors are not consistent' in str(err.value)


def test_fmt_type_too_long():
    """The combination of type and format has a maximum character length of 24"""
    lines = copy.copy(SIMPLE_LINES)
    lines[6] = "field[b] = float8:.0000000000000000000000000000001f"
    with pytest.raises(TdatFormatError) as err:
        Table.read("\n".join(lines), format="ascii.tdat")
    assert "The type:fmt specifier" in str(err.value)


def assert_objects_equal(obj1, obj2, attrs, compare_class=True):
    if compare_class:
        assert obj1.__class__ is obj2.__class__

    assert obj1.shape == obj2.shape

    info_attrs = [
        "info.name",
        "info.format",
        "info.unit",
        "info.description",
        "info.dtype",
    ]
    for attr in attrs + info_attrs:
        a1 = obj1
        a2 = obj2
        for subattr in attr.split("."):
            try:
                a1 = getattr(a1, subattr)
                a2 = getattr(a2, subattr)
            except AttributeError:
                a1 = a1[subattr]
                a2 = a2[subattr]

        if isinstance(a1, np.ndarray) and a1.dtype.kind == "f":
            assert quantity_allclose(a1, a2, rtol=1e-10)
        else:
            assert np.all(a1 == a2)

    # For no attrs that means we just compare directly.
    if not attrs:
        if isinstance(obj1, np.ndarray) and obj1.dtype.kind == "f":
            assert quantity_allclose(obj1, obj2, rtol=1e-15)
        else:
            assert np.all(obj1 == obj2)


def test_round_trip_masked_table_default(tmp_path):
    """Test (mostly) round-trip of MaskedColumn through tdat using default serialization
    that uses an empty string "" to mark NULL values.  Note:

    >>> simple_table(masked=True)
    <Table masked=True length=3>
      a      b     c
    int64 float64 str1
    ----- ------- ----
       --     1.0    c
        2     2.0   --
        3      --    e
    """
    filename = tmp_path / "x.tdat"

    t = simple_table(masked=True)  # int, float, and str cols with one masked element
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=TdatFormatWarning)
        t.write(filename)

        t2 = Table.read(filename)
        assert t2.masked is False
        assert t2.colnames == t.colnames
        for name in t2.colnames:
            # From formal perspective the round-trip columns are the "same"
            assert np.all(t2[name].mask == t[name].mask)
            assert np.all(t2[name] == t[name])

            # But peeking under the mask shows that the underlying data are changed
            # because by default ECSV uses "" to represent masked elements.
            t[name].mask = False
            t2[name].mask = False
            assert not np.all(t2[name] == t[name])  # Expected diff


def test_deprecated_keyword():
    """Deprecated and obsolete keywords should raise warnings"""
    test_data = test_dat.copy()
    data_lines = test_data[35:-1]
    data_lines = "\123".join(data_lines)  # now a single string with a delimiter
    data_lines = data_lines.replace("|", "!")
    del test_data[35:-1]
    test_data.insert(35, data_lines)
    test_data.insert(8, 'record_delimiter = "\123"')
    test_data.insert(8, 'field_delimiter = "!"')
    with pytest.warns(TdatFormatWarning, match="keyword is deprecated"):
        t = Table.read(test_data, format="ascii.tdat")

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=TdatFormatWarning)
        t = Table.read(test_data, format="ascii.tdat")
        assert all(t == test_table)

    test_data = test_dat.copy()
    test_data.insert(8, "relate[ra] = dec")
    with pytest.warns(TdatFormatWarning, match="keyword is obsolete"):
        t = Table.read(test_data, format="ascii.tdat")


def test_tdat_format_error():
    """Test the basic error message"""
    _STD_MSG = (
        "See details in https://heasarc.gsfc.nasa.gov/docs/software/dbdocs/tdat.html"
    )
    error_msg = "Test error message"
    with pytest.raises(TdatFormatError) as exc_info:
        raise TdatFormatError(error_msg)

    assert str(exc_info.value).startswith(error_msg)

    # Test that the standard message is appended
    assert _STD_MSG in str(exc_info.value)

    # Test with empty error message
    with pytest.raises(TdatFormatError) as exc_info:
        raise TdatFormatError()

    assert str(exc_info.value).startswith("\nSee details")

    # Test with a different error message
    another_msg = "Another test message"
    with pytest.raises(TdatFormatError) as exc_info:
        raise TdatFormatError(another_msg)

    assert str(exc_info.value).startswith(another_msg)
    assert _STD_MSG in str(exc_info.value)


def test_delimiter_in_data():
    """Test escaped delimiters in data"""
    lines = copy.deepcopy(test_dat)
    lines[38] = "4|20|\\|||||"
    t = Table.read(lines, format="ascii.tdat")
    assert t["name"][3] == "|"

    lines.insert(8, 'field_delimiter = "|$"')
    lines[39] = "4$20|\\|\\$||$|"
    with pytest.warns(TdatFormatWarning, match="keyword is deprecated"):
        t = Table.read(lines, format="ascii.tdat")
        assert t["name"][3] == "|$"
        out = StringIO()
        t.write(out, format="ascii.tdat")
        assert out.getvalue().split("\n")[36] == "4|20|\\|$||||"
