import contextlib
import datetime
import decimal
import io
import textwrap
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from astropy.io.misc.pyarrow.csv import strip_comment_lines
from astropy.table import Table
from astropy.utils.compat.optional_deps import HAS_PYARROW

if HAS_PYARROW:
    import pyarrow as pa
else:
    pytest.skip("pyarrow is not available", allow_module_level=True)


def check_tables_equal(t1, t2):
    """Check that two tables are equal.

    This function checks that the two tables have the same column names,
    lengths, and that the data in each column is equal. It also checks
    that the masks (if any) are equal.
    """
    assert t1.colnames == t2.colnames
    assert len(t1) == len(t2)
    for col1, col2 in zip(t1.itercols(), t2.itercols()):
        assert col1.name == col2.name
        assert col1.dtype == col2.dtype
        assert_array_equal(col1, col2)
        has_mask1 = hasattr(col1, "mask")
        has_mask2 = hasattr(col2, "mask")
        assert has_mask1 == has_mask2
        if has_mask1:
            assert_array_equal(col1.mask, col2.mask)


def convert_table_to_text(tbl, delimiter=",", **kwargs):
    """Convert a table to text using the io.ascii CSV format.

    This function is used to convert a table to text for testing purposes.
    It uses the `ascii.csv` format to write the table to a string buffer.
    """
    text_io = io.StringIO()
    tbl.write(text_io, format="ascii.csv", delimiter=delimiter, **kwargs)
    return text_io.getvalue()


@pytest.fixture(scope="module")
def tbl():
    """Masked table with bool, int, float, and non-ASCII string types."""
    return Table(
        rows=[
            [np.ma.masked, 1, np.ma.masked, np.ma.masked],
            [False, np.ma.masked, 2.5, 'bår q"ux'],
            [True, 2, 3.5, "bazœo"],
        ],
        names=["a", "b", "c 2", "d"],
    )


@pytest.fixture(scope="module")
def tbl_text(tbl):
    """Default CSV text representation for simple table above.

    This is equivalent to::

    '''
    a,b,ç 2,d
    ,,,
    False,1,2.5,"bår q""ux"
    True,2,3.5,bazœo
    '''
    """
    return convert_table_to_text(tbl)


@contextlib.contextmanager
def get_input_file(
    text: str,
    input_type: str,
    encoding: str | None = None,
    tmp_path: Path | None = None,
):
    """
    Generate an input file or stream based on the specified type and encoding.

    Parameters
    ----------
    text : str
        The text content to be written to the input file or stream.
    input_type : str
        The type of input to generate. Supported values are: - "str": A temporary file
        path as a string. - "path": A temporary file path as a `Path` object. -
        "bytesio": A `BytesIO` stream.
    encoding : str | None, optional
        The encoding to use when writing the text. If None, no encoding is applied.
    tmp_path : Path | None, optional
        The temporary path to use for the input file for "str" or "path" input type.

    Yields
    ------
    str | Path | io.BytesIO
        The generated input file or stream based on the specified `input_type`.

    Notes
    -----
    - For "str" and "path" input types, a temporary file in ``tmp_path`` is created and
      its content is written with the specified encoding (if provided).
    - For "bytesio", the text is encoded into bytes and returned as a `BytesIO` stream.
    - The temporary file is automatically cleaned up after use.
    """
    encoding = encoding or "utf-8"
    if input_type in ("str", "path"):
        path = tmp_path / "test.csv"
        path.write_text(text, encoding=encoding)
        yield str(path) if input_type == "str" else path
    elif input_type == "bytesio":
        yield io.BytesIO(text.encode(encoding))
    else:
        raise ValueError(f"Unknown input_type: {input_type}")


def table_read_csv(
    text: str,
    input_type: str = "bytesio",
    encoding: str | None = None,
    tmp_path: Path | None = None,
    **kwargs,
):
    """Read ``text`` using ``Table.read`` with format="pyarrow.csv".

    The ``input_type`` parameter determines how the text is passed to ``Table.read``.
    For "str" or "path", a named temporary file is created and that name ("str") or path
    ("path") is passed to ``Table.read``. For "bytesio", a ``BytesIO`` object is created
    and passed to ``Table.read``. The text is encoded using the specified encoding (if
    any) before being passed to ``BytesIO``. The default is "bytesio".

    The ``encoding`` parameter is used only if ``input_type`` is "bytesio". If
    ``encoding`` is None, the text is encoded using the default encoding.

    Parameters
    ----------
    text : str
        The text to read.
    input_type : str
        The type of input for Table.read(). One of "str", "path", or "bytesio".
    encoding : str | None
        The encoding to use when reading the text. If None, no encoding is applied.
    tmp_path : Path | None
        The temporary path to use for the input file for input_type="str" or "path".
    **kwargs : dict
        Additional keyword arguments to pass to Table.read().
    """
    if encoding is not None:
        kwargs["encoding"] = encoding

    with get_input_file(text, input_type, encoding, tmp_path) as input_file:
        out = Table.read(input_file, format="pyarrow.csv", **kwargs)
    return out


@pytest.mark.parametrize("input_type", ["str", "path", "bytesio"])
@pytest.mark.parametrize("encoding", [None, "utf-8", "utf-16"])
def test_read_tbl_simple_input_type_encoding(
    input_type, encoding, tbl, tbl_text, tmp_path
):
    """Test reading a simple CSV file with different input types.

    This tests:
    - input_file : str, PathLike, or binary file-like object
    - encoding : None, "utf-8", or "utf-16"
    """
    out = table_read_csv(tbl_text, input_type, encoding, tmp_path=tmp_path)
    check_tables_equal(tbl, out)


@pytest.mark.parametrize("delimiter", ["|", "\t", " "])
def test_read_delimiter(tbl, delimiter):
    """Test reading a simple CSV file with different delimiters.

    This tests:
    - delimiter : "|", "\t", and " "
    """
    tbl_text = convert_table_to_text(tbl, delimiter=delimiter)
    out = table_read_csv(tbl_text, delimiter=delimiter)
    check_tables_equal(tbl, out)


def test_read_dtypes(tbl_text):
    """Test reading a simple CSV file with different input types.

    This tests:
    - dtypes: dict of types, with a sampling of types including numpy and pyarrow
    """
    tbl_text = textwrap.dedent(
        """
        a,b,c,d,e,f,time
        0,1,2,3,0,5,12:34:56.123456
        ,,,,,,
        """
    )

    dtypes = {
        "a": "str",
        "b": pa.decimal128(precision=5, scale=2),
        "c": np.uint8,
        "d": pa.float32(),
        "e": "bool",
        "f": pa.utf8(),
        "time": pa.time64("us"),
    }
    out = table_read_csv(tbl_text, dtypes=dtypes)
    assert out["a"].dtype == "U1"
    assert out["b"].dtype == "object"
    assert isinstance(out["b"][0], decimal.Decimal)
    assert out["c"].dtype == "uint8"
    assert out["d"].dtype == "float32"
    assert out["e"].dtype == "bool"
    assert out["f"].dtype == "U1"
    assert out["time"].dtype == "object"
    assert isinstance(out["time"][0], datetime.time)

    # Unmask the table and check the zero values
    for col in out.itercols():
        col.mask[:] = False
    assert out["a"][1] == ""
    assert out["b"][1] == decimal.Decimal(0)
    assert out["c"][1] == 0
    assert out["d"][1] == 0.0
    assert out["e"][1] == False  # noqa: E712
    assert out["f"][1] == ""
    assert out["time"][1] == datetime.time(0, 0, 0)


def test_read_dtypes_invalid_conversion(tbl_text):
    """Test reading with an dtype that is inconsistent with the file data.

    This tests:
    - dtypes: dict of types, with invalid type for conversion provided
    """
    dtypes = {"a": "int"}  # True/False are not convertible to int
    with pytest.raises(pa.lib.ArrowInvalid, match="In CSV column #0"):
        table_read_csv(tbl_text, dtypes=dtypes)


def test_read_dtypes_invalid_dtype(tbl_text):
    """Test reading with an invalid dtype.

    This tests:
    - dtypes: dict of types, with invalid type for conversion provided
    """
    dtypes = {"a": "asdf"}
    with pytest.raises(TypeError, match="data type 'asdf' not understood"):
        table_read_csv(tbl_text, dtypes=dtypes)


def test_read_quotechar():
    """Test reading a simple CSV file with a single quote quotechar.

    This tests:
    - quotechar : single quote and False
    """
    tbl_text = textwrap.dedent("""
    a,b,c
    0,'''',2.5
    ","2",'3.5'
    """)

    out = table_read_csv(tbl_text, quotechar="'")
    exp = [
        " a    b      c   ",
        "str1 str3 float64",
        "---- ---- -------",
        "   0    '     2.5",  # '''' turns into a single quote
        '   "  "2"     3.5',  # '3.5' is still float
    ]
    assert out.pformat(show_dtype=True) == exp

    out = table_read_csv(tbl_text, quotechar=False)
    # No quoting conversions with quotechar=False
    exp = [
        " a    b     c  ",
        "str1 str4  str5",
        "---- ---- -----",
        "   0 ''''   2.5",
        '   "  "2" \'3.5\'',
    ]
    assert out.pformat(show_dtype=True) == exp


def test_read_doublequote():
    """Test reading a simple CSV file with doublequote.

    This tests:
    - doublequote : False (default True is tested above)
    """
    tbl_text = textwrap.dedent('''
    a,b,c
    0,"""",2.5
    ''')
    out = table_read_csv(tbl_text, doublequote=False)
    exp = [
        "  a    b      c   ",
        "int64 str2 float64",
        "----- ---- -------",
        '    0   ""     2.5',
    ]
    assert out.pformat(show_dtype=True) == exp


@pytest.mark.parametrize("ec", ["\\", "/"])
def test_read_escapechar(ec):
    r"""Test reading a simple CSV file with escapechar.

    This tests:
    - escapechar : "\" and "/"
    - Blank lines are ignored (default behavior, tested elsewhere as well)
    """
    tbl_text = textwrap.dedent(f"""
    a,b,c
    0,{ec},,{ec}"{ec}'
    """)
    out = table_read_csv(tbl_text, escapechar=ec)
    exp = [
        "  a    b    c  ",
        "int64 str1 str2",
        "----- ---- ----",
        "    0    ,   \"'",
    ]
    assert out.pformat(show_dtype=True) == exp


def test_read_header_start_1():
    """Test reading a simple CSV file with header_start.

    This tests:
    - header_start : 1
    """
    tbl_text = "a,b,c\n0,1,2\n3,4,5"
    out = table_read_csv(tbl_text, header_start=1)
    exp = [
        "  0     1     2  ",
        "int64 int64 int64",
        "----- ----- -----",
        "    3     4     5",
    ]
    assert out.pformat(show_dtype=True) == exp


def test_read_header_start_none():
    """Test reading a simple CSV file with header_start.

    This tests:
    - header_start : None
    """
    tbl_text = "a,b,c\n0,1,2\n3,4,5"
    out = table_read_csv(tbl_text, header_start=None)
    exp = [
        " f0   f1   f2 ",
        "str1 str1 str1",
        "---- ---- ----",
        "   a    b    c",
        "   0    1    2",
        "   3    4    5",
    ]
    assert out.pformat(show_dtype=True) == exp


def test_read_header_start_none_data_start_1():
    """Test reading a simple CSV file with data_start.

    This tests:
    - header_start : None
    - data_start : 1
    - names : provided list
    """
    tbl_text = "a,b,c\n0,1,2\n3,4,5\n6,7,8"
    names = ["x", "y", "z"]
    out = table_read_csv(tbl_text, header_start=None, data_start=1, names=names)
    exp = [
        "  x     y     z  ",
        "int64 int64 int64",
        "----- ----- -----",
        "    0     1     2",
        "    3     4     5",
        "    6     7     8",
    ]
    assert out.pformat(show_dtype=True) == exp


def test_read_data_start_2():
    """Test reading a simple CSV file with data_start.

    This tests:
    - data_start : 2
    """
    tbl_text = "a,b,c\n0,1,2\n3,4,5\n6,7,8"
    out = table_read_csv(tbl_text, data_start=2)
    exp = [
        "  a     b     c  ",
        "int64 int64 int64",
        "----- ----- -----",
        "    3     4     5",
        "    6     7     8",
    ]
    assert out.pformat(show_dtype=True) == exp


def test_read_names_exception(tbl_text):
    """Test that providing names without header_start=None raises an exception.

    This tests:
    - names : provided list
    """
    with pytest.raises(
        ValueError, match="cannot specify `names` unless `header_start=None`"
    ):
        table_read_csv(tbl_text, names=["x", "y"])


def test_read_include_names():
    tbl_text = "a,b,c\n0,1,2\n3,4,5"
    out = table_read_csv(tbl_text, include_names=["a", "b"])
    exp = [
        "  a     b  ",
        "int64 int64",
        "----- -----",
        "    0     1",
        "    3     4",
    ]
    assert out.pformat(show_dtype=True) == exp


def test_read_include_names_invalid(tbl_text):
    with pytest.raises(
        pa.lib.ArrowKeyError,
        match="Column 'x' in include_columns does not exist in CSV file",
    ):
        table_read_csv(tbl_text, include_names=["x", "y"])


@pytest.mark.parametrize("on_top", [True, False])
def test_strip_comment_lines(on_top):
    """Low level test of stripping comment lines."""
    tbl_text = textwrap.dedent("""\
    # Comment 1
     # Comment 2 (leading whitespace)
    # Comment 3
    a,b,c
    0,1,2
    3,4,5
    """)
    if not on_top:
        tbl_text += "# Comment 4\n# Comment 5\n"

    exp = textwrap.dedent("""\
    a,b,c
    0,1,2
    3,4,5
    """).encode("utf-8")

    input_file = io.BytesIO(tbl_text.encode("utf-8"))

    # For the default header_start and data_start values, the behavior depends on
    # whether there comment lines are only at the top of the file or not.
    output_file, header_start = strip_comment_lines(
        input_file, comment="#", header_start=0, data_start=None, encoding="utf-8"
    )
    if on_top:
        assert output_file is None
        assert header_start == 3
    else:
        assert output_file.getvalue() == exp
        assert header_start is None

    # In both of the following cases, strip_comment_lines should return a new BytesIO
    # buffer with the comment lines stripped, and the header_start should be None.
    input_file.seek(0)
    output_file, header_start = strip_comment_lines(
        input_file, comment="#", header_start=1, data_start=2, encoding="utf-8"
    )
    assert output_file.getvalue() == exp
    assert header_start is None

    input_file.seek(0)
    output_file, header_start = strip_comment_lines(
        input_file, comment="#", header_start=None, data_start=1, encoding="utf-8"
    )
    assert output_file.getvalue() == exp
    assert header_start is None


@pytest.mark.parametrize("on_top", [True, False])
def test_read_comments(on_top):
    """Test reading a simple CSV file with comments.

    This tests:
    - comment : "#"
    - All comment lines are at the top of file
    """
    tbl_text = textwrap.dedent("""\
    # Comment 1
     # Comment 2 (leading whitespace)
    # Comment 3
    a,b,c
    0,1,2
    3,4,5
    """)
    if not on_top:
        tbl_text += "\n# Comment 4\n# Comment 5\n"

    out = table_read_csv(tbl_text, comment="#")
    exp = [
        "  a     b     c  ",
        "int64 int64 int64",
        "----- ----- -----",
        "    0     1     2",
        "    3     4     5",
    ]
    assert out.pformat(show_dtype=True) == exp

    out = table_read_csv(tbl_text, header_start=None, data_start=1, comment="#")
    exp = [
        "  f0    f1    f2 ",
        "int64 int64 int64",
        "----- ----- -----",
        "    0     1     2",
        "    3     4     5",
    ]
    assert out.pformat(show_dtype=True) == exp

    out = table_read_csv(tbl_text, header_start=1, data_start=2, comment="#")
    exp = [
        "  0     1     2  ",
        "int64 int64 int64",
        "----- ----- -----",
        "    3     4     5",
    ]
    assert out.pformat(show_dtype=True) == exp


def test_read_null_values():
    """Test reading a simple CSV file with null values.

    This tests:
    - null_values : ["", "NLL"], ["NLL"], []
    """
    tbl_text = textwrap.dedent("""\
    a,b,c
    NLL,2,3
    4,,6
    """)
    out = table_read_csv(tbl_text, null_values=["", "NLL"])
    exp = Table(
        rows=[
            [np.ma.masked, 2, 3],
            [4, np.ma.masked, 6],
        ],
        names=["a", "b", "c"],
    )
    check_tables_equal(exp, out)

    out = table_read_csv(tbl_text, null_values=["NLL"])
    exp = Table(
        rows=[
            [np.ma.masked, "2", 3],
            [4, "", 6],
        ],
        names=["a", "b", "c"],
    )
    check_tables_equal(exp, out)

    out = table_read_csv(tbl_text, null_values=[])
    exp = Table(
        rows=[
            ["NLL", "2", 3],
            ["4", "", 6],
        ],
        names=["a", "b", "c"],
    )
    check_tables_equal(exp, out)


def test_read_null_values_default():
    """Test that the default null value is [""].

    PyArrow by default includes other null values like "NaN" or "nan".
    """
    tbl_text = textwrap.dedent("""\
    a,b,c
    ,nan,NaN
    4,,#N/A
    """)
    out = table_read_csv(tbl_text)
    exp = [
        "  a      b     c  ",
        "int64 float64 str4",
        "----- ------- ----",
        "   --     nan  NaN",
        "    4      -- #N/A",
    ]
    assert out.pformat(show_dtype=True) == exp


def test_read_newlines_in_values():
    """Test reading a simple CSV file with newlines in values.

    This tests:
    - newlines_in_values : True (default is False)
    """
    tbl_text = textwrap.dedent("""\
    a,b,c
    0,"1
    2",3
    4,5,"6
    7"
    """)
    out = table_read_csv(tbl_text, newlines_in_values=True)
    exp = [
        "  a    b    c  ",
        "int64 str3 str3",
        "----- ---- ----",
        "    0 1\\n2    3",
        "    4    5 6\\n7",
    ]
    assert out.pformat(show_dtype=True) == exp


def test_read_dates_times():
    """Test reading CSV file with times, dates, and timestamps."""
    tbl_text = textwrap.dedent("""
    date,time_of_day,timestamp_s,timestamp_ns
    2023-12-25,12:34:56,2023-12-25T12:34:56,2023-12-25T12:34:56.123456
    2024-01-01,23:59:59,2024-01-01T23:59:59,2024-01-01T23:59:59.987654
    """)
    out = table_read_csv(tbl_text)
    exp = [
        "     date     time_of_day     timestamp_s              timestamp_ns        ",
        "datetime64[D]    object      datetime64[s]            datetime64[ns]       ",
        "------------- ----------- ------------------- -----------------------------",
        "   2023-12-25    12:34:56 2023-12-25T12:34:56 2023-12-25T12:34:56.123456000",
        "   2024-01-01    23:59:59 2024-01-01T23:59:59 2024-01-01T23:59:59.987654000",
    ]
    assert out.pformat(show_dtype=True) == exp


def test_read_dates_times_masked():
    """Test reading CSV file with missing times, dates, and timestamps."""
    tbl_text = textwrap.dedent("""
    date,time_of_day,timestamp_s,timestamp_ns
    ,,,
    2024-01-01,23:59:59,2024-01-01T23:59:59,2024-01-01T23:59:59.987654
    """)
    out = table_read_csv(tbl_text)
    exp = [
        "     date     time_of_day     timestamp_s              timestamp_ns        ",
        "datetime64[D]    object      datetime64[s]            datetime64[ns]       ",
        "------------- ----------- ------------------- -----------------------------",
        "           --          --                  --                            --",
        "   2024-01-01    23:59:59 2024-01-01T23:59:59 2024-01-01T23:59:59.987654000",
    ]
    assert out.pformat(show_dtype=True) == exp

    # Confirm the values under the mask
    for col in out.itercols():
        col.mask = False
    exp = [
        "     date     time_of_day     timestamp_s              timestamp_ns        ",
        "datetime64[D]    object      datetime64[s]            datetime64[ns]       ",
        "------------- ----------- ------------------- -----------------------------",
        "   1970-01-01    00:00:00 1970-01-01T00:00:00 1970-01-01T00:00:00.000000000",
        "   2024-01-01    23:59:59 2024-01-01T23:59:59 2024-01-01T23:59:59.987654000",
    ]
    assert out.pformat(show_dtype=True) == exp


def test_read_dates_times_custom():
    """Test reading CSV file with custom timestamp formats."""
    tbl_text = textwrap.dedent("""
    date,timestamp_s,timestamp_custom
    2023-12-25,2023-12-25T12:34:56,01/31/2024 14:23:55
    2024-01-01,2024-01-01T23:59:59,12/25/2023 09:00:00
    """)
    out = table_read_csv(
        tbl_text, timestamp_parsers=[pa.csv.ISO8601, "%m/%d/%Y %H:%M:%S"]
    )
    exp = [
        "     date         timestamp_s       timestamp_custom ",
        "datetime64[D]    datetime64[s]       datetime64[s]   ",
        "------------- ------------------- -------------------",
        "   2023-12-25 2023-12-25T12:34:56 2024-01-31T14:23:55",
        "   2024-01-01 2024-01-01T23:59:59 2023-12-25T09:00:00",
    ]
    assert out.pformat(show_dtype=True) == exp


def test_read_whitespace_handling():
    """Test reading CSV file with leading and trailing whitespace.

    All whitespace in header and string columns are significant.
    Whitespace in a numeric column is ignored.
    """
    tbl_text = """
   a  , b , c
    0, 1.0 , x y
    3, 4.0 , w z """
    out = table_read_csv(tbl_text)
    assert out.colnames == ["   a  ", " b ", " c"]
    exp = [
        "   a      b       c ",
        "int64  float64  str5",
        "------ ------- -----",
        "     0     1.0   x y",
        "     3     4.0  w z ",
    ]
    assert out.pformat(show_dtype=True) == exp
