# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io
import locale
import pathlib
import platform
import re
from collections import OrderedDict
from io import BytesIO, StringIO

import numpy as np
import pytest

from astropy import table
from astropy.io import ascii
from astropy.io.ascii import core
from astropy.io.ascii.core import convert_numpy
from astropy.io.ascii.ui import _probably_html, get_read_trace
from astropy.table import MaskedColumn, Table
from astropy.table.table_helpers import simple_table
from astropy.units import Unit

# NOTE: Python can be built without bz2.
from astropy.utils.compat.optional_deps import HAS_BZ2
from astropy.utils.data import get_pkg_data_path
from astropy.utils.exceptions import AstropyDeprecationWarning, AstropyWarning

# setup/teardown function to have the tests run in the correct directory
from .common import (
    assert_almost_equal,
    assert_equal,
    assert_true,
    setup_function,  # noqa: F401
    teardown_function,  # noqa: F401
)


def asciiIO(x):
    return BytesIO(x.encode("ascii"))


@pytest.fixture
def home_is_data(monkeypatch, request):
    """
    Pytest fixture to run a test case with tilde-prefixed paths.

    In the tilde-path case, environment variables are temporarily
    modified so that '~' resolves to the data directory.
    """
    path = get_pkg_data_path("data")
    # For Unix
    monkeypatch.setenv("HOME", path)
    # For Windows
    monkeypatch.setenv("USERPROFILE", path)


@pytest.mark.parametrize(
    "fast_reader",
    [True, False, {"use_fast_converter": False}, {"use_fast_converter": True}, "force"],
)
def test_convert_overflow(fast_reader):
    """
    Test reading an extremely large integer, which falls through to
    string due to an overflow error (#2234). The C parsers used to
    return inf (kind 'f') for this.
    """
    expected_kind = "U"
    with pytest.warns(
        AstropyWarning, match="OverflowError converting to IntType in column a"
    ):
        dat = ascii.read(
            ["a", "1" * 10000], format="basic", fast_reader=fast_reader, guess=False
        )
    assert dat["a"].dtype.kind == expected_kind


def test_read_specify_converters_with_names():
    """
    Exact example from #9701: When using ascii.read with both the names and
    converters arguments, the converters dictionary ignores the user-supplied
    names and requires that you know the guessed names.
    """
    csv_text = ["a,b,c", "1,2,3", "4,5,6"]
    names = ["A", "B", "C"]

    converters = {
        "A": [ascii.convert_numpy(float)],
        "B": [ascii.convert_numpy(int)],
        "C": [ascii.convert_numpy(str)],
    }
    t = ascii.read(csv_text, format="csv", names=names, converters=converters)
    assert t["A"].dtype.kind == "f"
    assert t["B"].dtype.kind == "i"
    assert t["C"].dtype.kind == "U"


def test_read_remove_and_rename_columns():
    csv_text = ["a,b,c", "1,2,3", "4,5,6"]
    reader = ascii.get_reader(reader_cls=ascii.Csv)
    reader.read(csv_text)
    header = reader.header
    with pytest.raises(KeyError, match="Column NOT-EXIST does not exist"):
        header.remove_columns(["NOT-EXIST"])

    header.remove_columns(["c"])
    assert header.colnames == ("a", "b")

    header.rename_column("a", "aa")
    assert header.colnames == ("aa", "b")

    with pytest.raises(KeyError, match="Column NOT-EXIST does not exist"):
        header.rename_column("NOT-EXIST", "aa")


def test_guess_with_names_arg():
    """
    Make sure reading a table with guess=True gives the expected result when
    the names arg is specified.
    """
    # This is a NoHeader format table and so `names` should replace
    # the default col0, col1 names.  It fails as a Basic format
    # table when guessing because the column names would be '1', '2'.
    dat = ascii.read(["1,2", "3,4"], names=("a", "b"))
    assert len(dat) == 2
    assert dat.colnames == ["a", "b"]

    # This is a Basic format table and the first row
    # gives the column names 'c', 'd', which get replaced by 'a', 'b'
    dat = ascii.read(["c,d", "3,4"], names=("a", "b"))
    assert len(dat) == 1
    assert dat.colnames == ["a", "b"]

    # This is also a Basic format table and the first row
    # gives the column names 'c', 'd', which get replaced by 'a', 'b'
    dat = ascii.read(["c d", "e f"], names=("a", "b"))
    assert len(dat) == 1
    assert dat.colnames == ["a", "b"]


def test_guess_with_format_arg():
    """
    When the format is explicitly given then disable the strict column name checking in
    guessing.
    """
    dat = ascii.read(["1,2", "3,4"], format="basic")
    assert len(dat) == 1
    assert dat.colnames == ["1", "2"]

    dat = ascii.read(["1,2", "3,4"], names=("a", "b"), format="basic")
    assert len(dat) == 1
    assert dat.colnames == ["a", "b"]

    # For good measure check the same in the unified I/O interface
    dat = Table.read(["1,2", "3,4"], format="ascii.basic")
    assert len(dat) == 1
    assert dat.colnames == ["1", "2"]

    dat = Table.read(["1,2", "3,4"], format="ascii.basic", names=("a", "b"))
    assert len(dat) == 1
    assert dat.colnames == ["a", "b"]


def test_guess_with_delimiter_arg():
    """
    When the delimiter is explicitly given then do not try others in guessing.
    """
    fields = ["10.1E+19", "3.14", "2048", "-23"]
    values = [1.01e20, 3.14, 2048, -23]

    # Default guess should recognise CSV with optional spaces
    t0 = ascii.read(asciiIO(", ".join(fields)), guess=True)
    for n, v in zip(t0.colnames, values):
        assert t0[n][0] == v

    # Forcing space as delimiter produces type str columns ('10.1E+19,')
    t1 = ascii.read(asciiIO(", ".join(fields)), guess=True, delimiter=" ")
    for n, v in zip(t1.colnames[:-1], fields[:-1]):
        assert t1[n][0] == v + ","


def test_reading_mixed_delimiter_tabs_spaces():
    # Regression test for https://github.com/astropy/astropy/issues/6770
    dat = ascii.read("1 2\t3\n1 2\t3", format="no_header", names=list("abc"))
    assert len(dat) == 2

    Table.read(["1 2\t3", "1 2\t3"], format="ascii.no_header", names=["a", "b", "c"])
    assert len(dat) == 2


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
def test_read_with_names_arg(fast_reader):
    """
    Test that a bad value of `names` raises an exception.
    """
    # CParser only uses columns in `names` and thus reports mismatch in num_col
    with pytest.raises(ascii.InconsistentTableError):
        ascii.read(["c d", "e f"], names=("a",), guess=False, fast_reader=fast_reader)


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
@pytest.mark.parametrize("path_format", ["plain", "tilde-str", "tilde-pathlib"])
def test_read_all_files(fast_reader, path_format, home_is_data):
    for testfile in get_testfiles():
        if testfile.get("skip"):
            print(f"\n\n******** SKIPPING {testfile['name']}")
            continue
        if "tilde" in path_format:
            if "str" in path_format:
                testfile["name"] = "~/" + testfile["name"][5:]
            else:
                testfile["name"] = pathlib.Path("~/", testfile["name"][5:])
        print(f"\n\n******** READING {testfile['name']}")
        for guess in (True, False):
            test_opts = testfile["opts"].copy()
            if "guess" not in test_opts:
                test_opts["guess"] = guess
            if (
                "format" in test_opts
                and f"fast_{test_opts['format']}" in core.FAST_CLASSES
            ):  # has fast version
                if "inputter_cls" not in test_opts:  # fast reader doesn't allow this
                    test_opts["fast_reader"] = fast_reader
            table = ascii.read(testfile["name"], **test_opts)
            assert_equal(table.dtype.names, testfile["cols"])
            for colname in table.dtype.names:
                assert_equal(len(table[colname]), testfile["nrows"])


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
@pytest.mark.parametrize("path_format", ["plain", "tilde-str", "tilde-pathlib"])
def test_read_all_files_via_table(fast_reader, path_format, home_is_data):
    for testfile in get_testfiles():
        if testfile.get("skip"):
            print(f"\n\n******** SKIPPING {testfile['name']}")
            continue
        if "tilde" in path_format:
            if "str" in path_format:
                testfile["name"] = "~/" + testfile["name"][5:]
            else:
                testfile["name"] = pathlib.Path("~/", testfile["name"][5:])
        print(f"\n\n******** READING {testfile['name']}")
        for guess in (True, False):
            test_opts = testfile["opts"].copy()
            if "guess" not in test_opts:
                test_opts["guess"] = guess
            if "format" in test_opts:
                format = f"ascii.{test_opts['format']}"
                del test_opts["format"]
            else:
                format = "ascii"
            if f"fast_{format}" in core.FAST_CLASSES:
                test_opts["fast_reader"] = fast_reader
            table = Table.read(testfile["name"], format=format, **test_opts)
            assert_equal(table.dtype.names, testfile["cols"])
            for colname in table.dtype.names:
                assert_equal(len(table[colname]), testfile["nrows"])


def test_guess_all_files():
    for testfile in get_testfiles():
        if testfile.get("skip"):
            print(f"\n\n******** SKIPPING {testfile['name']}")
            continue
        if not testfile["opts"].get("guess", True):
            continue
        print(f"\n\n******** READING {testfile['name']}")
        for filter_read_opts in (["format", "delimiter", "quotechar"], []):
            # Copy read options except for those in filter_read_opts
            guess_opts = {
                k: v for k, v in testfile["opts"].items() if k not in filter_read_opts
            }
            table = ascii.read(testfile["name"], guess=True, **guess_opts)
            assert_equal(table.dtype.names, testfile["cols"])
            for colname in table.dtype.names:
                assert_equal(len(table[colname]), testfile["nrows"])


def test_validate_read_kwargs():
    lines = ["a b", "1 2", "3 4"]
    # Check that numpy integers are allowed
    out = ascii.read(lines, data_start=np.int16(2))
    assert np.all(out["a"] == [3])

    with pytest.raises(
        TypeError,
        match=r"read\(\) argument 'data_end' must be a "
        r"<class 'int'> object, "
        r"got <class 'str'> instead",
    ):
        ascii.read(lines, data_end="needs integer")

    with pytest.raises(
        TypeError,
        match=r"read\(\) argument 'fill_include_names' must "
        r"be a list-like object, got <class 'str'> instead",
    ):
        ascii.read(lines, fill_include_names="ID")


def test_daophot_indef():
    """Test that INDEF is correctly interpreted as a missing value"""
    table = ascii.read("data/daophot2.dat", format="daophot")
    for col in table.itercols():
        # Four columns have all INDEF values and are masked, rest are normal Column
        if col.name in ("OTIME", "MAG", "MERR", "XAIRMASS"):
            assert np.all(col.mask)
        else:
            assert not hasattr(col, "mask")


def test_daophot_types():
    """
    Test specific data types which are different from what would be
    inferred automatically based only data values.  DAOphot reader uses
    the header information to assign types.
    """
    table = ascii.read("data/daophot2.dat", format="daophot")
    assert table["LID"].dtype.char in "fd"  # float or double
    assert table["MAG"].dtype.char in "fd"  # even without any data values
    assert (
        table["PIER"].dtype.char in "US"
    )  # string (data values are consistent with int)
    assert table["ID"].dtype.kind == "i"  # int types: int, long, int64


def test_daophot_header_keywords():
    table = ascii.read("data/daophot.dat", format="daophot")
    expected_keywords = (
        ("NSTARFILE", "test.nst.1", "filename", "%-23s"),
        ("REJFILE", '"hello world"', "filename", "%-23s"),
        ("SCALE", "1.", "units/pix", "%-23.7g"),
    )

    keywords = table.meta["keywords"]  # Ordered dict of keyword structures
    for name, value, units, format_ in expected_keywords:
        keyword = keywords[name]
        assert_equal(keyword["value"], value)
        assert_equal(keyword["units"], units)
        assert_equal(keyword["format"], format_)


def test_daophot_multiple_aperture():
    table = ascii.read("data/daophot3.dat", format="daophot")
    assert "MAG5" in table.colnames  # MAG5 is one of the newly created column names
    assert table["MAG5"][4] == 22.13  # A sample entry in daophot3.dat file
    assert table["MERR2"][0] == 1.171
    assert np.all(
        table["RAPERT5"] == 23.3
    )  # assert all the 5th apertures are same 23.3


def test_daophot_multiple_aperture2():
    table = ascii.read("data/daophot4.dat", format="daophot")
    assert "MAG15" in table.colnames  # MAG15 is one of the newly created column name
    assert table["MAG15"][1] == -7.573  # A sample entry in daophot4.dat file
    assert table["MERR2"][0] == 0.049
    assert np.all(table["RAPERT5"] == 5.0)  # assert all the 5th apertures are same 5.0


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
def test_empty_table_no_header(fast_reader):
    with pytest.raises(ascii.InconsistentTableError):
        ascii.read(
            "data/no_data_without_header.dat",
            format="no_header",
            guess=False,
            fast_reader=fast_reader,
        )


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
def test_wrong_quote(fast_reader):
    with pytest.raises(ascii.InconsistentTableError):
        ascii.read("data/simple.txt", guess=False, fast_reader=fast_reader)


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
def test_extra_data_col(fast_reader):
    with pytest.raises(ascii.InconsistentTableError):
        ascii.read("data/bad.txt", fast_reader=fast_reader)


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
def test_extra_data_col2(fast_reader):
    with pytest.raises(ascii.InconsistentTableError):
        ascii.read("data/simple5.txt", delimiter="|", fast_reader=fast_reader)


def test_missing_file():
    with pytest.raises(OSError):
        ascii.read("does_not_exist")


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
def test_set_names(fast_reader):
    names = ("c1", "c2", "c3", "c4", "c5", "c6")
    data = ascii.read(
        "data/simple3.txt", names=names, delimiter="|", fast_reader=fast_reader
    )
    assert_equal(data.dtype.names, names)


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
def test_set_include_names(fast_reader):
    names = ("c1", "c2", "c3", "c4", "c5", "c6")
    include_names = ("c1", "c3")
    data = ascii.read(
        "data/simple3.txt",
        names=names,
        include_names=include_names,
        delimiter="|",
        fast_reader=fast_reader,
    )
    assert_equal(data.dtype.names, include_names)


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
def test_set_exclude_names(fast_reader):
    exclude_names = ("Y", "object")
    data = ascii.read(
        "data/simple3.txt",
        exclude_names=exclude_names,
        delimiter="|",
        fast_reader=fast_reader,
    )
    assert_equal(data.dtype.names, ("obsid", "redshift", "X", "rad"))


def test_include_names_daophot():
    include_names = ("ID", "MAG", "PIER")
    data = ascii.read("data/daophot.dat", include_names=include_names)
    assert_equal(data.dtype.names, include_names)


def test_exclude_names_daophot():
    exclude_names = ("ID", "YCENTER", "MERR", "NITER", "CHI", "PERROR")
    data = ascii.read("data/daophot.dat", exclude_names=exclude_names)
    assert_equal(data.dtype.names, ("XCENTER", "MAG", "MSKY", "SHARPNESS", "PIER"))


def test_custom_process_lines():
    def process_lines(lines):
        bars_at_ends = re.compile(r"^\| | \|$", re.VERBOSE)
        striplines = (x.strip() for x in lines)
        return [bars_at_ends.sub("", x) for x in striplines if len(x) > 0]

    reader = ascii.get_reader(delimiter="|")
    reader.inputter.process_lines = process_lines
    data = reader.read("data/bars_at_ends.txt")
    assert_equal(data.dtype.names, ("obsid", "redshift", "X", "Y", "object", "rad"))
    assert_equal(len(data), 3)


def test_custom_process_line():
    def process_line(line):
        line_out = re.sub(r"^\|\s*", "", line.strip())
        return line_out

    reader = ascii.get_reader(data_start=2, delimiter="|")
    reader.header.splitter.process_line = process_line
    reader.data.splitter.process_line = process_line
    data = reader.read("data/nls1_stackinfo.dbout")
    cols = get_testfiles("data/nls1_stackinfo.dbout")["cols"]
    assert_equal(data.dtype.names, cols[1:])


def test_custom_splitters():
    reader = ascii.get_reader()
    reader.header.splitter = ascii.BaseSplitter()
    reader.data.splitter = ascii.BaseSplitter()
    f = "data/test4.dat"
    data = reader.read(f)
    testfile = get_testfiles(f)
    assert_equal(data.dtype.names, testfile["cols"])
    assert_equal(len(data), testfile["nrows"])
    assert_almost_equal(data.field("zabs1.nh")[2], 0.0839710433091)
    assert_almost_equal(data.field("p1.gamma")[2], 1.25997502704)
    assert_almost_equal(data.field("p1.ampl")[2], 0.000696444029148)
    assert_equal(data.field("statname")[2], "chi2modvar")
    assert_almost_equal(data.field("statval")[2], 497.56468441)


def test_start_end():
    data = ascii.read("data/test5.dat", header_start=1, data_start=3, data_end=-5)
    assert_equal(len(data), 13)
    assert_equal(data.field("statname")[0], "chi2xspecvar")
    assert_equal(data.field("statname")[-1], "chi2gehrels")


def test_set_converters():
    converters = {
        "zabs1.nh": [ascii.convert_numpy("int32"), ascii.convert_numpy("float32")],
        "p1.gamma": [ascii.convert_numpy("str")],
    }
    data = ascii.read("data/test4.dat", converters=converters)
    assert_equal(str(data["zabs1.nh"].dtype), "float32")
    assert_equal(data["p1.gamma"][0], "1.26764500000")


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
def test_from_string(fast_reader):
    f = "data/simple.txt"
    with open(f) as fd:
        table = fd.read()
    testfile = get_testfiles(f)[0]
    data = ascii.read(table, fast_reader=fast_reader, **testfile["opts"])
    assert_equal(data.dtype.names, testfile["cols"])
    assert_equal(len(data), testfile["nrows"])


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
def test_from_filelike(fast_reader):
    f = "data/simple.txt"
    testfile = get_testfiles(f)[0]
    with open(f, "rb") as fd:
        data = ascii.read(fd, fast_reader=fast_reader, **testfile["opts"])
    assert_equal(data.dtype.names, testfile["cols"])
    assert_equal(len(data), testfile["nrows"])


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
def test_from_lines(fast_reader):
    f = "data/simple.txt"
    with open(f) as fd:
        table = fd.readlines()
    testfile = get_testfiles(f)[0]
    data = ascii.read(table, fast_reader=fast_reader, **testfile["opts"])
    assert_equal(data.dtype.names, testfile["cols"])
    assert_equal(len(data), testfile["nrows"])


def test_comment_lines():
    table = ascii.get_reader(reader_cls=ascii.Rdb)
    data = table.read("data/apostrophe.rdb")
    assert_equal(table.comment_lines, ["# first comment", "  # second comment"])
    assert_equal(data.meta["comments"], ["first comment", "second comment"])


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
def test_fill_values(fast_reader):
    f = "data/fill_values.txt"
    testfile = get_testfiles(f)
    data = ascii.read(
        f, fill_values=("a", "1"), fast_reader=fast_reader, **testfile["opts"]
    )
    assert_true((data["a"].mask == [False, True]).all())
    assert_true((data["a"] == [1, 1]).all())
    assert_true((data["b"].mask == [False, True]).all())
    assert_true((data["b"] == [2, 1]).all())


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
def test_fill_values_col(fast_reader):
    f = "data/fill_values.txt"
    testfile = get_testfiles(f)
    data = ascii.read(
        f, fill_values=("a", "1", "b"), fast_reader=fast_reader, **testfile["opts"]
    )
    check_fill_values(data)


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
def test_fill_values_include_names(fast_reader):
    f = "data/fill_values.txt"
    testfile = get_testfiles(f)
    data = ascii.read(
        f,
        fill_values=("a", "1"),
        fast_reader=fast_reader,
        fill_include_names=["b"],
        **testfile["opts"],
    )
    check_fill_values(data)


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
def test_fill_values_exclude_names(fast_reader):
    f = "data/fill_values.txt"
    testfile = get_testfiles(f)
    data = ascii.read(
        f,
        fill_values=("a", "1"),
        fast_reader=fast_reader,
        fill_exclude_names=["a"],
        **testfile["opts"],
    )
    check_fill_values(data)


def check_fill_values(data):
    """compare array column by column with expectation"""
    assert not hasattr(data["a"], "mask")
    assert_true((data["a"] == ["1", "a"]).all())
    assert_true((data["b"].mask == [False, True]).all())
    # Check that masked value is "do not care" in comparison
    assert_true((data["b"] == [2, -999]).all())
    data["b"].mask = False  # explicitly unmask for comparison
    assert_true((data["b"] == [2, 1]).all())


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
def test_fill_values_list(fast_reader):
    f = "data/fill_values.txt"
    testfile = get_testfiles(f)
    data = ascii.read(
        f,
        fill_values=[("a", "42"), ("1", "42", "a")],
        fast_reader=fast_reader,
        **testfile["opts"],
    )
    data["a"].mask = False  # explicitly unmask for comparison
    assert_true((data["a"] == [42, 42]).all())


def test_masking_Cds_Mrt():
    f = "data/cds.dat"  # Tested for CDS and MRT
    for testfile in get_testfiles(f):
        data = ascii.read(f, **testfile["opts"])
        assert_true(data["AK"].mask[0])
        assert not hasattr(data["Fit"], "mask")


def test_null_Ipac():
    f = "data/ipac.dat"
    testfile = get_testfiles(f)[0]
    data = ascii.read(f, **testfile["opts"])
    mask = np.array(
        [(True, False, True, False, True), (False, False, False, False, False)],
        dtype=[
            ("ra", "|b1"),
            ("dec", "|b1"),
            ("sai", "|b1"),
            ("v2", "|b1"),
            ("sptype", "|b1"),
        ],
    )
    assert np.all(data.mask == mask)


def test_Ipac_meta():
    keywords = OrderedDict(
        (
            ("intval", 1),
            ("floatval", 2.3e3),
            ("date", "Wed Sp 20 09:48:36 1995"),
            ("key_continue", "IPAC keywords can continue across lines"),
        )
    )
    comments = ["This is an example of a valid comment"]
    f = "data/ipac.dat"
    testfile = get_testfiles(f)[0]
    data = ascii.read(f, **testfile["opts"])
    assert data.meta["keywords"].keys() == keywords.keys()
    for data_kv, kv in zip(data.meta["keywords"].values(), keywords.values()):
        assert data_kv["value"] == kv
    assert data.meta["comments"] == comments


def test_set_guess_kwarg():
    """Read a file using guess with one of the typical guess_kwargs explicitly set."""
    data = ascii.read("data/space_delim_no_header.dat", delimiter=",", guess=True)
    assert data.dtype.names == ("1 3.4 hello",)
    assert len(data) == 1


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
def test_read_rdb_wrong_type(fast_reader):
    """Read RDB data with inconsistent data type (except failure)"""
    table = """col1\tcol2
N\tN
1\tHello"""
    with pytest.raises(ValueError):
        ascii.read(table, format="rdb", fast_reader=fast_reader)


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
def test_default_missing(fast_reader):
    """
    Read a table with empty values and ensure that corresponding entries are masked
    """
    table = "\n".join(  # noqa: FLY002
        [
            "a,b,c,d",
            "1,3,,",
            "2, , 4.0 , ss ",
        ]
    )
    dat = ascii.read(table, fast_reader=fast_reader)
    assert dat.masked is False
    assert dat.pformat() == [
        " a   b   c   d ",
        "--- --- --- ---",
        "  1   3  --  --",
        "  2  -- 4.0  ss",
    ]

    # Single row table with a single missing element
    table = """ a \n "" """
    dat = ascii.read(table, fast_reader=fast_reader)
    assert dat.pformat() == [" a ", "---", " --"]
    assert dat["a"].dtype.kind == "i"

    # Same test with a fixed width reader
    table = "\n".join(  # noqa: FLY002
        [
            " a   b   c   d ",
            "--- --- --- ---",
            "  1   3        ",
            "  2     4.0  ss",
        ]
    )
    dat = ascii.read(table, format="fixed_width_two_line")
    assert dat.masked is False
    assert dat.pformat() == [
        " a   b   c   d ",
        "--- --- --- ---",
        "  1   3  --  --",
        "  2  -- 4.0  ss",
    ]

    dat = ascii.read(table, format="fixed_width_two_line", fill_values=None)
    assert dat.masked is False
    assert dat.pformat() == [
        " a   b   c   d ",
        "--- --- --- ---",
        "  1   3        ",
        "  2     4.0  ss",
    ]

    dat = ascii.read(table, format="fixed_width_two_line", fill_values=[])
    assert dat.masked is False
    assert dat.pformat() == [
        " a   b   c   d ",
        "--- --- --- ---",
        "  1   3        ",
        "  2     4.0  ss",
    ]


def get_testfiles(name=None):
    """Set up information about the columns, number of rows, and reader params to
    read a bunch of test files and verify columns and number of rows."""

    testfiles = [
        {
            "cols": ("agasc_id", "n_noids", "n_obs"),
            "name": "data/apostrophe.rdb",
            "nrows": 2,
            "opts": {"format": "rdb"},
        },
        {
            "cols": ("agasc_id", "n_noids", "n_obs"),
            "name": "data/apostrophe.tab",
            "nrows": 2,
            "opts": {"format": "tab"},
        },
        {
            "cols": (
                "Index",
                "RAh",
                "RAm",
                "RAs",
                "DE-",
                "DEd",
                "DEm",
                "DEs",
                "Match",
                "Class",
                "AK",
                "Fit",
            ),
            "name": "data/cds.dat",
            "nrows": 1,
            "opts": {"format": "cds"},
        },
        {
            "cols": (
                "Index",
                "RAh",
                "RAm",
                "RAs",
                "DE-",
                "DEd",
                "DEm",
                "DEs",
                "Match",
                "Class",
                "AK",
                "Fit",
            ),
            "name": "data/cds.dat",
            "nrows": 1,
            "opts": {"format": "mrt"},
        },
        # Test CDS and MRT files with dashes in column name
        {
            "cols": (
                "DefaultName",
                "#CompsOnThisRow",
            ),
            "name": "data/cds_mrt_dashes.txt",
            "nrows": 8,
            "opts": {"format": "mrt"},
        },
        {
            "cols": (
                "DefaultName",
                "#CompsOnThisRow",
            ),
            "name": "data/cds_mrt_dashes.txt",
            "nrows": 8,
            "opts": {"format": "cds"},
        },
        # Test malformed CDS file (issues #2241 #467)
        {
            "cols": (
                "Index",
                "RAh",
                "RAm",
                "RAs",
                "DE-",
                "DEd",
                "DEm",
                "DEs",
                "Match",
                "Class",
                "AK",
                "Fit",
            ),
            "name": "data/cds_malformed.dat",
            "nrows": 1,
            "opts": {"format": "cds", "data_start": "guess"},
        },
        {
            "cols": ("a", "b", "c"),
            "name": "data/commented_header.dat",
            "nrows": 2,
            "opts": {"format": "commented_header"},
        },
        {
            "cols": ("a", "b", "c"),
            "name": "data/commented_header2.dat",
            "nrows": 2,
            "opts": {"format": "commented_header", "header_start": -1},
        },
        {
            "cols": ("col1", "col2", "col3", "col4", "col5"),
            "name": "data/continuation.dat",
            "nrows": 2,
            "opts": {
                "inputter_cls": ascii.ContinuationLinesInputter,
                "format": "no_header",
            },
        },
        {
            "cols": (
                "ID",
                "XCENTER",
                "YCENTER",
                "MAG",
                "MERR",
                "MSKY",
                "NITER",
                "SHARPNESS",
                "CHI",
                "PIER",
                "PERROR",
            ),
            "name": "data/daophot.dat",
            "nrows": 2,
            "opts": {"format": "daophot"},
        },
        {
            "cols": (
                "NUMBER",
                "FLUX_ISO",
                "FLUXERR_ISO",
                "VALU-ES",
                "VALU-ES_1",
                "FLAG",
            ),
            "name": "data/sextractor.dat",
            "nrows": 3,
            "opts": {"format": "sextractor"},
        },
        {
            "cols": ("ra", "dec", "sai", "v2", "sptype"),
            "name": "data/ipac.dat",
            "nrows": 2,
            "opts": {"format": "ipac"},
        },
        {
            "cols": (
                "col0",
                "objID",
                "osrcid",
                "xsrcid",
                "SpecObjID",
                "ra",
                "dec",
                "obsid",
                "ccdid",
                "z",
                "modelMag_i",
                "modelMagErr_i",
                "modelMag_r",
                "modelMagErr_r",
                "expo",
                "theta",
                "rad_ecf_39",
                "detlim90",
                "fBlim90",
            ),
            "name": "data/nls1_stackinfo.dbout",
            "nrows": 58,
            "opts": {"data_start": 2, "delimiter": "|", "guess": False},
        },
        {
            "cols": (
                "Index",
                "RAh",
                "RAm",
                "RAs",
                "DE-",
                "DEd",
                "DEm",
                "DEs",
                "Match",
                "Class",
                "AK",
                "Fit",
            ),
            "name": "data/no_data_cds.dat",
            "nrows": 0,
            "opts": {"format": "cds"},
        },
        {
            "cols": (
                "Index",
                "RAh",
                "RAm",
                "RAs",
                "DE-",
                "DEd",
                "DEm",
                "DEs",
                "Match",
                "Class",
                "AK",
                "Fit",
            ),
            "name": "data/no_data_cds.dat",
            "nrows": 0,
            "opts": {"format": "mrt"},
        },
        {
            "cols": (
                "ID",
                "XCENTER",
                "YCENTER",
                "MAG",
                "MERR",
                "MSKY",
                "NITER",
                "SHARPNESS",
                "CHI",
                "PIER",
                "PERROR",
            ),
            "name": "data/no_data_daophot.dat",
            "nrows": 0,
            "opts": {"format": "daophot"},
        },
        {
            "cols": ("NUMBER", "FLUX_ISO", "FLUXERR_ISO", "VALUES", "VALUES_1", "FLAG"),
            "name": "data/no_data_sextractor.dat",
            "nrows": 0,
            "opts": {"format": "sextractor"},
        },
        {
            "cols": ("ra", "dec", "sai", "v2", "sptype"),
            "name": "data/no_data_ipac.dat",
            "nrows": 0,
            "opts": {"format": "ipac"},
        },
        {
            "cols": ("ra", "v2"),
            "name": "data/ipac.dat",
            "nrows": 2,
            "opts": {"format": "ipac", "include_names": ["ra", "v2"]},
        },
        {
            "cols": ("a", "b", "c"),
            "name": "data/no_data_with_header.dat",
            "nrows": 0,
            "opts": {},
        },
        {
            "cols": ("agasc_id", "n_noids", "n_obs"),
            "name": "data/short.rdb",
            "nrows": 7,
            "opts": {"format": "rdb"},
        },
        {
            "cols": ("agasc_id", "n_noids", "n_obs"),
            "name": "data/short.tab",
            "nrows": 7,
            "opts": {"format": "tab"},
        },
        {
            "cols": ("test 1a", "test2", "test3", "test4"),
            "name": "data/simple.txt",
            "nrows": 2,
            "opts": {"quotechar": "'"},
        },
        {
            "cols": ("top1", "top2", "top3", "top4"),
            "name": "data/simple.txt",
            "nrows": 1,
            "opts": {"quotechar": "'", "header_start": 1, "data_start": 2},
        },
        {
            "cols": ("top1", "top2", "top3", "top4"),
            "name": "data/simple.txt",
            "nrows": 1,
            "opts": {"quotechar": "'", "header_start": 1},
        },
        {
            "cols": ("top1", "top2", "top3", "top4"),
            "name": "data/simple.txt",
            "nrows": 2,
            "opts": {"quotechar": "'", "header_start": 1, "data_start": 1},
        },
        {
            "cols": ("obsid", "redshift", "X", "Y", "object", "rad"),
            "name": "data/simple2.txt",
            "nrows": 3,
            "opts": {"delimiter": "|"},
        },
        {
            "cols": ("obsid", "redshift", "X", "Y", "object", "rad"),
            "name": "data/simple3.txt",
            "nrows": 2,
            "opts": {"delimiter": "|"},
        },
        {
            "cols": ("col1", "col2", "col3", "col4", "col5", "col6"),
            "name": "data/simple4.txt",
            "nrows": 3,
            "opts": {"format": "no_header", "delimiter": "|"},
        },
        {
            "cols": ("col1", "col2", "col3"),
            "name": "data/space_delim_no_header.dat",
            "nrows": 2,
            "opts": {"format": "no_header"},
        },
        {
            "cols": ("col1", "col2", "col3"),
            "name": "data/space_delim_no_header.dat",
            "nrows": 2,
            "opts": {"format": "no_header", "header_start": None},
        },
        {
            "cols": ("obsid", "offset", "x", "y", "name", "oaa"),
            "name": "data/space_delim_blank_lines.txt",
            "nrows": 3,
            "opts": {},
        },
        {
            "cols": ("zabs1.nh", "p1.gamma", "p1.ampl", "statname", "statval"),
            "name": "data/test4.dat",
            "nrows": 9,
            "opts": {},
        },
        {
            "cols": ("a", "b", "c"),
            "name": "data/fill_values.txt",
            "nrows": 2,
            "opts": {"delimiter": ","},
        },
        {
            "name": "data/whitespace.dat",
            "cols": ("quoted colname with tab\tinside", "col2", "col3"),
            "nrows": 2,
            "opts": {"delimiter": r"\s"},
        },
        {
            "name": "data/simple_csv.csv",
            "cols": ("a", "b", "c"),
            "nrows": 2,
            "opts": {"format": "csv"},
        },
        {
            "name": "data/simple_csv_missing.csv",
            "cols": ("a", "b", "c"),
            "nrows": 2,
            "skip": True,
            "opts": {"format": "csv"},
        },
        {
            "cols": ("cola", "colb", "colc"),
            "name": "data/latex1.tex",
            "nrows": 2,
            "opts": {"format": "latex"},
        },
        {
            "cols": ("Facility", "Id", "exposure", "date"),
            "name": "data/latex2.tex",
            "nrows": 3,
            "opts": {"format": "aastex"},
        },
        {
            "cols": ("cola", "colb", "colc"),
            "name": "data/latex3.tex",
            "nrows": 2,
            "opts": {"format": "latex"},
        },
        {
            "cols": ("Col1", "Col2", "Col3", "Col4"),
            "name": "data/fixed_width_2_line.txt",
            "nrows": 2,
            "opts": {"format": "fixed_width_two_line"},
        },
    ]

    try:
        import bs4  # noqa: F401

        testfiles.append(
            {
                "cols": ("Column 1", "Column 2", "Column 3"),
                "name": "data/html.html",
                "nrows": 3,
                "opts": {"format": "html"},
            }
        )
    except ImportError:
        pass

    if name is not None:
        # If there are multiple matches then return a list, else return just
        # the one match.
        out = [x for x in testfiles if x["name"] == name]
        if len(out) == 1:
            out = out[0]
    else:
        out = testfiles

    return out


def test_header_start_exception():
    """Check certain Readers throw an exception if ``header_start`` is set

    For certain Readers it does not make sense to set the ``header_start``, they
    throw an exception if you try.
    This was implemented in response to issue #885.
    """
    for readerclass in [
        ascii.NoHeader,
        ascii.SExtractor,
        ascii.Ipac,
        ascii.BaseReader,
        ascii.FixedWidthNoHeader,
        ascii.Cds,
        ascii.Mrt,
        ascii.Daophot,
    ]:
        with pytest.raises(ValueError):
            ascii.core._get_reader(readerclass, header_start=5)


def test_csv_table_read():
    """
    Check for a regression introduced by #1935.  Pseudo-CSV file with
    commented header line.
    """
    lines = ["# a, b", "1, 2", "3, 4"]
    t = ascii.read(lines)
    assert t.colnames == ["a", "b"]


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
def test_overlapping_names(fast_reader):
    """
    Check that the names argument list can overlap with the existing column names.
    This tests the issue in #1991.
    """
    t = ascii.read(["a b", "1 2"], names=["b", "a"], fast_reader=fast_reader)
    assert t.colnames == ["b", "a"]


def test_sextractor_units():
    """
    Make sure that the SExtractor reader correctly inputs descriptions and units.
    """
    table = ascii.read("data/sextractor2.dat", format="sextractor", guess=False)
    expected_units = [
        None,
        Unit("pix"),
        Unit("pix"),
        Unit("mag"),
        Unit("mag"),
        None,
        Unit("pix**2"),
        Unit("m**(-6)"),
        Unit("mag * arcsec**(-2)"),
    ]
    expected_descrs = [
        "Running object number",
        "Windowed position estimate along x",
        "Windowed position estimate along y",
        "Kron-like elliptical aperture magnitude",
        "RMS error for AUTO magnitude",
        "Extraction flags",
        None,
        "Barycenter position along MAMA x axis",
        "Peak surface brightness above background",
    ]
    for i, colname in enumerate(table.colnames):
        assert table[colname].unit == expected_units[i]
        assert table[colname].description == expected_descrs[i]


def test_sextractor_last_column_array():
    """
    Make sure that the SExtractor reader handles the last column correctly when it is array-like.
    """
    table = ascii.read("data/sextractor3.dat", format="sextractor", guess=False)
    expected_columns = [
        "X_IMAGE",
        "Y_IMAGE",
        "ALPHA_J2000",
        "DELTA_J2000",
        "MAG_AUTO",
        "MAGERR_AUTO",
        "MAG_APER",
        "MAG_APER_1",
        "MAG_APER_2",
        "MAG_APER_3",
        "MAG_APER_4",
        "MAG_APER_5",
        "MAG_APER_6",
        "MAGERR_APER",
        "MAGERR_APER_1",
        "MAGERR_APER_2",
        "MAGERR_APER_3",
        "MAGERR_APER_4",
        "MAGERR_APER_5",
        "MAGERR_APER_6",
    ]
    expected_units = [
        Unit("pix"),
        Unit("pix"),
        Unit("deg"),
        Unit("deg"),
        Unit("mag"),
        Unit("mag"),
        Unit("mag"),
        Unit("mag"),
        Unit("mag"),
        Unit("mag"),
        Unit("mag"),
        Unit("mag"),
        Unit("mag"),
        Unit("mag"),
        Unit("mag"),
        Unit("mag"),
        Unit("mag"),
        Unit("mag"),
        Unit("mag"),
        Unit("mag"),
    ]
    expected_descrs = (
        [
            "Object position along x",
            None,
            "Right ascension of barycenter (J2000)",
            "Declination of barycenter (J2000)",
            "Kron-like elliptical aperture magnitude",
            "RMS error for AUTO magnitude",
        ]
        + ["Fixed aperture magnitude vector"] * 7
        + ["RMS error vector for fixed aperture mag."] * 7
    )
    for i, colname in enumerate(table.colnames):
        assert table[colname].name == expected_columns[i]
        assert table[colname].unit == expected_units[i]
        assert table[colname].description == expected_descrs[i]


def test_list_with_newlines():
    """
    Check that lists of strings where some strings consist of just a newline
    ("\n") are parsed correctly.
    """
    t = ascii.read(["abc", "123\n", "456\n", "\n", "\n"])
    assert t.colnames == ["abc"]
    assert len(t) == 2
    assert t[0][0] == 123
    assert t[1][0] == 456


def test_commented_csv():
    """
    Check that Csv reader does not have ignore lines with the # comment
    character which is defined for most Basic readers.
    """
    t = ascii.read(["#a,b", "1,2", "#3,4"], format="csv")
    assert t.colnames == ["#a", "b"]
    assert len(t) == 2
    assert t["#a"][1] == "#3"


def test_meta_comments():
    """
    Make sure that line comments are included in the ``meta`` attribute
    of the output Table.
    """
    t = ascii.read(["#comment1", "#   comment2 \t", "a,b,c", "1,2,3"])
    assert t.colnames == ["a", "b", "c"]
    assert t.meta["comments"] == ["comment1", "comment2"]


def test_guess_fail():
    """
    Check the error message when guess fails
    """
    with pytest.raises(ascii.InconsistentTableError) as err:
        ascii.read("asfdasdf\n1 2 3", format="basic")
    assert "** To figure out why the table did not read, use guess=False and" in str(
        err.value
    )

    # Test the case with guessing enabled but for a format that has no free params
    with pytest.raises(ValueError) as err:
        ascii.read("asfdasdf\n1 2 3", format="ipac")
    assert (
        "At least one header line beginning and ending with delimiter required"
        in str(err.value)
    )

    # Test the case with guessing enabled but with all params specified
    with pytest.raises(ValueError) as err:
        ascii.read(
            "asfdasdf\n1 2 3",
            format="basic",
            quotechar='"',
            delimiter=" ",
            fast_reader=False,
        )
    assert "Number of header columns (1) inconsistent with data columns (3)" in str(
        err.value
    )


@pytest.mark.xfail(not HAS_BZ2, reason="requires bz2")
def test_guessing_file_object():
    """
    Test guessing a file object.  Fixes #3013 and similar issue noted in #3019.
    """
    with open("data/ipac.dat.bz2", "rb") as fd:
        t = ascii.read(fd)
    assert t.colnames == ["ra", "dec", "sai", "v2", "sptype"]


def test_pformat_roundtrip():
    """Check that the screen output of ``print tab`` can be read. See #3025."""
    """Read a table with empty values and ensure that corresponding entries are masked"""
    table = "\n".join(  # noqa: FLY002
        [
            "a,b,c,d",
            "1,3,1.11,1",
            "2, 2, 4.0 , ss ",
        ]
    )
    dat = ascii.read(table)
    out = ascii.read(dat.pformat())
    assert len(dat) == len(out)
    assert dat.colnames == out.colnames
    for c in dat.colnames:
        assert np.all(dat[c] == out[c])


def test_ipac_abbrev():
    lines = [
        "| c1 | c2 | c3   |   c4 | c5| c6 | c7  | c8 | c9|c10|c11|c12|",
        "| r  | rE | rea  | real | D | do | dou | f  | i | l | da| c |",
        "  1    2    3       4     5   6    7     8    9   10  11  12 ",
    ]
    dat = ascii.read(lines, format="ipac")
    for name in dat.columns[0:8]:
        assert dat[name].dtype.kind == "f"
    for name in dat.columns[8:10]:
        assert dat[name].dtype.kind == "i"
    for name in dat.columns[10:12]:
        assert dat[name].dtype.kind in ("U", "S")


def test_almost_but_not_quite_daophot():
    """Regression test for #3319.
    This tables looks so close to a daophot table, that the daophot reader gets
    quite far before it fails with an AttributeError.

    Note that this table will actually be read as Commented Header table with
    the columns ['some', 'header', 'info'].
    """
    lines = [
        "# some header info",
        "#F header info beginning with 'F'",
        "1 2 3",
        "4 5 6",
        "7 8 9",
    ]
    dat = ascii.read(lines)
    assert len(dat) == 3


@pytest.mark.parametrize("fast", [False, "force"])
def test_commented_header_comments(fast):
    """
    Test that comments in commented_header are as expected with header_start
    at different positions, and that the table round-trips.
    """
    comments = ["comment 1", "comment 2", "comment 3"]
    lines = ["# a b", "# comment 1", "# comment 2", "# comment 3", "1 2", "3 4"]
    dat = ascii.read(lines, format="commented_header", fast_reader=fast)
    assert dat.meta["comments"] == comments
    assert dat.colnames == ["a", "b"]

    out = StringIO()
    ascii.write(dat, out, format="commented_header", fast_writer=fast)
    assert out.getvalue().splitlines() == lines

    lines.insert(1, lines.pop(0))
    dat = ascii.read(lines, format="commented_header", header_start=1, fast_reader=fast)
    assert dat.meta["comments"] == comments
    assert dat.colnames == ["a", "b"]

    lines.insert(2, lines.pop(1))
    dat = ascii.read(lines, format="commented_header", header_start=2, fast_reader=fast)
    assert dat.meta["comments"] == comments
    assert dat.colnames == ["a", "b"]
    dat = ascii.read(
        lines, format="commented_header", header_start=-2, fast_reader=fast
    )
    assert dat.meta["comments"] == comments
    assert dat.colnames == ["a", "b"]

    lines.insert(3, lines.pop(2))
    dat = ascii.read(
        lines, format="commented_header", header_start=-1, fast_reader=fast
    )
    assert dat.meta["comments"] == comments
    assert dat.colnames == ["a", "b"]

    lines = ["# a b", "1 2", "3 4"]
    dat = ascii.read(lines, format="commented_header", fast_reader=fast)
    assert "comments" not in dat.meta
    assert dat.colnames == ["a", "b"]


def test_probably_html(home_is_data):
    """
    Test the routine for guessing if a table input to ascii.read is probably HTML
    """
    for tabl0 in (
        "data/html.html",
        "~/html.html",
        "http://blah.com/table.html",
        "https://blah.com/table.html",
        "file://blah/table.htm",
        "ftp://blah.com/table.html",
        "file://blah.com/table.htm",
        " <! doctype html > hello world",
        "junk < table baz> <tr foo > <td bar> </td> </tr> </table> junk",
        [
            "junk < table baz>",
            " <tr foo >",
            " <td bar> ",
            "</td> </tr>",
            "</table> junk",
        ],
        (" <! doctype html > ", " hello world"),
    ):
        assert _probably_html(tabl0) is True

    for tabl0 in (
        "data/html.htms",
        "Xhttp://blah.com/table.html",
        " https://blah.com/table.htm",
        "fole://blah/table.htm",
        " < doctype html > hello world",
        "junk < tble baz> <tr foo > <td bar> </td> </tr> </table> junk",
        [
            "junk < table baz>",
            " <t foo >",
            " <td bar> ",
            "</td> </tr>",
            "</table> junk",
        ],
        (" <! doctype htm > ", " hello world"),
        [[1, 2, 3]],
    ):
        assert _probably_html(tabl0) is False


@pytest.mark.parametrize("fast_reader", [True, False, "force"])
def test_data_header_start(fast_reader):
    tests = [
        (
            [
                "# comment",
                "",
                " ",
                "skip this line",  # line 0
                "a b",  # line 1
                "1 2",
            ],  # line 2
            [{"header_start": 1}, {"header_start": 1, "data_start": 2}],
        ),
        (
            [
                "# comment",
                "",
                " \t",
                "skip this line",  # line 0
                "a b",  # line 1
                "",
                " \t",
                "skip this line",  # line 2
                "1 2",
            ],  # line 3
            [{"header_start": 1, "data_start": 3}],
        ),
        (
            [
                "# comment",
                "",
                " ",
                "a b",  # line 0
                "",
                " ",
                "skip this line",  # line 1
                "1 2",
            ],  # line 2
            [{"header_start": 0, "data_start": 2}, {"data_start": 2}],
        ),
    ]

    for lines, kwargs_list in tests:
        for kwargs in kwargs_list:
            t = ascii.read(
                lines, format="basic", fast_reader=fast_reader, guess=True, **kwargs
            )
            assert t.colnames == ["a", "b"]
            assert len(t) == 1
            assert np.all(t["a"] == [1])
            # Sanity check that the expected reader_cls is being used
            assert get_read_trace()[-1]["kwargs"]["reader_cls"] is (
                ascii.Basic if (fast_reader is False) else ascii.FastBasic
            )


def test_table_with_no_newline():
    """
    Test that an input file which is completely empty fails in the expected way.
    Test that an input file with one line but no newline succeeds.
    """
    # With guessing
    table = BytesIO()
    with pytest.raises(ascii.InconsistentTableError):
        ascii.read(table)

    # Without guessing
    table = BytesIO()
    with pytest.raises(ValueError) as err:
        ascii.read(table, guess=False, fast_reader=False, format="basic")
    assert "No header line found" in str(err.value)

    table = BytesIO()
    t = ascii.read(table, guess=False, fast_reader=True, format="fast_basic")
    assert not t and t.as_array().size == 0

    # Put a single line of column names but with no newline
    for kwargs in [
        {},
        {"guess": False, "fast_reader": False, "format": "basic"},
        {"guess": False, "fast_reader": True, "format": "fast_basic"},
    ]:
        table = BytesIO()
        table.write(b"a b")
        t = ascii.read(table, **kwargs)
        assert t.colnames == ["a", "b"]
        assert len(t) == 0


def test_path_object():
    fpath = pathlib.Path("data/simple.txt")
    data = ascii.read(fpath)

    assert len(data) == 2
    assert sorted(data.columns) == ["test 1a", "test2", "test3", "test4"]
    assert data["test2"][1] == "hat2"


def test_column_conversion_error():
    """
    Test that context information (upstream exception message) from column
    conversion error is provided.
    """
    ipac = """\
| col0   |
| double |
 1  2
"""
    with pytest.raises(ValueError) as err:
        ascii.read(ipac, guess=False, format="ipac")
    assert "Column col0 failed to convert:" in str(err.value)

    with pytest.raises(ValueError) as err:
        ascii.read(["a b", "1 2"], guess=False, format="basic", converters={"a": []})
    assert "no converters" in str(err.value)


def test_non_C_locale_with_fast_reader():
    """Test code that forces "C" locale while calling fast reader (#4364)"""
    current = locale.setlocale(locale.LC_ALL)

    try:
        if platform.system() == "Darwin":
            locale.setlocale(locale.LC_ALL, "fr_FR")
        else:
            locale.setlocale(locale.LC_ALL, "fr_FR.utf8")

        for fast_reader in (
            True,
            False,
            {"use_fast_converter": False},
            {"use_fast_converter": True},
        ):
            t = ascii.read(
                ["a b", "1.5 2"], format="basic", guess=False, fast_reader=fast_reader
            )
            assert t["a"].dtype.kind == "f"
    except locale.Error as e:
        pytest.skip(f"Locale error: {e}")
    finally:
        locale.setlocale(locale.LC_ALL, current)


def test_no_units_for_char_columns():
    """Test that a char column of a Table is assigned no unit and not
    a dimensionless unit."""
    t1 = Table([["A"]], names="B")
    out = StringIO()
    ascii.write(t1, out, format="ipac")
    t2 = ascii.read(out.getvalue(), format="ipac", guess=False)
    assert t2["B"].unit is None


def test_initial_column_fill_values():
    """Regression test for #5336, #5338."""

    class TestHeader(ascii.BasicHeader):
        def _set_cols_from_names(self):
            self.cols = [ascii.Column(name=x) for x in self.names]
            # Set some initial fill values
            for col in self.cols:
                col.fill_values = {"--": "0"}

    class Tester(ascii.Basic):
        header_class = TestHeader

    reader = ascii.get_reader(reader_cls=Tester)

    assert (
        reader.read(
            """# Column definition is the first uncommented line
# Default delimiter is the space character.
a b c
# Data starts after the header column definition, blank lines ignored
-- 2 3
4 5 6 """
        )["a"][0]
        is np.ma.masked
    )


def test_latex_no_trailing_backslash():
    """
    Test that latex/aastex file with no trailing backslash can be read.
    """
    lines = r"""
\begin{table}
\begin{tabular}{ccc}
a & b & c \\
1 & 1.0 & c \\ % comment
3\% & 3.0 & e  % comment
\end{tabular}
\end{table}
"""
    dat = ascii.read(lines, format="latex")
    assert dat.colnames == ["a", "b", "c"]
    assert np.all(dat["a"] == ["1", r"3\%"])
    assert np.all(dat["c"] == ["c", "e"])


def text_aastex_no_trailing_backslash():
    lines = r"""
\begin{deluxetable}{ccc}
\tablehead{\colhead{a} & \colhead{b} & \colhead{c}}
\startdata
1 & 1.0 & c \\
2 & 2.0 & d \\ % comment
3\% & 3.0 & e  % comment
\enddata
\end{deluxetable}
"""
    dat = ascii.read(lines, format="aastex")
    assert dat.colnames == ["a", "b", "c"]
    assert np.all(dat["a"] == ["1", r"3\%"])
    assert np.all(dat["c"] == ["c", "e"])


@pytest.mark.parametrize("encoding", ["utf8", "latin1", "cp1252"])
def test_read_with_encoding(tmp_path, encoding):
    data = {"commented_header": "# à b è \n 1 2 héllo", "csv": "à,b,è\n1,2,héllo"}

    testfile = tmp_path / "test.txt"
    for fmt, content in data.items():
        with open(testfile, "w", encoding=encoding) as f:
            f.write(content)

        table = ascii.read(testfile, encoding=encoding)
        assert table.pformat() == [" à   b    è  ", "--- --- -----", "  1   2 héllo"]

        for guess in (True, False):
            table = ascii.read(
                testfile, format=fmt, fast_reader=False, encoding=encoding, guess=guess
            )
            assert table["è"].dtype.kind == "U"
            assert table.pformat() == [
                " à   b    è  ",
                "--- --- -----",
                "  1   2 héllo",
            ]


def test_unsupported_read_with_encoding():
    # Fast reader is not supported, make sure it raises an exception
    with pytest.raises(ascii.ParameterError):
        ascii.read(
            "data/simple3.txt",
            guess=False,
            fast_reader="force",
            encoding="latin1",
            format="fast_csv",
        )


def test_read_chunks_input_types():
    """
    Test chunked reading for different input types: file path, file object,
    and string input.
    """
    fpath = "data/test5.dat"
    t1 = ascii.read(
        fpath,
        header_start=1,
        data_start=3,
    )

    with open(fpath) as fd1, open(fpath) as fd2:
        for fp in (fpath, fd1, fd2.read()):
            t_gen = ascii.read(
                fp,
                header_start=1,
                data_start=3,
                guess=False,
                format="fast_basic",
                fast_reader={"chunk_size": 400, "chunk_generator": True},
            )
            ts = list(t_gen)
            for t in ts:
                for col, col1 in zip(t.columns.values(), t1.columns.values()):
                    assert col.name == col1.name
                    assert col.dtype.kind == col1.dtype.kind

            assert len(ts) == 4
            t2 = table.vstack(ts)
            assert np.all(t1 == t2)

    with open(fpath) as fd1, open(fpath) as fd2:
        for fp in (fpath, fd1, fd2.read()):
            # Now read the full table in chunks
            t3 = ascii.read(
                fp, header_start=1, data_start=3, fast_reader={"chunk_size": 300}
            )
            assert np.all(t1 == t3)


@pytest.mark.parametrize("masked", [True, False])
def test_read_chunks_formats(masked):
    """
    Test different supported formats for chunked reading.
    """
    t1 = simple_table(size=102, cols=10, kinds="fS", masked=masked)
    for i, name in enumerate(t1.colnames):
        t1.rename_column(name, f"col{i + 1}")

    # TO DO commented_header does not currently work due to the special-cased
    # implementation of header parsing.

    for format in "tab", "csv", "no_header", "rdb", "basic":
        out = StringIO()
        ascii.write(t1, out, format=format)
        t_gen = ascii.read(
            out.getvalue(),
            format=format,
            fast_reader={"chunk_size": 400, "chunk_generator": True},
        )
        ts = list(t_gen)
        for t in ts:
            for col, col1 in zip(t.columns.values(), t1.columns.values()):
                assert col.name == col1.name
                assert col.dtype.kind == col1.dtype.kind

        assert len(ts) > 4
        t2 = table.vstack(ts)
        assert np.all(t1 == t2)

        # Now read the full table in chunks
        t3 = ascii.read(out.getvalue(), format=format, fast_reader={"chunk_size": 400})
        assert np.all(t1 == t3)


def test_read_chunks_chunk_size_too_small():
    fpath = "data/test5.dat"
    with pytest.raises(ValueError) as err:
        ascii.read(fpath, header_start=1, data_start=3, fast_reader={"chunk_size": 10})
    assert "no newline found in chunk (chunk_size too small?)" in str(err.value)


def test_read_chunks_table_changes():
    """Column changes type or size between chunks.  This also tests the case with
    no final newline.
    """
    col = ["a b c"] + ["1.12334 xyz a"] * 50 + ["abcdefg 555 abc"] * 50
    table = "\n".join(col)
    t1 = ascii.read(table, guess=False)
    t2 = ascii.read(table, fast_reader={"chunk_size": 100})

    # This also confirms that the dtypes are exactly the same, i.e.
    # the string itemsizes are the same.
    assert np.all(t1 == t2)


def test_read_non_ascii():
    """Test that pure-Python reader is used in case the file contains non-ASCII characters
    in it.
    """
    table = Table.read(["col1, col2", "\u2119, \u01b4", "1, 2"], format="csv")
    assert np.all(table["col1"] == ["\u2119", "1"])
    assert np.all(table["col2"] == ["\u01b4", "2"])


@pytest.mark.parametrize("enable", [True, False, "force"])
def test_kwargs_dict_guess(enable):
    """Test that fast_reader dictionary is preserved through guessing sequence."""
    # Fails for enable=(True, 'force') - #5578
    ascii.read("a\tb\n 1\t2\n3\t 4.0", fast_reader={"enable": enable})
    assert get_read_trace()[-1]["kwargs"]["reader_cls"] is (
        ascii.Tab if (enable is False) else ascii.FastTab
    )
    for k in get_read_trace():
        if not k.get("status", "Disabled").startswith("Disabled"):
            assert k.get("kwargs").get("fast_reader").get("enable") is enable


def _get_lines(rdb):
    lines = ["a a_2 a_1 a a"]
    if rdb:
        lines += ["N N N N N"]
    lines += ["1 2 3 4 5", "10 20 30 40 50"]

    if rdb:
        lines = ["\t".join(line.split()) for line in lines]
    return lines


@pytest.mark.parametrize("rdb", [False, True])
@pytest.mark.parametrize("fast_reader", [False, "force"])
def test_deduplicate_names_basic(rdb, fast_reader):
    """Test that duplicate column names are successfully de-duplicated for the
    basic format.  Skip the case of rdb=True and fast_reader='force' when selecting
    include_names, since that fails and is tested below.
    """
    lines = _get_lines(rdb)

    dat = ascii.read(lines, fast_reader=fast_reader)
    assert dat.colnames == ["a", "a_2", "a_1", "a_3", "a_4"]
    assert len(dat) == 2

    dat = ascii.read(lines, fast_reader=fast_reader, include_names=["a", "a_2", "a_3"])
    assert len(dat) == 2
    assert dat.colnames == ["a", "a_2", "a_3"]
    assert np.all(dat["a"] == [1, 10])
    assert np.all(dat["a_2"] == [2, 20])
    assert np.all(dat["a_3"] == [4, 40])

    dat = ascii.read(
        lines,
        fast_reader=fast_reader,
        names=["b1", "b2", "b3", "b4", "b5"],
        include_names=["b1", "b2", "a_4", "b4"],
    )
    assert len(dat) == 2
    assert dat.colnames == ["b1", "b2", "b4"]
    assert np.all(dat["b1"] == [1, 10])
    assert np.all(dat["b2"] == [2, 20])
    assert np.all(dat["b4"] == [4, 40])

    dat = ascii.read(
        lines,
        fast_reader=fast_reader,
        names=["b1", "b2", "b3", "b4", "b5"],
        exclude_names=["b3", "b5", "a_3", "a_4"],
    )
    assert len(dat) == 2
    assert dat.colnames == ["b1", "b2", "b4"]
    assert np.all(dat["b1"] == [1, 10])
    assert np.all(dat["b2"] == [2, 20])
    assert np.all(dat["b4"] == [4, 40])


def test_include_names_rdb_fast():
    """Test that selecting column names via `include_names` works for the RDB format
    with fast reader. This is testing the fix for a bug identified in #9939.
    """
    lines = _get_lines(True)
    lines[0] = "a\ta_2\ta_1\ta_3\ta_4"
    dat = ascii.read(lines, fast_reader="force", include_names=["a", "a_2", "a_3"])
    assert len(dat) == 2
    assert dat["a"].dtype.kind == "i"
    assert dat["a_2"].dtype.kind == "i"


@pytest.mark.parametrize("fast_reader", [False, "force"])
def test_deduplicate_names_with_types(fast_reader):
    """Test that on selecting column names via `include_names` in the RDB format with
    different types and duplicate column names type assignment is correctly preserved.
    """
    lines = _get_lines(True)
    lines[1] = "N\tN\tN\tS\tS"

    dat = ascii.read(lines, fast_reader=fast_reader, include_names=["a", "a_2", "a_3"])
    assert len(dat) == 2
    assert dat["a_2"].dtype.kind == "i"
    assert dat["a_3"].dtype.kind == "U"

    dat = ascii.read(
        lines,
        fast_reader=fast_reader,
        names=["b1", "b2", "b3", "b4", "b5"],
        include_names=["a1", "a_2", "b1", "b2", "b4"],
    )
    assert len(dat) == 2
    assert dat.colnames == ["b1", "b2", "b4"]
    assert dat["b2"].dtype.kind == "i"
    assert dat["b4"].dtype.kind == "U"


@pytest.mark.parametrize("rdb", [False, True])
@pytest.mark.parametrize("fast_reader", [False, "force"])
def test_set_invalid_names(rdb, fast_reader):
    """
    Test exceptions for invalid (duplicate or `None`) names specified via argument.
    """
    lines = _get_lines(rdb)
    if rdb:
        fmt = "rdb"
    else:
        fmt = "basic"

    with pytest.raises(ValueError) as err:
        ascii.read(
            lines,
            fast_reader=fast_reader,
            format=fmt,
            guess=rdb,
            names=["b1", "b2", "b1", "b4", "b5"],
        )
    assert "Duplicate column names" in str(err.value)

    with pytest.raises(TypeError) as err:
        ascii.read(
            lines,
            fast_reader=fast_reader,
            format=fmt,
            guess=rdb,
            names=["b1", "b2", "b1", None, None],
        )
    assert "Cannot have None for column name" in str(err.value)


def test_read_masked_bool():
    txt = """\
col0 col1
1       1
0       2
True    3
""      4
False   5
"""
    # Reading without converters returns col0 as a string
    dat = ascii.read(txt, format="basic")
    col = dat["col0"]
    assert isinstance(col, MaskedColumn)
    assert col.dtype.kind == "U"
    assert col[0] == "1"

    # Force col0 to be read as bool
    converters = {"col0": [convert_numpy(bool)]}
    dat = ascii.read(txt, format="basic", converters=converters)
    col = dat["col0"]
    assert isinstance(col, MaskedColumn)
    assert col.dtype.kind == "b"
    assert np.all(col.mask == [False, False, False, True, False])
    assert np.all(col == [True, False, True, False, False])


def test_read_converters_wildcard():
    """Test converters where the column name is specified with
    a wildcard.
    """
    converters = {"F*": [ascii.convert_numpy(np.float32)]}
    t = ascii.read(["Fabc Iabc", "1 2"], converters=converters)
    assert np.issubdtype(t["Fabc"].dtype, np.float32)
    assert not np.issubdtype(t["Iabc"].dtype, np.float32)


def test_read_converters_simplified():
    """Test providing io.ascii read converters as type or dtypes instead of
    convert_numpy(type) outputs"""
    t = Table()
    t["a"] = [1, 2]
    t["b"] = [3.5, 4]
    t["c"] = ["True", "False"]
    t["d"] = ["true", "false"]  # Looks kindof like boolean but actually a string
    t["e"] = [5, 6]

    out = StringIO()
    t.write(out, format="ascii.basic")

    converters = {"a": str, "e": np.float32}
    t2 = Table.read(out.getvalue(), format="ascii.basic", converters=converters)
    assert t2.pformat(show_dtype=True) == [
        " a      b      c     d      e   ",
        "str1 float64  str5  str5 float32",
        "---- ------- ----- ----- -------",
        "   1     3.5  True  true     5.0",
        "   2     4.0 False false     6.0",
    ]

    converters = {"a": float, "*": [np.int64, float, bool, str]}
    t2 = Table.read(out.getvalue(), format="ascii.basic", converters=converters)
    assert t2.pformat_all(show_dtype=True) == [
        "   a       b      c     d     e  ",
        "float64 float64  bool  str5 int64",
        "------- ------- ----- ----- -----",
        "    1.0     3.5  True  true     5",
        "    2.0     4.0 False false     6",
    ]

    # Test failures
    for converters in (
        {"*": [int, 1, bool, str]},  # bad converter type
        # Tuple converter where 2nd element is not a subclass of NoType
        {"a": [(int, int)]},
        # Tuple converter with 3 elements not 2
        {"a": [(int, int, int)]},
    ):
        with pytest.raises(ValueError, match="Error: invalid format for converters"):
            t2 = Table.read(
                out.getvalue(), format="ascii.basic", converters=converters, guess=False
            )


def test_read_deprecations():
    def check_warns(func, *args):
        with pytest.warns(AstropyDeprecationWarning) as warns:
            out = func(
                *args,
                Reader=ascii.Basic,
                Inputter=ascii.BaseInputter,
                Outputter=ascii.TableOutputter,
                header_Splitter=ascii.DefaultSplitter,
                data_Splitter=ascii.DefaultSplitter,
            )
            assert len(warns) == 5
            for kwarg in (
                "Reader",
                "Inputter",
                "Outputter",
                "header_Splitter",
                "data_Splitter",
            ):
                msg = f'"{kwarg}" was deprecated'
                assert any(warn.message.args[0].startswith(msg) for warn in warns)
            return out

    tbl = check_warns(ascii.read, ["a b", "1 2"])
    assert tbl.pformat_all() == [" a   b ", "--- ---", "  1   2"]

    reader = check_warns(ascii.get_reader)
    tbl = reader.read(["a b", "1 2"])
    assert tbl.pformat_all() == [" a   b ", "--- ---", "  1   2"]


def test_write_deprecations():
    t = simple_table()
    out = io.StringIO()
    with pytest.warns(AstropyDeprecationWarning, match='"Writer" was deprecated'):
        ascii.write(t, out, Writer=ascii.Csv)
    assert out.getvalue().splitlines() == ["a,b,c", "1,1.0,c", "2,2.0,d", "3,3.0,e"]

    with pytest.warns(AstropyDeprecationWarning, match='"Writer" was deprecated'):
        writer = ascii.get_writer(Writer=ascii.Csv)
    out = io.StringIO()
    writer.write(t, out)
    assert out.getvalue().splitlines() == ["a,b,c", "1,1.0,c", "2,2.0,d", "3,3.0,e"]
