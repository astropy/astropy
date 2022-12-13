# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test the conversion to/from astropy.table
"""
import io
import os
import pathlib

import numpy as np
import pytest

from astropy.config import reload_config, set_temp_config
from astropy.io.votable import conf, from_table, is_votable, tree, validate
from astropy.io.votable.exceptions import E25, W39, VOWarning
from astropy.io.votable.table import parse, writeto
from astropy.table import Column, Table
from astropy.table.table_helpers import simple_table
from astropy.units import Unit
from astropy.utils.data import (
    get_pkg_data_filename,
    get_pkg_data_fileobj,
    get_pkg_data_path,
)
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.utils.misc import _NOT_OVERWRITING_MSG_MATCH


@pytest.fixture
def home_is_data(monkeypatch):
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


@pytest.fixture
def home_is_tmpdir(monkeypatch, tmp_path):
    """
    Pytest fixture to run a test case with tilde-prefixed paths.

    In the tilde-path case, environment variables are temporarily
    modified so that '~' resolves to the temp directory.
    """
    # For Unix
    monkeypatch.setenv("HOME", str(tmp_path))
    # For Windows
    monkeypatch.setenv("USERPROFILE", str(tmp_path))


def test_table(tmp_path):
    # Read the VOTABLE
    with np.errstate(over="ignore"):
        # https://github.com/astropy/astropy/issues/13341
        votable = parse(get_pkg_data_filename("data/regression.xml"))
    table = votable.get_first_table()
    astropy_table = table.to_table()

    for name in table.array.dtype.names:
        assert np.all(astropy_table.mask[name] == table.array.mask[name])

    votable2 = tree.VOTableFile.from_table(astropy_table)
    t = votable2.get_first_table()

    field_types = [
        ("string_test", {"datatype": "char", "arraysize": "*"}),
        ("string_test_2", {"datatype": "char", "arraysize": "10"}),
        ("unicode_test", {"datatype": "unicodeChar", "arraysize": "*"}),
        ("fixed_unicode_test", {"datatype": "unicodeChar", "arraysize": "10"}),
        ("string_array_test", {"datatype": "char", "arraysize": "4"}),
        ("unsignedByte", {"datatype": "unsignedByte"}),
        ("short", {"datatype": "short"}),
        ("int", {"datatype": "int"}),
        ("long", {"datatype": "long"}),
        ("double", {"datatype": "double"}),
        ("float", {"datatype": "float"}),
        ("array", {"datatype": "long", "arraysize": "2*"}),
        ("bit", {"datatype": "bit"}),
        ("bitarray", {"datatype": "bit", "arraysize": "3x2"}),
        ("bitvararray", {"datatype": "bit", "arraysize": "*"}),
        ("bitvararray2", {"datatype": "bit", "arraysize": "3x2*"}),
        ("floatComplex", {"datatype": "floatComplex"}),
        ("doubleComplex", {"datatype": "doubleComplex"}),
        ("doubleComplexArray", {"datatype": "doubleComplex", "arraysize": "*"}),
        ("doubleComplexArrayFixed", {"datatype": "doubleComplex", "arraysize": "2"}),
        ("boolean", {"datatype": "bit"}),
        ("booleanArray", {"datatype": "bit", "arraysize": "4"}),
        ("nulls", {"datatype": "int"}),
        ("nulls_array", {"datatype": "int", "arraysize": "2x2"}),
        ("precision1", {"datatype": "double"}),
        ("precision2", {"datatype": "double"}),
        ("doublearray", {"datatype": "double", "arraysize": "*"}),
        ("bitarray2", {"datatype": "bit", "arraysize": "16"}),
    ]

    for field, type in zip(t.fields, field_types):
        name, d = type
        assert field.ID == name
        assert (
            field.datatype == d["datatype"]
        ), f'{name} expected {d["datatype"]} but get {field.datatype}'
        if "arraysize" in d:
            assert field.arraysize == d["arraysize"]

    # W39: Bit values can not be masked
    with pytest.warns(W39):
        writeto(votable2, str(tmp_path / "through_table.xml"))


def test_read_from_tilde_path(home_is_data):
    # Just test that these run without error for tilde-paths
    path = os.path.join("~", "regression.xml")
    with np.errstate(over="ignore"):
        # https://github.com/astropy/astropy/issues/13341
        votable = parse(path)
        Table.read(path, format="votable", table_id="main_table")


def test_read_through_table_interface(tmp_path):
    with np.errstate(over="ignore"):
        # https://github.com/astropy/astropy/issues/13341
        with get_pkg_data_fileobj("data/regression.xml", encoding="binary") as fd:
            t = Table.read(fd, format="votable", table_id="main_table")

    assert len(t) == 5

    # Issue 8354
    assert t["float"].format is None

    fn = tmp_path / "table_interface.xml"

    # W39: Bit values can not be masked
    with pytest.warns(W39):
        t.write(fn, table_id="FOO", format="votable")

    with open(fn, "rb") as fd:
        t2 = Table.read(fd, format="votable", table_id="FOO")

    assert len(t2) == 5


def test_read_through_table_interface2():
    with np.errstate(over="ignore"):
        # https://github.com/astropy/astropy/issues/13341
        with get_pkg_data_fileobj("data/regression.xml", encoding="binary") as fd:
            t = Table.read(fd, format="votable", table_id="last_table")

    assert len(t) == 0


def test_pass_kwargs_through_table_interface():
    # Table.read() should pass on keyword arguments meant for parse()
    filename = get_pkg_data_filename("data/nonstandard_units.xml")
    t = Table.read(filename, format="votable", unit_format="generic")
    assert t["Flux1"].unit == Unit("erg / (Angstrom cm2 s)")


def test_names_over_ids():
    with get_pkg_data_fileobj("data/names.xml", encoding="binary") as fd:
        votable = parse(fd)

    table = votable.get_first_table().to_table(use_names_over_ids=True)

    assert table.colnames == [
        "Name",
        "GLON",
        "GLAT",
        "RAdeg",
        "DEdeg",
        "Jmag",
        "Hmag",
        "Kmag",
        "G3.6mag",
        "G4.5mag",
        "G5.8mag",
        "G8.0mag",
        "4.5mag",
        "8.0mag",
        "Emag",
        "24mag",
        "f_Name",
    ]


def test_explicit_ids():
    with get_pkg_data_fileobj("data/names.xml", encoding="binary") as fd:
        votable = parse(fd)

    table = votable.get_first_table().to_table(use_names_over_ids=False)

    assert table.colnames == [
        "col1",
        "col2",
        "col3",
        "col4",
        "col5",
        "col6",
        "col7",
        "col8",
        "col9",
        "col10",
        "col11",
        "col12",
        "col13",
        "col14",
        "col15",
        "col16",
        "col17",
    ]


def test_table_read_with_unnamed_tables():
    """
    Issue #927
    """
    with get_pkg_data_fileobj("data/names.xml", encoding="binary") as fd:
        t = Table.read(fd, format="votable")

    assert len(t) == 1


def test_votable_path_object():
    """
    Testing when votable is passed as pathlib.Path object #4412.
    """
    fpath = pathlib.Path(get_pkg_data_filename("data/names.xml"))
    table = parse(fpath).get_first_table().to_table()

    assert len(table) == 1
    assert int(table[0][3]) == 266


def test_from_table_without_mask():
    t = Table()
    c = Column(data=[1, 2, 3], name="a")
    t.add_column(c)
    output = io.BytesIO()
    t.write(output, format="votable")


def test_write_with_format():
    t = Table()
    c = Column(data=[1, 2, 3], name="a")
    t.add_column(c)

    output = io.BytesIO()
    t.write(output, format="votable", tabledata_format="binary")
    obuff = output.getvalue()
    assert b'VOTABLE version="1.4"' in obuff
    assert b"BINARY" in obuff
    assert b"TABLEDATA" not in obuff

    output = io.BytesIO()
    t.write(output, format="votable", tabledata_format="binary2")
    obuff = output.getvalue()
    assert b'VOTABLE version="1.4"' in obuff
    assert b"BINARY2" in obuff
    assert b"TABLEDATA" not in obuff


def test_write_overwrite(tmp_path):
    t = simple_table(3, 3)
    filename = tmp_path / "overwrite_test.vot"
    t.write(filename, format="votable")
    with pytest.raises(OSError, match=_NOT_OVERWRITING_MSG_MATCH):
        t.write(filename, format="votable")
    t.write(filename, format="votable", overwrite=True)


def test_write_tilde_path(home_is_tmpdir):
    fname = os.path.join("~", "output")

    t = Table()
    t["a"] = [1, 2, 3]
    t.write(fname, format="votable", tabledata_format="binary")

    # Ensure the tilde-prefixed path wasn't treated literally
    assert not os.path.exists(fname)

    with open(os.path.expanduser(fname)) as f:
        obuff = f.read()
    assert 'VOTABLE version="1.4"' in obuff
    assert "BINARY" in obuff
    assert "TABLEDATA" not in obuff


@pytest.mark.parametrize("path_format", ["plain", "tilde"])
def test_writeto(path_format, tmp_path, home_is_tmpdir):
    if path_format == "plain":
        # pathlib.Path objects are not accepted by votable.writeto, so convert
        # to a string
        fname = str(tmp_path / "writeto_test.vot")
    else:
        fname = os.path.join("~", "writeto_test.vot")

    t = Table()
    t["a"] = [1, 2, 3]
    vt = from_table(t)
    writeto(vt, fname)

    if path_format == "tilde":
        # Ensure the tilde-prefixed path wasn't treated literally
        assert not os.path.exists(fname)

    with open(os.path.expanduser(fname)) as f:
        obuff = f.read()
    assert 'VOTABLE version="1.4"' in obuff
    assert "BINARY" not in obuff
    assert "TABLEDATA" in obuff


def test_empty_table():
    votable = parse(get_pkg_data_filename("data/empty_table.xml"))
    table = votable.get_first_table()
    table.to_table()


def test_no_field_not_empty_table():
    votable = parse(get_pkg_data_filename("data/no_field_not_empty_table.xml"))
    table = votable.get_first_table()
    assert len(table.fields) == 0
    assert len(table.infos) == 1


def test_no_field_not_empty_table_exception():
    with pytest.raises(E25):
        parse(
            get_pkg_data_filename("data/no_field_not_empty_table.xml"),
            verify="exception",
        )


def test_binary2_masked_strings():
    """
    Issue #8995
    """
    # Read a VOTable which sets the null mask bit for each empty string value.
    votable = parse(get_pkg_data_filename("data/binary2_masked_strings.xml"))
    table = votable.get_first_table()
    astropy_table = table.to_table()

    # Ensure string columns have no masked values and can be written out
    assert not np.any(table.array.mask["epoch_photometry_url"])
    output = io.BytesIO()
    astropy_table.write(output, format="votable")


def test_validate_output_invalid():
    """
    Issue #12603. Test that we get the correct output from votable.validate with an invalid
    votable.
    """

    # A votable with errors
    invalid_votable_filepath = get_pkg_data_filename("data/regression.xml")

    # When output is None, check that validate returns validation output as a string
    validate_out = validate(invalid_votable_filepath, output=None)
    assert isinstance(validate_out, str)
    # Check for known error string
    assert "E02: Incorrect number of elements in array." in validate_out

    # When output is not set, check that validate returns a bool
    validate_out = validate(invalid_votable_filepath)
    assert isinstance(validate_out, bool)
    # Check that validation output is correct (votable is not valid)
    assert validate_out is False


def test_validate_output_valid():
    """
    Issue #12603. Test that we get the correct output from votable.validate with a valid
    votable
    """

    # A valid votable. (Example from the votable standard:
    # https://www.ivoa.net/documents/VOTable/20191021/REC-VOTable-1.4-20191021.html )
    valid_votable_filepath = get_pkg_data_filename("data/valid_votable.xml")

    # When output is None, check that validate returns validation output as a string
    validate_out = validate(valid_votable_filepath, output=None)
    assert isinstance(validate_out, str)
    # Check for known good output string
    assert "astropy.io.votable found no violations" in validate_out

    # When output is not set, check that validate returns a bool
    validate_out = validate(valid_votable_filepath)
    assert isinstance(validate_out, bool)
    # Check that validation output is correct (votable is valid)
    assert validate_out is True


def test_validate_tilde_path(home_is_data):
    validate(os.path.join("~", "valid_votable.xml"))


def test_is_votable_tilde_path(home_is_data):
    assert is_votable(os.path.join("~", "valid_votable.xml"))


class TestVerifyOptions:
    # Start off by checking the default (ignore)

    def test_default(self):
        parse(get_pkg_data_filename("data/gemini.xml"))

    # Then try the various explicit options

    def test_verify_ignore(self):
        parse(get_pkg_data_filename("data/gemini.xml"), verify="ignore")

    def test_verify_warn(self):
        with pytest.warns(VOWarning) as w:
            parse(get_pkg_data_filename("data/gemini.xml"), verify="warn")
        assert len(w) == 24

    def test_verify_exception(self):
        with pytest.raises(VOWarning):
            parse(get_pkg_data_filename("data/gemini.xml"), verify="exception")

    # Make sure the deprecated pedantic option still works for now

    def test_pedantic_false(self):
        with pytest.warns(VOWarning) as w:
            parse(get_pkg_data_filename("data/gemini.xml"), pedantic=False)
        assert len(w) == 25

    def test_pedantic_true(self):
        with pytest.warns(AstropyDeprecationWarning):
            with pytest.raises(VOWarning):
                parse(get_pkg_data_filename("data/gemini.xml"), pedantic=True)

    # Make sure that the default behavior can be set via configuration items

    def test_conf_verify_ignore(self):
        with conf.set_temp("verify", "ignore"):
            parse(get_pkg_data_filename("data/gemini.xml"))

    def test_conf_verify_warn(self):
        with conf.set_temp("verify", "warn"):
            with pytest.warns(VOWarning) as w:
                parse(get_pkg_data_filename("data/gemini.xml"))
            assert len(w) == 24

    def test_conf_verify_exception(self):
        with conf.set_temp("verify", "exception"):
            with pytest.raises(VOWarning):
                parse(get_pkg_data_filename("data/gemini.xml"))

    # And make sure the old configuration item will keep working

    def test_conf_pedantic_false(self, tmp_path):
        with set_temp_config(tmp_path):
            with open(tmp_path / "astropy" / "astropy.cfg", "w") as f:
                f.write("[io.votable]\npedantic = False")

            reload_config("astropy.io.votable")

            with pytest.warns(VOWarning) as w:
                parse(get_pkg_data_filename("data/gemini.xml"))
            assert len(w) == 25

    def test_conf_pedantic_true(self, tmp_path):
        with set_temp_config(tmp_path):
            with open(tmp_path / "astropy" / "astropy.cfg", "w") as f:
                f.write("[io.votable]\npedantic = True")

            reload_config("astropy.io.votable")

            with pytest.warns(AstropyDeprecationWarning):
                with pytest.raises(VOWarning):
                    parse(get_pkg_data_filename("data/gemini.xml"))
