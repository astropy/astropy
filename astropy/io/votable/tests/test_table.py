# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test the conversion to/from astropy.table.
"""

import io
import os
import pathlib
import struct
import warnings

import numpy as np
import pytest

from astropy import units as u
from astropy.io.votable import conf, from_table, is_votable, tree, validate
from astropy.io.votable.converters import Char, UnicodeChar, get_converter
from astropy.io.votable.exceptions import E01, E25, W39, W46, W50, VOWarning
from astropy.io.votable.table import parse, writeto
from astropy.io.votable.tree import Field, VOTableFile
from astropy.table import Column, Table
from astropy.table.table_helpers import simple_table
from astropy.utils.compat.optional_deps import HAS_PYARROW
from astropy.utils.data import (
    get_pkg_data_filename,
    get_pkg_data_fileobj,
    get_pkg_data_path,
)
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
        ("string_array_test", {"datatype": "char", "arraysize": "4*"}),
        ("unsignedByte", {"datatype": "unsignedByte"}),
        ("short", {"datatype": "short"}),
        ("int", {"datatype": "int"}),
        ("intNoNull", {"datatype": "int"}),
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

    for field, (name, d) in zip(t.fields, field_types):
        assert field.ID == name
        assert field.datatype == d["datatype"], (
            f"{name} expected {d['datatype']} but get {field.datatype}"
        )
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
    assert t["Flux1"].unit == u.Unit("erg / (Angstrom cm2 s)")


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
    Issue #927.
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


@pytest.mark.skipif(not HAS_PYARROW, reason="requires pyarrow")
@pytest.mark.parametrize("overwrite", [True, False])
def test_read_write_votable_parquet(tmp_path, overwrite):
    """
    Test to write and read VOTable with Parquet serialization
    """

    # Create some fake data
    number_of_objects = 10
    ids = [f"COSMOS_{ii:03g}" for ii in range(number_of_objects)]
    redshift = np.random.uniform(low=0, high=3, size=number_of_objects)
    mass = np.random.uniform(low=1e8, high=1e10, size=number_of_objects)
    sfr = np.random.uniform(low=1, high=100, size=number_of_objects)
    astropytab = Table([ids, redshift, mass, sfr], names=["id", "z", "mass", "sfr"])

    # Create Column metadata
    column_metadata = {
        "id": {"unit": "", "ucd": "meta.id", "utype": "none"},
        "z": {"unit": "", "ucd": "src.redshift", "utype": "none"},
        "mass": {"unit": "solMass", "ucd": "phys.mass", "utype": "none"},
        "sfr": {"unit": "solMass / yr", "ucd": "phys.SFR", "utype": "none"},
    }

    # Write VOTable with Parquet serialization
    filename = tmp_path / "test_votable_parquet.vot"
    astropytab.write(
        filename,
        column_metadata=column_metadata,
        overwrite=overwrite,
        format="votable.parquet",
    )

    # Check both files are written out
    assert set(os.listdir(tmp_path)) == {
        "test_votable_parquet.vot",
        "test_votable_parquet.vot.parquet",
    }

    # Open created VOTable with Parquet serialization
    with warnings.catch_warnings():
        warnings.simplefilter("always", ResourceWarning)
        votable = parse(filename)

    # Get table out
    votable_table = votable.resources[0].tables[0].array

    # compare tables
    assert (astropytab == votable_table).all()

    # Compare metadata
    # Note: VOTable parses empty units ("") as "---". This is
    # taken into account below by .replace("---","").
    saved_bool = []
    for kk, key in enumerate(column_metadata.keys()):
        for tag in column_metadata[key].keys():
            saved_bool.append(
                column_metadata[key][tag]
                == str(
                    eval(f"votable.resources[0].tables[0].fields[{kk}].{tag}")
                ).replace("---", "")
            )
    assert np.asarray(saved_bool).all()


@pytest.mark.skipif(not HAS_PYARROW, reason="requires pyarrow")
@pytest.mark.parametrize("format", [None, "votable.parquet"])
def test_stored_parquet_votable(format):
    # Ensures that parquet is found as relative to the votable and not the test file
    with warnings.catch_warnings():
        warnings.simplefilter("always", ResourceWarning)
        if format is None:
            stored_votable = Table.read(
                get_pkg_data_filename("data/parquet_binary.xml")
            )
        else:
            stored_votable = Table.read(
                get_pkg_data_filename("data/parquet_binary.xml"), format=format
            )

    assert len(stored_votable) == 10
    assert stored_votable.colnames == ["id", "z", "mass", "sfr"]
    assert stored_votable["sfr"].unit == u.solMass / u.year


def test_write_jybeam_unit(tmp_path, recwarn):
    t = Table(
        {
            "flux": [5 * (u.Jy / u.beam)],
            "foo": [0 * u.Unit("Crab", format="ogip")],
            "bar": [1 * u.def_unit("my_unit")],
        }
    )

    # Crab raises warning outside of VO standards, purely from units.
    assert len(recwarn) == 1
    assert issubclass(recwarn[0].category, u.UnitsWarning)
    assert "Crab" in str(recwarn[0].message)

    filename = tmp_path / "test.xml"
    t.write(filename, format="votable", overwrite=True)

    # Have to use recwarn instead of pytest.warns() because the second run in
    # the double run job does not see these warnings; perhaps something to do
    # with io.votable warning handling. The first run should produce 2 warnings.
    n_warns = len(recwarn)
    assert n_warns in (1, 3)
    if n_warns == 3:
        assert issubclass(recwarn[1].category, W50)
        assert "Crab" in str(recwarn[1].message)
        assert issubclass(recwarn[2].category, W50)
        assert "my_unit" in str(recwarn[2].message)

    t_rt = Table.read(filename, format="votable")

    # No new warnings are emitted on roundtrip read.
    assert len(recwarn) == n_warns
    assert t_rt["flux"].unit == t["flux"].unit

    # These are not VOUnit so while string would match, not same unit instance.
    assert t_rt["foo"].unit.to_string() == t["foo"].unit.to_string()
    assert t_rt["bar"].unit.to_string() == t["bar"].unit.to_string()


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
    Issue #8995.
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
    votable.
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


def test_binary2_single_bounded_char_array():
    """
    Test that variable length char arrays with a max length
    (arraysize="128*") are read correctly in BINARY2 format.
    """
    votable = parse(get_pkg_data_filename("data/binary2_variable_length_char.xml"))
    table = votable.get_first_table()
    astropy_table = table.to_table()

    assert len(astropy_table) == 1
    assert (
        astropy_table[0]["access_format"]
        == "application/x-votable+xml;content=datalink"
    )

    output = io.BytesIO()
    astropy_table.write(output, format="votable", tabledata_format="binary2")
    output.seek(0)

    votable2 = parse(output)
    table2 = votable2.get_first_table()
    astropy_table2 = table2.to_table()

    assert len(astropy_table2) == 1
    assert (
        astropy_table2[0]["access_format"]
        == "application/x-votable+xml;content=datalink"
    )


def test_binary2_bounded_variable_length_char():
    """
    Test that max length variable length char arrays (arraysize="128*")
    are correctly serialized and deserialized correctly
    """
    votable = tree.VOTableFile()
    resource = tree.Resource()
    votable.resources.append(resource)
    table = tree.TableElement(votable)
    resource.tables.append(table)

    table.fields.append(
        tree.Field(
            votable,
            name="bounded_char",
            datatype="char",
            arraysize="128*",
            ID="bounded_char",
        )
    )

    table.create_arrays(2)
    table.array[0]["bounded_char"] = "Short string"
    table.array[1]["bounded_char"] = "Longer string that still fits the 128 characters"

    for format_name in ["tabledata", "binary", "binary2"]:
        bio = io.BytesIO()

        if format_name == "binary2":
            votable.version = "1.3"
            table._config = table._config or {}
            table._config["version_1_3_or_later"] = True

        table.format = format_name
        votable.to_xml(bio)
        bio.seek(0)

        votable2 = parse(bio)
        table2 = votable2.get_first_table()

        assert table2.fields[0].arraysize == "128*", (
            f"Failed to preserve arraysize with format {format_name}"
        )
        assert table2.array[0]["bounded_char"] == "Short string", (
            f"Data mismatch with format {format_name}"
        )
        assert (
            table2.array[1]["bounded_char"]
            == "Longer string that still fits the 128 characters"
        ), f"Data mismatch with format {format_name}"


def test_binary2_bounded_variable_length_char_edge_cases():
    """
    Test edge cases for bounded variable-length char arrays in BINARY2 format
    """
    votable = tree.VOTableFile()
    resource = tree.Resource()
    votable.resources.append(resource)
    table = tree.TableElement(votable)
    resource.tables.append(table)

    table._config = table._config or {}

    table.fields.append(
        tree.Field(
            votable,
            name="bounded_char",
            datatype="char",
            arraysize="10*",
            ID="bounded_char",
        )
    )

    table.create_arrays(3)
    table.array[0]["bounded_char"] = ""
    table.array[1]["bounded_char"] = "1234567890"
    table.array[2]["bounded_char"] = "12345"

    bio = io.BytesIO()

    votable.version = "1.3"
    table._config["version_1_3_or_later"] = True
    table.format = "binary2"

    votable.to_xml(bio)
    bio.seek(0)

    votable2 = parse(bio)
    table2 = votable2.get_first_table()

    assert table2.fields[0].arraysize == "10*", "Failed to preserve arraysize"
    assert table2.array[0]["bounded_char"] == "", "Empty string not preserved"
    assert table2.array[1]["bounded_char"] == "1234567890", (
        "Maximum length string not preserved"
    )
    assert table2.array[2]["bounded_char"] == "12345", (
        "Partial length string not preserved"
    )


def test_binary2_char_fields_vizier_data():
    """
    Test parsing a VOTable from a Vizier query which includes bounded
    variable-length char arrays (e.g., arraysize="7*").

    Related astropy issue:
    https://github.com/astropy/astropy/issues/8737
    """
    sample_file = get_pkg_data_filename("data/vizier_b2_votable.xml")

    votable = parse(sample_file)
    table = votable.get_first_table()

    assert table.fields[7].ID == "CS"
    assert table.fields[7].arraysize == "7*"
    assert table.fields[8].ID == "Morph"
    assert table.fields[8].arraysize == "7*"

    assert table.array[0]["CS"] == ""

    bio = io.BytesIO()
    table.format = "binary2"
    votable.version = "1.3"
    table._config["version_1_3_or_later"] = True
    votable.to_xml(bio)
    bio.seek(0)

    votable2 = parse(bio)
    table2 = votable2.get_first_table()

    assert table2.fields[7].arraysize == "7*"
    assert table2.fields[8].arraysize == "7*"

    assert table2.array[0]["CS"] == ""
    assert table2.array[4]["CS"] == ""


def test_char_bounds_validation():
    """
    Test that char converter correctly validates the bounds.
    """
    votable = tree.VOTableFile()
    field = tree.Field(
        votable, name="bounded_char", datatype="char", arraysize="2*", ID="bounded_char"
    )

    converter = get_converter(field)

    long_string = "abcdefgh"
    with pytest.warns(W46) as record:
        value, mask = converter.parse(long_string, {}, None)

    assert any("char" in str(w.message) and "2" in str(w.message) for w in record)

    if hasattr(converter, "_binoutput_var"):
        with pytest.warns(W46) as record:
            result = converter._binoutput_var(long_string, False)

        assert any("char" in str(w.message) and "2" in str(w.message) for w in record)


def test_binary2_bounded_char_warnings():
    """
    Test that appropriate warnings are issued when strings exceed
    the maximum length in bounded variable-length char arrays.
    """
    votable = tree.VOTableFile()
    resource = tree.Resource()
    votable.resources.append(resource)
    table = tree.TableElement(votable)
    resource.tables.append(table)

    table.fields.append(
        tree.Field(
            votable,
            name="bounded_char",
            datatype="char",
            arraysize="5*",
            ID="bounded_char",
        )
    )

    table.create_arrays(1)
    table.array[0]["bounded_char"] = "abcdefghijklmnopqrstu"

    bio = io.BytesIO()
    votable.version = "1.3"
    table._config = table._config or {}
    table._config["version_1_3_or_later"] = True
    table.format = "binary2"
    votable.to_xml(bio)
    bio.seek(0)

    votable2 = parse(bio)
    with pytest.warns(W46) as record:
        field = votable2.get_first_table().fields[0]
        value, mask = field.converter.parse("Too long value", {}, None)

    assert any("arraysize" in str(w.message) or "5" in str(w.message) for w in record)


def test_binary2_bounded_unichar_valid():
    """
    Test that variable length unicodeChars with a maximum length are correctly
    serialized and deserialized in BINARY2 format.
    """
    votable = tree.VOTableFile()
    resource = tree.Resource()
    votable.resources.append(resource)
    table = tree.TableElement(votable)
    resource.tables.append(table)

    table.fields.append(
        tree.Field(
            votable,
            name="bounded_unichar",
            datatype="unicodeChar",
            arraysize="10*",
            ID="bounded_unichar",
        )
    )

    table.create_arrays(3)
    table.array[0]["bounded_unichar"] = ""
    table.array[1]["bounded_unichar"] = "áéíóú"
    table.array[2]["bounded_unichar"] = "áéíóúÁÉÍÓ"

    for format_name in ["tabledata", "binary", "binary2"]:
        bio = io.BytesIO()

        if format_name == "binary2":
            votable.version = "1.3"
            table._config = table._config or {}
            table._config["version_1_3_or_later"] = True

        table.format = format_name
        votable.to_xml(bio)
        bio.seek(0)

        votable2 = parse(bio)
        table2 = votable2.get_first_table()

        assert table2.fields[0].arraysize == "10*", (
            f"Arraysize not preserved in {format_name} format"
        )

        assert table2.array[0]["bounded_unichar"] == "", (
            f"Empty string not preserved in {format_name} format"
        )
        assert table2.array[1]["bounded_unichar"] == "áéíóú", (
            f"Unicode string not preserved in {format_name} format"
        )
        assert table2.array[2]["bounded_unichar"] == "áéíóúÁÉÍÓ", (
            f"Full-length string not preserved in {format_name} format"
        )


def test_unichar_bounds_validation():
    """
    Test that unicodeChar converter correctly validates bounds and issues warnings.
    """
    votable = tree.VOTableFile()
    field = tree.Field(
        votable,
        name="bounded_unichar",
        datatype="unicodeChar",
        arraysize="5*",
        ID="bounded_unichar",
    )

    converter = get_converter(field)

    long_string = "áéíóúÁÉÍÓÚabcdefgh"
    with pytest.warns(W46) as record:
        value, mask = converter.parse(long_string, {}, None)

    assert any(
        "unicodeChar" in str(w.message) and "5" in str(w.message) for w in record
    )

    if hasattr(converter, "_binoutput_var"):
        with pytest.warns(W46) as record:
            result = converter._binoutput_var(long_string, False)

        assert any(
            "unicodeChar" in str(w.message) and "5" in str(w.message) for w in record
        )


def test_char_invalid_arraysize_raises_e01():
    """
    Test that an invalid arraysize for a char field raises E01.
    """

    class DummyField:
        name = "dummy"
        ID = "dummy_id"
        arraysize = "invalid*"

    field = DummyField()
    config = {}

    with pytest.raises(E01) as excinfo:
        Char(field, config)

    err_msg = str(excinfo.value)
    assert "invalid" in err_msg
    assert "char" in err_msg
    assert "dummy_id" in err_msg


def test_unicodechar_invalid_arraysize_raises_e01():
    """
    Test that an invalid arraysize for a unicodechar field raises E01.
    """

    class DummyField:
        name = "dummy"
        ID = "dummy_id"
        arraysize = "invalid*"

    field = DummyField()
    config = {}

    with pytest.raises(E01) as excinfo:
        UnicodeChar(field, config)

    err_msg = str(excinfo.value)
    assert "invalid" in err_msg
    assert "unicode" in err_msg
    assert "dummy_id" in err_msg


def test_char_exceeding_arraysize():
    """
    Test that appropriate warnings are issued when strings exceed
    the maximum length.
    """
    votable = VOTableFile()

    field = Field(
        votable, name="field_name", datatype="char", arraysize="5", ID="field_id"
    )
    converter = Char(field)
    test_value = "Value that is too long for the field"

    with pytest.warns(W46) as record:
        result, mask = converter.parse(test_value)

    assert any("char" in str(w.message) and "5" in str(w.message) for w in record)

    def mock_read(size):
        if size == 4:
            return struct.pack(">I", 10)
        else:
            return b"0123456789"

    with pytest.warns(W46) as record:
        result, mask = converter._binparse_var(mock_read)

    assert any("char" in str(w.message) and "5" in str(w.message) for w in record)
    assert result == "0123456789"
    assert mask is False


def test_unicodechar_binparse_var_exceeds_arraysize():
    """
    Test that a warning is issued when reading a unicodeChar array in
    binary2 serialization that exceeds the specified max length.
    """
    votable = tree.VOTableFile()
    field = tree.Field(
        votable,
        name="field_name",
        datatype="unicodeChar",
        arraysize="5*",
        ID="field_id",
    )

    converter = get_converter(field)

    def mock_read(size):
        if size == 4:
            return struct.pack(">I", 10)
        else:
            return "0123456789".encode("utf_16_be")

    with pytest.warns(W46) as record:
        result, mask = converter.binparse(mock_read)

    assert any(
        "unicodeChar" in str(w.message) and "5" in str(w.message) for w in record
    )

    assert result == "0123456789"
    assert mask is False


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
