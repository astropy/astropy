# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io

# THIRD-PARTY
import numpy as np
import pytest
from numpy.testing import assert_array_equal

# LOCAL
from astropy.io.votable import converters, exceptions, tree
from astropy.io.votable.table import parse_single_table
from astropy.utils.data import get_pkg_data_filename


def test_invalid_arraysize():
    with pytest.raises(exceptions.E13):
        field = tree.Field(None, name="broken", datatype="char", arraysize="foo")
        converters.get_converter(field)


def test_oversize_char():
    config = {"verify": "exception"}
    with pytest.warns(exceptions.W47) as w:
        field = tree.Field(None, name="c", datatype="char", config=config)
        c = converters.get_converter(field, config=config)
    assert len(w) == 1

    with pytest.warns(exceptions.W46) as w:
        c.parse("XXX")
    assert len(w) == 1


def test_char_mask():
    config = {"verify": "exception"}
    field = tree.Field(None, name="c", arraysize="1", datatype="char", config=config)
    c = converters.get_converter(field, config=config)
    assert c.output("Foo", True) == ""


def test_oversize_unicode():
    config = {"verify": "exception"}
    with pytest.warns(exceptions.W46) as w:
        field = tree.Field(
            None, name="c2", datatype="unicodeChar", arraysize="1", config=config
        )
        c = converters.get_converter(field, config=config)
        c.parse("XXX")
    assert len(w) == 1


def test_unicode_mask():
    config = {"verify": "exception"}
    field = tree.Field(
        None, name="c", arraysize="1", datatype="unicodeChar", config=config
    )
    c = converters.get_converter(field, config=config)
    assert c.output("Foo", True) == ""


def test_unicode_as_char():
    config = {"verify": "exception"}
    field = tree.Field(
        None, name="unicode_in_char", datatype="char", arraysize="*", config=config
    )
    c = converters.get_converter(field, config=config)

    # Test parsing.
    c.parse("XYZ")  # ASCII succeeds
    with pytest.warns(
        exceptions.W55,
        match=(
            r'FIELD \(unicode_in_char\) has datatype="char" but contains non-ASCII'
            r" value"
        ),
    ):
        c.parse("zła")  # non-ASCII

    # Test output.
    c.output("XYZ", False)  # ASCII str succeeds
    c.output(b"XYZ", False)  # ASCII bytes succeeds
    value = "zła"
    value_bytes = value.encode("utf-8")
    with pytest.warns(exceptions.E24, match=r"E24: Attempt to write non-ASCII value"):
        c.output(value, False)  # non-ASCII str raises
    with pytest.warns(exceptions.E24, match=r"E24: Attempt to write non-ASCII value"):
        c.output(value_bytes, False)  # non-ASCII bytes raises


def test_unicode_as_char_binary():
    config = {"verify": "exception"}

    field = tree.Field(
        None, name="unicode_in_char", datatype="char", arraysize="*", config=config
    )
    c = converters.get_converter(field, config=config)
    c._binoutput_var("abc", False)  # ASCII succeeds
    with pytest.raises(exceptions.E24, match=r"E24: Attempt to write non-ASCII value"):
        c._binoutput_var("zła", False)

    field = tree.Field(
        None, name="unicode_in_char", datatype="char", arraysize="3", config=config
    )
    c = converters.get_converter(field, config=config)
    c._binoutput_fixed("xyz", False)
    with pytest.raises(exceptions.E24, match=r"E24: Attempt to write non-ASCII value"):
        c._binoutput_fixed("zła", False)


def test_wrong_number_of_elements():
    config = {"verify": "exception"}
    field = tree.Field(None, name="c", datatype="int", arraysize="2x3*", config=config)
    c = converters.get_converter(field, config=config)
    with pytest.raises(exceptions.E02):
        c.parse("2 3 4 5 6")


def test_float_mask():
    config = {"verify": "exception"}
    field = tree.Field(None, name="c", datatype="float", config=config)
    c = converters.get_converter(field, config=config)
    assert c.parse("") == (c.null, True)
    with pytest.raises(ValueError):
        c.parse("null")


def test_float_mask_permissive():
    config = {"verify": "ignore"}
    field = tree.Field(None, name="c", datatype="float", config=config)

    # config needs to be also passed into parse() to work.
    # https://github.com/astropy/astropy/issues/8775
    c = converters.get_converter(field, config=config)
    assert c.parse("null", config=config) == (c.null, True)


@pytest.mark.usefixtures("without_legacy_printoptions")
def test_double_array():
    config = {"verify": "exception", "version_1_3_or_later": True}
    field = tree.Field(None, name="c", datatype="double", arraysize="3", config=config)
    data = (1.0, 2.0, 3.0)
    c = converters.get_converter(field, config=config)
    assert c.output(1.0, False) == "1"
    assert c.output(1.0, [False, False]) == "1"
    assert c.output(data, False) == "1 2 3"
    assert c.output(data, [False, False, False]) == "1 2 3"
    assert c.output(data, [False, False, True]) == "1 2 NaN"
    assert c.output(data, [False, False]) == "1 2"

    a = c.parse("1 2 3", config=config)
    assert_array_equal(a[0], data)
    assert_array_equal(a[1], False)

    with pytest.raises(exceptions.E02):
        c.parse("1", config=config)

    with pytest.raises(AttributeError), pytest.warns(exceptions.E02):
        c.parse("1")

    with pytest.raises(exceptions.E02):
        c.parse("2 3 4 5 6", config=config)

    with pytest.warns(exceptions.E02):
        a = c.parse("2 3 4 5 6")

    assert_array_equal(a[0], [2, 3, 4])
    assert_array_equal(a[1], False)


def test_complex_array_vararray():
    config = {"verify": "exception"}
    field = tree.Field(
        None, name="c", datatype="floatComplex", arraysize="2x3*", config=config
    )
    c = converters.get_converter(field, config=config)
    with pytest.raises(exceptions.E02):
        c.parse("2 3 4 5 6")


def test_complex_array_vararray2():
    config = {"verify": "exception"}
    field = tree.Field(
        None, name="c", datatype="floatComplex", arraysize="2x3*", config=config
    )
    c = converters.get_converter(field, config=config)
    x = c.parse("")
    assert len(x[0]) == 0


def test_complex_array_vararray3():
    config = {"verify": "exception"}
    field = tree.Field(
        None, name="c", datatype="doubleComplex", arraysize="2x3*", config=config
    )
    c = converters.get_converter(field, config=config)
    x = c.parse("1 2 3 4 5 6 7 8 9 10 11 12")
    assert len(x) == 2
    assert np.all(x[0][0][0] == complex(1, 2))


def test_complex_vararray():
    config = {"verify": "exception"}
    field = tree.Field(
        None, name="c", datatype="doubleComplex", arraysize="*", config=config
    )
    c = converters.get_converter(field, config=config)
    x = c.parse("1 2 3 4")
    assert len(x) == 2
    assert x[0][0] == complex(1, 2)


def test_complex():
    config = {"verify": "exception"}
    field = tree.Field(None, name="c", datatype="doubleComplex", config=config)
    c = converters.get_converter(field, config=config)
    with pytest.raises(exceptions.E03):
        c.parse("1 2 3")


def test_bit():
    config = {"verify": "exception"}
    field = tree.Field(None, name="c", datatype="bit", config=config)
    c = converters.get_converter(field, config=config)
    with pytest.raises(exceptions.E04):
        c.parse("T")


def test_bit_mask():
    config = {"verify": "exception"}
    with pytest.warns(exceptions.W39) as w:
        field = tree.Field(None, name="c", datatype="bit", config=config)
        c = converters.get_converter(field, config=config)
        c.output(True, True)
    assert len(w) == 1


def test_boolean():
    config = {"verify": "exception"}
    field = tree.Field(None, name="c", datatype="boolean", config=config)
    c = converters.get_converter(field, config=config)
    with pytest.raises(exceptions.E05):
        c.parse("YES")


def test_boolean_array():
    config = {"verify": "exception"}
    field = tree.Field(None, name="c", datatype="boolean", arraysize="*", config=config)
    c = converters.get_converter(field, config=config)
    r, mask = c.parse("TRUE FALSE T F 0 1")
    assert_array_equal(r, [True, False, True, False, False, True])


def test_invalid_type():
    config = {"verify": "exception"}
    with pytest.raises(exceptions.E06):
        field = tree.Field(None, name="c", datatype="foobar", config=config)
        converters.get_converter(field, config=config)


def test_precision():
    config = {"verify": "exception"}

    field = tree.Field(None, name="c", datatype="float", precision="E4", config=config)
    c = converters.get_converter(field, config=config)
    assert c.output(266.248, False) == "266.2"

    field = tree.Field(None, name="c", datatype="float", precision="F4", config=config)
    c = converters.get_converter(field, config=config)
    assert c.output(266.248, False) == "266.2480"


def test_integer_overflow():
    config = {"verify": "exception"}

    field = tree.Field(None, name="c", datatype="int", config=config)
    c = converters.get_converter(field, config=config)
    with pytest.raises(exceptions.W51):
        c.parse("-2208988800", config=config)


@pytest.mark.usefixtures("without_legacy_printoptions")
def test_float_default_precision():
    config = {"verify": "exception"}

    field = tree.Field(None, name="c", datatype="float", arraysize="4", config=config)
    c = converters.get_converter(field, config=config)
    assert (
        c.output([1, 2, 3, 8.9990234375], [False, False, False, False])
        == "1 2 3 8.9990234375"
    )


def test_vararray():
    votable = tree.VOTableFile()
    resource = tree.Resource()
    votable.resources.append(resource)
    table = tree.TableElement(votable)
    resource.tables.append(table)

    tabarr = []
    heads = ["headA", "headB", "headC"]
    types = ["char", "double", "int"]

    vals = [["A", 1.0, 2], ["B", 2.0, 3], ["C", 3.0, 4]]
    for i in range(len(heads)):
        tabarr.append(
            tree.Field(votable, name=heads[i], datatype=types[i], arraysize="*")
        )

    table.fields.extend(tabarr)
    table.create_arrays(len(vals))
    for i in range(len(vals)):
        values = tuple(vals[i])
        table.array[i] = values
    buff = io.BytesIO()
    votable.to_xml(buff)


def test_gemini_v1_2():
    """
    see Pull Request 4782 or Issue 4781 for details.
    """
    table = parse_single_table(get_pkg_data_filename("data/gemini.xml"))
    assert table is not None

    tt = table.to_table()
    assert (
        tt["access_url"][0]
        == "http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/GEMINI/"
        "S20120515S0064?runid=bx9b1o8cvk1qesrt"
    )
