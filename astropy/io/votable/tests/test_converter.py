# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io
import warnings

# THIRD-PARTY
import numpy as np
import numpy.ma as ma
import pytest
from numpy.testing import assert_array_equal

# LOCAL
from astropy.io.votable import converters, exceptions, tree
from astropy.io.votable.table import parse, parse_single_table
from astropy.table import Table
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


def test_char_fix_array():
    config = {"verify": "exception"}
    field = tree.Field(None, name="c", arraysize="3x2", datatype="char", config=config)
    c = converters.get_converter(field, config=config)
    assert c.output([["Foo", "Bar"], ["Foo", "Bar"]], True) == "            "

    assert (
        c.output([["Foo", "Bar"], ["Foo", "Bar"]], [[False, False], [False, False]])
        == "FooBarFooBar"
    )

    assert (
        c.output([["Foo", "Bar"], ["Foo", "Bar"]], [[True, False], [True, False]])
        == "   Bar   Bar"
    )


def test_char_variable_array():
    config = {"verify": "exception"}
    field = tree.Field(None, name="c", arraysize="3x*", datatype="char", config=config)
    c = converters.get_converter(field, config=config)
    assert c.output([["Foo", "Bar"], ["Foo", "Bar"]], True) == "            "

    assert (
        c.output([["Foo", "Bar"], ["Foo", "Bar"]], [[False, False], [False, False]])
        == "FooBarFooBar"
    )
    assert (
        c.output([["Foo", "Bar"], ["Foo", "Bar"]], [[True, False], [True, False]])
        == "   Bar   Bar"
    )


def test_unicodeChar_fix_array():
    config = {"verify": "exception"}
    field = tree.Field(
        None, name="c", arraysize="3x2", datatype="unicodeChar", config=config
    )
    c = converters.get_converter(field, config=config)
    assert c.output([["áéí", "ÁÉÍ"], ["áéí", "ÁÉÍ"]], True) == "            "

    assert (
        c.output([["áéí", "ÁÉÍ"], ["áéí", "ÁÉÍ"]], [[False, False], [False, False]])
        == "áéíÁÉÍáéíÁÉÍ"
    )
    assert (
        c.output([["áéí", "ÁÉÍ"], ["áéí", "ÁÉÍ"]], [[True, False], [True, False]])
        == "   ÁÉÍ   ÁÉÍ"
    )


def test_unicodeChar_variable_array():
    config = {"verify": "exception"}

    field = tree.Field(
        None, name="c", arraysize="3x*", datatype="unicodeChar", config=config
    )

    c = converters.get_converter(field, config=config)
    assert c.output([["Foo", "Bar"], ["Foo", "Bar"]], True) == "            "

    assert (
        c.output([["Foo", "Bar"], ["Foo", "Bar"]], [[False, False], [False, False]])
        == "FooBarFooBar"
    )
    assert (
        c.output([["Foo", "Bar"], ["Foo", "Bar"]], [[True, False], [True, False]])
        == "   Bar   Bar"
    )

    c = converters.get_converter(field, config=config)
    assert c.output([["你好你", "αβγ"], ["你好你", "αβγ"]], True) == "            "

    assert (
        c.output(
            [["你好你", "αβγ"], ["你好你", "αβγ"]], [[False, False], [False, False]]
        )
        == "你好你αβγ你好你αβγ"
    )
    assert (
        c.output([["你好你", "αβγ"], ["你好你", "αβγ"]], [[True, False], [True, False]])
        == "   αβγ   αβγ"
    )


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


def test_bounded_variable_size_unicode():
    # regression test for #18075
    field = tree.Field(None, name="unicode", datatype="unicodeChar", arraysize="20*")
    unicode_converter = converters.get_converter(field)
    assert unicode_converter.parse("XXX") == ("XXX", False)


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


def test_unicode_as_char_fix_array():
    config = {"verify": "exception"}
    field = tree.Field(
        None, name="unicode_in_char", datatype="char", arraysize="3x2", config=config
    )
    c = converters.get_converter(field, config=config)

    # Test parsing.
    c.parse("XYZABC")  # ASCII succeeds
    with pytest.warns(
        exceptions.W55,
        match=(
            r'FIELD \(unicode_in_char\) has datatype="char" but contains non-ASCII'
            r" value"
        ),
    ):
        c.parse("złazła")  # non-ASCII

    # Test output.
    c.output("XYZABC", False)  # ASCII str succeeds
    c.output(b"XYZABC", False)  # ASCII bytes succeeds
    value = "zła zła"
    value_bytes = value.encode("utf-8")
    with pytest.warns(exceptions.E24, match=r"E24: Attempt to write non-ASCII value"):
        c.output(value, False)  # non-ASCII str raises
    with pytest.warns(exceptions.E24, match=r"E24: Attempt to write non-ASCII value"):
        c.output(value_bytes, False)  # non-ASCII bytes raises


def test_unicode_as_char_variable_array():
    config = {"verify": "exception"}
    field = tree.Field(
        None, name="unicode_in_char", datatype="char", arraysize="3x*", config=config
    )
    c = converters.get_converter(field, config=config)

    # Test parsing.
    c.parse("XYZABC")  # ASCII succeeds
    with pytest.warns(
        exceptions.W55,
        match=(
            r'FIELD \(unicode_in_char\) has datatype="char" but contains non-ASCII'
            r" value"
        ),
    ):
        c.parse("złazła")  # non-ASCII

    # Test output.
    c.output("XYZABC", False)  # ASCII str succeeds
    c.output(b"XYZABC", False)  # ASCII bytes succeeds
    value = "złazła"
    value_bytes = value.encode("utf-8")
    with pytest.warns(exceptions.E24, match=r"E24: Attempt to write non-ASCII value"):
        c.output(value, False)  # non-ASCII str raises
    with pytest.warns(exceptions.E24, match=r"E24: Attempt to write non-ASCII value"):
        c.output(value_bytes, False)  # non-ASCII bytes raises


def test_unicode_as_char_binary_fix_array():
    config = {"verify": "exception"}

    field = tree.Field(
        None, name="unicode_in_char", datatype="char", arraysize="3x2", config=config
    )
    converter = converters.get_converter(field, config=config)
    converter.binoutput(["abc", "def"], [False, False])  # ASCII succeeds
    with pytest.raises(exceptions.E24, match=r"E24: Attempt to write non-ASCII value"):
        converter.binoutput(["zła", "zła"], [False, False])

    field = tree.Field(
        None, name="unicode_in_char", datatype="char", arraysize="3", config=config
    )
    converter = converters.get_converter(field, config=config)
    converter.binoutput("xyz abc", False)
    with pytest.raises(exceptions.E24, match=r"E24: Attempt to write non-ASCII value"):
        converter.binoutput("zła", False)


def test_unicode_as_char_binary_variable_array():
    config = {"verify": "exception"}

    field = tree.Field(
        None, name="unicode_in_char", datatype="char", arraysize="3x*", config=config
    )
    converter = converters.get_converter(field, config=config)
    converter.binoutput(
        np.ma.array(["abc", "def"], mask=[False, False]), False
    )  # ASCII succeeds
    with pytest.raises(exceptions.E24, match=r"E24: Attempt to write non-ASCII value"):
        converter.binoutput(np.ma.array(["zła", "zła"], mask=[False, False]), False)

    field = tree.Field(
        None, name="unicode_in_char", datatype="char", arraysize="3", config=config
    )
    converter = converters.get_converter(field, config=config)
    converter.binoutput("xyz abc", False)
    with pytest.raises(exceptions.E24, match=r"E24: Attempt to write non-ASCII value"):
        converter.binoutput("zła", False)


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

    config = {"verify": "exception"}
    field = tree.Field(
        None, name="c", arraysize="2x2", datatype="double", config=config
    )
    c = converters.get_converter(field, config=config)
    assert c.output([[1, 2], [3, 4]], [[False, False], [False, False]]) == "1 2 3 4"

    config = {"verify": "exception"}
    field = tree.Field(
        None, name="c", arraysize="2x2", datatype="double", config=config
    )
    c = converters.get_converter(field, config=config)
    assert c.output([[1, 2], [3, 4]], [[True, False], [True, True]]) == "NaN 2 NaN NaN"


def test_complex_array_vararray():
    config = {"verify": "exception", "version_1_3_or_later": True}
    field = tree.Field(
        None, name="c", datatype="floatComplex", arraysize="2x3*", config=config
    )
    c = converters.get_converter(field, config=config)
    with pytest.raises(exceptions.E02):
        c.parse("2 3 4 5 6", config=config)


def test_int_array_vararray():
    config = {"verify": "exception", "version_1_3_or_later": True}

    field = tree.Field(None, name="c", datatype="int", arraysize="2x3", config=config)
    c = converters.get_converter(field, config=config)
    with pytest.raises(exceptions.E02):
        c.parse("1 2 3 4 5", config=config)


def test_complex_array_vararray2():
    config = {"verify": "exception", "version_1_3_or_later": True}
    field = tree.Field(
        None, name="c", datatype="floatComplex", arraysize="2x3*", config=config
    )
    c = converters.get_converter(field, config=config)
    x = c.parse("", config=config)
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


def test_overflow_binary_converter():
    """Regression for #20302: overflow for binary conversion
    We test that we can parse a binary votable that contains all numeric types with
    extreme values.
    """
    content = b"""<?xml version="1.0" encoding="utf-8"?>
    <VOTABLE version="1.4" xmlns="http://www.ivoa.net/xml/VOTable/v1.3" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.ivoa.net/xml/VOTable/v1.3 http://www.ivoa.net/xml/VOTable/VOTable-1.4.xsd">
    <RESOURCE type="results">
    <TABLE nrows="3">
    <FIELD ID="bit_column" datatype="bit" name="bit_column"/>
    <FIELD ID="unsignedByte_column" datatype="unsignedByte" name="unsignedByte_column"/>
    <FIELD ID="short_column" datatype="short" name="short_column">
        <VALUES null="-32768"/>
    </FIELD>
    <FIELD ID="int_column" datatype="int" name="int_column">
        <VALUES null="-2147483648"/>
    </FIELD>
    <FIELD ID="long_column" datatype="long" name="long_column">
        <VALUES null="-9223372036854775808"/>
    </FIELD>
    <FIELD ID="float_column" datatype="float" name="float_column">
        <VALUES null="-3.4028235e+38"/>
    </FIELD>
    <FIELD ID="double_column" datatype="double" name="double_column">
        <VALUES null="-1.7976931348623157e+308"/>
    </FIELD>
    <DATA>
        <BINARY>
        <STREAM encoding="base64">
    AACAAIAAAACAAAAAAAAAAP9/////7////////wj/f/9/////f/////////9/f///f+////////8IgAPoAAAAKgAAAAJUC+QAQEj1w0AFvwqLFFdp     </STREAM>
        </BINARY>
    </DATA>
    </TABLE>
    </RESOURCE>
    </VOTABLE>"""
    table = parse_single_table(io.BytesIO(content))
    assert list(table.to_table()["bit_column"]) == [False, True, True]


def test_oversize_char_array():

    # exception
    config = {"verify": "exception", "version_1_3_or_later": True}

    field = tree.Field(None, name="c", datatype="char", arraysize="3x2", config=config)
    c = converters.get_converter(field, config=config)
    with pytest.raises(exceptions.E02) as exc_info:
        x = c.parse("abc", config=config)

    assert (
        "Incorrect number of elements in array. Expected multiple of 2, got 1"
        in str(exc_info.value)
    )

    with pytest.raises(exceptions.E02) as exc_info:
        x = c.parse("ab", config=config)

    assert (
        "Incorrect number of elements in array. Expected multiple of 2, got 1"
        in str(exc_info.value)
    )

    # warn
    config = {"verify": "warn", "version_1_3_or_later": True}

    field = tree.Field(None, name="c", datatype="char", arraysize="3x2", config=config)
    c = converters.get_converter(field, config=config)
    with pytest.warns(exceptions.E02) as exc_info:
        x = c.parse("abc", config=config)

    assert len(x) == 2

    with pytest.raises(exceptions.E02) as exc_info:
        x = c.parse("ab", config=config)

    assert (
        "Incorrect number of elements in array. Expected multiple of 2, got 1"
        in str(exc_info.value)
    )
    # Ignore
    config = {"verify": "ignore", "version_1_3_or_later": True}

    field = tree.Field(None, name="c", datatype="char", arraysize="3x2", config=config)
    c = converters.get_converter(field, config=config)

    with warnings.catch_warnings(record=True) as exc_info:
        warnings.resetwarnings()
        warnings.simplefilter("always")
        x = c.parse("abc", config=config)

        assert len(x) == 2


def test_oversize_unicodeChar_array():
    # exception
    config = {"verify": "exception", "version_1_3_or_later": True}

    field = tree.Field(
        None, name="c", datatype="unicodeChar", arraysize="3x2", config=config
    )
    c = converters.get_converter(field, config=config)

    with pytest.raises(exceptions.E02) as exc_info:
        x = c.parse("αβγ", config=config)

    assert (
        "Incorrect number of elements in array. Expected multiple of 2, got 1"
        in str(exc_info.value)
    )

    with pytest.raises(exceptions.E02) as exc_info:
        x = c.parse("αβ", config=config)

    assert (
        "Incorrect number of elements in array. Expected multiple of 2, got 1"
        in str(exc_info.value)
    )

    # warn
    config = {"verify": "warn", "version_1_3_or_later": True}

    field = tree.Field(
        None, name="c", datatype="unicodeChar", arraysize="3x2", config=config
    )
    c = converters.get_converter(field, config=config)

    with pytest.warns(exceptions.E02) as exc_info:
        x = c.parse("αβγ", config=config)

    assert len(x) == 2

    with pytest.raises(exceptions.E02) as exc_info:
        x = c.parse("αβ", config=config)

    assert (
        "Incorrect number of elements in array. Expected multiple of 2, got 1"
        in str(exc_info.value)
    )

    # Ignore
    config = {"verify": "ignore", "version_1_3_or_later": True}

    field = tree.Field(
        None, name="c", datatype="unicodeChar", arraysize="3x2", config=config
    )
    c = converters.get_converter(field, config=config)

    with warnings.catch_warnings(record=True) as exc_info:
        warnings.simplefilter("always")
        x = c.parse("αβγ", config=config)

        assert len(x) == 2


@pytest.mark.parametrize(
    "dtype, expected",
    [
        (np.dtype("float64"), {"datatype": "double"}),
        (np.dtype("float32"), {"datatype": "float"}),
        (np.dtype("bool"), {"datatype": "bit"}),
        (np.dtype("uint8"), {"datatype": "unsignedByte"}),
        (np.dtype("int16"), {"datatype": "short"}),
        (np.dtype("int32"), {"datatype": "int"}),
        (np.dtype("int64"), {"datatype": "long"}),
        (np.dtype("complex64"), {"datatype": "floatComplex"}),
        (np.dtype("complex128"), {"datatype": "doubleComplex"}),
    ],
)
def test_numpy_to_votable_dtype_scalar(dtype, expected):
    result = converters.numpy_to_votable_dtype(dtype, ())
    assert result == expected


@pytest.mark.parametrize(
    "dtype, expected",
    [
        (
            np.dtype("S10"),
            {
                "datatype": "char",
                "arraysize": "10",
            },
        ),
        (
            np.dtype("U10"),
            {
                "datatype": "unicodeChar",
                "arraysize": "10",
            },
        ),
        (
            np.dtype("S1"),
            {
                "datatype": "char",
                "arraysize": "1",
            },
        ),
        (
            np.dtype("U1"),
            {
                "datatype": "unicodeChar",
                "arraysize": "1",
            },
        ),
    ],
)
def test_numpy_to_votable_dtype_strings(dtype, expected):
    result = converters.numpy_to_votable_dtype(dtype, ())
    assert result == expected


@pytest.mark.parametrize(
    "dtype, shape, expected",
    [
        (
            np.dtype("int32"),
            (3,),
            {"datatype": "int", "arraysize": "3"},
        ),
        (
            np.dtype("float64"),
            (2, 3),
            {"datatype": "double", "arraysize": "2x3"},
        ),
        (
            np.dtype("float32"),
            (2, 3, 4),
            {"datatype": "float", "arraysize": "2x3x4"},
        ),
    ],
)
def test_numpy_to_votable_dtype_arrays(dtype, shape, expected):
    result = converters.numpy_to_votable_dtype(dtype, shape)
    assert result == expected


def test_numpy_to_votable_dtype_unsupported():
    dtype = np.dtype("object")

    with pytest.raises(TypeError, match="can not be represented in VOTable"):
        converters.numpy_to_votable_dtype(dtype, ())


def test_numpy_to_votable_dtype_empty_shape():
    result = converters.numpy_to_votable_dtype(np.dtype("int32"), ())
    assert result == {"datatype": "int"}


@pytest.mark.parametrize(
    ("values", "dtype", "expected_datatype", "expected_values"),
    [
        pytest.param(
            [["abc", "de"]], None, "unicodeChar", ["abc", "de "], id="ascii-default"
        ),
        pytest.param(
            [["αβγ", "αβ"]], None, "unicodeChar", ["αβγ", "αβ "], id="unicode-default"
        ),
        pytest.param(
            [["abc", "de"]], "U", "unicodeChar", ["abc", "de "], id="ascii-unicode"
        ),
        pytest.param(
            [["αβγ", "αβ"]], "U", "unicodeChar", ["αβγ", "αβ "], id="unicode-unicode"
        ),
        pytest.param([["abc", "de"]], "S", "char", ["abc", "de "], id="ascii-bytes"),
        pytest.param([[b"abc", b"de"]], "S", "char", ["abc", "de "], id="bytes"),
    ],
)
def test_numpy_to_votable_dtype_string_array(
    values, dtype, expected_values, expected_datatype
):
    astropy_table = Table({"c": np.array(values, dtype=dtype)})
    col = astropy_table["c"]

    result = converters.numpy_to_votable_dtype(col.dtype, col.shape[1:])

    assert result == {
        "datatype": expected_datatype,
        "arraysize": "3x2",
    }

    output = io.BytesIO()
    astropy_table.write(output, format="votable")  # , tabledata_format="binary2")
    output.seek(0)

    votable2 = parse(output)
    table2 = votable2.get_first_table()
    astropy_table2 = table2.to_table()
    actual_row = astropy_table2[0]

    assert len(astropy_table2) == 1
    actual = astropy_table2[0]

    assert_array_equal(ma.getdata(actual["c"]), expected_values)
