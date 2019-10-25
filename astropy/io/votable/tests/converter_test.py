# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io

# THIRD-PARTY
import numpy as np
from numpy.testing import assert_array_equal

# LOCAL
from astropy.io.votable import converters
from astropy.io.votable import exceptions
from astropy.io.votable import tree

from astropy.io.votable.table import parse_single_table
from astropy.tests.helper import raises, catch_warnings
from astropy.utils.data import get_pkg_data_filename


@raises(exceptions.E13)
def test_invalid_arraysize():
    field = tree.Field(
        None, name='broken', datatype='char', arraysize='foo')
    converters.get_converter(field)


def test_oversize_char():
    config = {'verify': 'exception'}
    with catch_warnings(exceptions.W47) as w:
        field = tree.Field(
            None, name='c', datatype='char',
            config=config)
        c = converters.get_converter(field, config=config)
    assert len(w) == 1

    with catch_warnings(exceptions.W46) as w:
        c.parse("XXX")
    assert len(w) == 1


def test_char_mask():
    config = {'verify': 'exception'}
    field = tree.Field(None, name='c', arraysize='1', datatype='char',
                       config=config)
    c = converters.get_converter(field, config=config)
    assert c.output("Foo", True) == ''


def test_oversize_unicode():
    config = {'verify': 'exception'}
    with catch_warnings(exceptions.W46) as w:
        field = tree.Field(
            None, name='c2', datatype='unicodeChar',
            config=config)
        c = converters.get_converter(field, config=config)

        c.parse("XXX")
    assert len(w) == 1


def test_unicode_mask():
    config = {'verify': 'exception'}
    field = tree.Field(None, name='c', arraysize='1', datatype='unicodeChar',
                       config=config)
    c = converters.get_converter(field, config=config)
    assert c.output("Foo", True) == ''


@raises(exceptions.E02)
def test_wrong_number_of_elements():
    config = {'verify': 'exception'}
    field = tree.Field(
        None, name='c', datatype='int', arraysize='2x3*',
        config=config)
    c = converters.get_converter(field, config=config)
    c.parse("2 3 4 5 6")


@raises(ValueError)
def test_float_mask():
    config = {'verify': 'exception'}
    field = tree.Field(
        None, name='c', datatype='float',
        config=config)
    c = converters.get_converter(field, config=config)
    assert c.parse('') == (c.null, True)
    c.parse('null')


def test_float_mask_permissive():
    config = {'verify': 'ignore'}
    field = tree.Field(
        None, name='c', datatype='float',
        config=config)

    # config needs to be also passed into parse() to work.
    # https://github.com/astropy/astropy/issues/8775
    c = converters.get_converter(field, config=config)
    assert c.parse('null', config=config) == (c.null, True)


@raises(exceptions.E02)
def test_complex_array_vararray():
    config = {'verify': 'exception'}
    field = tree.Field(
        None, name='c', datatype='floatComplex', arraysize='2x3*',
        config=config)
    c = converters.get_converter(field, config=config)
    c.parse("2 3 4 5 6")


def test_complex_array_vararray2():
    config = {'verify': 'exception'}
    field = tree.Field(
        None, name='c', datatype='floatComplex', arraysize='2x3*',
        config=config)
    c = converters.get_converter(field, config=config)
    x = c.parse("")
    assert len(x[0]) == 0


def test_complex_array_vararray3():
    config = {'verify': 'exception'}
    field = tree.Field(
        None, name='c', datatype='doubleComplex', arraysize='2x3*',
        config=config)
    c = converters.get_converter(field, config=config)
    x = c.parse("1 2 3 4 5 6 7 8 9 10 11 12")
    assert len(x) == 2
    assert np.all(x[0][0][0] == complex(1, 2))


def test_complex_vararray():
    config = {'verify': 'exception'}
    field = tree.Field(
        None, name='c', datatype='doubleComplex', arraysize='*',
        config=config)
    c = converters.get_converter(field, config=config)
    x = c.parse("1 2 3 4")
    assert len(x) == 2
    assert x[0][0] == complex(1, 2)


@raises(exceptions.E03)
def test_complex():
    config = {'verify': 'exception'}
    field = tree.Field(
        None, name='c', datatype='doubleComplex',
        config=config)
    c = converters.get_converter(field, config=config)
    c.parse("1 2 3")


@raises(exceptions.E04)
def test_bit():
    config = {'verify': 'exception'}
    field = tree.Field(
        None, name='c', datatype='bit',
        config=config)
    c = converters.get_converter(field, config=config)
    c.parse("T")


def test_bit_mask():
    config = {'verify': 'exception'}
    with catch_warnings(exceptions.W39) as w:
        field = tree.Field(
            None, name='c', datatype='bit',
            config=config)
        c = converters.get_converter(field, config=config)
        c.output(True, True)
    assert len(w) == 1


@raises(exceptions.E05)
def test_boolean():
    config = {'verify': 'exception'}
    field = tree.Field(
        None, name='c', datatype='boolean',
        config=config)
    c = converters.get_converter(field, config=config)
    c.parse('YES')


def test_boolean_array():
    config = {'verify': 'exception'}
    field = tree.Field(
        None, name='c', datatype='boolean', arraysize='*',
        config=config)
    c = converters.get_converter(field, config=config)
    r, mask = c.parse('TRUE FALSE T F 0 1')
    assert_array_equal(r, [True, False, True, False, False, True])


@raises(exceptions.E06)
def test_invalid_type():
    config = {'verify': 'exception'}
    field = tree.Field(
        None, name='c', datatype='foobar',
        config=config)
    converters.get_converter(field, config=config)


def test_precision():
    config = {'verify': 'exception'}

    field = tree.Field(
        None, name='c', datatype='float', precision="E4",
        config=config)
    c = converters.get_converter(field, config=config)
    assert c.output(266.248, False) == '266.2'

    field = tree.Field(
        None, name='c', datatype='float', precision="F4",
        config=config)
    c = converters.get_converter(field, config=config)
    assert c.output(266.248, False) == '266.2480'


@raises(exceptions.W51)
def test_integer_overflow():
    config = {'verify': 'exception'}

    field = tree.Field(
        None, name='c', datatype='int', config=config)
    c = converters.get_converter(field, config=config)
    c.parse('-2208988800', config=config)


def test_float_default_precision():
    config = {'verify': 'exception'}

    field = tree.Field(
        None, name='c', datatype='float', arraysize="4",
        config=config)
    c = converters.get_converter(field, config=config)
    assert (c.output([1, 2, 3, 8.9990234375], [False, False, False, False]) ==
            '1 2 3 8.9990234375')


def test_vararray():
    votable = tree.VOTableFile()
    resource = tree.Resource()
    votable.resources.append(resource)
    table = tree.Table(votable)
    resource.tables.append(table)

    tabarr = []
    heads = ['headA', 'headB', 'headC']
    types = ["char", "double", "int"]

    vals = [["A", 1.0, 2],
            ["B", 2.0, 3],
            ["C", 3.0, 4]]
    for i in range(len(heads)):
        tabarr.append(tree.Field(
            votable, name=heads[i], datatype=types[i], arraysize="*"))

    table.fields.extend(tabarr)
    table.create_arrays(len(vals))
    for i in range(len(vals)):
        values = tuple(vals[i])
        table.array[i] = values
    buff = io.BytesIO()
    votable.to_xml(buff)


def test_gemini_v1_2():
    '''
    see Pull Request 4782 or Issue 4781 for details
    '''
    table = parse_single_table(get_pkg_data_filename('data/gemini.xml'))
    assert table is not None

    tt = table.to_table()
    assert tt['access_url'][0] == (
        b'http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/GEMINI/'
        b'S20120515S0064?runid=bx9b1o8cvk1qesrt')
