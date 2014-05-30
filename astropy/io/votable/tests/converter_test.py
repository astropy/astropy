# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

from __future__ import absolute_import, division, print_function, unicode_literals

import io

# THIRD-PARTY
import numpy as np
from numpy.testing import assert_array_equal

# LOCAL
from .. import converters
from .. import exceptions
from .. import tree
from ....tests.helper import raises, catch_warnings


@raises(exceptions.E13)
def test_invalid_arraysize():
    field = tree.Field(
        None, name='broken', datatype='char', arraysize='foo')
    converters.get_converter(field)


def test_oversize_char():
    config = {'pedantic': True}
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
    config = {'pedantic': True}
    field = tree.Field(
        None, name='c', datatype='char',
        config=config)
    c = converters.get_converter(field, config=config)
    assert c.output("Foo", True) == ''


def test_oversize_unicode():
    config = {'pedantic': True}
    with catch_warnings(exceptions.W46) as w:
        field = tree.Field(
            None, name='c2', datatype='unicodeChar',
            config=config)
        c = converters.get_converter(field, config=config)

        c.parse("XXX")
    assert len(w) == 1


def test_unicode_mask():
    config = {'pedantic': True}
    field = tree.Field(
        None, name='c', datatype='unicodeChar',
        config=config)
    c = converters.get_converter(field, config=config)
    assert c.output("Foo", True) == ''


@raises(exceptions.E02)
def test_wrong_number_of_elements():
    config = {'pedantic': True}
    field = tree.Field(
        None, name='c', datatype='int', arraysize='2x3*',
        config=config)
    c = converters.get_converter(field, config=config)
    c.parse("2 3 4 5 6")


@raises(ValueError)
def test_float_mask():
    config = {'pedantic': True}
    field = tree.Field(
        None, name='c', datatype='float',
        config=config)
    c = converters.get_converter(field, config=config)
    assert c.parse('') == (c.null, True)
    c.parse('null')


def test_float_mask_permissive():
    config = {'pedantic': False}
    field = tree.Field(
        None, name='c', datatype='float',
        config=config)
    c = converters.get_converter(field, config=config)
    assert c.parse('null') == (c.null, True)


@raises(exceptions.E02)
def test_complex_array_vararray():
    config = {'pedantic': True}
    field = tree.Field(
        None, name='c', datatype='floatComplex', arraysize='2x3*',
        config=config)
    c = converters.get_converter(field, config=config)
    c.parse("2 3 4 5 6")


def test_complex_array_vararray2():
    config = {'pedantic': True}
    field = tree.Field(
        None, name='c', datatype='floatComplex', arraysize='2x3*',
        config=config)
    c = converters.get_converter(field, config=config)
    x = c.parse("")
    assert len(x[0]) == 0


def test_complex_array_vararray3():
    config = {'pedantic': True}
    field = tree.Field(
        None, name='c', datatype='doubleComplex', arraysize='2x3*',
        config=config)
    c = converters.get_converter(field, config=config)
    x = c.parse("1 2 3 4 5 6 7 8 9 10 11 12")
    assert len(x) == 2
    assert np.all(x[0][0][0] == complex(1, 2))


def test_complex_vararray():
    config = {'pedantic': True}
    field = tree.Field(
        None, name='c', datatype='doubleComplex', arraysize='*',
        config=config)
    c = converters.get_converter(field, config=config)
    x = c.parse("1 2 3 4")
    assert len(x) == 2
    assert x[0][0] == complex(1, 2)


@raises(exceptions.E03)
def test_complex():
    config = {'pedantic': True}
    field = tree.Field(
        None, name='c', datatype='doubleComplex',
        config=config)
    c = converters.get_converter(field, config=config)
    x = c.parse("1 2 3")


@raises(exceptions.E04)
def test_bit():
    config = {'pedantic': True}
    field = tree.Field(
        None, name='c', datatype='bit',
        config=config)
    c = converters.get_converter(field, config=config)
    x = c.parse("T")


def test_bit_mask(recwarn):
    config = {'pedantic': True}
    with catch_warnings(exceptions.W39) as w:
        field = tree.Field(
            None, name='c', datatype='bit',
            config=config)
        c = converters.get_converter(field, config=config)
        c.output(True, True)
    assert len(w) == 1


@raises(exceptions.E05)
def test_boolean():
    config = {'pedantic': True}
    field = tree.Field(
        None, name='c', datatype='boolean',
        config=config)
    c = converters.get_converter(field, config=config)
    c.parse('YES')


def test_boolean_array():
    config = {'pedantic': True}
    field = tree.Field(
        None, name='c', datatype='boolean', arraysize='*',
        config=config)
    c = converters.get_converter(field, config=config)
    r, mask = c.parse('TRUE FALSE T F 0 1')
    assert_array_equal(r, [True, False, True, False, False, True])


@raises(exceptions.E06)
def test_invalid_type():
    config = {'pedantic': True}
    field = tree.Field(
        None, name='c', datatype='foobar',
        config=config)
    c = converters.get_converter(field, config=config)


def test_precision():
    config = {'pedantic': True}

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
    config = {'pedantic': True}

    field = tree.Field(
        None, name='c', datatype='int', config=config)
    c = converters.get_converter(field, config=config)
    c.parse('-2208988800', config=config)


def test_float_default_precision():
    config = {'pedantic': True}

    field = tree.Field(
        None, name='c', datatype='float', arraysize="4",
        config=config)
    c = converters.get_converter(field, config=config)
    assert (c.output([1, 2, 3, 8.999999], [False, False, False, False]) ==
            '1 2 3 8.9999990000000007')


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
