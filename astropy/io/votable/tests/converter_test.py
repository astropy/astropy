# Licensed under a 3-clause BSD style license - see LICENSE.rst
# THIRD-PARTY
import numpy as np
from numpy.testing import assert_array_equal

# LOCAL
from .. import converters
from .. import exceptions
from .. import tree
from ....tests.helper import raises


@raises(exceptions.E13)
def test_invalid_arraysize():
    field = tree.Field(
        None, name='broken', datatype='char', arraysize='foo')
    converters.get_converter(field)


def test_oversize_char(recwarn):
    config = {'pedantic': True}
    field = tree.Field(
        None, name='c', datatype='char',
        config=config)
    c = converters.get_converter(field, config=config)
    w = recwarn.pop(exceptions.W47)

    c.parse(u"XXX")
    w = recwarn.pop(exceptions.W46)


def test_char_mask():
    config = {'pedantic': True}
    field = tree.Field(
        None, name='c', datatype='char',
        config=config)
    c = converters.get_converter(field, config=config)
    assert c.output("Foo", True) == ''


def test_oversize_unicode(recwarn):
    config = {'pedantic': True}
    field = tree.Field(
        None, name='c2', datatype='unicodeChar',
        config=config)
    c = converters.get_converter(field, config=config)

    c.parse(u"XXX")
    w = recwarn.pop(exceptions.W46)


def test_unicode_mask():
    config = {'pedantic': True}
    field = tree.Field(
        None, name='c', datatype='unicodeChar',
        config=config)
    c = converters.get_converter(field, config=config)
    assert c.output(u"Foo", True) == u''


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


def test_complex_array_vararray():
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
    field = tree.Field(
        None, name='c', datatype='bit',
        config=config)
    c = converters.get_converter(field, config=config)
    c.output(True, True)
    recwarn.pop(exceptions.W39)


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
