from copy import copy

import pytest
import numpy as np

from ..io_registry import _readers, _writers, _identifiers
from .. import Table, register_reader, register_writer, \
               register_identifier, get_reader, get_writer
from ... import table

_READERS_ORIGINAL = copy(_readers)
_WRITERS_ORIGINAL = copy(_writers)
_IDENTIFIERS_ORIGINAL = copy(_identifiers)


def setup_function(function):
    for key in _readers.keys():
        _readers.pop(key)
    for key in _writers.keys():
        _writers.pop(key)
    for key in _identifiers.keys():
        _identifiers.pop(key)


def empty_reader(*args, **kwargs):
    return Table()


def empty_writer(table, *args, **kwargs):
    return Table()


def empty_identifier(args, kwargs):
    return True


def test_get_reader_invalid():
    with pytest.raises(Exception) as exc:
        get_reader('test')
    assert exc.value.args[0] == "No reader defined for format 'test'"


def test_get_writer_invalid():
    with pytest.raises(Exception) as exc:
        get_writer('test')
    assert exc.value.args[0] == "No writer defined for format 'test'"


def test_register_reader():
    register_reader('test1', empty_reader)
    register_reader('test2', empty_reader)
    assert get_reader('test1') == empty_reader
    assert get_reader('test2') == empty_reader


def test_register_writer():
    register_writer('test1', empty_writer)
    register_writer('test2', empty_writer)
    assert get_writer('test1') == empty_writer
    assert get_writer('test2') == empty_writer


def test_register_identifier():
    register_identifier('test1', empty_identifier)
    register_identifier('test2', empty_identifier)


def test_register_reader_invalid():
    register_reader('test', empty_reader)
    with pytest.raises(Exception) as exc:
        register_reader('test', empty_reader)
    assert exc.value.args[0] == "Reader for format test is already defined"


def test_register_writer_invalid():
    register_writer('test', empty_writer)
    with pytest.raises(Exception) as exc:
        register_writer('test', empty_writer)
    assert exc.value.args[0] == "Writer for format test is already defined"


def test_register_identifier_invalid():
    register_identifier('test', empty_identifier)
    with pytest.raises(Exception) as exc:
        register_identifier('test', empty_identifier)
    assert exc.value.args[0] == "Identifier for format test is already defined"


def test_register_reader_force():
    register_reader('test', empty_reader)
    register_reader('test', empty_reader, force=True)


def test_register_writer_force():
    register_writer('test', empty_writer)
    register_writer('test', empty_writer, force=True)


def test_register_identifier_force():
    register_identifier('test', empty_identifier)
    register_identifier('test', empty_identifier, force=True)


def test_read_noformat():
    with pytest.raises(Exception) as exc:
        table.read()
    assert exc.value.args[0] == "Format could not be identified"


def test_write_noformat():
    with pytest.raises(Exception) as exc:
        Table().write()
    assert exc.value.args[0] == "Format could not be identified"


def test_read_toomanyformats():
    register_identifier('test1', lambda x, y: True)
    register_identifier('test2', lambda x, y: True)
    with pytest.raises(Exception) as exc:
        table.read()
    assert exc.value.args[0] == "Format is ambiguous - options are: test1, test2"


def test_write_toomanyformats():
    register_identifier('test1', lambda x, y: True)
    register_identifier('test2', lambda x, y: True)
    with pytest.raises(Exception) as exc:
        Table().write()
    assert exc.value.args[0] == "Format is ambiguous - options are: test1, test2"


def test_read_format_noreader():
    with pytest.raises(Exception) as exc:
        table.read(format='test')
    assert exc.value.args[0] == "No reader defined for format 'test'"


def test_write_format_nowriter():
    with pytest.raises(Exception) as exc:
        Table().write(format='test')
    assert exc.value.args[0] == "No writer defined for format 'test'"


def test_read_identifier():

    register_identifier('test1', lambda x, y: x[0].startswith('a'))
    register_identifier('test2', lambda x, y: x[0].startswith('b'))

    # Now check that we got past the identifier and are trying to get
    # the reader. The get_reader will fail but the error message will
    # tell us if the identifier worked.

    with pytest.raises(Exception) as exc:
        table.read('abc')
    assert exc.value.args[0] == "No reader defined for format 'test1'"

    with pytest.raises(Exception) as exc:
        table.read('bac')
    assert exc.value.args[0] == "No reader defined for format 'test2'"


def test_write_identifier():

    register_identifier('test1', lambda x, y: x[0].startswith('a'))
    register_identifier('test2', lambda x, y: x[0].startswith('b'))

    # Now check that we got past the identifier and are trying to get
    # the reader. The get_writer will fail but the error message will
    # tell us if the identifier worked.

    with pytest.raises(Exception) as exc:
        Table().write('abc')
    assert exc.value.args[0] == "No writer defined for format 'test1'"

    with pytest.raises(Exception) as exc:
        Table().write('bac')
    assert exc.value.args[0] == "No writer defined for format 'test2'"


def test_read_valid_return():
    register_reader('test', lambda: Table())
    t = table.read(format='test')
    assert isinstance(t, Table)


def test_read_invalid_return():
    register_reader('test', lambda: 'spam')
    with pytest.raises(TypeError) as exc:
        table.read(format='test')
    assert exc.value.args[0] == "reader should return a Table instance"


def test_read_basic():
    data = np.array(zip([1, 2, 3], ['a', 'b', 'c']), dtype=[('A', int), ('B', '|S1')])
    register_reader('test', lambda x: Table(x))
    t = table.read(data, format='test')
    assert t.keys() == ['A', 'B']
    for i in range(3):
        assert t['A'][i] == data['A'][i]
        assert t['B'][i] == data['B'][i]


def teardown_function(function):
    for key in _READERS_ORIGINAL:
        _readers[key] = _READERS_ORIGINAL[key]
    for key in _WRITERS_ORIGINAL:
        _writers[key] = _WRITERS_ORIGINAL[key]
    for key in _IDENTIFIERS_ORIGINAL:
        _identifiers[key] = _IDENTIFIERS_ORIGINAL[key]
