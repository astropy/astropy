from copy import copy

import numpy as np

from ...tests.helper import pytest
from ...table.io_registry import _readers, _writers, _identifiers
from ...table import Table
from ...table import io_registry

_READERS_ORIGINAL = copy(_readers)
_WRITERS_ORIGINAL = copy(_writers)
_IDENTIFIERS_ORIGINAL = copy(_identifiers)


def setup_function(function):
    _readers.clear()
    _writers.clear()
    _identifiers.clear()


def empty_reader(*args, **kwargs):
    return Table()


def empty_writer(table, *args, **kwargs):
    pass


def empty_identifier(args, kwargs):
    return True


def test_get_reader_invalid():
    with pytest.raises(Exception) as exc:
        io_registry.get_reader('test')
    assert exc.value.args[0] == "No reader defined for format 'test'"


def test_get_writer_invalid():
    with pytest.raises(Exception) as exc:
        io_registry.get_writer('test')
    assert exc.value.args[0] == "No writer defined for format 'test'"


def test_register_reader():
    io_registry.register_reader('test1', empty_reader)
    io_registry.register_reader('test2', empty_reader)
    assert io_registry.get_reader('test1') == empty_reader
    assert io_registry.get_reader('test2') == empty_reader


def test_register_writer():
    io_registry.register_writer('test1', empty_writer)
    io_registry.register_writer('test2', empty_writer)
    assert io_registry.get_writer('test1') == empty_writer
    assert io_registry.get_writer('test2') == empty_writer


def test_register_identifier():
    io_registry.register_identifier('test1', empty_identifier)
    io_registry.register_identifier('test2', empty_identifier)


def test_register_reader_invalid():
    io_registry.register_reader('test', empty_reader)
    with pytest.raises(Exception) as exc:
        io_registry.register_reader('test', empty_reader)
    assert exc.value.args[0] == "Reader for format test is already defined"


def test_register_writer_invalid():
    io_registry.register_writer('test', empty_writer)
    with pytest.raises(Exception) as exc:
        io_registry.register_writer('test', empty_writer)
    assert exc.value.args[0] == "Writer for format test is already defined"


def test_register_identifier_invalid():
    io_registry.register_identifier('test', empty_identifier)
    with pytest.raises(Exception) as exc:
        io_registry.register_identifier('test', empty_identifier)
    assert exc.value.args[0] == "Identifier for format test is already defined"


def test_register_reader_force():
    io_registry.register_reader('test', empty_reader)
    io_registry.register_reader('test', empty_reader, force=True)


def test_register_writer_force():
    io_registry.register_writer('test', empty_writer)
    io_registry.register_writer('test', empty_writer, force=True)


def test_register_identifier_force():
    io_registry.register_identifier('test', empty_identifier)
    io_registry.register_identifier('test', empty_identifier, force=True)


def test_read_noformat():
    with pytest.raises(Exception) as exc:
        Table.read()
    assert exc.value.args[0] == "Format could not be identified"


def test_write_noformat():
    with pytest.raises(Exception) as exc:
        Table().write()
    assert exc.value.args[0] == "Format could not be identified"


def test_read_noformat_arbitrary():
    """Test that all identifier functions can accept arbitary input"""
    _identifiers.update(_IDENTIFIERS_ORIGINAL)
    with pytest.raises(Exception) as exc:
        Table.read(object())
    assert exc.value.args[0] == "Format could not be identified"


def test_write_noformat_arbitrary():
    """Test that all identifier functions can accept arbitary input"""
    _identifiers.update(_IDENTIFIERS_ORIGINAL)
    with pytest.raises(Exception) as exc:
        Table().write(object())
    assert exc.value.args[0] == "Format could not be identified"


def test_read_toomanyformats():
    io_registry.register_identifier('test1', lambda o, x, y: True)
    io_registry.register_identifier('test2', lambda o, x, y: True)
    with pytest.raises(Exception) as exc:
        Table.read()
    assert exc.value.args[0] == "Format is ambiguous - options are: test1, test2"


def test_write_toomanyformats():
    io_registry.register_identifier('test1', lambda o, x, y: True)
    io_registry.register_identifier('test2', lambda o, x, y: True)
    with pytest.raises(Exception) as exc:
        Table().write()
    assert exc.value.args[0] == "Format is ambiguous - options are: test1, test2"


def test_read_format_noreader():
    with pytest.raises(Exception) as exc:
        Table.read(format='test')
    assert exc.value.args[0] == "No reader defined for format 'test'"


def test_write_format_nowriter():
    with pytest.raises(Exception) as exc:
        Table().write(format='test')
    assert exc.value.args[0] == "No writer defined for format 'test'"


def test_read_identifier():

    io_registry.register_identifier('test1', lambda o, x, y: x[0].startswith('a'))
    io_registry.register_identifier('test2', lambda o, x, y: x[0].startswith('b'))

    # Now check that we got past the identifier and are trying to get
    # the reader. The io_registry.get_reader will fail but the error message will
    # tell us if the identifier worked.

    with pytest.raises(Exception) as exc:
        Table.read('abc')
    assert exc.value.args[0] == "No reader defined for format 'test1'"

    with pytest.raises(Exception) as exc:
        Table.read('bac')
    assert exc.value.args[0] == "No reader defined for format 'test2'"


def test_write_identifier():

    io_registry.register_identifier('test1', lambda o, x, y: x[0].startswith('a'))
    io_registry.register_identifier('test2', lambda o, x, y: x[0].startswith('b'))

    # Now check that we got past the identifier and are trying to get
    # the reader. The io_registry.get_writer will fail but the error message will
    # tell us if the identifier worked.

    with pytest.raises(Exception) as exc:
        Table().write('abc')
    assert exc.value.args[0] == "No writer defined for format 'test1'"

    with pytest.raises(Exception) as exc:
        Table().write('bac')
    assert exc.value.args[0] == "No writer defined for format 'test2'"


def test_identifier_origin():

    io_registry.register_identifier('test1', lambda o, x, y: o == 'read')
    io_registry.register_identifier('test2', lambda o, x, y: o == 'write')
    io_registry.register_reader('test1', empty_reader)
    io_registry.register_writer('test2', empty_writer)

    # There should not be too many formats defined
    Table.read()
    Table().write()

    with pytest.raises(Exception) as exc:
        Table.read(format='test2')
    assert exc.value.args[0] == "No reader defined for format 'test2'"

    with pytest.raises(Exception) as exc:
        Table().write(format='test1')
    assert exc.value.args[0] == "No writer defined for format 'test1'"


def test_read_valid_return():
    io_registry.register_reader('test', lambda: Table())
    t = Table.read(format='test')
    assert isinstance(t, Table)


def test_read_invalid_return():
    io_registry.register_reader('test', lambda: 'spam')
    with pytest.raises(TypeError) as exc:
        Table.read(format='test')
    assert exc.value.args[0] == "reader should return a Table instance"


def test_read_basic():
    data = np.array(zip([1, 2, 3], ['a', 'b', 'c']), dtype=[('A', int), ('B', '|S1')])
    io_registry.register_reader('test', lambda x: Table(x))
    t = Table.read(data, format='test')
    assert t.keys() == ['A', 'B']
    for i in range(3):
        assert t['A'][i] == data['A'][i]
        assert t['B'][i] == data['B'][i]


def teardown_function(function):
    _readers.update(_READERS_ORIGINAL)
    _writers.update(_WRITERS_ORIGINAL)
    _identifiers.update(_IDENTIFIERS_ORIGINAL)
