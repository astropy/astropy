# Licensed under a 3-clause BSD style license - see LICENSE.rst
from copy import copy

import numpy as np

from ...tests.helper import pytest
from ..registry import _readers, _writers, _identifiers, DataIO
from .. import registry as io_registry
from ...table import Table

_READERS_ORIGINAL = copy(_readers)
_WRITERS_ORIGINAL = copy(_writers)
_IDENTIFIERS_ORIGINAL = copy(_identifiers)


def setup_function(function):
    _readers.clear()
    _writers.clear()
    _identifiers.clear()


def empty_reader(*args, **kwargs):
    return DataIO()


def empty_writer(table, *args, **kwargs):
    pass


def empty_identifier(args, kwargs):
    return True


def test_get_reader_invalid():
    with pytest.raises(Exception) as exc:
        io_registry.get_reader('test', DataIO)
    assert exc.value.args[0] == "No reader defined for format 'test' and class 'DataIO'"


def test_get_writer_invalid():
    with pytest.raises(Exception) as exc:
        io_registry.get_writer('test', DataIO)
    assert exc.value.args[0] == "No writer defined for format 'test' and class 'DataIO'"


def test_register_reader():
    io_registry.register_reader('test1', DataIO, empty_reader)
    io_registry.register_reader('test2', DataIO, empty_reader)
    assert io_registry.get_reader('test1', DataIO) == empty_reader
    assert io_registry.get_reader('test2', DataIO) == empty_reader


def test_register_writer():
    io_registry.register_writer('test1', DataIO, empty_writer)
    io_registry.register_writer('test2', DataIO, empty_writer)
    assert io_registry.get_writer('test1', DataIO) == empty_writer
    assert io_registry.get_writer('test2', DataIO) == empty_writer


def test_register_identifier():
    io_registry.register_identifier('test1', DataIO, empty_identifier)
    io_registry.register_identifier('test2', DataIO, empty_identifier)


def test_register_reader_invalid():
    io_registry.register_reader('test', DataIO, empty_reader)
    with pytest.raises(Exception) as exc:
        io_registry.register_reader('test', DataIO, empty_reader)
    assert exc.value.args[0] == "Reader for format 'test' and class 'DataIO' is already defined"


def test_register_writer_invalid():
    io_registry.register_writer('test', DataIO, empty_writer)
    with pytest.raises(Exception) as exc:
        io_registry.register_writer('test', DataIO, empty_writer)
    assert exc.value.args[0] == "Writer for format 'test' and class 'DataIO' is already defined"


def test_register_identifier_invalid():
    io_registry.register_identifier('test', DataIO, empty_identifier)
    with pytest.raises(Exception) as exc:
        io_registry.register_identifier('test', DataIO, empty_identifier)
    assert exc.value.args[0] == "Identifier for format 'test' and class 'DataIO' is already defined"


def test_register_reader_force():
    io_registry.register_reader('test', DataIO, empty_reader)
    io_registry.register_reader('test', DataIO, empty_reader, force=True)


def test_register_writer_force():
    io_registry.register_writer('test', DataIO, empty_writer)
    io_registry.register_writer('test', DataIO, empty_writer, force=True)


def test_register_identifier_force():
    io_registry.register_identifier('test', DataIO, empty_identifier)
    io_registry.register_identifier('test', DataIO, empty_identifier, force=True)


def test_read_noformat():
    with pytest.raises(Exception) as exc:
        DataIO.read()
    assert exc.value.args[0] == "Format could not be identified"


def test_write_noformat():
    with pytest.raises(Exception) as exc:
        DataIO().write()
    assert exc.value.args[0] == "Format could not be identified"


def test_read_noformat_arbitrary():
    """Test that all identifier functions can accept arbitary input"""
    _identifiers.update(_IDENTIFIERS_ORIGINAL)
    with pytest.raises(Exception) as exc:
        DataIO.read(object())
    assert exc.value.args[0] == "Format could not be identified"


def test_write_noformat_arbitrary():
    """Test that all identifier functions can accept arbitary input"""
    _identifiers.update(_IDENTIFIERS_ORIGINAL)
    with pytest.raises(Exception) as exc:
        DataIO().write(object())
    assert exc.value.args[0] == "Format could not be identified"


def test_read_toomanyformats():
    io_registry.register_identifier('test1', DataIO, lambda o, x, y: True)
    io_registry.register_identifier('test2', DataIO, lambda o, x, y: True)
    with pytest.raises(Exception) as exc:
        DataIO.read()
    assert exc.value.args[0] == "Format is ambiguous - options are: test1, test2"


def test_write_toomanyformats():
    io_registry.register_identifier('test1', DataIO, lambda o, x, y: True)
    io_registry.register_identifier('test2', DataIO, lambda o, x, y: True)
    with pytest.raises(Exception) as exc:
        DataIO().write()
    assert exc.value.args[0] == "Format is ambiguous - options are: test1, test2"


def test_read_format_noreader():
    with pytest.raises(Exception) as exc:
        DataIO.read(format='test')
    assert exc.value.args[0] == "No reader defined for format 'test' and class 'DataIO'"


def test_write_format_nowriter():
    with pytest.raises(Exception) as exc:
        DataIO().write(format='test')
    assert exc.value.args[0] == "No writer defined for format 'test' and class 'DataIO'"


def test_read_identifier():

    io_registry.register_identifier('test1', DataIO, lambda o, x, y: x[0].startswith('a'))
    io_registry.register_identifier('test2', DataIO, lambda o, x, y: x[0].startswith('b'))

    # Now check that we got past the identifier and are trying to get
    # the reader. The io_registry.get_reader will fail but the error message will
    # tell us if the identifier worked.

    with pytest.raises(Exception) as exc:
        DataIO.read('abc')
    assert exc.value.args[0] == "No reader defined for format 'test1' and class 'DataIO'"

    with pytest.raises(Exception) as exc:
        DataIO.read('bac')
    assert exc.value.args[0] == "No reader defined for format 'test2' and class 'DataIO'"


def test_write_identifier():

    io_registry.register_identifier('test1', DataIO, lambda o, x, y: x[0].startswith('a'))
    io_registry.register_identifier('test2', DataIO, lambda o, x, y: x[0].startswith('b'))

    # Now check that we got past the identifier and are trying to get
    # the reader. The io_registry.get_writer will fail but the error message will
    # tell us if the identifier worked.

    with pytest.raises(Exception) as exc:
        DataIO().write('abc')
    assert exc.value.args[0] == "No writer defined for format 'test1' and class 'DataIO'"

    with pytest.raises(Exception) as exc:
        DataIO().write('bac')
    assert exc.value.args[0] == "No writer defined for format 'test2' and class 'DataIO'"


def test_identifier_origin():

    io_registry.register_identifier('test1', DataIO, lambda o, x, y: o == 'read')
    io_registry.register_identifier('test2', DataIO, lambda o, x, y: o == 'write')
    io_registry.register_reader('test1', DataIO, empty_reader)
    io_registry.register_writer('test2', DataIO, empty_writer)

    # There should not be too many formats defined
    DataIO.read()
    DataIO().write()

    with pytest.raises(Exception) as exc:
        DataIO.read(format='test2')
    assert exc.value.args[0] == "No reader defined for format 'test2' and class 'DataIO'"

    with pytest.raises(Exception) as exc:
        DataIO().write(format='test1')
    assert exc.value.args[0] == "No writer defined for format 'test1' and class 'DataIO'"


def test_read_valid_return():
    io_registry.register_reader('test', DataIO, lambda: DataIO())
    t = DataIO.read(format='test')
    assert isinstance(t, DataIO)


def test_read_invalid_return():
    io_registry.register_reader('test', DataIO, lambda: 'spam')
    with pytest.raises(TypeError) as exc:
        DataIO.read(format='test')
    assert exc.value.args[0] == "reader should return a DataIO instance"


def test_read_basic_table():
    data = np.array(zip([1, 2, 3], ['a', 'b', 'c']), dtype=[('A', int), ('B', '|S1')])
    io_registry.register_reader('test', Table, lambda x: Table(x))
    t = Table.read(data, format='test')
    assert t.keys() == ['A', 'B']
    for i in range(3):
        assert t['A'][i] == data['A'][i]
        assert t['B'][i] == data['B'][i]


def teardown_function(function):
    _readers.update(_READERS_ORIGINAL)
    _writers.update(_WRITERS_ORIGINAL)
    _identifiers.update(_IDENTIFIERS_ORIGINAL)
