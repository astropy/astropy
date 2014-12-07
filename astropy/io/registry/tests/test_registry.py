# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os

from copy import copy

import numpy as np

from ....tests.helper import pytest
from ....table import Table
from ....extern.six.moves import zip
from ....extern.six import StringIO

from ..registry import _registry
from .. import registry as io_registry


class TestData(object):
    read = classmethod(io_registry.read)
    write = io_registry.write


def setup_function(function):
    # Reset the registry to a clean state
    io_registry._registry.clear()
    io_registry._load_builtins()


def teardown_function(function):
    # Reset the registry to a clean state
    io_registry._registry.clear()
    io_registry._load_builtins()


def empty_reader(*args, **kwargs):
    return TestData()


def empty_writer(table, *args, **kwargs):
    pass


def empty_identifier(*args, **kwargs):
    return True


def test_get_reader_invalid():
    with pytest.raises(Exception) as exc:
        io_registry.get_reader('test', TestData)
    assert exc.value.args[0].startswith(
        "No reader defined for format 'test' and class 'TestData'")


def test_get_writer_invalid():
    with pytest.raises(Exception) as exc:
        io_registry.get_writer('test', TestData)
    assert exc.value.args[0].startswith(
        "No writer defined for format 'test' and class 'TestData'")


def test_register_reader():
    io_registry.register_reader('test1', TestData, empty_reader)
    io_registry.register_reader('test2', TestData, empty_reader)
    assert io_registry.get_reader('test1', TestData) == empty_reader
    assert io_registry.get_reader('test2', TestData) == empty_reader


def test_register_writer():
    io_registry.register_writer('test1', TestData, empty_writer)
    io_registry.register_writer('test2', TestData, empty_writer)
    assert io_registry.get_writer('test1', TestData) == empty_writer
    assert io_registry.get_writer('test2', TestData) == empty_writer


def test_register_identifier():
    io_registry.register_identifier('test1', TestData, empty_identifier)
    io_registry.register_identifier('test2', TestData, empty_identifier)


def test_register_reader_invalid():
    io_registry.register_reader('test', TestData, empty_reader)
    with pytest.raises(Exception) as exc:
        io_registry.register_reader('test', TestData, empty_reader)
    assert exc.value.args[0] == "Reader for format 'test' and class 'TestData' is already defined"


def test_register_writer_invalid():
    io_registry.register_writer('test', TestData, empty_writer)
    with pytest.raises(Exception) as exc:
        io_registry.register_writer('test', TestData, empty_writer)
    assert exc.value.args[0] == "Writer for format 'test' and class 'TestData' is already defined"


def test_register_identifier_invalid():
    io_registry.register_identifier('test', TestData, empty_identifier)
    with pytest.raises(Exception) as exc:
        io_registry.register_identifier('test', TestData, empty_identifier)
    assert exc.value.args[0] == "Identifier for format 'test' and class 'TestData' is already defined"


def test_register_reader_force():
    io_registry.register_reader('test', TestData, empty_reader)
    io_registry.register_reader('test', TestData, empty_reader, force=True)


def test_register_writer_force():
    io_registry.register_writer('test', TestData, empty_writer)
    io_registry.register_writer('test', TestData, empty_writer, force=True)


def test_register_identifier_force():
    io_registry.register_identifier('test', TestData, empty_identifier)
    io_registry.register_identifier('test', TestData, empty_identifier, force=True)


def test_read_noformat():
    with pytest.raises(Exception) as exc:
        TestData.read()
    assert exc.value.args[0].startswith("Format could not be identified.")


def test_write_noformat():
    with pytest.raises(Exception) as exc:
        TestData().write()
    assert exc.value.args[0].startswith("Format could not be identified.")


def test_read_noformat_arbitrary():
    """Test that all built-in identifier functions can accept arbitary input"""
    io_registry._builtin_registered = False
    with pytest.raises(Exception) as exc:
        TestData.read(object())
    assert exc.value.args[0].startswith("Format could not be identified.")


def test_read_noformat_arbitrary_file(tmpdir):
    """Tests that all built-in reader functions can accept arbitrary files"""

    testfile = str(tmpdir.join('foo.example'))
    with open(testfile, 'w') as f:
        f.write("Hello world")

    with pytest.raises(Exception) as exc:
        Table.read(testfile)
    assert exc.value.args[0].startswith("Format could not be identified.")


def test_write_noformat_arbitrary():
    """Test that all built-in writer functions can accept arbitary input"""
    with pytest.raises(Exception) as exc:
        TestData().write(object())
    assert exc.value.args[0].startswith("Format could not be identified.")


def test_write_noformat_arbitrary_file(tmpdir):
    """Tests that all built-in writer functions can accept arbitrary files"""
    testfile = str(tmpdir.join('foo.example'))

    with pytest.raises(Exception) as exc:
        Table().write(testfile)
    assert exc.value.args[0].startswith("Format could not be identified.")


def test_read_toomanyformats():
    io_registry.register_identifier('test1', TestData, lambda o, *x, **y: True)
    io_registry.register_identifier('test2', TestData, lambda o, *x, **y: True)
    with pytest.raises(Exception) as exc:
        TestData.read()
    assert exc.value.args[0] == "Format is ambiguous - options are: test1, test2"


def test_write_toomanyformats():
    io_registry.register_identifier('test1', TestData, lambda o, *x, **y: True)
    io_registry.register_identifier('test2', TestData, lambda o, *x, **y: True)
    with pytest.raises(Exception) as exc:
        TestData().write()
    assert exc.value.args[0] == "Format is ambiguous - options are: test1, test2"


def test_read_format_noreader():
    with pytest.raises(Exception) as exc:
        TestData.read(format='test')
    assert exc.value.args[0].startswith(
        "No reader defined for format 'test' and class 'TestData'")


def test_write_format_nowriter():
    with pytest.raises(Exception) as exc:
        TestData().write(format='test')
    assert exc.value.args[0].startswith(
        "No writer defined for format 'test' and class 'TestData'")


def test_read_identifier():

    io_registry.register_identifier(
        'test1', TestData,
        lambda o, path, fileobj, *x, **y: path.startswith('a'))
    io_registry.register_identifier(
        'test2', TestData,
        lambda o, path, fileobj, *x, **y: path.startswith('b'))

    # Now check that we got past the identifier and are trying to get
    # the reader. The io_registry.get_reader will fail but the error message will
    # tell us if the identifier worked.

    with pytest.raises(Exception) as exc:
        TestData.read('abc')
    assert exc.value.args[0].startswith(
        "No reader defined for format 'test1' and class 'TestData'")

    with pytest.raises(Exception) as exc:
        TestData.read('bac')
    assert exc.value.args[0].startswith(
        "No reader defined for format 'test2' and class 'TestData'")


def test_write_identifier():

    io_registry.register_identifier('test1', TestData, lambda o, *x, **y: x[0].startswith('a'))
    io_registry.register_identifier('test2', TestData, lambda o, *x, **y: x[0].startswith('b'))

    # Now check that we got past the identifier and are trying to get
    # the reader. The io_registry.get_writer will fail but the error message will
    # tell us if the identifier worked.

    with pytest.raises(Exception) as exc:
        TestData().write('abc')
    assert exc.value.args[0].startswith(
        "No writer defined for format 'test1' and class 'TestData'")

    with pytest.raises(Exception) as exc:
        TestData().write('bac')
    assert exc.value.args[0].startswith(
        "No writer defined for format 'test2' and class 'TestData'")


def test_identifier_origin():

    io_registry.register_identifier('test1', TestData, lambda o, *x, **y: o == 'read')
    io_registry.register_identifier('test2', TestData, lambda o, *x, **y: o == 'write')
    io_registry.register_reader('test1', TestData, empty_reader)
    io_registry.register_writer('test2', TestData, empty_writer)

    # There should not be too many formats defined
    TestData.read()
    TestData().write()

    with pytest.raises(Exception) as exc:
        TestData.read(format='test2')
    assert exc.value.args[0].startswith(
        "No reader defined for format 'test2' and class 'TestData'")

    with pytest.raises(Exception) as exc:
        TestData().write(format='test1')
    assert exc.value.args[0].startswith(
        "No writer defined for format 'test1' and class 'TestData'")


def test_read_valid_return():
    io_registry.register_reader('test', TestData, lambda: TestData())
    t = TestData.read(format='test')
    assert isinstance(t, TestData)


def test_read_invalid_return():
    io_registry.register_reader('test', TestData, lambda: 'spam')
    with pytest.raises(TypeError) as exc:
        TestData.read(format='test')
    assert exc.value.args[0] == "reader should return a TestData instance"


def test_read_basic_table():
    data = np.array(list(zip([1, 2, 3], ['a', 'b', 'c'])),
                    dtype=[(str('A'), int), (str('B'), '|S1')])
    io_registry.register_reader('test', Table, lambda x: Table(x))
    t = Table.read(data, format='test')
    assert t.keys() == ['A', 'B']
    for i in range(3):
        assert t['A'][i] == data['A'][i]
        assert t['B'][i] == data['B'][i]


def test_register_readers_with_same_name_on_different_classes():
    # No errors should be generated if the same name is registered for
    # different objects...but this failed under python3
    io_registry.register_reader('test', TestData, lambda: TestData())
    io_registry.register_reader('test', Table, lambda: Table())
    t = TestData.read(format='test')
    assert isinstance(t, TestData)
    tbl = Table.read(format='test')
    assert isinstance(tbl, Table)


class TestSubclass:
    """
    Test using registry with a Table sub-class
    """
    def test_read_table_subclass(self):
        class MyTable(Table):
            pass
        data = ['a b', '1 2']
        mt = MyTable.read(data, format='ascii')
        t = Table.read(data, format='ascii')
        assert np.all(mt == t)
        assert mt.colnames == t.colnames
        assert type(mt) is MyTable

    def test_write_table_subclass(self):
        buffer = StringIO()
        class MyTable(Table):
            pass
        mt = MyTable([[1], [2]], names=['a', 'b'])
        mt.write(buffer, format='ascii')
        assert buffer.getvalue() == os.linesep.join(['a b', '1 2', ''])
