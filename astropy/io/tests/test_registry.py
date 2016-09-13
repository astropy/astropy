# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
from copy import copy

import numpy as np

from ...tests.helper import pytest
from ..registry import _readers, _writers, _identifiers
from .. import registry as io_registry
from ...table import Table
from ...extern.six.moves import zip, range
from ...extern.six import StringIO

_READERS_ORIGINAL = copy(_readers)
_WRITERS_ORIGINAL = copy(_writers)
_IDENTIFIERS_ORIGINAL = copy(_identifiers)


class TestData(object):
    read = classmethod(io_registry.read)
    write = io_registry.write


def setup_function(function):
    _readers.clear()
    _writers.clear()
    _identifiers.clear()


def empty_reader(*args, **kwargs):
    return TestData()


def empty_writer(table, *args, **kwargs):
    pass


def empty_identifier(*args, **kwargs):
    return True


def test_get_reader_invalid():
    with pytest.raises(io_registry.IORegistryError) as exc:
        io_registry.get_reader('test', TestData)
    assert str(exc.value).startswith(
        "No reader defined for format 'test' and class 'TestData'")


def test_get_writer_invalid():
    with pytest.raises(io_registry.IORegistryError) as exc:
        io_registry.get_writer('test', TestData)
    assert str(exc.value).startswith(
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
    with pytest.raises(io_registry.IORegistryError) as exc:
        io_registry.register_reader('test', TestData, empty_reader)
    assert (str(exc.value) == "Reader for format 'test' and class 'TestData' "
                              "is already defined")


def test_register_writer_invalid():
    io_registry.register_writer('test', TestData, empty_writer)
    with pytest.raises(io_registry.IORegistryError) as exc:
        io_registry.register_writer('test', TestData, empty_writer)
    assert (str(exc.value) == "Writer for format 'test' and class 'TestData' "
                              "is already defined")


def test_register_identifier_invalid():
    io_registry.register_identifier('test', TestData, empty_identifier)
    with pytest.raises(io_registry.IORegistryError) as exc:
        io_registry.register_identifier('test', TestData, empty_identifier)
    assert (str(exc.value) == "Identifier for format 'test' and class "
                              "'TestData' is already defined")


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
    with pytest.raises(io_registry.IORegistryError) as exc:
        TestData.read()
    assert str(exc.value).startswith("Format could not be identified.")


def test_write_noformat():
    with pytest.raises(io_registry.IORegistryError) as exc:
        TestData().write()
    assert str(exc.value).startswith("Format could not be identified.")


def test_read_noformat_arbitrary():
    """Test that all identifier functions can accept arbitrary input"""
    _identifiers.update(_IDENTIFIERS_ORIGINAL)
    with pytest.raises(io_registry.IORegistryError) as exc:
        TestData.read(object())
    assert str(exc.value).startswith("Format could not be identified.")


def test_read_noformat_arbitrary_file(tmpdir):
    """Tests that all identifier functions can accept arbitrary files"""
    _readers.update(_READERS_ORIGINAL)
    testfile = str(tmpdir.join('foo.example'))
    with open(testfile, 'w') as f:
        f.write("Hello world")

    with pytest.raises(io_registry.IORegistryError) as exc:
        Table.read(testfile)
    assert str(exc.value).startswith("Format could not be identified.")


def test_write_noformat_arbitrary():
    """Test that all identifier functions can accept arbitrary input"""
    _identifiers.update(_IDENTIFIERS_ORIGINAL)
    with pytest.raises(io_registry.IORegistryError) as exc:
        TestData().write(object())
    assert str(exc.value).startswith("Format could not be identified.")


def test_write_noformat_arbitrary_file(tmpdir):
    """Tests that all identifier functions can accept arbitrary files"""
    _writers.update(_WRITERS_ORIGINAL)
    testfile = str(tmpdir.join('foo.example'))

    with pytest.raises(io_registry.IORegistryError) as exc:
        Table().write(testfile)
    assert str(exc.value).startswith("Format could not be identified.")


def test_read_toomanyformats():
    io_registry.register_identifier('test1', TestData, lambda o, *x, **y: True)
    io_registry.register_identifier('test2', TestData, lambda o, *x, **y: True)
    with pytest.raises(io_registry.IORegistryError) as exc:
        TestData.read()
    assert str(exc.value) == "Format is ambiguous - options are: test1, test2"


def test_write_toomanyformats():
    io_registry.register_identifier('test1', TestData, lambda o, *x, **y: True)
    io_registry.register_identifier('test2', TestData, lambda o, *x, **y: True)
    with pytest.raises(io_registry.IORegistryError) as exc:
        TestData().write()
    assert str(exc.value) == "Format is ambiguous - options are: test1, test2"


def test_read_format_noreader():
    with pytest.raises(io_registry.IORegistryError) as exc:
        TestData.read(format='test')
    assert str(exc.value).startswith(
        "No reader defined for format 'test' and class 'TestData'")


def test_write_format_nowriter():
    with pytest.raises(io_registry.IORegistryError) as exc:
        TestData().write(format='test')
    assert str(exc.value).startswith(
        "No writer defined for format 'test' and class 'TestData'")


def test_read_identifier(tmpdir):

    io_registry.register_identifier(
        'test1', TestData,
        lambda o, path, fileobj, *x, **y: path.endswith('a'))
    io_registry.register_identifier(
        'test2', TestData,
        lambda o, path, fileobj, *x, **y: path.endswith('b'))

    # Now check that we got past the identifier and are trying to get
    # the reader. The io_registry.get_reader will fail but the error message
    # will tell us if the identifier worked.

    filename = tmpdir.join("testfile.a").strpath
    open(filename, 'w').close()
    with pytest.raises(io_registry.IORegistryError) as exc:
        TestData.read(filename)
    assert str(exc.value).startswith(
        "No reader defined for format 'test1' and class 'TestData'")

    filename = tmpdir.join("testfile.b").strpath
    open(filename, 'w').close()
    with pytest.raises(io_registry.IORegistryError) as exc:
        TestData.read(filename)
    assert str(exc.value).startswith(
        "No reader defined for format 'test2' and class 'TestData'")


def test_write_identifier():

    io_registry.register_identifier('test1', TestData, lambda o, *x, **y: x[0].startswith('a'))
    io_registry.register_identifier('test2', TestData, lambda o, *x, **y: x[0].startswith('b'))

    # Now check that we got past the identifier and are trying to get
    # the reader. The io_registry.get_writer will fail but the error message
    # will tell us if the identifier worked.

    with pytest.raises(io_registry.IORegistryError) as exc:
        TestData().write('abc')
    assert str(exc.value).startswith(
        "No writer defined for format 'test1' and class 'TestData'")

    with pytest.raises(io_registry.IORegistryError) as exc:
        TestData().write('bac')
    assert str(exc.value).startswith(
        "No writer defined for format 'test2' and class 'TestData'")


def test_identifier_origin():

    io_registry.register_identifier('test1', TestData, lambda o, *x, **y: o == 'read')
    io_registry.register_identifier('test2', TestData, lambda o, *x, **y: o == 'write')
    io_registry.register_reader('test1', TestData, empty_reader)
    io_registry.register_writer('test2', TestData, empty_writer)

    # There should not be too many formats defined
    TestData.read()
    TestData().write()

    with pytest.raises(io_registry.IORegistryError) as exc:
        TestData.read(format='test2')
    assert str(exc.value).startswith(
        "No reader defined for format 'test2' and class 'TestData'")

    with pytest.raises(io_registry.IORegistryError) as exc:
        TestData().write(format='test1')
    assert str(exc.value).startswith(
        "No writer defined for format 'test1' and class 'TestData'")


def test_read_valid_return():
    io_registry.register_reader('test', TestData, lambda: TestData())
    t = TestData.read(format='test')
    assert isinstance(t, TestData)


def test_read_invalid_return():
    io_registry.register_reader('test', TestData, lambda: 'spam')
    with pytest.raises(TypeError) as exc:
        TestData.read(format='test')
    assert str(exc.value) == "reader should return a TestData instance"


def test_non_existing_unknown_ext():
    """Raise the correct error when attempting to read a non-existing
    file with an unknown extension."""
    with pytest.raises(IOError):
        data = Table.read('non-existing-file-with-unknown.ext')


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


def teardown_function(function):
    _readers.update(_READERS_ORIGINAL)
    _writers.update(_WRITERS_ORIGINAL)
    _identifiers.update(_IDENTIFIERS_ORIGINAL)


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
