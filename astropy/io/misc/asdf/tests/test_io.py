# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import pytest

import asdf

from astropy.table import Table


def make_table():

    a = [1, 4, 5]
    b = [2.0, 5.0, 8.2]
    c = ['x', 'y', 'z']
    return Table([a, b, c], names=('a', 'b', 'c'), meta={'name': 'first table'})


def test_table_io(tmpdir):

    tmpfile = str(tmpdir.join('table.asdf'))

    table = make_table()

    table.write(tmpfile)

    # Simple sanity check using ASDF directly
    with asdf.open(tmpfile) as af:
        assert 'data' in af.keys()
        assert isinstance(af['data'], Table)
        assert all(af['data'] == table)

    # Now test using the table reader
    new_t = Table.read(tmpfile)
    assert all(new_t == table)


def test_table_io_custom_key(tmpdir):

    tmpfile = str(tmpdir.join('table.asdf'))

    table = make_table()

    table.write(tmpfile, data_key='something')

    # Simple sanity check using ASDF directly
    with asdf.open(tmpfile) as af:
        assert 'something' in af.keys()
        assert 'data' not in af.keys()
        assert isinstance(af['something'], Table)
        assert all(af['something'] == table)

    # Now test using the table reader
    with pytest.raises(KeyError):
        new_t = Table.read(tmpfile)

    new_t = Table.read(tmpfile, data_key='something')
    assert all(new_t == table)
