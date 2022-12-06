# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

asdf = pytest.importorskip("asdf")

from astropy.table import Table


def make_table():
    a = [1, 4, 5]
    b = [2.0, 5.0, 8.2]
    c = ["x", "y", "z"]
    return Table([a, b, c], names=("a", "b", "c"), meta={"name": "first table"})


def test_table_io(tmpdir):
    tmpfile = str(tmpdir.join("table.asdf"))

    table = make_table()

    table.write(tmpfile)

    # Simple sanity check using ASDF directly
    with asdf.open(tmpfile) as af:
        assert "data" in af.keys()
        assert isinstance(af["data"], Table)
        assert all(af["data"] == table)

    # Now test using the table reader
    new_t = Table.read(tmpfile)
    assert all(new_t == table)


def test_table_io_custom_key(tmpdir):
    tmpfile = str(tmpdir.join("table.asdf"))

    table = make_table()

    table.write(tmpfile, data_key="something")

    # Simple sanity check using ASDF directly
    with asdf.open(tmpfile) as af:
        assert "something" in af.keys()
        assert "data" not in af.keys()
        assert isinstance(af["something"], Table)
        assert all(af["something"] == table)

    # Now test using the table reader
    with pytest.raises(KeyError):
        new_t = Table.read(tmpfile)

    new_t = Table.read(tmpfile, data_key="something")
    assert all(new_t == table)


def test_table_io_custom_tree(tmpdir):
    tmpfile = str(tmpdir.join("table.asdf"))

    table = make_table()

    def make_custom_tree(tab):
        return dict(foo=dict(bar=tab))

    table.write(tmpfile, make_tree=make_custom_tree)

    # Simple sanity check using ASDF directly
    with asdf.open(tmpfile) as af:
        assert "foo" in af.keys()
        assert "bar" in af["foo"]
        assert "data" not in af.keys()
        assert all(af["foo"]["bar"] == table)

    # Now test using table reader
    with pytest.raises(KeyError):
        new_t = Table.read(tmpfile)

    def find_table(asdffile):
        return asdffile["foo"]["bar"]

    new_t = Table.read(tmpfile, find_table=find_table)
    assert all(new_t == table)
