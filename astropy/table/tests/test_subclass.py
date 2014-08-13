# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

from ... import table
from .. import pprint

class MyRow(table.Row):
    def __str__(self):
        return str(self.data)


class MyColumn(table.Column):
    pass


class MyMaskedColumn(table.Column):
    pass


class MyTableColumns(table.TableColumns):
    pass


class MyTableFormatter(pprint.TableFormatter):
    pass


class MyTable(table.Table):
    Row = MyRow
    Column = MyColumn
    MaskedColumn = MyMaskedColumn
    TableColumns = MyTableColumns
    TableFormatter = MyTableFormatter

def test_simple_subclass():
    t = MyTable([[1, 2], [3, 4]])
    row = t[0]
    assert isinstance(row, MyRow)
    assert isinstance(t['col0'], MyColumn)
    assert isinstance(t.columns, MyTableColumns)
    assert isinstance(t.formatter, MyTableFormatter)

    t2 = MyTable(t)
    row = t2[0]
    assert isinstance(row, MyRow)
    assert str(row) == '(1, 3)'

    t3 = table.Table(t)
    row = t3[0]
    assert not isinstance(row, MyRow)
    assert str(row) != '(1, 3)'

    t = MyTable([[1, 2], [3, 4]], masked=True)
    row = t[0]
    assert isinstance(row, MyRow)
    assert str(row) == '(1, 3)'
    assert isinstance(t['col0'], MyMaskedColumn)
    assert isinstance(t.formatter, MyTableFormatter)

