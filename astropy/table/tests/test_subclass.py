# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

from ... import table


class MyRow(table.Row):
    def __str__(self):
        return str(self.data)


class MyColumn(table.Column):
    def cool(self):
        return 'Cool!'


class MyMaskedColumn(table.Column):
    def cool(self):
        return 'MaskedCool!'


class MyTableColumns(table.TableColumns):
    def cool(self):
        return 'CoolTableColumns!'


class MyTable(table.Table):
    _Row = MyRow
    _Column = MyColumn
    _MaskedColumn = MyMaskedColumn
    _TableColumns = MyTableColumns


def test_simple_subclass():
    t = MyTable([[1, 2], [3, 4]])
    row = t[0]
    assert isinstance(row, MyRow)
    assert str(row) == '(1, 3)'
    assert t['col0'].cool() == 'Cool!'
    assert t.columns.cool() == 'CoolTableColumns!'

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
    assert t['col0'].cool() == 'MaskedCool!'
    
