# Licensed under a 3-clause BSD style license - see LICENSE.rst


from astropy import table
from astropy.table import pprint


class MyRow(table.Row):
    def __str__(self):
        return str(self.as_void())


class MyColumn(table.Column):
    pass


class MyMaskedColumn(table.MaskedColumn):
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
    assert isinstance(t["col0"], MyColumn)
    assert isinstance(t.columns, MyTableColumns)
    assert isinstance(t.formatter, MyTableFormatter)

    t2 = MyTable(t)
    row = t2[0]
    assert isinstance(row, MyRow)
    assert str(row) == "(1, 3)"

    t3 = table.Table(t)
    row = t3[0]
    assert not isinstance(row, MyRow)
    assert str(row) != "(1, 3)"

    t = MyTable([[1, 2], [3, 4]], masked=True)
    row = t[0]
    assert isinstance(row, MyRow)
    assert str(row) == "(1, 3)"
    assert isinstance(t["col0"], MyMaskedColumn)
    assert isinstance(t.formatter, MyTableFormatter)


class ParamsRow(table.Row):
    """
    Row class that allows access to an arbitrary dict of parameters
    stored as a dict object in the ``params`` column.
    """

    def __getitem__(self, item):
        if item not in self.colnames:
            return super().__getitem__("params")[item]
        else:
            return super().__getitem__(item)

    def keys(self):
        out = [name for name in self.colnames if name != "params"]
        params = [key.lower() for key in sorted(self["params"])]
        return out + params

    def values(self):
        return [self[key] for key in self.keys()]


class ParamsTable(table.Table):
    Row = ParamsRow


def test_params_table():
    t = ParamsTable(names=["a", "b", "params"], dtype=["i", "f", "O"])
    t.add_row((1, 2.0, {"x": 1.5, "y": 2.5}))
    t.add_row((2, 3.0, {"z": "hello", "id": 123123}))
    assert t["params"][0] == {"x": 1.5, "y": 2.5}
    assert t[0]["params"] == {"x": 1.5, "y": 2.5}
    assert t[0]["y"] == 2.5
    assert t[1]["id"] == 123123
    assert list(t[1].keys()) == ["a", "b", "id", "z"]
    assert list(t[1].values()) == [2, 3.0, 123123, "hello"]
