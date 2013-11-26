import numpy as np

from ...table import Table, Column
from ...time import Time

a = Column([3, 4, 5], name='a')
time_obj = Time([1, 2, 3], format='unix').replicate(format='iso')
time_obj.name = 'b'


def test_time_mixin_init():
    """
    Basic test of initializing a table with a mixin column.
    """
    t1 = Table([a])
    t1['b'] = time_obj

    t2 = Table([a, time_obj])

    t3 = Table()
    t3['b'] = time_obj
    t3.add_column(a, index=0)

    for t in (t1, t2, t3):
        assert t.colnames == ['a', 'b']
        assert t._data.dtype.names == ('a', 'b__jd1', 'b__jd2')
        assert type(t['b']) is Time
        assert t.pformat() == [' a             b           ',
                               '--- -----------------------',
                               '  3 1970-01-01 00:00:01.000',
                               '  4 1970-01-01 00:00:02.000',
                               '  5 1970-01-01 00:00:03.000']


def test_time_mixin_slice():
    """
    Test slicing with mixin columns.
    """
    t = Table([a, time_obj])
    t = t[:2]
    assert t.colnames == ['a', 'b']
    assert t._data.dtype.names == ('a', 'b__jd1', 'b__jd2')
    assert type(t['b']) is Time
    assert t.pformat() == [' a             b           ',
                           '--- -----------------------',
                           '  3 1970-01-01 00:00:01.000',
                           '  4 1970-01-01 00:00:02.000']


def test_time_mixin_item_access():
    """
    Test slicing with mixin columns.
    """
    t = Table([a, time_obj])
    assert t[1]['b'].iso == '1970-01-01 00:00:02.000'
    assert np.all(t['b'].iso == ['1970-01-01 00:00:01.000',
                                 '1970-01-01 00:00:02.000',
                                 '1970-01-01 00:00:03.000'])
