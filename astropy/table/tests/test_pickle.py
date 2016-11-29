from ...extern.six.moves import cPickle as pickle

import numpy as np

from ...table import Table, Column, MaskedColumn, QTable
from ...table.table_helpers import simple_table
from ...units import Quantity, deg
from ...time import Time
from ...coordinates import SkyCoord


def test_pickle_column(protocol):
    c = Column(data=[1, 2], name='a', format='%05d', description='col a', unit='cm', meta={'a': 1})
    cs = pickle.dumps(c)
    cp = pickle.loads(cs)
    assert np.all(cp == c)
    assert cp.attrs_equal(c)
    assert cp._parent_table is None
    assert repr(c) == repr(cp)


def test_pickle_masked_column(protocol):
    c = MaskedColumn(data=[1, 2], name='a', format='%05d', description='col a', unit='cm',
                     meta={'a': 1})
    c.mask[1] = True
    c.fill_value = -99

    cs = pickle.dumps(c)
    cp = pickle.loads(cs)

    assert np.all(cp._data == c._data)
    assert np.all(cp.mask == c.mask)
    assert cp.attrs_equal(c)
    assert cp.fill_value == -99
    assert cp._parent_table is None
    assert repr(c) == repr(cp)


def test_pickle_multidimensional_column(protocol):
    """Regression test for https://github.com/astropy/astropy/issues/4098"""

    a = np.zeros((3, 2))
    c = Column(a, name='a')
    cs = pickle.dumps(c)
    cp = pickle.loads(cs)

    assert np.all(c == cp)
    assert c.shape == cp.shape
    assert cp.attrs_equal(c)
    assert repr(c) == repr(cp)


def test_pickle_table(protocol):
    a = Column(data=[1, 2], name='a', format='%05d', description='col a', unit='cm', meta={'a': 1})
    b = Column(data=[3.0, 4.0], name='b', format='%05d', description='col b', unit='cm',
               meta={'b': 1})

    for table_class in Table, QTable:
        t = table_class([a, b], meta={'a': 1, 'b': Quantity(10, unit='s')})
        t['c'] = Quantity([1, 2], unit='m')
        t['d'] = Time(['2001-01-02T12:34:56', '2001-02-03T00:01:02'])
        t['e'] = SkyCoord([125.0,180.0]*deg, [-45.0,36.5]*deg)

        ts = pickle.dumps(t)
        tp = pickle.loads(ts)

        assert tp.__class__ is table_class
        assert np.all(tp['a'] == t['a'])
        assert np.all(tp['b'] == t['b'])

        # test mixin columns
        assert np.all(tp['c'] == t['c'])
        assert np.all(tp['d'] == t['d'])
        assert np.all(tp['e'].ra == t['e'].ra)
        assert np.all(tp['e'].dec == t['e'].dec)
        assert type(tp['c']) is type(t['c'])  # nopep8
        assert type(tp['d']) is type(t['d'])  # nopep8
        assert type(tp['e']) is type(t['e'])  # nopep8
        assert tp.meta == t.meta
        assert type(tp) is type(t)

        assert isinstance(tp['c'], Quantity if (table_class is QTable) else Column)

def test_pickle_masked_table(protocol):
    a = Column(data=[1, 2], name='a', format='%05d', description='col a', unit='cm', meta={'a': 1})
    b = Column(data=[3.0, 4.0], name='b', format='%05d', description='col b', unit='cm',
               meta={'b': 1})
    t = Table([a, b], meta={'a': 1}, masked=True)
    t['a'].mask[1] = True
    t['a'].fill_value = -99

    ts = pickle.dumps(t)
    tp = pickle.loads(ts)

    for colname in ('a', 'b'):
        for attr in ('_data', 'mask', 'fill_value'):
            assert np.all(getattr(tp[colname], attr) == getattr(tp[colname], attr))

    assert tp['a'].attrs_equal(t['a'])
    assert tp['b'].attrs_equal(t['b'])
    assert tp.meta == t.meta


def test_pickle_indexed_table(protocol):
    """
    Ensure that any indices that have been added will survive pickling.
    """
    t = simple_table()
    t.add_index('a')
    t.add_index(['a', 'b'])
    ts = pickle.dumps(t)
    tp = pickle.loads(ts)

    assert len(t.indices) == len(tp.indices)
    for index, indexp in zip(t.indices, tp.indices):
        assert np.all(index.data.data == indexp.data.data)
        assert index.data.data.colnames == indexp.data.data.colnames
