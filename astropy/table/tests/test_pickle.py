from ...extern.six.moves import cPickle as pickle

import numpy as np
import pytest

from ...table import Table, Column, MaskedColumn


def test_pickle_column(protocol):
    c = Column(data=[1, 2], name='a', format='%05d', description='col a', unit='cm', meta={'a': 1})
    cs = pickle.dumps(c)
    cp = pickle.loads(cs)
    assert np.all(cp == c)
    assert cp.attrs_equal(c)


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


def test_pickle_table(protocol):
    a = Column(data=[1, 2], name='a', format='%05d', description='col a', unit='cm', meta={'a': 1})
    b = Column(data=[3.0, 4.0], name='b', format='%05d', description='col b', unit='cm',
               meta={'b': 1})
    t = Table([a, b], meta={'a': 1})
    ts = pickle.dumps(t)
    tp = pickle.loads(ts)

    assert np.all(tp['a'] == t['a'])
    assert np.all(tp['b'] == t['b'])
    assert tp['a'].attrs_equal(t['a'])
    assert tp['b'].attrs_equal(t['b'])
    assert tp.meta == t.meta


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
