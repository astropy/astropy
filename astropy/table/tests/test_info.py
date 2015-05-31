# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

import numpy as np

from ...extern import six
from ... import units as u
from ... import time
from ...utils.column_info import column_info_factory

STRING8 = 'string8' if six.PY2 else 'bytes8'

def test_info_attributes(table_types):
    """
    Test the info() method of printing a summary of table column attributes
    """
    a = np.array([1, 2, 3], dtype='int32')
    b = np.array([1, 2, 3], dtype='float32')
    c = np.array(['a', 'c', 'e'], dtype='|S1')
    t = table_types.Table([a, b, c], names=['a', 'b', 'c'])

    # Minimal output for a typical table
    tinfo = t.info(out=None)
    assert tinfo.colnames == ['name', 'dtype']
    assert np.all(tinfo['name'] == ['a', 'b', 'c'])
    assert np.all(tinfo['dtype'] == ['int32', 'float32', STRING8])

    # All output fields including a mixin column
    t['d'] = [1,2,3] * u.m
    t['d'].description = 'description'
    t['a'].format = '%02d'
    t['e'] = time.Time([1,2,3], format='mjd')
    tinfo = t.info(out=None)
    assert tinfo.colnames == 'name  dtype  unit format description class'.split()
    assert np.all(tinfo['name'] == 'a b c d e'.split())
    assert np.all(tinfo['dtype'] == ['int32', 'float32', STRING8, 'float64', 'object'])
    assert np.all(tinfo['unit'] == ['', '', '', 'm', ''])
    assert np.all(tinfo['format'] == ['%02d', '', '', '', ''])
    assert np.all(tinfo['description'] == ['', '', '', 'description', ''])
    assert np.all(tinfo['class'] == ['', '', '', '', 'Time'])

def test_info_stats(table_types):
    """
    Test the info() method of printing a summary of table column statistics
    """
    a = np.array([1, 2, 1, 2], dtype='int32')
    b = np.array([1, 2, 1, 2], dtype='float32')
    c = np.array(['a', 'c', 'e', 'f'], dtype='|S1')
    d = time.Time([1, 2, 1, 2], format='mjd')
    t = table_types.Table([a, b, c, d], names=['a', 'b', 'c', 'd'])

    # option = 'stats'
    masked = 'masked=True ' if t.masked else ''
    out = six.moves.cStringIO()
    t.info('stats', out=out)
    table_header_line = '<{0} {1}length=4>'.format(t.__class__.__name__, masked)
    exp = [table_header_line,
           'name mean std min max',
           '---- ---- --- --- ---',
           '   a  1.5 0.5   1   2',
           '   b  1.5 0.5 1.0 2.0',
           '   c   --  --  --  --',
           '   d   --  -- 1.0 2.0']
    assert out.getvalue().splitlines() == exp

    # option = ['attributes', 'stats']
    tinfo = t.info(['attributes', 'stats'], out=None)
    assert tinfo.colnames == 'name  dtype  class mean std min max'.split()
    assert np.all(tinfo['mean'] == ['1.5', '1.5', '--', '--'])
    assert np.all(tinfo['std'] == ['0.5', '0.5', '--', '--'])
    assert np.all(tinfo['min'] == ['1', '1.0', '--', '1.0'])
    assert np.all(tinfo['max'] == ['2', '2.0', '--', '2.0'])

    # option = ['attributes', custom]
    custom = column_info_factory(names=['sum', 'first'],
                                 funcs=[np.sum, lambda col: col[0]])
    out = six.moves.cStringIO()
    tinfo = t.info(['attributes', custom], out=None)
    assert tinfo.colnames == 'name dtype class sum first'.split()
    assert np.all(tinfo['name'] == ['a', 'b', 'c', 'd'])
    assert np.all(tinfo['dtype'] == ['int32', 'float32', STRING8, 'object'])
    assert np.all(tinfo['class'] == ['', '', '', 'Time'])
    assert np.all(tinfo['sum'] == ['6', '6.0', '--', '--'])
    assert np.all(tinfo['first'] == ['1', '1.0', 'a' if six.PY2 else "b'a'", '1.0'])
