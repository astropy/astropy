# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst


import warnings
from io import StringIO
from collections import OrderedDict
from copy import deepcopy

import numpy as np
import pytest

from ... import units as u
from ... import time
from ... import coordinates
from ... import table
from ..info import serialize_method_as
from ...utils.data_info import data_info_factory, dtype_info_name
from ..table_helpers import simple_table


def test_table_info_attributes(table_types):
    """
    Test the info() method of printing a summary of table column attributes
    """
    a = np.array([1, 2, 3], dtype='int32')
    b = np.array([1, 2, 3], dtype='float32')
    c = np.array(['a', 'c', 'e'], dtype='|S1')
    t = table_types.Table([a, b, c], names=['a', 'b', 'c'])

    # Minimal output for a typical table
    tinfo = t.info(out=None)
    subcls = ['class'] if table_types.Table.__name__ == 'MyTable' else []
    assert tinfo.colnames == ['name', 'dtype', 'shape', 'unit', 'format',
                              'description', 'class', 'n_bad', 'length']
    assert np.all(tinfo['name'] == ['a', 'b', 'c'])
    assert np.all(tinfo['dtype'] == ['int32', 'float32', dtype_info_name('S1')])
    if subcls:
        assert np.all(tinfo['class'] == ['MyColumn'] * 3)

    # All output fields including a mixin column
    t['d'] = [1, 2, 3] * u.m
    t['d'].description = 'quantity'
    t['a'].format = '%02d'
    t['e'] = time.Time([1, 2, 3], format='mjd')
    t['e'].info.description = 'time'
    t['f'] = coordinates.SkyCoord([1, 2, 3], [1, 2, 3], unit='deg')
    t['f'].info.description = 'skycoord'

    tinfo = t.info(out=None)
    assert np.all(tinfo['name'] == 'a b c d e f'.split())
    assert np.all(tinfo['dtype'] == ['int32', 'float32', dtype_info_name('S1'), 'float64',
                                     'object', 'object'])
    assert np.all(tinfo['unit'] == ['', '', '', 'm', '', 'deg,deg'])
    assert np.all(tinfo['format'] == ['%02d', '', '', '', '', ''])
    assert np.all(tinfo['description'] == ['', '', '', 'quantity', 'time', 'skycoord'])
    cls = t.ColumnClass.__name__
    assert np.all(tinfo['class'] == [cls, cls, cls, cls, 'Time', 'SkyCoord'])

    # Test that repr(t.info) is same as t.info()
    out = StringIO()
    t.info(out=out)
    assert repr(t.info) == out.getvalue()


def test_table_info_stats(table_types):
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
    out = StringIO()
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
    assert tinfo.colnames == ['name', 'dtype', 'shape', 'unit', 'format', 'description',
                              'class', 'mean', 'std', 'min', 'max', 'n_bad', 'length']
    assert np.all(tinfo['mean'] == ['1.5', '1.5', '--', '--'])
    assert np.all(tinfo['std'] == ['0.5', '0.5', '--', '--'])
    assert np.all(tinfo['min'] == ['1', '1.0', '--', '1.0'])
    assert np.all(tinfo['max'] == ['2', '2.0', '--', '2.0'])

    out = StringIO()
    t.info('stats', out=out)
    exp = [table_header_line,
           'name mean std min max',
           '---- ---- --- --- ---',
           '   a  1.5 0.5   1   2',
           '   b  1.5 0.5 1.0 2.0',
           '   c   --  --  --  --',
           '   d   --  -- 1.0 2.0']
    assert out.getvalue().splitlines() == exp

    # option = ['attributes', custom]
    custom = data_info_factory(names=['sum', 'first'],
                               funcs=[np.sum, lambda col: col[0]])
    out = StringIO()
    tinfo = t.info(['attributes', custom], out=None)
    assert tinfo.colnames == ['name', 'dtype', 'shape', 'unit', 'format', 'description',
                              'class', 'sum', 'first', 'n_bad', 'length']
    assert np.all(tinfo['name'] == ['a', 'b', 'c', 'd'])
    assert np.all(tinfo['dtype'] == ['int32', 'float32', dtype_info_name('S1'), 'object'])
    assert np.all(tinfo['sum'] == ['6', '6.0', '--', '--'])
    assert np.all(tinfo['first'] == ['1', '1.0', 'a', '1.0'])


def test_data_info():
    """
    Test getting info for just a column.
    """
    cols = [table.Column([1.0, 2.0, np.nan], name='name',
                         description='description', unit='m/s'),
            table.MaskedColumn([1.0, 2.0, 3.0], name='name',
                               description='description', unit='m/s',
                               mask=[False, False, True])]
    for c in cols:
        # Test getting the full ordered dict
        cinfo = c.info(out=None)
        assert cinfo == OrderedDict([('name', 'name'),
                                     ('dtype', 'float64'),
                                     ('shape', ''),
                                     ('unit', 'm / s'),
                                     ('format', ''),
                                     ('description', 'description'),
                                     ('class', type(c).__name__),
                                     ('n_bad', 1),
                                     ('length', 3)])

        # Test the console (string) version which omits trivial values
        out = StringIO()
        c.info(out=out)
        exp = ['name = name',
               'dtype = float64',
               'unit = m / s',
               'description = description',
               'class = {0}'.format(type(c).__name__),
               'n_bad = 1',
               'length = 3']
        assert out.getvalue().splitlines() == exp

        # repr(c.info) gives the same as c.info()
        assert repr(c.info) == out.getvalue()

        # Test stats info
        cinfo = c.info('stats', out=None)
        assert cinfo == OrderedDict([('name', 'name'),
                                     ('mean', '1.5'),
                                     ('std', '0.5'),
                                     ('min', '1.0'),
                                     ('max', '2.0'),
                                     ('n_bad', 1),
                                     ('length', 3)])


def test_data_info_subclass():
    class Column(table.Column):
        """
        Confusingly named Column on purpose, but that is legal.
        """
        pass
    for data in ([], [1, 2]):
        c = Column(data, dtype='int64')
        cinfo = c.info(out=None)
        assert cinfo == OrderedDict([('dtype', 'int64'),
                                     ('shape', ''),
                                     ('unit', ''),
                                     ('format', ''),
                                     ('description', ''),
                                     ('class', 'Column'),
                                     ('n_bad', 0),
                                     ('length', len(data))])


def test_scalar_info():
    """
    Make sure info works with scalar values
    """
    c = time.Time('2000:001')
    cinfo = c.info(out=None)
    assert cinfo['n_bad'] == 0
    assert 'length' not in cinfo


def test_empty_table():
    t = table.Table()
    out = StringIO()
    t.info(out=out)
    exp = ['<Table length=0>', '<No columns>']
    assert out.getvalue().splitlines() == exp


def test_class_attribute():
    """
    Test that class info column is suppressed only for identical non-mixin
    columns.
    """
    vals = [[1] * u.m, [2] * u.m]

    texp = ['<Table length=1>',
            'name  dtype  unit',
            '---- ------- ----',
            'col0 float64    m',
            'col1 float64    m']

    qexp = ['<QTable length=1>',
            'name  dtype  unit  class  ',
            '---- ------- ---- --------',
            'col0 float64    m Quantity',
            'col1 float64    m Quantity']

    for table_cls, exp in ((table.Table, texp),
                           (table.QTable, qexp)):
        t = table_cls(vals)
        out = StringIO()
        t.info(out=out)
        assert out.getvalue().splitlines() == exp


def test_ignore_warnings():
    t = table.Table([[np.nan, np.nan]])
    with warnings.catch_warnings(record=True) as warns:
        t.info('stats', out=None)
        assert len(warns) == 0


def test_no_deprecation_warning():
    # regression test for #5459, where numpy deprecation warnings were
    # emitted unnecessarily.
    t = simple_table()
    with warnings.catch_warnings(record=True) as warns:
        t.info()
        assert len(warns) == 0


def test_lost_parent_error():
    c = table.Column([1, 2, 3], name='a')
    with pytest.raises(AttributeError) as err:
        c[:].info.name
    assert 'failed access "info" attribute' in str(err)


def test_info_serialize_method():
    """
    Unit test of context manager to set info.serialize_method.  Normally just
    used to set this for writing a Table to file (FITS, ECSV, HDF5).
    """
    t = table.Table({'tm': time.Time([1, 2], format='cxcsec'),
                     'sc': coordinates.SkyCoord([1, 2], [1, 2], unit='deg'),
                     'mc': table.MaskedColumn([1, 2], mask=[True, False]),
                     'mc2': table.MaskedColumn([1, 2], mask=[True, False])}
                    )

    origs = {}
    for name in ('tm', 'mc', 'mc2'):
        origs[name] = deepcopy(t[name].info.serialize_method)

    # Test setting by name and getting back to originals
    with serialize_method_as(t, {'tm': 'test_tm', 'mc': 'test_mc'}):
        for name in ('tm', 'mc'):
            assert all(t[name].info.serialize_method[key] == 'test_' + name
                       for key in t[name].info.serialize_method)
        assert t['mc2'].info.serialize_method == origs['mc2']
        assert not hasattr(t['sc'].info, 'serialize_method')

    for name in ('tm', 'mc', 'mc2'):
        assert t[name].info.serialize_method == origs[name]  # dict compare
    assert not hasattr(t['sc'].info, 'serialize_method')

    # Test setting by name and class, where name takes precedence.  Also
    # test that it works for subclasses.
    with serialize_method_as(t, {'tm': 'test_tm', 'mc': 'test_mc',
                                 table.Column: 'test_mc2'}):
        for name in ('tm', 'mc', 'mc2'):
            assert all(t[name].info.serialize_method[key] == 'test_' + name
                       for key in t[name].info.serialize_method)
        assert not hasattr(t['sc'].info, 'serialize_method')

    for name in ('tm', 'mc', 'mc2'):
        assert t[name].info.serialize_method == origs[name]  # dict compare
    assert not hasattr(t['sc'].info, 'serialize_method')

    # Test supplying a single string that all applies to all columns with
    # a serialize_method.
    with serialize_method_as(t, 'test'):
        for name in ('tm', 'mc', 'mc2'):
            assert all(t[name].info.serialize_method[key] == 'test'
                       for key in t[name].info.serialize_method)
        assert not hasattr(t['sc'].info, 'serialize_method')

    for name in ('tm', 'mc', 'mc2'):
        assert t[name].info.serialize_method == origs[name]  # dict compare
    assert not hasattr(t['sc'].info, 'serialize_method')


def test_info_serialize_method_exception():
    """
    Unit test of context manager to set info.serialize_method.  Normally just
    used to set this for writing a Table to file (FITS, ECSV, HDF5).
    """
    t = simple_table(masked=True)
    origs = deepcopy(t['a'].info.serialize_method)
    try:
        with serialize_method_as(t, 'test'):
            assert all(t['a'].info.serialize_method[key] == 'test'
                       for key in t['a'].info.serialize_method)
            raise ZeroDivisionError()
    except ZeroDivisionError:
        pass

    assert t['a'].info.serialize_method == origs  # dict compare
