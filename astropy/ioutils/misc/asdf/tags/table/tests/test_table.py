# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import pytest
import numpy as np

from astropy import table

asdf = pytest.importorskip('asdf', minversion='2.0.0')
from asdf.tests import helpers


def test_table(tmpdir):
    data_rows = [(1, 2.0, 'x'),
                 (4, 5.0, 'y'),
                 (5, 8.2, 'z')]
    t = table.Table(rows=data_rows, names=('a', 'b', 'c'),
                    dtype=('i4', 'f8', 'S1'))
    t.columns['a'].description = 'RA'
    t.columns['a'].unit = 'degree'
    t.columns['a'].meta = {'foo': 'bar'}
    t.columns['c'].description = 'Some description of some sort'

    def check(ff):
        assert len(ff.blocks) == 3

    helpers.assert_roundtrip_tree({'table': t}, tmpdir, asdf_check_func=check)


def test_array_columns(tmpdir):
    a = np.array([([[1, 2], [3, 4]], 2.0, 'x'),
                 ([[5, 6], [7, 8]], 5.0, 'y'),
                  ([[9, 10], [11, 12]], 8.2, 'z')],
                 dtype=[(str('a'), str('<i4'), (2, 2)),
                        (str('b'), str('<f8')),
                        (str('c'), str('|S1'))])

    t = table.Table(a, copy=False)
    assert t.columns['a'].shape == (3, 2, 2)

    def check(ff):
        assert len(ff.blocks) == 1

    helpers.assert_roundtrip_tree({'table': t}, tmpdir, asdf_check_func=check)


def test_structured_array_columns(tmpdir):
    a = np.array([((1, 'a'), 2.0, 'x'),
                  ((4, 'b'), 5.0, 'y'),
                  ((5, 'c'), 8.2, 'z')],
                 dtype=[(str('a'), [(str('a0'), str('<i4')),
                                    (str('a1'), str('|S1'))]),
                        (str('b'), str('<f8')),
                        (str('c'), str('|S1'))])

    t = table.Table(a, copy=False)

    def check(ff):
        assert len(ff.blocks) == 1

    helpers.assert_roundtrip_tree({'table': t}, tmpdir, asdf_check_func=check)


def test_table_row_order(tmpdir):
    a = np.array([(1, 2.0, 'x'),
                  (4, 5.0, 'y'),
                  (5, 8.2, 'z')],
                 dtype=[(str('a'), str('<i4')),
                        (str('b'), str('<f8')),
                        (str('c'), str('|S1'))])

    t = table.Table(a, copy=False)
    t.columns['a'].description = 'RA'
    t.columns['a'].unit = 'degree'
    t.columns['a'].meta = {'foo': 'bar'}
    t.columns['c'].description = 'Some description of some sort'

    def check(ff):
        assert len(ff.blocks) == 1

    helpers.assert_roundtrip_tree({'table': t}, tmpdir, asdf_check_func=check)


def test_table_inline(tmpdir):
    data_rows = [(1, 2.0, 'x'),
                 (4, 5.0, 'y'),
                 (5, 8.2, 'z')]
    t = table.Table(rows=data_rows, names=('a', 'b', 'c'),
                    dtype=('i4', 'f8', 'S1'))
    t.columns['a'].description = 'RA'
    t.columns['a'].unit = 'degree'
    t.columns['a'].meta = {'foo': 'bar'}
    t.columns['c'].description = 'Some description of some sort'

    def check(ff):
        assert len(list(ff.blocks.internal_blocks)) == 0

    helpers.assert_roundtrip_tree({'table': t}, tmpdir, asdf_check_func=check,
                                  write_options={'auto_inline': 64})


def test_mismatched_columns():
    yaml = """
table: !core/table
  columns:
  - !core/column
    data: !core/ndarray
      data: [0, 1, 2]
    name: a
  - !core/column
    data: !core/ndarray
      data: [0, 1, 2, 3]
    name: b
    """

    buff = helpers.yaml_to_asdf(yaml)

    with pytest.raises(ValueError):
        with asdf.AsdfFile.open(buff) as ff:
            pass


def test_masked_table(tmpdir):
    data_rows = [(1, 2.0, 'x'),
                 (4, 5.0, 'y'),
                 (5, 8.2, 'z')]
    t = table.Table(rows=data_rows, names=('a', 'b', 'c'),
                    dtype=('i4', 'f8', 'S1'), masked=True)
    t.columns['a'].description = 'RA'
    t.columns['a'].unit = 'degree'
    t.columns['a'].meta = {'foo': 'bar'}
    t.columns['a'].mask = [True, False, True]
    t.columns['c'].description = 'Some description of some sort'

    def check(ff):
        assert len(ff.blocks) == 4

    helpers.assert_roundtrip_tree({'table': t}, tmpdir, asdf_check_func=check)
