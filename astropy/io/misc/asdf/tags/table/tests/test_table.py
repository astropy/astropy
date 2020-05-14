# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import pytest
import numpy as np

import astropy.units as u
from astropy import table
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.io.misc.asdf.tags.helpers import skycoord_equal

from astropy import __minimum_asdf_version__
asdf = pytest.importorskip('asdf', minversion=__minimum_asdf_version__)
from asdf.tests import helpers
from asdf.tags.core.ndarray import NDArrayType

from astropy.io.misc.asdf.tags.tests.helpers import run_schema_example_test


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
                 dtype=[('a', '<i4', (2, 2)),
                        ('b', '<f8'),
                        ('c', '|S1')])

    t = table.Table(a, copy=False)
    assert t.columns['a'].shape == (3, 2, 2)

    def check(ff):
        assert len(ff.blocks) == 1

    helpers.assert_roundtrip_tree({'table': t}, tmpdir, asdf_check_func=check)


def test_structured_array_columns(tmpdir):
    a = np.array([((1, 'a'), 2.0, 'x'),
                  ((4, 'b'), 5.0, 'y'),
                  ((5, 'c'), 8.2, 'z')],
                 dtype=[('a', [('a0', '<i4'), ('a1', '|S1')]),
                        ('b', '<f8'), ('c', '|S1')])

    t = table.Table(a, copy=False)

    def check(ff):
        assert len(ff.blocks) == 1

    helpers.assert_roundtrip_tree({'table': t}, tmpdir, asdf_check_func=check)


def test_table_row_order(tmpdir):
    a = np.array([(1, 2.0, 'x'), (4, 5.0, 'y'),  (5, 8.2, 'z')],
                 dtype=[('a', '<i4'), ('b', '<f8'), ('c', '|S1')])

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
table: !<tag:astropy.org:astropy/table/table-1.0.0>
  columns:
  - !core/column-1.0.0
    data: !core/ndarray-1.0.0
      data: [0, 1, 2]
    name: a
  - !core/column-1.0.0
    data: !core/ndarray-1.0.0
      data: [0, 1, 2, 3]
    name: b
  colnames: [a, b]
    """

    buff = helpers.yaml_to_asdf(yaml)

    with pytest.raises(ValueError) as err:
        with asdf.open(buff) as ff:
            pass
    assert 'Inconsistent data column lengths' in str(err.value)


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


def test_quantity_mixin(tmpdir):

    t = table.QTable()
    t['a'] = [1, 2, 3]
    t['b'] = ['x', 'y', 'z']
    t['c'] = [2.0, 5.0, 8.2] * u.m

    def check(ff):
        assert isinstance(ff['table']['c'], u.Quantity)

    helpers.assert_roundtrip_tree({'table': t}, tmpdir, asdf_check_func=check)


def test_time_mixin(tmpdir):

    t = table.Table()
    t['a'] = [1, 2]
    t['b'] = ['x', 'y']
    t['c'] = Time(['2001-01-02T12:34:56', '2001-02-03T00:01:02'])

    def check(ff):
        assert isinstance(ff['table']['c'], Time)

    helpers.assert_roundtrip_tree({'table': t}, tmpdir, asdf_check_func=check)


def test_timedelta_mixin(tmpdir):

    t = table.Table()
    t['a'] = [1, 2]
    t['b'] = ['x', 'y']
    t['c'] = TimeDelta([1, 2] * u.day)

    def check(ff):
        assert isinstance(ff['table']['c'], TimeDelta)

    helpers.assert_roundtrip_tree({'table': t}, tmpdir, asdf_check_func=check)


def test_skycoord_mixin(tmpdir):

    t = table.Table()
    t['a'] = [1, 2]
    t['b'] = ['x', 'y']
    t['c'] = SkyCoord([1, 2], [3, 4], unit='deg,deg', frame='fk4', obstime='J1990.5')

    def check(ff):
        assert isinstance(ff['table']['c'], SkyCoord)

    def tree_match(old, new):
        NDArrayType.assert_equal(new['a'], old['a'])
        NDArrayType.assert_equal(new['b'], old['b'])
        assert skycoord_equal(new['c'], old['c'])

    helpers.assert_roundtrip_tree({'table': t}, tmpdir, asdf_check_func=check,
                                  tree_match_func=tree_match)


def test_earthlocation_mixin(tmpdir):

    t = table.Table()
    t['a'] = [1, 2]
    t['b'] = ['x', 'y']
    t['c'] = EarthLocation(x=[1, 2] * u.km, y=[3, 4] * u.km, z=[5, 6] * u.km)

    def check(ff):
        assert isinstance(ff['table']['c'], EarthLocation)

    helpers.assert_roundtrip_tree({'table': t}, tmpdir, asdf_check_func=check)


def test_ndarray_mixin(tmpdir):

    t = table.Table()
    t['a'] = [1, 2]
    t['b'] = ['x', 'y']
    t['c'] = table.NdarrayMixin([5, 6])

    helpers.assert_roundtrip_tree({'table': t}, tmpdir)


def test_backwards_compat():
    """
    Make sure that we can continue to read tables that use the schema from
    the ASDF Standard.

    This test uses the examples in the table schema from the ASDF Standard,
    since these make no reference to Astropy's own table definition.
    """

    def check(asdffile):
        assert isinstance(asdffile['example'], table.Table)

    run_schema_example_test('stsci.edu', 'asdf', 'core/table', '1.0.0', check)
