# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import os
import urllib.parse

import yaml

import pytest
import numpy as np

import astropy.units as u
from astropy import table
from astropy import __minimum_asdf_version__

asdf = pytest.importorskip('asdf', minversion=__minimum_asdf_version__)

from asdf import treeutil
from asdf.tests import helpers
from asdf.types import format_tag
from asdf.resolver import default_resolver


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
    assert 'Inconsistent data column lengths' in str(err)


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


def test_backwards_compat(tmpdir):
    """
    Make sure that we can continue to read tables that use the schema from
    the ASDF Standard.

    This test uses the examples in the table schema from the ASDF Standard,
    since these make no reference to Astropy's own table definition.
    """

    tag = format_tag('stsci.edu', 'asdf', '1.0.0', 'core/table')
    schema_path = urllib.parse.urlparse(default_resolver(tag)).path

    with open(schema_path, 'rb') as ff:
        schema = yaml.load(ff)

    examples = []
    for node in treeutil.iter_tree(schema):
        if (isinstance(node, dict) and
            'examples' in node and
            isinstance(node['examples'], list)):
            for desc, example in node['examples']:
                examples.append(example)

    for example in examples:
        buff = helpers.yaml_to_asdf('example: ' + example.strip())
        ff = asdf.AsdfFile(uri=schema_path)
        # Add some dummy blocks so that the ndarray examples work
        for i in range(3):
            b = asdf.block.Block(np.zeros((1024*1024*8), dtype=np.uint8))
            b._used = True
            ff.blocks.add(b)
        ff._open_impl(ff, buff, mode='r')
        assert isinstance(ff['example'], table.Table)
