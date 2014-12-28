# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module tests some of the methods related to the ``ECSV``
reader/writer.

Requires `pyyaml <http://pyyaml.org/>`_ to be installed.
"""
import os

import numpy as np

from ....table import Table, Column
from ....table.table_helpers import simple_table

from ....tests.helper import pytest
from ....extern.six.moves import StringIO

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

DTYPES = ['bool', 'int8', 'int16', 'int32', 'int64', 'uint8', 'uint16', 'uint32',
          'uint64', 'float16', 'float32', 'float64', 'float128',
          'str']
if os.name == 'nt':
    DTYPES.remove('float128')

T_DTYPES = Table()

for dtype in DTYPES:
    if dtype == 'bool':
        data = np.array([False, True, False])
    elif dtype == 'str':
        data = np.array(['ab 0', 'ab, 1', 'ab2'])
    else:
        data = np.arange(3, dtype=dtype)
    c = Column(data, unit='m / s', description='descr_' + dtype,
               meta={'meta ' + dtype: 1})
    T_DTYPES[dtype] = c

T_DTYPES.meta['comments'] = ['comment1', 'comment2']

@pytest.mark.skipif('not HAS_YAML')
def test_write_simple():
    """
    Write a simple table with common types
    """
    t = simple_table()

    lines = ['# %ECSV 1.0',
             '# ---',
             '# columns:',
             '# - {name: a, type: int32}',
             '# - {name: b, type: float32}',
             '# - {name: c, type: string}',
             'a b c',
             '1 1.0 c',
             '2 2.0 d',
             '3 3.0 e']

    out = StringIO()
    t.write(out, format='ascii.ecsv')
    assert out.getvalue().splitlines() == lines

@pytest.mark.skipif('not HAS_YAML')
def test_write_full():
    """
    Write a full-featured table with common types and explicitly checkout output
    """
    t = T_DTYPES['bool', 'int64', 'float64', 'str']
    lines = ['# %ECSV 1.0',
             '# ---',
             '# columns:',
             '# - name: bool',
             '#   unit: m / s',
             '#   type: bool',
             '#   description: descr_bool',
             '#   meta: {meta bool: 1}',
             '# - name: int64',
             '#   unit: m / s',
             '#   type: int64',
             '#   description: descr_int64',
             '#   meta: {meta int64: 1}',
             '# - name: float64',
             '#   unit: m / s',
             '#   type: float64',
             '#   description: descr_float64',
             '#   meta: {meta float64: 1}',
             '# - name: str',
             '#   unit: m / s',
             '#   type: string',
             '#   description: descr_str',
             '#   meta: {meta str: 1}',
             '# table_meta: !!omap',
             '# - comments: [comment1, comment2]',
             'bool int64 float64 str',
             'False 0 0.0 "ab 0"',
             'True 1 1.0 "ab, 1"',
             'False 2 2.0 ab2']

    out = StringIO()
    t.write(out, format='ascii.ecsv')
    assert out.getvalue().splitlines() == lines


@pytest.mark.skipif('not HAS_YAML')
def test_write_read_roundtrip():
    """
    Write a full-featured table with all types and see that it round-trips on
    readback.
    """
    t = T_DTYPES
    out = StringIO()
    t.write(out, format='ascii.ecsv')

    t2 = Table.read(out.getvalue(), format='ascii.ecsv')
    assert t.meta == t2.meta
    for name in t.colnames:
        assert t[name].attrs_equal(t2[name])
        assert np.all(t[name] == t2[name])


@pytest.mark.skipif('not HAS_YAML')
def test_bad_delimiter():
    """
    Passing a delimiter other than space or comma gives an exception
    """
    pass

@pytest.mark.skipif('not HAS_YAML')
def test_bad_header_start():
    """
    Bad header without initial # %ECSV 1.0
    """
    pass
