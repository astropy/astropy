# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module tests some of the methods related to the ``ECSV``
reader/writer.

Requires `pyyaml <https://pyyaml.org/>`_ to be installed.
"""
from astropy.table.column import MaskedColumn
import os
import copy
import sys
from io import StringIO

import pytest
import numpy as np

from astropy.table import Table, Column
from astropy.table.table_helpers import simple_table
from astropy.utils.exceptions import AstropyWarning

from astropy.io.ascii.ecsv import DELIMITERS
from astropy.io import ascii
from astropy import units as u
from astropy.utils.compat.optional_deps import HAS_YAML  # noqa

DTYPES = ['bool', 'int8', 'int16', 'int32', 'int64', 'uint8', 'uint16', 'uint32',
          'uint64', 'float16', 'float32', 'float64', 'float128',
          'str']
if os.name == 'nt' or sys.maxsize <= 2**32:
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

# Corresponds to simple_table()
SIMPLE_LINES = ['# %ECSV 0.9',
                '# ---',
                '# datatype:',
                '# - {name: a, datatype: int64}',
                '# - {name: b, datatype: float64}',
                '# - {name: c, datatype: string}',
                '# schema: astropy-2.0',
                'a b c',
                '1 1.0 c',
                '2 2.0 d',
                '3 3.0 e']


@pytest.mark.skipif('not HAS_YAML')
def test_write_simple():
    """
    Write a simple table with common types.  This shows the compact version
    of serialization with one line per column.
    """
    t = simple_table()

    out = StringIO()
    t.write(out, format='ascii.ecsv')
    assert out.getvalue().splitlines() == SIMPLE_LINES


@pytest.mark.skipif('not HAS_YAML')
def test_write_full():
    """
    Write a full-featured table with common types and explicitly checkout output
    """
    t = T_DTYPES['bool', 'int64', 'float64', 'str']
    lines = ['# %ECSV 0.9',
             '# ---',
             '# datatype:',
             '# - name: bool',
             '#   unit: m / s',
             '#   datatype: bool',
             '#   description: descr_bool',
             '#   meta: {meta bool: 1}',
             '# - name: int64',
             '#   unit: m / s',
             '#   datatype: int64',
             '#   description: descr_int64',
             '#   meta: {meta int64: 1}',
             '# - name: float64',
             '#   unit: m / s',
             '#   datatype: float64',
             '#   description: descr_float64',
             '#   meta: {meta float64: 1}',
             '# - name: str',
             '#   unit: m / s',
             '#   datatype: string',
             '#   description: descr_str',
             '#   meta: {meta str: 1}',
             '# meta: !!omap',
             '# - comments: [comment1, comment2]',
             '# schema: astropy-2.0',
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
    readback.  Use both space and comma delimiters.
    """
    t = T_DTYPES
    for delimiter in DELIMITERS:
        out = StringIO()
        t.write(out, format='ascii.ecsv', delimiter=delimiter)

        t2s = [Table.read(out.getvalue(), format='ascii.ecsv'),
               Table.read(out.getvalue(), format='ascii'),
               ascii.read(out.getvalue()),
               ascii.read(out.getvalue(), format='ecsv', guess=False),
               ascii.read(out.getvalue(), format='ecsv')]
        for t2 in t2s:
            assert t.meta == t2.meta
            for name in t.colnames:
                assert t[name].attrs_equal(t2[name])
                assert np.all(t[name] == t2[name])


@pytest.mark.skipif('not HAS_YAML')
def test_bad_delimiter():
    """
    Passing a delimiter other than space or comma gives an exception
    """
    out = StringIO()
    with pytest.raises(ValueError) as err:
        T_DTYPES.write(out, format='ascii.ecsv', delimiter='|')
    assert 'only space and comma are allowed' in str(err.value)


@pytest.mark.skipif('not HAS_YAML')
def test_bad_header_start():
    """
    Bad header without initial # %ECSV x.x
    """
    lines = copy.copy(SIMPLE_LINES)
    lines[0] = '# %ECV 0.9'
    with pytest.raises(ascii.InconsistentTableError):
        Table.read('\n'.join(lines), format='ascii.ecsv', guess=False)


@pytest.mark.skipif('not HAS_YAML')
def test_bad_delimiter_input():
    """
    Illegal delimiter in input
    """
    lines = copy.copy(SIMPLE_LINES)
    lines.insert(2, '# delimiter: |')
    with pytest.raises(ValueError) as err:
        Table.read('\n'.join(lines), format='ascii.ecsv', guess=False)
    assert 'only space and comma are allowed' in str(err.value)


@pytest.mark.skipif('not HAS_YAML')
def test_multidim_input():
    """
    Multi-dimensional column in input
    """
    t = Table()
    t['a'] = np.arange(24).reshape(2, 3, 4)
    t['a'].info.description = 'description'
    t['a'].info.meta = {1: 2}
    t['b'] = [1, 2]

    out = StringIO()
    t.write(out, format='ascii.ecsv')
    t2 = Table.read(out.getvalue(), format='ascii.ecsv')

    assert np.all(t2['a'] == t['a'])
    assert t2['a'].shape == t['a'].shape
    assert t2['a'].dtype == t['a'].dtype
    assert t2['a'].info.description == t['a'].info.description
    assert t2['a'].info.meta == t['a'].info.meta

    assert np.all(t2['b'] == t['b'])


@pytest.mark.skipif('not HAS_YAML')
def test_round_trip_empty_table():
    """Test fix in #5010 for issue #5009 (ECSV fails for empty type with bool type)"""
    t = Table(dtype=[bool, 'i', 'f'], names=['a', 'b', 'c'])
    out = StringIO()
    t.write(out, format='ascii.ecsv')
    t2 = Table.read(out.getvalue(), format='ascii.ecsv')
    assert t.dtype == t2.dtype
    assert len(t2) == 0


@pytest.mark.skipif('not HAS_YAML')
def test_csv_ecsv_colnames_mismatch():
    """
    Test that mismatch in column names from normal CSV header vs.
    ECSV YAML header raises the expected exception.
    """
    lines = copy.copy(SIMPLE_LINES)
    header_index = lines.index('a b c')
    lines[header_index] = 'a b d'
    with pytest.raises(ValueError) as err:
        ascii.read(lines, format='ecsv')
    assert "column names from ECSV header ['a', 'b', 'c']" in str(err.value)


@pytest.mark.skipif('not HAS_YAML')
def test_regression_5604():
    """
    See https://github.com/astropy/astropy/issues/5604 for more.
    """
    t = Table()
    t.meta = {"foo": 5 * u.km, "foo2": u.s}
    t["bar"] = [7] * u.km

    out = StringIO()
    t.write(out, format="ascii.ecsv")

    assert '!astropy.units.Unit' in out.getvalue()
    assert '!astropy.units.Quantity' in out.getvalue()


@pytest.mark.skipif('HAS_YAML')
def test_ecsv_but_no_yaml_warning():
    """
    Test that trying to read an ECSV without PyYAML installed when guessing
    emits a warning, but reading with guess=False gives an exception.
    """
    with pytest.warns(AstropyWarning, match=r'file looks like ECSV format but '
                      'PyYAML is not installed') as w:
        ascii.read(SIMPLE_LINES)
    assert len(w) == 1

    with pytest.raises(ascii.InconsistentTableError, match='unable to parse yaml'), \
            pytest.warns(AstropyWarning, match=r'PyYAML is not installed'):
        ascii.read(SIMPLE_LINES, format='ecsv')


@pytest.mark.skipif('not HAS_YAML')
def test_round_trip_masked_table_default(tmpdir):
    """Test (mostly) round-trip of MaskedColumn through ECSV using default serialization
    that uses an empty string "" to mark NULL values.  Note:

    >>> simple_table(masked=True)
    <Table masked=True length=3>
      a      b     c
    int64 float64 str1
    ----- ------- ----
       --     1.0    c
        2     2.0   --
        3      --    e
    """
    filename = str(tmpdir.join('test.ecsv'))

    t = simple_table(masked=True)  # int, float, and str cols with one masked element
    t.write(filename)

    t2 = Table.read(filename)
    assert t2.masked is False
    assert t2.colnames == t.colnames
    for name in t2.colnames:
        # From formal perspective the round-trip columns are the "same"
        assert np.all(t2[name].mask == t[name].mask)
        assert np.all(t2[name] == t[name])

        # But peeking under the mask shows that the underlying data are changed
        # because by default ECSV uses "" to represent masked elements.
        t[name].mask = False
        t2[name].mask = False
        assert not np.all(t2[name] == t[name])  # Expected diff


@pytest.mark.skipif('not HAS_YAML')
def test_read_masked_bool():
    txt = """\
# %ECSV 0.9
# ---
# datatype:
# - {name: col0, datatype: bool}
# schema: astropy-2.0
col0
1
0
True
""
False
"""
    dat = ascii.read(txt, format='ecsv')
    col = dat['col0']
    assert isinstance(col, MaskedColumn)
    assert np.all(col.mask == [False, False, False, True, False])
    assert np.all(col == [True, False, True, False, False])


@pytest.mark.skipif('not HAS_YAML')
@pytest.mark.parametrize('serialize_method', ['null_value', 'data_mask'])
@pytest.mark.parametrize('dtype', [np.int64, np.float64, np.bool, np.str])
def test_roundtrip_multidim_masked_array(serialize_method, dtype):
    # TODO also test empty string with null value
    t = Table()
    col = MaskedColumn(np.arange(12).reshape(2, 3, 2), dtype=dtype)
    col.mask[0, 0, 0] = True
    col.mask[1, 1, 1] = True
    t['a'] = col
    out = StringIO()
    t.write(out, format='ascii.ecsv', serialize_method=serialize_method)
    t2 = Table.read(out.getvalue(), format='ascii.ecsv')

    assert t2.masked is False
    assert t2.colnames == t.colnames
    for name in t2.colnames:
        assert np.all(t2[name].mask == t[name].mask)
        assert np.all(t2[name] == t[name])
