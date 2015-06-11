# Licensed under a 3-clause BSD style license - see LICENSE.rst

try:
    from cStringIO import StringIO
except ImportError: # cStringIO doesn't exist in Python 3
    from io import BytesIO
    StringIO = lambda x: BytesIO(x.encode('ascii'))
import os
import functools

from textwrap import dedent

import numpy as np
from numpy import ma

from ....table import Table, MaskedColumn
from ... import ascii
from ...ascii.core import ParameterError, FastOptionsError
from ...ascii.cparser import CParserError
from ..fastbasic import FastBasic, FastCsv, FastTab, FastCommentedHeader, \
    FastRdb, FastNoHeader
from .common import assert_equal, assert_almost_equal, assert_true
from ....tests.helper import pytest
from ....extern import six

TRAVIS = os.environ.get('TRAVIS', False)

def assert_table_equal(t1, t2, check_meta=False):
    assert_equal(len(t1), len(t2))
    assert_equal(t1.colnames, t2.colnames)
    if check_meta:
        assert_equal(t1.meta, t2.meta)
    for name in t1.colnames:
        if len(t1) != 0:
            assert_equal(t1[name].dtype.kind, t2[name].dtype.kind)
        if not isinstance(t1[name], MaskedColumn):
            for i, el in enumerate(t1[name]):
                try:
                    if not isinstance(el, six.string_types) and np.isnan(el):
                        assert_true(not isinstance(t2[name][i], six.string_types) and np.isnan(t2[name][i]))
                    elif isinstance(el, six.string_types):
                        assert_equal(el, t2[name][i])
                    else:
                        assert_almost_equal(el, t2[name][i])
                except (TypeError, NotImplementedError):
                    pass # ignore for now

# Use this counter to create a unique filename for each file created in a test
# if this function is called more than once in a single test
_filename_counter = 0

def _read(tmpdir, table, Reader=None, format=None, parallel=False, check_meta=False, **kwargs):
    # make sure we have a newline so table can't be misinterpreted as a filename
    global _filename_counter

    table += '\n'
    reader = Reader(**kwargs)
    t1 = reader.read(table)
    t2 = reader.read(StringIO(table))
    t3 = reader.read(table.splitlines())
    t4 = ascii.read(table, format=format, guess=False, **kwargs)
    t5 = ascii.read(table, format=format, guess=False, fast_reader=False, **kwargs)
    assert_table_equal(t1, t2, check_meta=check_meta)
    assert_table_equal(t2, t3, check_meta=check_meta)
    assert_table_equal(t3, t4, check_meta=check_meta)
    assert_table_equal(t4, t5, check_meta=check_meta)

    if parallel:
        if TRAVIS:
            pytest.xfail("Multiprocessing can sometimes fail on Travis CI")
        elif os.name == 'nt':
            pytest.xfail("Multiprocessing is currently unsupported on Windows")
        t6 = ascii.read(table, format=format, guess=False, fast_reader={
            'parallel': True}, **kwargs)
        assert_table_equal(t1, t6, check_meta=check_meta)

    filename = str(tmpdir.join('table{0}.txt'.format(_filename_counter)))
    _filename_counter += 1

    with open(filename, 'wb') as f:
        f.write(table.encode('ascii'))
        f.flush()

    t7 = ascii.read(filename, format=format, guess=False, **kwargs)
    if parallel:
        t8 = ascii.read(filename, format=format, guess=False, fast_reader={
            'parallel': True}, **kwargs)

    assert_table_equal(t1, t7, check_meta=check_meta)
    if parallel:
        assert_table_equal(t1, t8, check_meta=check_meta)
    return t1


@pytest.fixture(scope='function')
def read_basic(tmpdir, request):
    return functools.partial(_read, tmpdir, Reader=FastBasic, format='basic')


@pytest.fixture(scope='function')
def read_csv(tmpdir, request):
    return functools.partial(_read, tmpdir, Reader=FastCsv, format='csv')


@pytest.fixture(scope='function')
def read_tab(tmpdir, request):
    return functools.partial(_read, tmpdir, Reader=FastTab, format='tab')


@pytest.fixture(scope='function')
def read_commented_header(tmpdir, request):
    return functools.partial(_read, tmpdir, Reader=FastCommentedHeader,
                             format='commented_header')


@pytest.fixture(scope='function')
def read_rdb(tmpdir, request):
    return functools.partial(_read, tmpdir, Reader=FastRdb, format='rdb')


@pytest.fixture(scope='function')
def read_no_header(tmpdir, request):
    return functools.partial(_read, tmpdir, Reader=FastNoHeader,
                             format='no_header')


@pytest.mark.parametrize("parallel", [True, False])
def test_simple_data(parallel, read_basic):
    """
    Make sure the fast reader works with basic input data.
    """
    table = read_basic("A B C\n1 2 3\n4 5 6", parallel=parallel)
    expected = Table([[1, 4], [2, 5], [3, 6]], names=('A', 'B', 'C'))
    assert_table_equal(table, expected)


def test_read_types():
    """
    Make sure that the read() function takes filenames,
    strings, and lists of strings in addition to file-like objects.
    """
    t1 = ascii.read("a b c\n1 2 3\n4 5 6", format='fast_basic', guess=False)
    #TODO: also read from file
    t2 = ascii.read(StringIO("a b c\n1 2 3\n4 5 6"), format='fast_basic', guess=False)
    t3 = ascii.read(["a b c", "1 2 3", "4 5 6"], format='fast_basic', guess=False)
    assert_table_equal(t1, t2)
    assert_table_equal(t2, t3)


@pytest.mark.parametrize("parallel", [True, False])
def test_supplied_names(parallel, read_basic):
    """
    If passed as a parameter, names should replace any
    column names found in the header.
    """
    table = read_basic("A B C\n1 2 3\n4 5 6", names=('X', 'Y', 'Z'), parallel=parallel)
    expected = Table([[1, 4], [2, 5], [3, 6]], names=('X', 'Y', 'Z'))
    assert_table_equal(table, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_no_header(parallel, read_basic, read_no_header):
    """
    The header should not be read when header_start=None. Unless names is
    passed, the column names should be auto-generated.
    """
    t1 = read_basic("A B C\n1 2 3\n4 5 6", header_start=None, data_start=0, parallel=parallel)
    t2 = read_no_header("A B C\n1 2 3\n4 5 6", parallel=parallel)
    expected = Table([['A', '1', '4'], ['B', '2', '5'], ['C', '3', '6']], names=('col1', 'col2', 'col3'))
    assert_table_equal(t1, expected)
    assert_table_equal(t2, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_no_header_supplied_names(parallel, read_basic):
    """
    If header_start=None and names is passed as a parameter, header
    data should not be read and names should be used instead.
    """
    table = read_basic("A B C\n1 2 3\n4 5 6", header_start=None, data_start=0,
                       names=('X', 'Y', 'Z'), parallel=parallel)
    expected = Table([['A', '1', '4'], ['B', '2', '5'], ['C', '3', '6']], names=('X', 'Y', 'Z'))
    assert_table_equal(table, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_comment(parallel, read_basic):
    """
    Make sure that line comments are ignored by the C reader.
    """
    table = read_basic("# comment\nA B C\n # another comment\n1 2 3\n4 5 6", parallel=parallel)
    expected = Table([[1, 4], [2, 5], [3, 6]], names=('A', 'B', 'C'))
    assert_table_equal(table, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_empty_lines(parallel, read_basic):
    """
    Make sure that empty lines are ignored by the C reader.
    """
    table = read_basic("\n\nA B C\n1 2 3\n\n\n4 5 6\n\n\n\n", parallel=parallel)
    expected = Table([[1, 4], [2, 5], [3, 6]], names=('A', 'B', 'C'))
    assert_table_equal(table, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_lstrip_whitespace(parallel, read_basic):
    """
    Test to make sure the reader ignores whitespace at the beginning of fields.
    """
    text = """
     1,  2,   \t3
 A,\t\t B,  C
  a, b,   c
""" + '  \n'

    table = read_basic(text, delimiter=',', parallel=parallel)
    expected = Table([['A', 'a'], ['B', 'b'], ['C', 'c']], names=('1', '2', '3'))
    assert_table_equal(table, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_rstrip_whitespace(parallel, read_basic):
    """
    Test to make sure the reader ignores whitespace at the end of fields.
    """
    text = ' 1 ,2 \t,3  \nA\t,B ,C\t \t \n  \ta ,b , c \n'
    table = read_basic(text, delimiter=',', parallel=parallel)
    expected = Table([['A', 'a'], ['B', 'b'], ['C', 'c']], names=('1', '2', '3'))
    assert_table_equal(table, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_conversion(parallel, read_basic):
    """
    The reader should try to convert each column to ints. If this fails, the
    reader should try to convert to floats. Failing this, it should fall back
    to strings.
    """
    text = """
A B C D E
1 a 3 4 5
2. 1 9 10 -5.3e4
4 2 -12 .4 six
"""
    table = read_basic(text, parallel=parallel)
    assert_equal(table['A'].dtype.kind, 'f')
    assert table['B'].dtype.kind in ('S', 'U')
    assert_equal(table['C'].dtype.kind, 'i')
    assert_equal(table['D'].dtype.kind, 'f')
    assert table['E'].dtype.kind in ('S', 'U')


@pytest.mark.parametrize("parallel", [True, False])
def test_delimiter(parallel, read_basic):
    """
    Make sure that different delimiters work as expected.
    """
    text = """
COL1 COL2 COL3
1 A -1
2 B -2
"""
    expected = Table([[1, 2], ['A', 'B'], [-1, -2]], names=('COL1', 'COL2', 'COL3'))

    for sep in ' ,\t#;':
        table = read_basic(text.replace(' ', sep), delimiter=sep, parallel=parallel)
        assert_table_equal(table, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_include_names(parallel, read_basic):
    """
    If include_names is not None, the parser should read only those columns in include_names.
    """
    table = read_basic("A B C D\n1 2 3 4\n5 6 7 8", include_names=['A', 'D'], parallel=parallel)
    expected = Table([[1, 5], [4, 8]], names=('A', 'D'))
    assert_table_equal(table, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_exclude_names(parallel, read_basic):
    """
    If exclude_names is not None, the parser should exclude the columns in exclude_names.
    """
    table = read_basic("A B C D\n1 2 3 4\n5 6 7 8", exclude_names=['A', 'D'], parallel=parallel)
    expected = Table([[2, 6], [3, 7]], names=('B', 'C'))
    assert_table_equal(table, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_include_exclude_names(parallel, read_basic):
    """
    Make sure that include_names is applied before exclude_names if both are specified.
    """
    text = """
A B C D E F G H
1 2 3 4 5 6 7 8
9 10 11 12 13 14 15 16
"""
    table = read_basic(text, include_names=['A', 'B', 'D', 'F', 'H'],
                       exclude_names=['B', 'F'], parallel=parallel)
    expected = Table([[1, 9], [4, 12], [8, 16]], names=('A', 'D', 'H'))
    assert_table_equal(table, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_quoted_fields(parallel, read_basic):
    """
    The character quotechar (default '"') should denote the start of a field which can
    contain the field delimiter and newlines.
    """
    if parallel:
        pytest.xfail("Multiprocessing can fail with quoted fields")
    text = """
"A B" C D
1.5 2.1 -37.1
a b "   c
 d"
"""
    table = read_basic(text, parallel=parallel)
    expected = Table([['1.5', 'a'], ['2.1', 'b'], ['-37.1', 'cd']], names=('A B', 'C', 'D'))
    assert_table_equal(table, expected)
    table = read_basic(text.replace('"', "'"), quotechar="'", parallel=parallel)
    assert_table_equal(table, expected)


def test_invalid_parameters():
    """
    Make sure the C reader raises an error if passed parameters it can't handle.
    """
    int_converter = ascii.convert_numpy(np.uint)
    converters = dict((i + 1, ascii.convert_numpy(np.uint)) for i in range(3))
    invalid_params = {'delimiter': ',,', # multi-char delimiter
                      'comment': '##', # multi-char comment
                      'data_start': None, # data_start=None
                      'quotechar': '##', # multi-char quote signifier
                      'data_start': -1, # negative data_start
                      'header_start': -1, # negative header_start
                      'converters': converters, # passing converters
                      'Inputter': ascii.ContinuationLinesInputter, # passing Inputter
                      'header_Splitter': ascii.DefaultSplitter, # passing Splitter
                      'data_Splitter': ascii.DefaultSplitter
                      }
    for key, val in invalid_params.items():
        with pytest.raises(ParameterError):
            print('Trying {0}={1} using constructor'.format(key, val))
            table = FastBasic(**{key: val}).read('1 2 3\n4 5 6')
        with pytest.raises(ParameterError):
            print('Trying {0}={1} using ascii.read'.format(key, val))
            table = ascii.read('1 2 3\n4 5 6', format='fast_basic', guess=False, **{key: val})

    with pytest.raises(TypeError):
        table = FastBasic(foo=7).read('1 2 3\n4 5 6') # unexpected argument
    with pytest.raises(FastOptionsError): # don't fall back on the slow reader
        table = ascii.read('1 2 3\n4 5 6', format='basic', fast_reader={'foo': 7})
    with pytest.raises(ParameterError):
        # Outputter cannot be specified in constructor
        table = FastBasic(Outputter=ascii.TableOutputter).read('1 2 3\n4 5 6')


def test_too_many_cols():
    """
    If a row contains too many columns, the C reader should raise an error.
    """
    text = """
A B C
1 2 3
4 5 6
7 8 9 10
11 12 13
"""
    with pytest.raises(CParserError) as e:
        table = FastBasic().read(text)
    assert 'CParserError: an error occurred while parsing table data: too many ' \
        'columns found in line 3 of data' in str(e)


@pytest.mark.parametrize("parallel", [True, False])
def test_not_enough_cols(parallel, read_csv):
    """
    If a row does not have enough columns, the FastCsv reader should add empty
    fields while the FastBasic reader should raise an error.
    """
    text = """
A,B,C
1,2,3
4,5
6,7,8
"""
    table = read_csv(text, parallel=parallel)
    assert table['B'][1] is not ma.masked
    assert table['C'][1] is ma.masked

    with pytest.raises(CParserError) as e:
        table = FastBasic(delimiter=',').read(text)


@pytest.mark.parametrize("parallel", [True, False])
def test_data_end(parallel, read_basic, read_rdb):
    """
    The parameter data_end should specify where data reading ends.
    """
    text = """
A B C
1 2 3
4 5 6
7 8 9
10 11 12
"""
    table = read_basic(text, data_end=3, parallel=parallel)
    expected = Table([[1, 4], [2, 5], [3, 6]], names=('A', 'B', 'C'))
    assert_table_equal(table, expected)

    # data_end supports negative indexing
    table = read_basic(text, data_end=-2, parallel=parallel)
    assert_table_equal(table, expected)

    text = """
A\tB\tC
N\tN\tS
1\t2\ta
3\t4\tb
5\t6\tc
"""
    # make sure data_end works with RDB
    table = read_rdb(text, data_end=-1, parallel=parallel)
    expected = Table([[1, 3], [2, 4], ['a', 'b']], names=('A', 'B', 'C'))
    assert_table_equal(table, expected)

    # positive index
    table = read_rdb(text, data_end=3, parallel=parallel)
    expected = Table([[1], [2], ['a']], names=('A', 'B', 'C'))
    assert_table_equal(table, expected)

    # empty table if data_end is too small
    table = read_rdb(text, data_end=1, parallel=parallel)
    expected = Table([[], [], []], names=('A', 'B', 'C'))
    assert_table_equal(table, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_inf_nan(parallel, read_basic):
    """
    Test that inf and nan-like values are correctly parsed on all platforms.

    Regression test for https://github.com/astropy/astropy/pull/3525
    """

    text = dedent("""\
        A
        nan
        +nan
        -nan
        inf
        infinity
        +inf
        +infinity
        -inf
        -infinity
    """)

    expected = Table({'A': [np.nan, np.nan, np.nan,
                            np.inf, np.inf, np.inf, np.inf,
                            -np.inf, -np.inf]})

    table = read_basic(text, parallel=parallel)
    assert table['A'].dtype.kind == 'f'
    assert_table_equal(table, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_fill_values(parallel, read_basic):
    """
    Make sure that the parameter fill_values works as intended. If fill_values
    is not specified, the default behavior should be to convert '' to 0.
    """
    text = """
A, B, C
, 2, nan
a, -999, -3.4
nan, 5, -9999
8, nan, 7.6e12
"""
    table = read_basic(text, delimiter=',', parallel=parallel)
    # The empty value in row A should become a masked '0'
    assert isinstance(table['A'], MaskedColumn)
    assert table['A'][0] is ma.masked
    # '0' rather than 0 because there is a string in the column
    assert_equal(table['A'].data.data[0], '0')
    assert table['A'][1] is not ma.masked

    table = read_basic(text, delimiter=',', fill_values=('-999', '0'), parallel=parallel)
    assert isinstance(table['B'], MaskedColumn)
    assert table['A'][0] is not ma.masked # empty value unaffected
    assert table['C'][2] is not ma.masked # -9999 is not an exact match
    assert table['B'][1] is ma.masked
    # Numeric because the rest of the column contains numeric data
    assert_equal(table['B'].data.data[1], 0.0)
    assert table['B'][0] is not ma.masked

    table = read_basic(text, delimiter=',', fill_values=[], parallel=parallel)
    # None of the columns should be masked
    for name in 'ABC':
        assert not isinstance(table[name], MaskedColumn)

    table = read_basic(text, delimiter=',', fill_values=[('', '0', 'A'),
                                ('nan', '999', 'A', 'C')], parallel=parallel)
    assert np.isnan(table['B'][3]) # nan filling skips column B
    assert table['B'][3] is not ma.masked # should skip masking as well as replacing nan
    assert table['A'][0] is ma.masked
    assert table['A'][2] is ma.masked
    assert_equal(table['A'].data.data[0], '0')
    assert_equal(table['A'].data.data[2], '999')
    assert table['C'][0] is ma.masked
    assert_almost_equal(table['C'].data.data[0], 999.0)
    assert_almost_equal(table['C'][1], -3.4) # column is still of type float


@pytest.mark.parametrize("parallel", [True, False])
def test_fill_include_exclude_names(parallel, read_csv):
    """
    fill_include_names and fill_exclude_names should filter missing/empty value handling
    in the same way that include_names and exclude_names filter output columns.
    """
    text = """
A, B, C
, 1, 2
3, , 4
5, 5,
"""
    table = read_csv(text, fill_include_names=['A', 'B'], parallel=parallel)
    assert table['A'][0] is ma.masked
    assert table['B'][1] is ma.masked
    assert table['C'][2] is not ma.masked # C not in fill_include_names

    table = read_csv(text, fill_exclude_names=['A', 'B'], parallel=parallel)
    assert table['C'][2] is ma.masked
    assert table['A'][0] is not ma.masked
    assert table['B'][1] is not ma.masked # A and B excluded from fill handling

    table = read_csv(text, fill_include_names=['A', 'B'], fill_exclude_names=['B'], parallel=parallel)
    assert table['A'][0] is ma.masked
    assert table['B'][1] is not ma.masked # fill_exclude_names applies after fill_include_names
    assert table['C'][2] is not ma.masked


@pytest.mark.parametrize("parallel", [True, False])
def test_many_rows(parallel, read_basic):
    """
    Make sure memory reallocation works okay when the number of rows
    is large (so that each column string is longer than INITIAL_COL_SIZE).
    """
    text = 'A B C\n'
    for i in range(500): # create 500 rows
        text += ' '.join([str(i) for i in range(3)])
        text += '\n'

    table = read_basic(text, parallel=parallel)
    expected = Table([[0] * 500, [1] * 500, [2] * 500], names=('A', 'B', 'C'))
    assert_table_equal(table, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_many_columns(parallel, read_basic):
    """
    Make sure memory reallocation works okay when the number of columns
    is large (so that each header string is longer than INITIAL_HEADER_SIZE).
    """
    # create a string with 500 columns and two data rows
    text = ' '.join([str(i) for i in range(500)])
    text += ('\n' + text + '\n' + text)
    table = read_basic(text, parallel=parallel)
    expected = Table([[i, i] for i in range(500)], names=[str(i) for i in range(500)])
    assert_table_equal(table, expected)


def test_fast_reader():
    """
    Make sure that ascii.read() works as expected by default and with
    fast_reader specified.
    """
    text = 'a b c\n1 2 3\n4 5 6'
    with pytest.raises(ParameterError): # C reader can't handle regex comment
        ascii.read(text, format='fast_basic', guess=False, comment='##')

    # Enable multiprocessing and the fast converter
    try:
        ascii.read(text, format='basic', guess=False,
                   fast_reader={'parallel': True, 'use_fast_converter': True})
    except NotImplementedError:
        # Might get this on Windows, try without parallel...
        if os.name == 'nt':
            ascii.read(text, format='basic', guess=False,
                       fast_reader={'parallel': False,
                                    'use_fast_converter': True})
        else:
            raise

    # Should raise an error if fast_reader has an invalid key
    with pytest.raises(FastOptionsError):
        ascii.read(text, format='fast_basic', guess=False, fast_reader={'foo': True})

    # Use the slow reader instead
    ascii.read(text, format='basic', guess=False, comment='##', fast_reader=False)
    # Will try the slow reader afterwards by default
    ascii.read(text, format='basic', guess=False, comment='##')


@pytest.mark.parametrize("parallel", [True, False])
def test_read_tab(parallel, read_tab):
    """
    The fast reader for tab-separated values should not strip whitespace, unlike
    the basic reader.
    """
    if parallel:
        pytest.xfail("Multiprocessing can fail with quoted fields")
    text = '1\t2\t3\n  a\t b \t\n c\t" d\n e"\t  '
    table = read_tab(text, parallel=parallel)
    assert_equal(table['1'][0], '  a')   # preserve line whitespace
    assert_equal(table['2'][0], ' b ')   # preserve field whitespace
    assert table['3'][0] is ma.masked    # empty value should be masked
    assert_equal(table['2'][1], ' d e')  # preserve whitespace in quoted fields
    assert_equal(table['3'][1], '  ')    # preserve end-of-line whitespace


@pytest.mark.parametrize("parallel", [True, False])
def test_default_data_start(parallel, read_basic):
    """
    If data_start is not explicitly passed to read(), data processing should
    beginning right after the header.
    """
    text = 'ignore this line\na b c\n1 2 3\n4 5 6'
    table = read_basic(text, header_start=1, parallel=parallel)
    expected = Table([[1, 4], [2, 5], [3, 6]], names=('a', 'b', 'c'))
    assert_table_equal(table, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_commented_header(parallel, read_commented_header):
    """
    The FastCommentedHeader reader should mimic the behavior of the
    CommentedHeader by overriding the default header behavior of FastBasic.
    """
    text = """
 # A B C
 1 2 3
 4 5 6
"""
    t1 = read_commented_header(text, parallel=parallel)
    expected = Table([[1, 4], [2, 5], [3, 6]], names=('A', 'B', 'C'))
    assert_table_equal(t1, expected)

    text = '# first commented line\n # second commented line\n\n' + text
    t2 = read_commented_header(text, header_start=2, data_start=0, parallel=parallel)
    assert_table_equal(t2, expected)
    t3 = read_commented_header(text, header_start=-1, data_start=0, parallel=parallel) # negative indexing allowed
    assert_table_equal(t3, expected)

    text += '7 8 9'
    t4 = read_commented_header(text, header_start=2, data_start=2, parallel=parallel)
    expected = Table([[7], [8], [9]], names=('A', 'B', 'C'))
    assert_table_equal(t4, expected)

    with pytest.raises(ParameterError):
        read_commented_header(text, header_start=-1, data_start=-1, parallel=parallel) # data_start cannot be negative


@pytest.mark.parametrize("parallel", [True, False])
def test_rdb(parallel, read_rdb):
    """
    Make sure the FastRdb reader works as expected.
    """
    text = """

A\tB\tC
1n\tS\t4N
1\t 9\t4.3
"""
    table = read_rdb(text, parallel=parallel)
    expected = Table([[1], [' 9'], [4.3]], names=('A', 'B', 'C'))
    assert_table_equal(table, expected)
    assert_equal(table['A'].dtype.kind, 'i')
    assert table['B'].dtype.kind in ('S', 'U')
    assert_equal(table['C'].dtype.kind, 'f')

    with pytest.raises(ValueError) as e:
        text = 'A\tB\tC\nN\tS\tN\n4\tb\ta' # C column contains non-numeric data
        read_rdb(text, parallel=parallel)
    assert 'Column C failed to convert' in str(e)

    with pytest.raises(ValueError) as e:
        text = 'A\tB\tC\nN\tN\n1\t2\t3' # not enough types specified
        read_rdb(text, parallel=parallel)
    assert 'mismatch between number of column names and column types' in str(e)

    with pytest.raises(ValueError) as e:
        text = 'A\tB\tC\nN\tN\t5\n1\t2\t3' # invalid type for column C
        read_rdb(text, parallel=parallel)
    assert 'type definitions do not all match [num](N|S)' in str(e)


@pytest.mark.parametrize("parallel", [True, False])
def test_data_start(parallel, read_basic):
    """
    Make sure that data parsing begins at data_start (ignoring empty and
    commented lines but not taking quoted values into account).
    """
    if parallel:
        pytest.xfail("Multiprocessing can fail with quoted fields")
    text = """
A B C
1 2 3
4 5 6

7 8 "9
 \t1"
# comment
10 11 12
"""
    table = read_basic(text, data_start=2, parallel=parallel)
    expected = Table([[4, 7, 10], [5, 8, 11], [6, 91, 12]], names=('A', 'B', 'C'))
    assert_table_equal(table, expected)

    table = read_basic(text, data_start=3, parallel=parallel)
    # ignore empty line
    expected = Table([[7, 10], [8, 11], [91, 12]], names=('A', 'B', 'C'))
    assert_table_equal(table, expected)

    with pytest.raises(CParserError) as e:
        # tries to begin in the middle of quoted field
        read_basic(text, data_start=4, parallel=parallel)
    assert 'not enough columns found in line 1 of data' in str(e)

    table = read_basic(text, data_start=5, parallel=parallel)
    # ignore commented line
    expected = Table([[10], [11], [12]], names=('A', 'B', 'C'))
    assert_table_equal(table, expected)

    text = """
A B C
1 2 3
4 5 6

7 8 9
# comment
10 11 12
"""
    # make sure reading works as expected in parallel
    table = read_basic(text, data_start=2, parallel=parallel)
    expected = Table([[4, 7, 10], [5, 8, 11], [6, 9, 12]], names=('A', 'B', 'C'))
    assert_table_equal(table, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_quoted_empty_values(parallel, read_basic):
    """
    Quoted empty values spanning multiple lines should be treated correctly.
    """
    if parallel:
        pytest.xfail("Multiprocessing can fail with quoted fields")
    text = 'a b c\n1 2 " \n "'
    table = read_basic(text, parallel=parallel)
    assert table['c'][0] is ma.masked # empty value masked by default


@pytest.mark.parametrize("parallel", [True, False])
def test_csv_comment_default(parallel, read_csv):
    """
    Unless the comment parameter is specified, the CSV reader should
    not treat any lines as comments.
    """
    text = 'a,b,c\n#1,2,3\n4,5,6'
    table = read_csv(text, parallel=parallel)
    expected = Table([['#1', '4'], [2, 5], [3, 6]], names=('a', 'b', 'c'))
    assert_table_equal(table, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_whitespace_before_comment(parallel, read_tab):
    """
    Readers that don't strip whitespace from data (Tab, RDB)
    should still treat lines with leading whitespace and then
    the comment char as comment lines.
    """
    text = 'a\tb\tc\n # comment line\n1\t2\t3'
    table = read_tab(text, parallel=parallel)
    expected = Table([[1], [2], [3]], names=('a', 'b', 'c'))
    assert_table_equal(table, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_strip_line_trailing_whitespace(parallel, read_basic):
    """
    Readers that strip whitespace from lines should ignore
    trailing whitespace after the last data value of each
    row.
    """
    text = 'a b c\n1 2 \n3 4 5'
    with pytest.raises(CParserError) as e:
        ascii.read(StringIO(text), format='fast_basic', guess=False)
    assert 'not enough columns found in line 1' in str(e)

    text = 'a b c\n 1 2 3   \t \n 4 5 6 '
    table = read_basic(text, parallel=parallel)
    expected = Table([[1, 4], [2, 5], [3, 6]], names=('a', 'b', 'c'))
    assert_table_equal(table, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_no_data(parallel, read_basic):
    """
    As long as column names are supplied, the C reader
    should return an empty table in the absence of data.
    """
    table = read_basic('a b c', parallel=parallel)
    expected = Table([[], [], []], names=('a', 'b', 'c'))
    assert_table_equal(table, expected)

    table = read_basic('a b c\n1 2 3', data_start=2, parallel=parallel)
    assert_table_equal(table, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_line_endings(parallel, read_basic, read_commented_header, read_rdb):
    """
    Make sure the fast reader accepts CR and CR+LF
    as newlines.
    """
    text = 'a b c\n1 2 3\n4 5 6\n7 8 9\n'
    expected = Table([[1, 4, 7], [2, 5, 8], [3, 6, 9]], names=('a', 'b', 'c'))

    for newline in ('\r\n', '\r'):
        table = read_basic(text.replace('\n', newline), parallel=parallel)
        assert_table_equal(table, expected)

    # Make sure the splitlines() method of FileString
    # works with CR/CR+LF line endings
    text = '#' + text
    for newline in ('\r\n', '\r'):
        table = read_commented_header(text.replace('\n', newline), parallel=parallel)
        assert_table_equal(table, expected)
    text = 'a\tb\tc\nN\tN\tN\n1\t2\t3\n4\t5\t6\n7\t8\t9\n'
    for newline in ('\r\n', '\r'):
        table = read_rdb(text.replace('\n', newline), parallel=parallel)
        assert_table_equal(table, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_store_comments(parallel, read_basic):
    """
    Make sure that the output Table produced by the fast
    reader stores any comment lines in its meta attribute.
    """
    text = """
# header comment
a b c
# comment 2
# comment 3
1 2 3
4 5 6
"""
    table = read_basic(text, parallel=parallel, check_meta=True)
    assert_equal(table.meta['comments'],
                 ['header comment', 'comment 2', 'comment 3'])

@pytest.mark.parametrize("parallel", [True, False])
def test_empty_quotes(parallel, read_basic):
    """
    Make sure the C reader doesn't segfault when the
    input data contains empty quotes. [#3407]
    """
    table = read_basic('a b\n1 ""\n2 ""', parallel=parallel)
    expected = Table([[1, 2], [0, 0]], names=('a', 'b'))
    assert_table_equal(table, expected)

@pytest.mark.parametrize("parallel", [True, False])
def test_fast_tab_with_names(parallel, read_tab):
    """
    Make sure the C reader doesn't segfault when the header for the
    first column is missing [#3545]
    """
    content = """#
\tdecDeg\tRate_pn_offAxis\tRate_mos2_offAxis\tObsID\tSourceID\tRADeg\tversion\tCounts_pn\tRate_pn\trun\tRate_mos1\tRate_mos2\tInserted_pn\tInserted_mos2\tbeta\tRate_mos1_offAxis\trcArcsec\tname\tInserted\tCounts_mos1\tInserted_mos1\tCounts_mos2\ty\tx\tCounts\toffAxis\tRot
-3.007559\t0.0000\t0.0010\t0013140201\t0\t213.462574\t0\t2\t0.0002\t0\t0.0001\t0.0001\t0\t1\t0.66\t0.0217\t3.0\tfakeXMMXCS J1413.8-0300\t3\t1\t2\t1\t398.000\t127.000\t5\t13.9\t72.3\t"""
    head = ['A{0}'.format(i) for i in range(28)]
    table = read_tab(content, data_start=1,
                     parallel=parallel, names=head)
