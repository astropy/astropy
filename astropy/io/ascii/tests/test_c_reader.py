# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import functools

from io import BytesIO
from textwrap import dedent

import pytest
import numpy as np
from numpy import ma

from ....table import Table, MaskedColumn
from ... import ascii
from ...ascii.core import ParameterError, FastOptionsError, InconsistentTableError
from ...ascii.cparser import CParserError
from ..fastbasic import (
    FastBasic, FastCsv, FastTab, FastCommentedHeader, FastRdb, FastNoHeader)
from .common import assert_equal, assert_almost_equal, assert_true


StringIO = lambda x: BytesIO(x.encode('ascii'))
TRAVIS = os.environ.get('TRAVIS', False)


def assert_table_equal(t1, t2, check_meta=False, rtol=1.e-15, atol=1.e-300):
    """
    Test equality of all columns in a table, with stricter tolerances for
    float columns than the np.allclose default.
    """
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
                    if not isinstance(el, str) and np.isnan(el):
                        assert_true(not isinstance(t2[name][i], str) and np.isnan(t2[name][i]))
                    elif isinstance(el, str):
                        assert_equal(el, t2[name][i])
                    else:
                        assert_almost_equal(el, t2[name][i], rtol=rtol, atol=atol)
                except (TypeError, NotImplementedError):
                    pass  # ignore for now


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
    # TODO: also read from file
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
    # Cannot set header_start=None for basic format
    with pytest.raises(ValueError):
        read_basic("A B C\n1 2 3\n4 5 6", header_start=None, data_start=0, parallel=parallel)

    t2 = read_no_header("A B C\n1 2 3\n4 5 6", parallel=parallel)
    expected = Table([['A', '1', '4'], ['B', '2', '5'], ['C', '3', '6']], names=('col1', 'col2', 'col3'))
    assert_table_equal(t2, expected)


@pytest.mark.parametrize("parallel", [True, False])
def test_no_header_supplied_names(parallel, read_basic, read_no_header):
    """
    If header_start=None and names is passed as a parameter, header
    data should not be read and names should be used instead.
    """
    table = read_no_header("A B C\n1 2 3\n4 5 6",
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


@pytest.mark.parametrize("key,val", [
    ('delimiter', ',,'),  # multi-char delimiter
    ('comment', '##'),  # multi-char comment
    ('data_start', None),  # data_start=None
    ('data_start', -1),  # data_start negative
    ('quotechar', '##'),  # multi-char quote signifier
    ('header_start', -1),  # negative header_start
    ('converters', dict((i + 1, ascii.convert_numpy(np.uint)) for i in range(3))),  # passing converters
    ('Inputter', ascii.ContinuationLinesInputter),  # passing Inputter
    ('header_Splitter', ascii.DefaultSplitter),  # passing Splitter
    ('data_Splitter', ascii.DefaultSplitter)])
def test_invalid_parameters(key, val):
    """
    Make sure the C reader raises an error if passed parameters it can't handle.
    """
    with pytest.raises(ParameterError):
        FastBasic(**{key: val}).read('1 2 3\n4 5 6')
    with pytest.raises(ParameterError):
        ascii.read('1 2 3\n4 5 6',
                   format='fast_basic', guess=False, **{key: val})


def test_invalid_parameters_other():
    with pytest.raises(TypeError):
        FastBasic(foo=7).read('1 2 3\n4 5 6')  # unexpected argument
    with pytest.raises(FastOptionsError):  # don't fall back on the slow reader
        ascii.read('1 2 3\n4 5 6', format='basic', fast_reader={'foo': 7})
    with pytest.raises(ParameterError):
        # Outputter cannot be specified in constructor
        FastBasic(Outputter=ascii.TableOutputter).read('1 2 3\n4 5 6')


def test_too_many_cols1():
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
    with pytest.raises(InconsistentTableError) as e:
        table = FastBasic().read(text)
    assert 'InconsistentTableError: Number of header columns (3) ' \
           'inconsistent with data columns in data line 2' in str(e)


def test_too_many_cols2():
    text = """\
aaa,bbb
1,2,
3,4,
"""
    with pytest.raises(InconsistentTableError) as e:
        table = FastCsv().read(text)
    assert 'InconsistentTableError: Number of header columns (2) ' \
           'inconsistent with data columns in data line 0' in str(e)


def test_too_many_cols3():
    text = """\
aaa,bbb
1,2,,
3,4,
"""
    with pytest.raises(InconsistentTableError) as e:
        table = FastCsv().read(text)
    assert 'InconsistentTableError: Number of header columns (2) ' \
           'inconsistent with data columns in data line 0' in str(e)


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

    with pytest.raises(InconsistentTableError) as e:
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
    assert table['A'][0] is not ma.masked  # empty value unaffected
    assert table['C'][2] is not ma.masked  # -9999 is not an exact match
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
    assert np.isnan(table['B'][3])  # nan filling skips column B
    assert table['B'][3] is not ma.masked  # should skip masking as well as replacing nan
    assert table['A'][0] is ma.masked
    assert table['A'][2] is ma.masked
    assert_equal(table['A'].data.data[0], '0')
    assert_equal(table['A'].data.data[2], '999')
    assert table['C'][0] is ma.masked
    assert_almost_equal(table['C'].data.data[0], 999.0)
    assert_almost_equal(table['C'][1], -3.4)  # column is still of type float


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
    assert table['C'][2] is not ma.masked  # C not in fill_include_names

    table = read_csv(text, fill_exclude_names=['A', 'B'], parallel=parallel)
    assert table['C'][2] is ma.masked
    assert table['A'][0] is not ma.masked
    assert table['B'][1] is not ma.masked  # A and B excluded from fill handling

    table = read_csv(text, fill_include_names=['A', 'B'], fill_exclude_names=['B'], parallel=parallel)
    assert table['A'][0] is ma.masked
    assert table['B'][1] is not ma.masked  # fill_exclude_names applies after fill_include_names
    assert table['C'][2] is not ma.masked


@pytest.mark.parametrize("parallel", [True, False])
def test_many_rows(parallel, read_basic):
    """
    Make sure memory reallocation works okay when the number of rows
    is large (so that each column string is longer than INITIAL_COL_SIZE).
    """
    text = 'A B C\n'
    for i in range(500):  # create 500 rows
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
    with pytest.raises(ParameterError):  # C reader can't handle regex comment
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
    t3 = read_commented_header(text, header_start=-1, data_start=0, parallel=parallel)  # negative indexing allowed
    assert_table_equal(t3, expected)

    text += '7 8 9'
    t4 = read_commented_header(text, header_start=2, data_start=2, parallel=parallel)
    expected = Table([[7], [8], [9]], names=('A', 'B', 'C'))
    assert_table_equal(t4, expected)

    with pytest.raises(ParameterError):
        read_commented_header(text, header_start=-1, data_start=-1, parallel=parallel)  # data_start cannot be negative


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
        text = 'A\tB\tC\nN\tS\tN\n4\tb\ta'  # C column contains non-numeric data
        read_rdb(text, parallel=parallel)
    assert 'Column C failed to convert' in str(e)

    with pytest.raises(ValueError) as e:
        text = 'A\tB\tC\nN\tN\n1\t2\t3'  # not enough types specified
        read_rdb(text, parallel=parallel)
    assert 'mismatch between number of column names and column types' in str(e)

    with pytest.raises(ValueError) as e:
        text = 'A\tB\tC\nN\tN\t5\n1\t2\t3'  # invalid type for column C
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

    with pytest.raises(InconsistentTableError) as e:
        # tries to begin in the middle of quoted field
        read_basic(text, data_start=4, parallel=parallel)
    assert 'header columns (3) inconsistent with data columns in data line 0' \
        in str(e)

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
    assert table['c'][0] is ma.masked  # empty value masked by default


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
    with pytest.raises(InconsistentTableError) as e:
        ascii.read(StringIO(text), format='fast_basic', guess=False)
    assert 'header columns (3) inconsistent with data columns in data line 0' \
        in str(e)

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

    expected = Table([[1, 4, 7], [2, 5, 8], [3, 6, 9]], names=('a', 'b', 'c'), masked=True)
    expected['a'][0] = np.ma.masked
    expected['c'][0] = np.ma.masked
    text = 'a\tb\tc\nN\tN\tN\n\t2\t\n4\t5\t6\n7\t8\t9\n'
    for newline in ('\r\n', '\r'):
        table = read_rdb(text.replace('\n', newline), parallel=parallel)
        assert_table_equal(table, expected)
        assert np.all(table == expected)


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


@pytest.mark.skipif(not os.getenv('TEST_READ_HUGE_FILE'),
                    reason='Environment variable TEST_READ_HUGE_FILE must be '
                    'defined to run this test')
def test_read_big_table(tmpdir):
    """Test reading of a huge file.

    This test generates a huge CSV file (~2.3Gb) before reading it (see
    https://github.com/astropy/astropy/pull/5319). The test is run only if the
    environment variable ``TEST_READ_HUGE_FILE`` is defined. Note that running
    the test requires quite a lot of memory (~18Gb when reading the file) !!

    """
    NB_ROWS = 250000
    NB_COLS = 500
    filename = str(tmpdir.join("big_table.csv"))

    print("Creating a {} rows table ({} columns).".format(NB_ROWS, NB_COLS))
    data = np.random.random(NB_ROWS)
    t = Table(data=[data]*NB_COLS, names=[str(i) for i in range(NB_COLS)])
    data = None

    print("Saving the table to {}".format(filename))
    t.write(filename, format='ascii.csv', overwrite=True)
    t = None

    print("Counting the number of lines in the csv, it should be {}"
          " + 1 (header).".format(NB_ROWS))
    assert sum(1 for line in open(filename)) == NB_ROWS + 1

    print("Reading the file with astropy.")
    t = Table.read(filename, format='ascii.csv', fast_reader=True)
    assert len(t) == NB_ROWS


@pytest.mark.skipif(not os.getenv('TEST_READ_HUGE_FILE'),
                    reason='Environment variable TEST_READ_HUGE_FILE must be '
                    'defined to run this test')
def test_read_big_table2(tmpdir):
    """Test reading of a file with a huge column.
    """
    # (2**32 // 2) : max value for int
    # // 10 : we use a value for rows that have 10 chars (1e9)
    # + 5 : add a few lines so the length cannot be stored by an int
    NB_ROWS = (2**32 // 2) // 10 + 5
    filename = str(tmpdir.join("big_table.csv"))

    print("Creating a {} rows table.".format(NB_ROWS))
    data = np.full(2**32 // 2 // 10 + 5, int(1e9), dtype=np.int32)
    t = Table(data=[data], names=['a'], copy=False)

    print("Saving the table to {}".format(filename))
    t.write(filename, format='ascii.csv', overwrite=True)
    t = None

    print("Counting the number of lines in the csv, it should be {}"
          " + 1 (header).".format(NB_ROWS))
    assert sum(1 for line in open(filename)) == NB_ROWS + 1

    print("Reading the file with astropy.")
    t = Table.read(filename, format='ascii.csv', fast_reader=True)
    assert len(t) == NB_ROWS


# Test these both with guessing turned on and off
@pytest.mark.parametrize("guess", [True, False])
# fast_reader configurations: False| 'use_fast_converter'=False|True
@pytest.mark.parametrize('fast_reader', [False, dict(use_fast_converter=False),
                                         dict(use_fast_converter=True)])
# Catch Windows environment since we cannot use _read() with custom fast_reader
@pytest.mark.parametrize("parallel", [False,
    pytest.param(True, marks=pytest.mark.xfail(os.name == 'nt', reason="Multiprocessing is currently unsupported on Windows"))])
def test_data_out_of_range(parallel, fast_reader, guess):
    """
    Numbers with exponents beyond float64 range (|~4.94e-324 to 1.7977e+308|)
    shall be returned as 0 and +-inf respectively by the C parser, just like
    the Python parser.
    Test fast converter only to nominal accuracy.
    """
    # Python reader and strtod() are expected to return precise results
    rtol = 1.e-30

    # Update fast_reader dict
    if fast_reader:
        fast_reader['parallel'] = parallel
        if fast_reader.get('use_fast_converter'):
            rtol = 1.e-15
        elif np.iinfo(np.int).dtype == np.dtype(np.int32):
            # On 32bit the standard C parser (strtod) returns strings for these
            pytest.xfail("C parser cannot handle float64 on 32bit systems")

    if parallel:
        if not fast_reader:
            pytest.skip("Multiprocessing only available in fast reader")
        elif TRAVIS:
            pytest.xfail("Multiprocessing can sometimes fail on Travis CI")

    fields = ['10.1E+199', '3.14e+313', '2048e+306', '0.6E-325', '-2.e345']
    values = np.array([1.01e200, np.inf, np.inf, 0.0, -np.inf])
    t = ascii.read(StringIO(' '.join(fields)), format='no_header',
                   guess=guess, fast_reader=fast_reader)
    read_values = np.array([col[0] for col in t.itercols()])
    assert_almost_equal(read_values, values, rtol=rtol, atol=1.e-324)

    # Test some additional corner cases
    fields = ['.0101E202', '0.000000314E+314', '1777E+305', '-1799E+305',
              '0.2e-323', '5200e-327', ' 0.0000000000000000000001024E+330']
    values = np.array([1.01e200, 3.14e307, 1.777e308, -np.inf, 0.0, 4.94e-324, 1.024e308])
    t = ascii.read(StringIO(' '.join(fields)), format='no_header',
                   guess=guess, fast_reader=fast_reader)
    read_values = np.array([col[0] for col in t.itercols()])
    assert_almost_equal(read_values, values, rtol=rtol, atol=1.e-324)

    # Test corner cases again with non-standard exponent_style (auto-detection)
    if fast_reader and fast_reader.get('use_fast_converter'):
        fast_reader.update({'exponent_style': 'A'})
    else:
        pytest.skip("Fortran exponent style only available in fast converter")

    fields = ['.0101D202', '0.000000314d+314', '1777+305', '-1799E+305',
              '0.2e-323', '2500-327', ' 0.0000000000000000000001024Q+330']
    t = ascii.read(StringIO(' '.join(fields)), format='no_header',
                   guess=guess, fast_reader=fast_reader)
    read_values = np.array([col[0] for col in t.itercols()])
    assert_almost_equal(read_values, values, rtol=rtol, atol=1.e-324)


@pytest.mark.parametrize("guess", [True, False])
# catch Windows environment since we cannot use _read() with custom fast_reader
@pytest.mark.parametrize("parallel", [
    pytest.param(True, marks=pytest.mark.xfail(os.name == 'nt', reason="Multiprocessing is currently unsupported on Windows")),
    False])

def test_int_out_of_range(parallel, guess):
    """
    Integer numbers outside int range shall be returned as string columns
    consistent with the standard (Python) parser (no 'upcasting' to float).
    """
    imin = np.iinfo(int).min+1
    imax = np.iinfo(int).max-1
    huge = '{:d}'.format(imax+2)

    text = 'P M S\n {:d} {:d} {:s}'.format(imax, imin, huge)
    expected = Table([[imax], [imin], [huge]], names=('P', 'M', 'S'))
    table = ascii.read(text, format='basic', guess=guess,
                       fast_reader={'parallel': parallel})
    assert_table_equal(table, expected)

    # check with leading zeroes to make sure strtol does not read them as octal
    text = 'P M S\n000{:d} -0{:d} 00{:s}'.format(imax, -imin, huge)
    expected = Table([[imax], [imin], ['00'+huge]], names=('P', 'M', 'S'))
    table = ascii.read(text, format='basic', guess=guess,
                       fast_reader={'parallel': parallel})
    assert_table_equal(table, expected)

    # Mixed columns should be returned as float, but if the out-of-range integer
    # shows up first, it will produce a string column - with both readers
    pytest.xfail("Integer fallback depends on order of rows")
    text = 'A B\n 12.3 {0:d}9\n {0:d}9 45.6e7'.format(imax)
    expected = Table([[12.3, 10.*imax], [10.*imax, 4.56e8]],
                     names=('A', 'B'))

    table = ascii.read(text, format='basic', guess=guess,
                       fast_reader={'parallel': parallel})
    assert_table_equal(table, expected)
    table = ascii.read(text, format='basic', guess=guess, fast_reader=False)
    assert_table_equal(table, expected)


@pytest.mark.parametrize("guess", [True, False])
@pytest.mark.parametrize("parallel", [
    pytest.param(True, marks=pytest.mark.xfail(os.name == 'nt', reason="Multiprocessing is currently unsupported on Windows")),
    False])

def test_fortran_reader(parallel, guess):
    """
    Make sure that ascii.read() can read Fortran-style exponential notation
    using the fast_reader.
    """

    # check for nominal np.float64 precision
    rtol = 1.e-15
    atol = 0.0
    text = 'A B C D\n100.01{:s}99       2.0  2.0{:s}-103 3\n' + \
           ' 4.2{:s}-1 5.0{:s}-1     0.6{:s}4 .017{:s}+309'
    expc = Table([[1.0001e101, 0.42], [2, 0.5], [2.e-103, 6.e3], [3, 1.7e307]],
                 names=('A', 'B', 'C', 'D'))

    expstyles = {'e': 6*('E'),
                 'D': ('D', 'd', 'd', 'D', 'd', 'D'),
                 'Q': 3*('q', 'Q'),
                  'Fortran': ('E', '0', 'D', 'Q', 'd', '0')}

    # C strtod (not-fast converter) can't handle Fortran exp
    with pytest.raises(FastOptionsError) as e:
        ascii.read(text.format(*(6*('D'))), format='basic', guess=guess,
                   fast_reader={'use_fast_converter': False,
                                'parallel': parallel, 'exponent_style': 'D'})
    assert 'fast_reader: exponent_style requires use_fast_converter' in str(e)

    # Enable multiprocessing and the fast converter iterate over
    # all style-exponent combinations, with auto-detection
    for s, c in expstyles.items():
        table = ascii.read(text.format(*c), guess=guess,
                           fast_reader={'parallel': parallel, 'exponent_style': s})
        assert_table_equal(table, expc, rtol=rtol, atol=atol)

    # Additional corner-case checks including triple-exponents without
    # any character and mixed whitespace separators
    text = 'A B\t\t C D\n1.0001+101 2.0+000\t 0.0002-099 3\n ' + \
           '0.42-000 \t 0.5 6.+003   0.000000000000000000000017+330'
    table = ascii.read(text, guess=guess,
                       fast_reader={'parallel': parallel, 'exponent_style': 'A'})
    assert_table_equal(table, expc, rtol=rtol, atol=atol)


@pytest.mark.parametrize("guess", [True, False])
@pytest.mark.parametrize("parallel", [
    pytest.param(True, marks=pytest.mark.xfail(os.name == 'nt', reason="Multiprocessing is currently unsupported on Windows")),
    False])
def test_fortran_invalid_exp(parallel, guess):
    """
    Test Fortran-style exponential notation in the fast_reader with invalid
    exponent-like patterns (no triple-digits) to make sure they are returned
    as strings instead, as with the standard C parser.
    """
    if parallel and TRAVIS:
        pytest.xfail("Multiprocessing can sometimes fail on Travis CI")

    formats = {'basic': ' ', 'tab': '\t', 'csv': ','}
    header = ['S1', 'F2', 'S2', 'F3', 'S3', 'F4', 'F5', 'S4', 'I1', 'F6', 'F7']
    # Tested entries and expected returns, first for auto-detect,
    # then for different specified exponents
    fields = ['1.0001+1', '.42d1', '2.3+10', '0.5', '3+1001', '3000.',
              '2', '4.56e-2.3', '8000', '4.2-022', '.00000145e314']
    vals_e = ['1.0001+1', '.42d1', '2.3+10',   0.5, '3+1001',  3.e3,
              2, '4.56e-2.3',    8000,  '4.2-022', 1.45e308]
    vals_d = ['1.0001+1',     4.2, '2.3+10',   0.5, '3+1001',  3.e3,
              2, '4.56e-2.3',    8000,  '4.2-022', '.00000145e314']
    vals_a = ['1.0001+1',     4.2, '2.3+10',   0.5, '3+1001',  3.e3,
              2, '4.56e-2.3',    8000,   4.2e-22,  1.45e308]
    vals_v = ['1.0001+1', 4.2, '2.3+10',   0.5, '3+1001',  3.e3,
               2, '4.56e-2.3',    8000,  '4.2-022', 1.45e308]

    # Iterate over supported format types and separators
    for f, s in formats.items():
        t1 = ascii.read(StringIO(s.join(header)+'\n'+s.join(fields)),
                        format=f, guess=guess,
                        fast_reader={'parallel': parallel, 'exponent_style': 'A'})
        assert_table_equal(t1, Table([[col] for col in vals_a], names=header))

    # Non-basic separators require guessing enabled to be detected
    if guess:
        formats['bar'] = '|'
    else:
        formats = {'basic': ' '}

    for s in formats.values():
        t2 = ascii.read(StringIO(s.join(header)+'\n'+s.join(fields)), guess=guess,
                fast_reader={'parallel': parallel, 'exponent_style': 'a'})

        assert_table_equal(t2, Table([[col] for col in vals_a], names=header))

    # Iterate for (default) expchar 'E'
    for s in formats.values():
        t3 = ascii.read(StringIO(s.join(header)+'\n'+s.join(fields)), guess=guess,
                fast_reader={'parallel': parallel, 'use_fast_converter': True})

        assert_table_equal(t3, Table([[col] for col in vals_e], names=header))

    # Iterate for expchar 'D'
    for s in formats.values():
        t4 = ascii.read(StringIO(s.join(header)+'\n'+s.join(fields)), guess=guess,
                fast_reader={'parallel': parallel, 'exponent_style': 'D'})

        assert_table_equal(t4, Table([[col] for col in vals_d], names=header))

    # Iterate for regular converter (strtod)
    for s in formats.values():
        t5 = ascii.read(StringIO(s.join(header)+'\n'+s.join(fields)), guess=guess,
                fast_reader={'parallel': parallel, 'use_fast_converter': False})

        read_values = [col[0] for col in t5.itercols()]
        if os.name == 'nt':
            # Apparently C strtod() on (some?) MSVC recognizes 'd' exponents!
            assert read_values == vals_v or read_values == vals_e
        else:
            assert read_values == vals_e


def test_fortran_reader_notbasic():
    """
    Check if readers without a fast option raise a value error when a
    fast_reader is asked for (implies the default 'guess=True').
    """

    tabstr = dedent("""
    a b
    1 1.23D4
    2 5.67D-8
    """)[1:-1]

    t1 = ascii.read(tabstr.split('\n'), fast_reader=dict(exponent_style='D'))

    assert t1['b'].dtype.kind == 'f'

    tabrdb = dedent("""
    a\tb
    # A simple RDB table
    N\tN
    1\t 1.23D4
    2\t 5.67-008
    """)[1:-1]

    t2 = ascii.read(tabrdb.split('\n'), format='rdb',
                    fast_reader=dict(exponent_style='fortran'))

    assert t2['b'].dtype.kind == 'f'

    tabrst = dedent("""
    = =======
    a b
    = =======
    1 1.23E4
    2 5.67E-8
    = =======
    """)[1:-1]

    t3 = ascii.read(tabrst.split('\n'), format='rst')

    assert t3['b'].dtype.kind == 'f'

    t4 = ascii.read(tabrst.split('\n'), guess=True)

    assert t4['b'].dtype.kind == 'f'

    # In the special case of fast_converter=True (the default),
    # incompatibility is ignored
    t5 = ascii.read(tabrst.split('\n'), format='rst', fast_reader=True)

    assert t5['b'].dtype.kind == 'f'

    with pytest.raises(ParameterError):
        t6 = ascii.read(tabrst.split('\n'), format='rst', guess=False,
                        fast_reader='force')

    with pytest.raises(ParameterError):
        t7 = ascii.read(tabrst.split('\n'), format='rst', guess=False,
                        fast_reader=dict(use_fast_converter=False))

    tabrst = tabrst.replace('E', 'D')

    with pytest.raises(ParameterError):
        t8 = ascii.read(tabrst.split('\n'), format='rst', guess=False,
                        fast_reader=dict(exponent_style='D'))


@pytest.mark.parametrize("guess", [True, False])
@pytest.mark.parametrize('fast_reader', [dict(exponent_style='D'),
                                         dict(exponent_style='A')])

def test_dict_kwarg_integrity(fast_reader, guess):
    """
    Check if dictionaries passed as kwargs (fast_reader in this test) are
    left intact by ascii.read()
    """
    expstyle = fast_reader.get('exponent_style', 'E')
    fields = ['10.1D+199', '3.14d+313', '2048d+306', '0.6D-325', '-2.d345']

    t = ascii.read(StringIO(' '.join(fields)), guess=guess,
                   fast_reader=fast_reader)
    assert fast_reader.get('exponent_style', None) == expstyle
