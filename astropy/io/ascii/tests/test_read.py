# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

import os
import re

import numpy as np

from ....utils import OrderedDict
from ....tests.helper import pytest
from ... import ascii
from ....table import Table
from ....units import Unit

from .common import (raises, assert_equal, assert_almost_equal,
                     assert_true, setup_function, teardown_function)
from .. import core



@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_convert_overflow(fast_reader):
    """
    Test reading an extremely large integer, which falls through to
    string due to an overflow error (#2234).
    """
    expected_kind = ('S', 'U')
    dat = ascii.read(['a', '1' * 10000], format='basic',
                     fast_reader=fast_reader, guess=False)
    assert dat['a'].dtype.kind in expected_kind


def test_guess_with_names_arg():
    """
    Make sure reading a table with guess=True gives the expected result when
    the names arg is specified.
    """
    # This is a NoHeader format table and so `names` should replace
    # the default col0, col1 names.  It fails as a Basic format
    # table when guessing because the column names would be '1', '2'.
    dat = ascii.read(['1,2', '3,4'], names=('a', 'b'))
    assert len(dat) == 2
    assert dat.colnames == ['a', 'b']

    # This is a Basic format table and the first row
    # gives the column names 'c', 'd', which get replaced by 'a', 'b'
    dat = ascii.read(['c,d', '3,4'], names=('a', 'b'))
    assert len(dat) == 1
    assert dat.colnames == ['a', 'b']

    # This is also a Basic format table and the first row
    # gives the column names 'c', 'd', which get replaced by 'a', 'b'
    dat = ascii.read(['c d', 'e f'], names=('a', 'b'))
    assert len(dat) == 1
    assert dat.colnames == ['a', 'b']


def test_guess_with_format_arg():
    """
    When the format or Reader is explicitly given then disable the
    strict column name checking in guessing.
    """
    dat = ascii.read(['1,2', '3,4'], format='basic')
    assert len(dat) == 1
    assert dat.colnames == ['1', '2']

    dat = ascii.read(['1,2', '3,4'], names=('a', 'b'), format='basic')
    assert len(dat) == 1
    assert dat.colnames == ['a', 'b']

    dat = ascii.read(['1,2', '3,4'], Reader=ascii.Basic)
    assert len(dat) == 1
    assert dat.colnames == ['1', '2']

    dat = ascii.read(['1,2', '3,4'], names=('a', 'b'), Reader=ascii.Basic)
    assert len(dat) == 1
    assert dat.colnames == ['a', 'b']

    # For good measure check the same in the unified I/O interface
    dat = Table.read(['1,2', '3,4'], format='ascii.basic')
    assert len(dat) == 1
    assert dat.colnames == ['1', '2']

    dat = Table.read(['1,2', '3,4'], format='ascii.basic', names=('a', 'b'))
    assert len(dat) == 1
    assert dat.colnames == ['a', 'b']


@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_read_with_names_arg(fast_reader):
    """
    Test that a bad value of `names` raises an exception.
    """
    with pytest.raises(ValueError):
        dat = ascii.read(['c d', 'e f'], names=('a', ), guess=False, fast_reader=fast_reader)


@pytest.mark.skipif(os.environ.get('APPVEYOR'), reason="fails on AppVeyor")
@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_read_all_files(fast_reader):
    for testfile in get_testfiles():
        if testfile.get('skip'):
            print('\n\n******** SKIPPING %s' % testfile['name'])
            continue
        print('\n\n******** READING %s' % testfile['name'])
        for guess in (True, False):
            test_opts = testfile['opts'].copy()
            if 'guess' not in test_opts:
                test_opts['guess'] = guess
            if 'Reader' in test_opts and 'fast_{0}'.format(test_opts['Reader']._format_name) \
                in core.FAST_CLASSES: # has fast version
                if 'Inputter' not in test_opts: # fast reader doesn't allow this
                    test_opts['fast_reader'] = fast_reader
            table = ascii.read(testfile['name'], **test_opts)
            assert_equal(table.dtype.names, testfile['cols'])
            for colname in table.dtype.names:
                assert_equal(len(table[colname]), testfile['nrows'])


@pytest.mark.skipif(os.environ.get('APPVEYOR'), reason="fails on AppVeyor")
@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_read_all_files_via_table(fast_reader):
    for testfile in get_testfiles():
        if testfile.get('skip'):
            print('\n\n******** SKIPPING %s' % testfile['name'])
            continue
        print('\n\n******** READING %s' % testfile['name'])
        for guess in (True, False):
            test_opts = testfile['opts'].copy()
            if 'guess' not in test_opts:
                test_opts['guess'] = guess
            if 'Reader' in test_opts:
                format = 'ascii.{0}'.format(test_opts['Reader']._format_name)
                del test_opts['Reader']
            else:
                format = 'ascii'
            if 'fast_{0}'.format(format) in core.FAST_CLASSES:
                test_opts['fast_reader'] = fast_reader
            table = Table.read(testfile['name'], format=format, **test_opts)
            assert_equal(table.dtype.names, testfile['cols'])
            for colname in table.dtype.names:
                assert_equal(len(table[colname]), testfile['nrows'])


def test_guess_all_files():
    for testfile in get_testfiles():
        if testfile.get('skip'):
            print('\n\n******** SKIPPING %s' % testfile['name'])
            continue
        if not testfile['opts'].get('guess', True):
            continue
        print('\n\n******** READING %s' % testfile['name'])
        for filter_read_opts in (['Reader', 'delimiter', 'quotechar'], []):
            # Copy read options except for those in filter_read_opts
            guess_opts = dict((k, v) for k, v in testfile['opts'].items()
                              if k not in filter_read_opts)
            table = ascii.read(testfile['name'], guess=True, **guess_opts)
            assert_equal(table.dtype.names, testfile['cols'])
            for colname in table.dtype.names:
                assert_equal(len(table[colname]), testfile['nrows'])


def test_daophot_indef():
    """Test that INDEF is correctly interpreted as a missing value"""
    table = ascii.read('t/daophot2.dat', Reader=ascii.Daophot)
    for colname in table.colnames:
        # Three columns have all INDEF values and are masked
        mask_value = colname in ('OTIME', 'MAG', 'MERR', 'XAIRMASS')
        assert np.all(table[colname].mask == mask_value)


def test_daophot_types():
    """
    Test specific data types which are different from what would be
    inferred automatically based only data values.  DAOphot reader uses
    the header information to assign types.
    """
    table = ascii.read('t/daophot2.dat', Reader=ascii.Daophot)
    assert table['LID'].dtype.char in 'fd'  # float or double
    assert table['MAG'].dtype.char in 'fd'  # even without any data values
    assert table['PIER'].dtype.char in 'US'  # string (data values are consistent with int)
    assert table['ID'].dtype.char in 'il'  # int or long


def test_daophot_header_keywords():
    table = ascii.read('t/daophot.dat', Reader=ascii.Daophot)
    expected_keywords = (('NSTARFILE', 'test.nst.1', 'filename', '%-23s'),
                         ('REJFILE', '"hello world"', 'filename', '%-23s'),
                         ('SCALE', '1.',  'units/pix', '%-23.7g'),)

    keywords = table.meta['keywords']  # Ordered dict of keyword structures
    for name, value, units, format_ in expected_keywords:
        keyword = keywords[name]
        assert_equal(keyword['value'], value)
        assert_equal(keyword['units'], units)
        assert_equal(keyword['format'], format_)


def test_daophot_multiple_aperture():
    table = ascii.read('t/daophot3.dat', Reader=ascii.Daophot)
    assert 'MAG5' in table.colnames  # MAG5 is one of the newly created column name
    assert table['MAG5'][4] == 22.13  # A sample entry in daophot3.dat file
    assert table['MERR2'][0] == 1.171
    assert np.all(table['RAPERT5'] == 23.3)  # assert all the 5th apertures are same 23.3


@pytest.mark.skipif(os.environ.get('APPVEYOR'), reason="fails on AppVeyor")
@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_empty_table_no_header(fast_reader):
    with pytest.raises(ascii.InconsistentTableError):
        table = ascii.read('t/no_data_without_header.dat', Reader=ascii.NoHeader,
                            guess=False, fast_reader=fast_reader)


@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_wrong_quote(fast_reader):
    with pytest.raises(ascii.InconsistentTableError):
        table = ascii.read('t/simple.txt', guess=False, fast_reader=fast_reader)


@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_extra_data_col(fast_reader):
    with pytest.raises(ascii.InconsistentTableError):
        table = ascii.read('t/bad.txt', fast_reader=fast_reader)


@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_extra_data_col2(fast_reader):
    with pytest.raises(ascii.InconsistentTableError):
        table = ascii.read('t/simple5.txt', delimiter='|', fast_reader=fast_reader)


@raises(IOError)
def test_missing_file():
    table = ascii.read('does_not_exist')


@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_set_names(fast_reader):
    names = ('c1', 'c2', 'c3', 'c4', 'c5', 'c6')
    data = ascii.read('t/simple3.txt', names=names, delimiter='|',
                      fast_reader=fast_reader)
    assert_equal(data.dtype.names, names)


@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_set_include_names(fast_reader):
    names = ('c1', 'c2', 'c3', 'c4', 'c5', 'c6')
    include_names = ('c1', 'c3')
    data = ascii.read('t/simple3.txt', names=names, include_names=include_names,
                           delimiter='|', fast_reader=fast_reader)
    assert_equal(data.dtype.names, include_names)


@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_set_exclude_names(fast_reader):
    exclude_names = ('Y', 'object')
    data = ascii.read('t/simple3.txt', exclude_names=exclude_names, delimiter='|',
                      fast_reader=fast_reader)
    assert_equal(data.dtype.names, ('obsid', 'redshift', 'X', 'rad'))


def test_include_names_daophot():
    include_names = ('ID', 'MAG', 'PIER')
    data = ascii.read('t/daophot.dat', include_names=include_names)
    assert_equal(data.dtype.names, include_names)


def test_exclude_names_daophot():
    exclude_names = ('ID', 'YCENTER', 'MERR', 'NITER', 'CHI', 'PERROR')
    data = ascii.read('t/daophot.dat', exclude_names=exclude_names)
    assert_equal(data.dtype.names, ('XCENTER', 'MAG', 'MSKY', 'SHARPNESS', 'PIER'))


def test_custom_process_lines():
    def process_lines(lines):
        bars_at_ends = re.compile(r'^\| | \|$', re.VERBOSE)
        striplines = (x.strip() for x in lines)
        return [bars_at_ends.sub('', x) for x in striplines if len(x) > 0]
    reader = ascii.get_reader(delimiter='|')
    reader.inputter.process_lines = process_lines
    data = reader.read('t/bars_at_ends.txt')
    assert_equal(data.dtype.names, ('obsid', 'redshift', 'X', 'Y', 'object', 'rad'))
    assert_equal(len(data), 3)


def test_custom_process_line():
    def process_line(line):
        line_out = re.sub(r'^\|\s*', '', line.strip())
        return line_out
    reader = ascii.get_reader(data_start=2, delimiter='|')
    reader.header.splitter.process_line = process_line
    reader.data.splitter.process_line = process_line
    data = reader.read('t/nls1_stackinfo.dbout')
    cols = get_testfiles('t/nls1_stackinfo.dbout')['cols']
    assert_equal(data.dtype.names, cols[1:])


def test_custom_splitters():
    reader = ascii.get_reader()
    reader.header.splitter = ascii.BaseSplitter()
    reader.data.splitter = ascii.BaseSplitter()
    f = 't/test4.dat'
    data = reader.read(f)
    testfile = get_testfiles(f)
    assert_equal(data.dtype.names, testfile['cols'])
    assert_equal(len(data), testfile['nrows'])
    assert_almost_equal(data.field('zabs1.nh')[2], 0.0839710433091)
    assert_almost_equal(data.field('p1.gamma')[2], 1.25997502704)
    assert_almost_equal(data.field('p1.ampl')[2], 0.000696444029148)
    assert_equal(data.field('statname')[2], 'chi2modvar')
    assert_almost_equal(data.field('statval')[2], 497.56468441)


def test_start_end():
    data = ascii.read('t/test5.dat', header_start=1, data_start=3, data_end=-5)
    assert_equal(len(data), 13)
    assert_equal(data.field('statname')[0], 'chi2xspecvar')
    assert_equal(data.field('statname')[-1], 'chi2gehrels')


def test_set_converters():
    converters = {'zabs1.nh': [ascii.convert_numpy('int32'),
                               ascii.convert_numpy('float32')],
                  'p1.gamma': [ascii.convert_numpy('str')]
                  }
    data = ascii.read('t/test4.dat', converters=converters)
    assert_equal(str(data['zabs1.nh'].dtype), 'float32')
    assert_equal(data['p1.gamma'][0], '1.26764544642')


@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_from_string(fast_reader):
    f = 't/simple.txt'
    with open(f) as fd:
        table = fd.read()
    testfile = get_testfiles(f)
    data = ascii.read(table, fast_reader=fast_reader, **testfile['opts'])
    assert_equal(data.dtype.names, testfile['cols'])
    assert_equal(len(data), testfile['nrows'])


@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_from_filelike(fast_reader):
    f = 't/simple.txt'
    testfile = get_testfiles(f)
    with open(f, 'rb') as fd:
        data = ascii.read(fd, fast_reader=fast_reader, **testfile['opts'])
    assert_equal(data.dtype.names, testfile['cols'])
    assert_equal(len(data), testfile['nrows'])


@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_from_lines(fast_reader):
    f = 't/simple.txt'
    with open(f) as fd:
        table = fd.readlines()
    testfile = get_testfiles(f)
    data = ascii.read(table, fast_reader=fast_reader, **testfile['opts'])
    assert_equal(data.dtype.names, testfile['cols'])
    assert_equal(len(data), testfile['nrows'])


def test_comment_lines():
    table = ascii.get_reader(Reader=ascii.Rdb)
    data = table.read('t/apostrophe.rdb')
    assert_equal(table.comment_lines, ['# first comment', '  # second comment'])
    assert_equal(data.meta['comments'], ['first comment', 'second comment'])


@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_fill_values(fast_reader):
    f = 't/fill_values.txt'
    testfile = get_testfiles(f)
    data = ascii.read(f, fill_values=('a', '1'), fast_reader=fast_reader,
                      **testfile['opts'])
    assert_true((data['a'].mask == [False, True]).all())
    assert_true((data['a'] == [1, 1]).all())
    assert_true((data['b'].mask == [False, True]).all())
    assert_true((data['b'] == [2, 1]).all())


@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_fill_values_col(fast_reader):
    f = 't/fill_values.txt'
    testfile = get_testfiles(f)
    data = ascii.read(f, fill_values=('a', '1', 'b'), fast_reader=fast_reader,
                      **testfile['opts'])
    check_fill_values(data)


@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_fill_values_include_names(fast_reader):
    f = 't/fill_values.txt'
    testfile = get_testfiles(f)
    data = ascii.read(f, fill_values=('a', '1'), fast_reader=fast_reader,
                           fill_include_names = ['b'], **testfile['opts'])
    check_fill_values(data)


@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_fill_values_exclude_names(fast_reader):
    f = 't/fill_values.txt'
    testfile = get_testfiles(f)
    data = ascii.read(f, fill_values=('a', '1'), fast_reader=fast_reader,
                           fill_exclude_names = ['a'], **testfile['opts'])
    check_fill_values(data)


def check_fill_values(data):
    """compare array column by column with expectation """
    assert_true((data['a'].mask == [False, False]).all())
    assert_true((data['a'] == ['1', 'a']).all())
    assert_true((data['b'].mask == [False, True]).all())
    # Check that masked value is "do not care" in comparison
    assert_true((data['b'] == [2, -999]).all())
    data['b'].mask = False  # explicitly unmask for comparison
    assert_true((data['b'] == [2, 1]).all())


@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_fill_values_list(fast_reader):
    f = 't/fill_values.txt'
    testfile = get_testfiles(f)
    data = ascii.read(f, fill_values=[('a', '42'), ('1', '42', 'a')],
                      fast_reader=fast_reader, **testfile['opts'])
    data['a'].mask = False  # explicitly unmask for comparison
    assert_true((data['a'] == [42, 42]).all())


def test_masking_Cds():
    f = 't/cds.dat'
    testfile = get_testfiles(f)
    data = ascii.read(f,
                           **testfile['opts'])
    assert_true(data['AK'].mask[0])
    assert_true(not data['Fit'].mask[0])


def test_null_Ipac():
    f = 't/ipac.dat'
    testfile = get_testfiles(f)
    data = ascii.read(f, **testfile['opts'])
    mask = np.array([(True, False, True, False, True),
                     (False, False, False, False, False)],
                    dtype=[(str('ra'), '|b1'),
                           (str('dec'), '|b1'),
                           (str('sai'), '|b1'),
                           (str('v2'), '|b1'),
                           (str('sptype'), '|b1')])
    assert np.all(data.mask == mask)


def test_Ipac_meta():
    keywords = OrderedDict((('intval', 1),
                            ('floatval', 2.3e3),
                            ('date', "Wed Sp 20 09:48:36 1995"),
                            ('key_continue', 'IPAC keywords can continue across lines')))
    comments = ['This is an example of a valid comment']
    f = 't/ipac.dat'
    testfile = get_testfiles(f)
    data = ascii.read(f, **testfile['opts'])
    assert data.meta['keywords'].keys() == keywords.keys()
    for data_kv, kv in zip(data.meta['keywords'].values(), keywords.values()):
        assert data_kv['value'] == kv
    assert data.meta['comments'] == comments


def test_set_guess_kwarg():
    """Read a file using guess with one of the typical guess_kwargs explicitly set."""
    data = ascii.read('t/space_delim_no_header.dat',
                           delimiter=',', guess=True)
    assert(data.dtype.names == ('1 3.4 hello',))
    assert(len(data) == 1)


@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_read_rdb_wrong_type(fast_reader):
    """Read RDB data with inconstent data type (except failure)"""
    table = """col1\tcol2
N\tN
1\tHello"""
    err_type = ValueError if not fast_reader else ascii.InconsistentTableError
    with pytest.raises(err_type):
        ascii.read(table, Reader=ascii.Rdb, fast_reader=fast_reader)


@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_default_missing(fast_reader):
    """Read a table with empty values and ensure that corresponding entries are masked"""
    table = '\n'.join(['a,b,c,d',
                       '1,3,,',
                       '2, , 4.0 , ss '])
    dat = ascii.read(table, fast_reader=fast_reader)
    assert dat.masked is True
    assert dat.pformat() == [' a   b   c   d ',
                             '--- --- --- ---',
                             '  1   3  --  --',
                             '  2  -- 4.0  ss']

    # Single row table with a single missing element
    table = """ a \n "" """
    dat = ascii.read(table, fast_reader=fast_reader)
    assert dat.pformat() == [' a ',
                             '---',
                             ' --']
    assert dat['a'].dtype.kind == 'i'

    # Same test with a fixed width reader
    table = '\n'.join([' a   b   c   d ',
                       '--- --- --- ---',
                       '  1   3        ',
                       '  2     4.0  ss'])
    dat = ascii.read(table, Reader=ascii.FixedWidthTwoLine)
    assert dat.masked is True
    assert dat.pformat() == [' a   b   c   d ',
                             '--- --- --- ---',
                             '  1   3  --  --',
                             '  2  -- 4.0  ss']

    dat = ascii.read(table, Reader=ascii.FixedWidthTwoLine, fill_values=None)
    assert dat.masked is False
    assert dat.pformat() == [' a   b   c   d ',
                             '--- --- --- ---',
                             '  1   3        ',
                             '  2     4.0  ss']

    dat = ascii.read(table, Reader=ascii.FixedWidthTwoLine, fill_values=[])
    assert dat.masked is False
    assert dat.pformat() == [' a   b   c   d ',
                             '--- --- --- ---',
                             '  1   3        ',
                             '  2     4.0  ss']


def get_testfiles(name=None):
    """Set up information about the columns, number of rows, and reader params to
    read a bunch of test files and verify columns and number of rows."""

    testfiles = [
        {'cols': ('agasc_id', 'n_noids', 'n_obs'),
         'name': 't/apostrophe.rdb',
         'nrows': 2,
         'opts': {'Reader': ascii.Rdb}},
        {'cols': ('agasc_id', 'n_noids', 'n_obs'),
         'name': 't/apostrophe.tab',
         'nrows': 2,
         'opts': {'Reader': ascii.Tab}},
        {'cols': ('Index',
                  'RAh',
                  'RAm',
                  'RAs',
                  'DE-',
                  'DEd',
                  'DEm',
                  'DEs',
                  'Match',
                  'Class',
                  'AK',
                  'Fit'),
         'name': 't/cds.dat',
         'nrows': 1,
         'opts': {'Reader': ascii.Cds}},
        # Test malformed CDS file (issues #2241 #467)
        {'cols': ('Index',
                  'RAh',
                  'RAm',
                  'RAs',
                  'DE-',
                  'DEd',
                  'DEm',
                  'DEs',
                  'Match',
                  'Class',
                  'AK',
                  'Fit'),
         'name': 't/cds_malformed.dat',
         'nrows': 1,
         'opts': {'Reader': ascii.Cds, 'data_start': 'guess'}},
        {'cols': ('a', 'b', 'c'),
         'name': 't/commented_header.dat',
         'nrows': 2,
         'opts': {'Reader': ascii.CommentedHeader}},
        {'cols': ('a', 'b', 'c'),
         'name': 't/commented_header2.dat',
         'nrows': 2,
         'opts': {'Reader': ascii.CommentedHeader, 'header_start': -1}},
        {'cols': ('col1', 'col2', 'col3', 'col4', 'col5'),
         'name': 't/continuation.dat',
         'nrows': 2,
         'opts': {'Inputter': ascii.ContinuationLinesInputter,
                  'Reader': ascii.NoHeader}},
        {'cols': ('ID',
                  'XCENTER',
                  'YCENTER',
                  'MAG',
                  'MERR',
                  'MSKY',
                  'NITER',
                  'SHARPNESS',
                  'CHI',
                  'PIER',
                  'PERROR'),
         'name': 't/daophot.dat',
         'nrows': 2,
         'opts': {'Reader': ascii.Daophot}},
        {'cols': ('NUMBER',
                  'FLUX_ISO',
                  'FLUXERR_ISO',
                  'VALUES',
                  'VALUES_1',
                  'FLAG'),
         'name': 't/sextractor.dat',
         'nrows': 3,
         'opts': {'Reader': ascii.SExtractor}},
        {'cols': ('ra', 'dec', 'sai', 'v2', 'sptype'),
         'name': 't/ipac.dat',
         'nrows': 2,
         'opts': {'Reader': ascii.Ipac}},
        {'cols': ('col0',
                  'objID',
                  'osrcid',
                  'xsrcid',
                  'SpecObjID',
                  'ra',
                  'dec',
                  'obsid',
                  'ccdid',
                  'z',
                  'modelMag_i',
                  'modelMagErr_i',
                  'modelMag_r',
                  'modelMagErr_r',
                  'expo',
                  'theta',
                  'rad_ecf_39',
                  'detlim90',
                  'fBlim90'),
         'name': 't/nls1_stackinfo.dbout',
         'nrows': 58,
         'opts': {'data_start': 2, 'delimiter': '|', 'guess': False}},
        {'cols': ('Index',
                  'RAh',
                  'RAm',
                  'RAs',
                  'DE-',
                  'DEd',
                  'DEm',
                  'DEs',
                  'Match',
                  'Class',
                  'AK',
                  'Fit'),
         'name': 't/no_data_cds.dat',
         'nrows': 0,
         'opts': {'Reader': ascii.Cds}},
        {'cols': ('ID',
                  'XCENTER',
                  'YCENTER',
                  'MAG',
                  'MERR',
                  'MSKY',
                  'NITER',
                  'SHARPNESS',
                  'CHI',
                  'PIER',
                  'PERROR'),
         'name': 't/no_data_daophot.dat',
         'nrows': 0,
         'opts': {'Reader': ascii.Daophot}},
        {'cols': ('NUMBER',
                  'FLUX_ISO',
                  'FLUXERR_ISO',
                  'VALUES',
                  'VALUES_1',
                  'FLAG'),
         'name': 't/no_data_sextractor.dat',
         'nrows': 0,
         'opts': {'Reader': ascii.SExtractor}},
        {'cols': ('ra', 'dec', 'sai', 'v2', 'sptype'),
         'name': 't/no_data_ipac.dat',
         'nrows': 0,
         'opts': {'Reader': ascii.Ipac}},
        {'cols': ('ra', 'v2'),
         'name': 't/ipac.dat',
         'nrows': 2,
         'opts': {'Reader': ascii.Ipac, 'include_names': ['ra', 'v2']}},
        {'cols': ('a', 'b', 'c'),
         'name': 't/no_data_with_header.dat',
         'nrows': 0,
         'opts': {}},
        {'cols': ('agasc_id', 'n_noids', 'n_obs'),
         'name': 't/short.rdb',
         'nrows': 7,
         'opts': {'Reader': ascii.Rdb}},
        {'cols': ('agasc_id', 'n_noids', 'n_obs'),
         'name': 't/short.tab',
         'nrows': 7,
         'opts': {'Reader': ascii.Tab}},
        {'cols': ('test 1a', 'test2', 'test3', 'test4'),
         'name': 't/simple.txt',
         'nrows': 2,
         'opts': {'quotechar': "'"}},
        {'cols': ('top1', 'top2', 'top3', 'top4'),
         'name': 't/simple.txt',
         'nrows': 1,
         'opts': {'quotechar': "'", 'header_start': 1, 'data_start': 2}},
        {'cols': ('top1', 'top2', 'top3', 'top4'),
         'name': 't/simple.txt',
         'nrows': 1,
         'opts': {'quotechar': "'", 'header_start': 1}},
        {'cols': ('top1', 'top2', 'top3', 'top4'),
         'name': 't/simple.txt',
         'nrows': 2,
         'opts': {'quotechar': "'", 'header_start': 1, 'data_start': 1}},
        {'cols': ('obsid', 'redshift', 'X', 'Y', 'object', 'rad'),
         'name': 't/simple2.txt',
         'nrows': 3,
         'opts': {'delimiter': '|'}},
        {'cols': ('obsid', 'redshift', 'X', 'Y', 'object', 'rad'),
         'name': 't/simple3.txt',
         'nrows': 2,
         'opts': {'delimiter': '|'}},
        {'cols': ('col1', 'col2', 'col3', 'col4', 'col5', 'col6'),
         'name': 't/simple4.txt',
         'nrows': 3,
         'opts': {'Reader': ascii.NoHeader, 'delimiter': '|'}},
        {'cols': ('col1', 'col2', 'col3'),
         'name': 't/space_delim_no_header.dat',
         'nrows': 2,
         'opts': {'Reader': ascii.NoHeader}},
        {'cols': ('col1', 'col2', 'col3'),
         'name': 't/space_delim_no_header.dat',
         'nrows': 2,
         'opts': {'Reader': ascii.NoHeader, 'header_start': None}},
        {'cols': ('obsid', 'offset', 'x', 'y', 'name', 'oaa'),
         'name': 't/space_delim_blank_lines.txt',
         'nrows': 3,
         'opts': {}},
        {'cols': ('zabs1.nh', 'p1.gamma', 'p1.ampl', 'statname', 'statval'),
         'name': 't/test4.dat',
         'nrows': 9,
         'opts': {}},
        {'cols': ('a', 'b', 'c'),
         'name': 't/fill_values.txt',
         'nrows': 2,
         'opts': {'delimiter': ','}},
        {'name': 't/whitespace.dat',
         'cols': ('quoted colname with tab\tinside', 'col2', 'col3'),
         'nrows': 2,
         'opts': {'delimiter': '\s'}},
        {'name': 't/simple_csv.csv',
         'cols': ('a','b','c'),
         'nrows': 2,
         'opts': {'Reader': ascii.Csv}},
        {'name': 't/simple_csv_missing.csv',
         'cols': ('a','b','c'),
         'nrows': 2,
         'skip': True,
         'opts': {'Reader': ascii.Csv}},
        {'cols': ('cola', 'colb', 'colc'),
         'name': 't/latex1.tex',
         'nrows': 2,
         'opts': {'Reader': ascii.Latex}},
        {'cols': ('Facility', 'Id', 'exposure', 'date'),
         'name': 't/latex2.tex',
         'nrows': 3,
         'opts': {'Reader': ascii.AASTex}},
        {'cols': ('Col1', 'Col2', 'Col3', 'Col4'),
         'name': 't/fixed_width_2_line.txt',
         'nrows': 2,
         'opts': {'Reader': ascii.FixedWidthTwoLine}},
    ]

    try:
        import bs4
        testfiles.append({'cols': ('Column 1', 'Column 2', 'Column 3'),
                          'name': 't/html.html',
                          'nrows': 3,
                          'opts': {'Reader': ascii.HTML}})
    except ImportError:
        pass

    if name is not None:
        return [x for x in testfiles if x['name'] == name][0]
    else:
        return testfiles

def test_header_start_exception():
    '''Check certain Readers throw an exception if ``header_start`` is set

    For certain Readers it does not make sense to set the ``header_start``, they
    throw an exception if you try.
    This was implemented in response to issue #885.
    '''
    for readerclass in [ascii.NoHeader, ascii.SExtractor, ascii.Ipac,
                   ascii.BaseReader, ascii.FixedWidthNoHeader, ascii.Cds, ascii.Daophot]:
        with pytest.raises(ValueError):
            reader = ascii.core._get_reader(readerclass, header_start=5)


def test_csv_table_read():
    """
    Check for a regression introduced by #1935.  Pseudo-CSV file with
    commented header line.
    """
    lines = ['# a, b',
             '1, 2',
             '3, 4']
    t = ascii.read(lines)
    assert t.colnames == ['a', 'b']

@pytest.mark.parametrize('fast_reader', [True, False, 'force'])
def test_overlapping_names(fast_reader):
    """
    Check that the names argument list can overlap with the existing column names.
    This tests the issue in #1991.
    """
    t = ascii.read(['a b', '1 2'], names=['b', 'a'], fast_reader=fast_reader)
    assert t.colnames == ['b', 'a']

def test_sextractor_units():
    """
    Make sure that the SExtractor reader correctly inputs descriptions and units.
    """
    table = ascii.read('t/sextractor2.dat', Reader=ascii.SExtractor, guess=False)
    expected_units = [None, Unit('pix'), Unit('pix'), Unit('mag'),
                Unit('mag'), None, Unit('pix**2'), Unit('m**(-6)'),
                Unit('mag * arcsec**(-2)')]
    expected_descrs = ['Running object number',
                       'Windowed position estimate along x',
                       'Windowed position estimate along y',
                       'Kron-like elliptical aperture magnitude',
                       'RMS error for AUTO magnitude',
                       'Extraction flags',
                       None,
                       'Barycenter position along MAMA x axis',
                       'Peak surface brightness above background']
    for i, colname in enumerate(table.colnames):
        assert table[colname].unit == expected_units[i]
        assert table[colname].description == expected_descrs[i]

def test_list_with_newlines():
    """
    Check that lists of strings where some strings consist of just a newline
    ("\n") are parsed correctly.
    """
    t = ascii.read(["abc", "123\n", "456\n", "\n", "\n"])
    assert t.colnames == ['abc']
    assert len(t) == 2
    assert t[0][0] == 123
    assert t[1][0] == 456

def test_commented_csv():
    """
    Check that Csv reader does not have ignore lines with the # comment
    character which is defined for most Basic readers.
    """
    t = ascii.read(['#a,b', '1,2', '#3,4'], format='csv')
    assert t.colnames == ['#a', 'b']
    assert len(t) == 2
    assert t['#a'][1] == '#3'

def test_meta_comments():
    """
    Make sure that line comments are included in the ``meta`` attribute
    of the output Table.
    """
    t = ascii.read(['#comment1', '#   comment2 \t', 'a,b,c', '1,2,3'])
    assert t.colnames == ['a', 'b', 'c']
    assert t.meta['comments'] == ['comment1', 'comment2']

def test_guess_fail():
    """
    Check the error message when guess fails
    """
    with pytest.raises(ascii.InconsistentTableError) as err:
        ascii.read('asfdasdf\n1 2 3', format='basic')
    assert "** To figure out why the table did not read, use guess=False and" in str(err.value)

    # Test the case with guessing enabled but for a format that has no free params
    with pytest.raises(ValueError) as err:
        ascii.read('asfdasdf\n1 2 3', format='ipac')
    assert 'At least one header line beginning and ending with delimiter required' in str(err.value)

    # Test the case with guessing enabled but with all params specified
    with pytest.raises(ValueError) as err:
        ascii.read('asfdasdf\n1 2 3', format='basic', quotechar='"', delimiter=' ', fast_reader=False)
    assert 'Number of header columns (1) inconsistent with data columns (3)' in str(err.value)


def test_guessing_file_object():
    """
    Test guessing a file object.  Fixes #3013 and similar issue noted in #3019.
    """
    t = ascii.read(open('t/ipac.dat.bz2', 'rb'))
    assert t.colnames == ['ra','dec','sai','v2','sptype']


def test_pformat_roundtrip():
    """Check that the screen output of ``print tab`` can be read. See #3025."""
    """Read a table with empty values and ensure that corresponding entries are masked"""
    table = '\n'.join(['a,b,c,d',
                       '1,3,1.11,1',
                       '2, 2, 4.0 , ss '])
    dat = ascii.read(table)
    out = ascii.read(dat.pformat())
    assert len(dat) == len(out)
    assert dat.colnames == out.colnames
    for c in dat.colnames:
        assert np.all(dat[c] == out[c])


def test_ipac_abbrev():
    lines = ['| c1 | c2 | c3   |   c4 | c5| c6 | c7  | c8 | c9|c10|c11|c12|',
             '| r  | rE | rea  | real | D | do | dou | f  | i | l | da| c |',
             '  1    2    3       4     5   6    7     8    9   10  11  12 ']
    dat = ascii.read(lines, format='ipac')
    for name in dat.columns[0:8]:
        assert dat[name].dtype.kind == 'f'
    for name in dat.columns[8:10]:
        assert dat[name].dtype.kind == 'i'
    for name in dat.columns[10:12]:
        assert dat[name].dtype.kind in ('U', 'S')


def test_almost_but_not_quite_daophot():
    '''Regression test for #3319.
    This tables looks so close to a daophot table, that the daophot reader gets
    quite far before it fails with an AttributeError.

    Note that this table will actually be read as Commented Header table with
    the columns ['some', 'header', 'info'].
    '''
    lines = ["# some header info",
             "#F header info beginning with 'F'",
             "1 2 3",
             "4 5 6",
             "7 8 9"]
    dat = ascii.read(lines)
    assert len(dat) == 3
