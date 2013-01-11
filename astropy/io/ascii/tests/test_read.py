# Licensed under a 3-clause BSD style license - see LICENSE.rst
import re
import glob
import math

import numpy as np

from ....tests.helper import pytest
from ... import ascii as asciitable

from .common import (raises, numpy_lt_1p5,
                     assert_equal, assert_almost_equal, assert_true,
                     setup_function, teardown_function, has_isnan)


def test_read_all_files():
    for testfile in get_testfiles():
        if testfile.get('skip'):
            print('\n\n******** SKIPPING %s' % testfile['name'])
            continue
        print('\n\n******** READING %s' % testfile['name'])
        for guess in (True, False):
            test_opts = testfile['opts'].copy()
            if 'guess' not in test_opts:
                test_opts['guess'] = guess
            table = asciitable.read(testfile['name'], **testfile['opts'])
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
            table = asciitable.read(testfile['name'], guess=True, **guess_opts)
            assert_equal(table.dtype.names, testfile['cols'])
            for colname in table.dtype.names:
                assert_equal(len(table[colname]), testfile['nrows'])


def test_daophot_header_keywords():
    table = asciitable.read('t/daophot.dat', Reader=asciitable.Daophot)
    expected_keywords = (('NSTARFILE', 'test.nst.1', 'filename', '%-23s'),
                         ('REJFILE', 'hello world', 'filename', '%-23s'),
                         ('SCALE', '1.',  'units/pix', '%-23.7g'),)

    keywords = table.meta['keywords']  # Ordered dict of keyword structures
    for name, value, units, format_ in expected_keywords:
        keyword = keywords[name]
        assert_equal(keyword['value'], value)
        assert_equal(keyword['units'], units)
        assert_equal(keyword['format'], format_)


@raises(asciitable.InconsistentTableError)
def test_empty_table_no_header():
    table = asciitable.read('t/no_data_without_header.dat', Reader=asciitable.NoHeader,
                            guess=False)


@raises(asciitable.InconsistentTableError)
def test_wrong_quote():
    table = asciitable.read('t/simple.txt', guess=False)


@raises(asciitable.InconsistentTableError)
def test_extra_data_col():
    table = asciitable.read('t/bad.txt')


@raises(asciitable.InconsistentTableError)
def test_extra_data_col2():
    table = asciitable.read('t/simple5.txt', delimiter='|')


@raises(IOError)
def test_missing_file():
    table = asciitable.read('does_not_exist')


def test_set_names():
    names = ('c1','c2','c3', 'c4', 'c5', 'c6')
    data = asciitable.read('t/simple3.txt', names=names, delimiter='|')
    assert_equal(data.dtype.names, names)


def test_set_include_names():
    names = ('c1','c2','c3', 'c4', 'c5', 'c6')
    include_names = ('c1', 'c3')
    data = asciitable.read('t/simple3.txt', names=names, include_names=include_names,
                           delimiter='|')
    assert_equal(data.dtype.names, include_names)


def test_set_exclude_names():
    exclude_names = ('Y', 'object')
    data = asciitable.read('t/simple3.txt', exclude_names=exclude_names, delimiter='|')
    assert_equal(data.dtype.names, ('obsid', 'redshift', 'X', 'rad'))


def test_include_names_daophot():
    include_names = ('ID', 'MAG', 'PIER')
    data = asciitable.read('t/daophot.dat', include_names=include_names)
    assert_equal(data.dtype.names, include_names)


def test_exclude_names_daophot():
    exclude_names = ('ID', 'YCENTER', 'MERR', 'NITER', 'CHI', 'PERROR')
    data = asciitable.read('t/daophot.dat', exclude_names=exclude_names)
    assert_equal(data.dtype.names, ('XCENTER', 'MAG', 'MSKY', 'SHARPNESS', 'PIER'))


def test_custom_process_lines():
    def process_lines(lines):
        bars_at_ends = re.compile(r'^\| | \|$', re.VERBOSE)
        striplines = (x.strip() for x in lines)
        return [bars_at_ends.sub('', x) for x in striplines if len(x) > 0]
    reader = asciitable.get_reader(delimiter='|')
    reader.inputter.process_lines = process_lines
    data = reader.read('t/bars_at_ends.txt')
    assert_equal(data.dtype.names, ('obsid', 'redshift', 'X', 'Y', 'object', 'rad'))
    assert_equal(len(data), 3)


@pytest.mark.xfail('numpy_lt_1p5')
def test_custom_process_line():
    def process_line(line):
        line_out = re.sub(r'^\|\s*', '', line.strip())
        return line_out
    reader = asciitable.get_reader(data_start=2, delimiter='|')
    reader.header.splitter.process_line = process_line
    reader.data.splitter.process_line = process_line
    data = reader.read('t/nls1_stackinfo.dbout')
    cols = get_testfiles('t/nls1_stackinfo.dbout')['cols']
    assert_equal(data.dtype.names, cols[1:])


def test_custom_splitters():
    reader = asciitable.get_reader()
    reader.header.splitter = asciitable.BaseSplitter()
    reader.data.splitter = asciitable.BaseSplitter()
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
    data = asciitable.read('t/test5.dat', header_start=1, data_start=3, data_end=-5)
    assert_equal(len(data), 13)
    assert_equal(data.field('statname')[0], 'chi2xspecvar')
    assert_equal(data.field('statname')[-1], 'chi2gehrels')


def test_set_converters():
    converters = {'zabs1.nh': [asciitable.convert_numpy('int32'),
                               asciitable.convert_numpy('float32')],
                  'p1.gamma': [asciitable.convert_numpy('str')]
                  }
    data = asciitable.read('t/test4.dat', converters=converters)
    assert_equal(str(data['zabs1.nh'].dtype), 'float32')
    assert_equal(data['p1.gamma'][0], '1.26764544642')


def test_from_string():
    f = 't/simple.txt'
    table = open(f).read()
    testfile = get_testfiles(f)
    data = asciitable.read(table, **testfile['opts'])
    assert_equal(data.dtype.names, testfile['cols'])
    assert_equal(len(data), testfile['nrows'])


def test_from_filelike():
    f = 't/simple.txt'
    table = open(f, 'rb')
    testfile = get_testfiles(f)
    data = asciitable.read(table, **testfile['opts'])
    assert_equal(data.dtype.names, testfile['cols'])
    assert_equal(len(data), testfile['nrows'])


def test_from_lines():
    f = 't/simple.txt'
    table = open(f).readlines()
    testfile = get_testfiles(f)
    data = asciitable.read(table, **testfile['opts'])
    assert_equal(data.dtype.names, testfile['cols'])
    assert_equal(len(data), testfile['nrows'])


def test_comment_lines():
    table = asciitable.get_reader(Reader=asciitable.Rdb)
    data = table.read('t/apostrophe.rdb')
    assert_equal(table.comment_lines, ['# first comment', '  # second comment'])


@pytest.mark.xfail('numpy_lt_1p5')
def test_fill_values():
    f = 't/fill_values.txt'
    testfile = get_testfiles(f)
    data = asciitable.read(f, fill_values=('a','1'), **testfile['opts'])
    assert_true((data['a'].mask==[False,True]).all())
    assert_true((data['a']==[1,1]).all())
    assert_true((data['b'].mask==[False,True]).all())
    assert_true((data['b']==[2,1]).all())


@pytest.mark.xfail('numpy_lt_1p5')
def test_fill_values_col():
    f = 't/fill_values.txt'
    testfile = get_testfiles(f)
    data = asciitable.read(f, fill_values=('a','1', 'b'), **testfile['opts'])
    check_fill_values(data)


@pytest.mark.xfail('numpy_lt_1p5')
def test_fill_values_include_names():
    f = 't/fill_values.txt'
    testfile = get_testfiles(f)
    data = asciitable.read(f, fill_values=('a','1'),
                           fill_include_names = ['b'], **testfile['opts'])
    check_fill_values(data)


@pytest.mark.xfail('numpy_lt_1p5')
def test_fill_values_exclude_names():
    f = 't/fill_values.txt'
    testfile = get_testfiles(f)
    data = asciitable.read(f, fill_values=('a','1'),
                           fill_exclude_names = ['a'], **testfile['opts'])
    check_fill_values(data)


@pytest.mark.xfail('numpy_lt_1p5')
def check_fill_values(data):
    """compare array column by column with expectation """
    assert_true((data['a'].mask==[False,False]).all())
    assert_true((data['a']==['1','a']).all())
    assert_true((data['b'].mask==[False,True]).all())
    # Check that masked value is "do not care" in comparison
    assert_true((data['b']==[2, -999]).all())
    data['b'].mask = False  # explicitly unmask for comparison
    assert_true((data['b']==[2,1]).all())


@pytest.mark.xfail('numpy_lt_1p5')
def test_fill_values_list():
    f = 't/fill_values.txt'
    testfile = get_testfiles(f)
    data = asciitable.read(f, fill_values=[('a','42'),('1','42','a')],
                           **testfile['opts'])
    data['a'].mask = False  # explicitly unmask for comparison
    assert_true((data['a']==[42,42]).all())


@pytest.mark.xfail('numpy_lt_1p5')
def test_masking_Cds():
    f = 't/cds.dat'
    testfile = get_testfiles(f)
    data = asciitable.read(f, 
                           **testfile['opts'])
    assert_true(data['AK'].mask[0])
    assert_true(not data['Fit'].mask[0])


@pytest.mark.xfail('numpy_lt_1p5')
def test_null_Ipac():
    f = 't/ipac.dat'
    testfile = get_testfiles(f)
    data = asciitable.read(f, **testfile['opts'])
    mask = np.array([(True, False, True, False, True),
                     (False, False, False, False, False)],
                    dtype=[('ra', '|b1'), ('dec', '|b1'), ('sai', '|b1'),
                           ('v2', '|b1'), ('sptype', '|b1')])
    assert np.all(data.mask == mask)


def test_set_guess_kwarg():
    """Read a file using guess with one of the typical guess_kwargs explicitly set."""
    data = asciitable.read('t/space_delim_no_header.dat',
                           delimiter=',', guess=True)
    assert(data.dtype.names == ('1 3.4 hello',))
    assert(len(data) == 1)


@raises(asciitable.InconsistentTableError)
def test_read_rdb_wrong_type():
    """Read RDB data with inconstent data type (except failure)"""
    table = """col1\tcol2
N\tN
1\tHello"""
    dat = asciitable.read(table, Reader=asciitable.Rdb)

def get_testfiles(name=None):
    """Set up information about the columns, number of rows, and reader params to
    read a bunch of test files and verify columns and number of rows."""

    testfiles = [
        {'cols': ('agasc_id', 'n_noids', 'n_obs'),
         'name': 't/apostrophe.rdb',
         'nrows': 2,
         'opts': {'Reader': asciitable.Rdb}},
        {'cols': ('agasc_id', 'n_noids', 'n_obs'),
         'name': 't/apostrophe.tab',
         'nrows': 2,
         'opts': {'Reader': asciitable.Tab}},
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
         'skip': numpy_lt_1p5,
         'opts': {'Reader': asciitable.Cds}},
        {'cols': ('a', 'b', 'c'),
         'name': 't/commented_header.dat',
         'nrows': 2,
         'opts': {'Reader': asciitable.CommentedHeader}},
        {'cols': ('a', 'b', 'c'),
         'name': 't/commented_header2.dat',
         'nrows': 2,
         'opts': {'Reader': asciitable.CommentedHeader, 'header_start': -1}},
        {'cols': ('col1', 'col2', 'col3', 'col4', 'col5'),
         'name': 't/continuation.dat',
         'nrows': 2,
         'opts': {'Inputter': asciitable.ContinuationLinesInputter,
                  'Reader': asciitable.NoHeader}},
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
         'requires_numpy': True,
         'opts': {'Reader': asciitable.Daophot}},
        {'cols': ('IMAGE',
                  'XINIT',
                  'YINIT',
                  'ID',
                  'COORDS',
                  'LID',
                  'XCENTER',
                  'YCENTER',
                  'XSHIFT',
                  'YSHIFT',
                  'XERR',
                  'YERR',
                  'CIER',
                  'CERROR',
                  'MSKY',
                  'STDEV',
                  'SSKEW',
                  'NSKY',
                  'NSREJ',
                  'SIER',
                  'SERROR',
                  'ITIME',
                  'XAIRMASS',
                  'IFILTER',
                  'OTIME',
                  'RAPERT',
                  'SUM',
                  'AREA',
                  'FLUX',
                  'MAG',
                  'MERR',
                  'PIER',
                  'PERROR'),
         'name': 't/daophot2.dat',
         'nrows': 2,
         'requires_numpy': True,
         'opts': {'Reader': asciitable.Daophot}},
        {'cols': ('NUMBER',
                  'FLUX_ISO',
                  'FLUXERR_ISO',
                  'VALUES',
                  'VALUES_1',
                  'FLAG'),
         'name': 't/sextractor.dat',
         'nrows': 3,
         'requires_numpy': True,
         'opts': {'Reader': asciitable.SExtractor}},
        {'cols': ('ra', 'dec', 'sai', 'v2', 'sptype'),
         'name': 't/ipac.dat',
         'nrows': 2,
         'skip': numpy_lt_1p5,
         'opts': {'Reader': asciitable.Ipac}},
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
         'skip': numpy_lt_1p5,
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
         'skip': numpy_lt_1p5,
         'nrows': 0,
         'opts': {'Reader': asciitable.Cds}},
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
         'requires_numpy': True,
         'opts': {'Reader': asciitable.Daophot}},
        {'cols': ('NUMBER',
                  'FLUX_ISO',
                  'FLUXERR_ISO',
                  'VALUES',
                  'VALUES_1',
                  'FLAG'),
         'name': 't/no_data_sextractor.dat',
         'nrows': 0,
         'requires_numpy': True,
         'opts': {'Reader': asciitable.SExtractor}},
        {'cols': ('ra', 'dec', 'sai', 'v2', 'sptype'),
         'name': 't/no_data_ipac.dat',
         'nrows': 0,
         'skip': numpy_lt_1p5,
         'opts': {'Reader': asciitable.Ipac}},
        {'cols': ('a', 'b', 'c'),
         'name': 't/no_data_with_header.dat',
         'nrows': 0,
         'opts': {}},
        {'cols': ('agasc_id', 'n_noids', 'n_obs'),
         'name': 't/short.rdb',
         'nrows': 7,
         'opts': {'Reader': asciitable.Rdb}},
        {'cols': ('agasc_id', 'n_noids', 'n_obs'),
         'name': 't/short.tab',
         'nrows': 7,
         'opts': {'Reader': asciitable.Tab}},
        {'cols': ('test 1a', 'test2', 'test3', 'test4'),
         'name': 't/simple.txt',
         'nrows': 2,
         'opts': {'quotechar': "'"}},
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
         'opts': {'Reader': asciitable.NoHeader, 'delimiter': '|'}},
        {'cols': ('col1', 'col2', 'col3'),
         'name': 't/space_delim_no_header.dat',
         'nrows': 2,
         'opts': {}},
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
        {'cols': ('cola', 'colb', 'colc'),
         'name': 't/latex1.tex',
         'nrows': 2,
         'opts': {'Reader': asciitable.Latex}},
        {'cols': ('Facility', 'Id', 'exposure', 'date'),
         'name': 't/latex2.tex',
         'nrows': 3,
         'opts': {'Reader': asciitable.AASTex}},
         ]

    if name is not None:
        return [x for x in testfiles if x['name'] == name][0]
    else:
        return testfiles
