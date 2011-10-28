import re
import glob
import math

try:
    from .. import ascii as asciitable
except ImportError:
    from .. import asciitable

if asciitable.has_numpy:
    import numpy as np

from .common import *

try:
    from math import isnan
except ImportError:
    try:
        from numpy import isnan
    except ImportError:
        print('Tests requiring isnan will fail')

@has_numpy_and_not_has_numpy
def test_read_all_files(numpy):
    for testfile in get_testfiles():
        print('\n\n******** READING %s' % testfile['name'])
        if testfile.get('requires_numpy') and not asciitable.has_numpy:
            return
        for guess in (True, False):
            test_opts = testfile['opts'].copy()
            if 'guess' not in test_opts:
                test_opts['guess'] = guess
            table = asciitable.read(testfile['name'], numpy=numpy, **testfile['opts'])
            assert_equal(table.dtype.names, testfile['cols'])
            for colname in table.dtype.names:
                assert_equal(len(table[colname]), testfile['nrows'])

@has_numpy_and_not_has_numpy
def test_guess_all_files(numpy):
    for testfile in get_testfiles():
        if not testfile['opts'].get('guess', True):
            continue
        print('\n\n******** READING %s' % testfile['name'])
        if testfile.get('requires_numpy') and not asciitable.has_numpy:
            return
        for filter_read_opts in (['Reader', 'delimiter', 'quotechar'], []):
            # Copy read options except for those in filter_read_opts
            guess_opts = dict((k, v) for k, v in testfile['opts'].items()
                              if k not in filter_read_opts)
            table = asciitable.read(testfile['name'], numpy=numpy, guess=True, **guess_opts)
            assert_equal(table.dtype.names, testfile['cols'])
            for colname in table.dtype.names:
                assert_equal(len(table[colname]), testfile['nrows'])

@has_numpy
def test_daophot_header_keywords(numpy):
    reader = asciitable.get_reader(Reader=asciitable.DaophotReader, numpy=numpy)
    table = reader.read('t/daophot.dat')
    expected_keywords = (('NSTARFILE', 'test.nst.1', 'filename', '%-23s'),
                         ('REJFILE', 'hello world', 'filename', '%-23s'),
                         ('SCALE', '1.',  'units/pix', '%-23.7g'),)

    for name, value, units, format_ in expected_keywords:
        for keyword in reader.keywords:
            if keyword.name == name:
                assert_equal(keyword.value, value)
                assert_equal(keyword.units, units)
                assert_equal(keyword.format, format_)
                break
        else:
            raise ValueError('Keyword not found')


@has_numpy_and_not_has_numpy
@raises(asciitable.InconsistentTableError)
def test_empty_table_no_header(numpy):
    table = asciitable.read('t/no_data_without_header.dat', Reader=asciitable.NoHeader,
                            numpy=numpy, guess=False)

@has_numpy_and_not_has_numpy
@raises(asciitable.InconsistentTableError)
def test_wrong_quote(numpy):
    table = asciitable.read('t/simple.txt', numpy=numpy, guess=False)

@has_numpy_and_not_has_numpy
@raises(asciitable.InconsistentTableError)
def test_extra_data_col(numpy):
    table = asciitable.read('t/bad.txt', numpy=numpy)

@has_numpy_and_not_has_numpy
@raises(asciitable.InconsistentTableError)
def test_extra_data_col2(numpy):
    table = asciitable.read('t/simple5.txt', delimiter='|', numpy=numpy)

@has_numpy_and_not_has_numpy
@raises(IOError)
def test_missing_file(numpy):
    table = asciitable.read('does_not_exist', numpy=numpy)

@has_numpy_and_not_has_numpy
def test_set_names(numpy):
    names = ('c1','c2','c3', 'c4', 'c5', 'c6')
    include_names = ('c1', 'c3')
    exclude_names = ('c4', 'c5', 'c6')
    data = asciitable.read('t/simple3.txt', names=names, delimiter='|', numpy=numpy)
    assert_equal(data.dtype.names, names)

@has_numpy_and_not_has_numpy
def test_set_include_names(numpy):
    names = ('c1','c2','c3', 'c4', 'c5', 'c6')
    include_names = ('c1', 'c3')
    data = asciitable.read('t/simple3.txt', names=names, include_names=include_names,
                           delimiter='|', numpy=numpy)
    assert_equal(data.dtype.names, include_names)

@has_numpy_and_not_has_numpy
def test_set_exclude_names(numpy):
    exclude_names = ('Y', 'object')
    data = asciitable.read('t/simple3.txt', exclude_names=exclude_names, delimiter='|', numpy=numpy)
    assert_equal(data.dtype.names, ('obsid', 'redshift', 'X', 'rad'))

@has_numpy_and_not_has_numpy
def test_custom_process_lines(numpy):
    def process_lines(lines):
        bars_at_ends = re.compile(r'^\| | \|$', re.VERBOSE)
        striplines = (x.strip() for x in lines)
        return [bars_at_ends.sub('', x) for x in striplines if len(x) > 0]
    reader = asciitable.get_reader(delimiter='|', numpy=numpy)
    reader.inputter.process_lines = process_lines
    data = reader.read('t/bars_at_ends.txt')
    assert_equal(data.dtype.names, ('obsid', 'redshift', 'X', 'Y', 'object', 'rad'))
    assert_equal(len(data), 3)

@has_numpy_and_not_has_numpy
def test_custom_process_line(numpy):
    def process_line(line):
        line_out = re.sub(r'^\|\s*', '', line.strip())
        return line_out
    reader = asciitable.get_reader(data_start=2, delimiter='|', numpy=numpy)
    reader.header.splitter.process_line = process_line
    reader.data.splitter.process_line = process_line
    data = reader.read('t/nls1_stackinfo.dbout')
    cols = get_testfiles('t/nls1_stackinfo.dbout')['cols']
    assert_equal(data.dtype.names, cols[1:])

@has_numpy_and_not_has_numpy
def test_custom_splitters(numpy):
    reader = asciitable.get_reader(numpy=numpy)
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
    
@has_numpy_and_not_has_numpy
def test_start_end(numpy):
    data = asciitable.read('t/test5.dat', header_start=1, data_start=3, data_end=-5, numpy=numpy)
    assert_equal(len(data), 13)
    assert_equal(data.field('statname')[0], 'chi2xspecvar')
    assert_equal(data.field('statname')[-1], 'chi2gehrels')

@has_numpy
def test_set_converters(numpy):
    converters = {'zabs1.nh': [asciitable.convert_numpy('int32'),
                               asciitable.convert_numpy('float32')],
                  'p1.gamma': [asciitable.convert_numpy('str')]
                  }
    data = asciitable.read('t/test4.dat', converters=converters, numpy=numpy)
    assert_equal(str(data['zabs1.nh'].dtype), 'float32')
    assert_equal(data['p1.gamma'][0], '1.26764544642')
    
@has_numpy_and_not_has_numpy
def test_from_string(numpy):
    f = 't/simple.txt'
    table = open(f).read()
    testfile = get_testfiles(f)
    data = asciitable.read(table, numpy=numpy, **testfile['opts'])
    assert_equal(data.dtype.names, testfile['cols'])
    assert_equal(len(data), testfile['nrows'])
    
@has_numpy_and_not_has_numpy
def test_from_filelike(numpy):
    f = 't/simple.txt'
    table = open(f)
    testfile = get_testfiles(f)
    data = asciitable.read(table, numpy=numpy, **testfile['opts'])
    assert_equal(data.dtype.names, testfile['cols'])
    assert_equal(len(data), testfile['nrows'])
    
@has_numpy_and_not_has_numpy
def test_from_lines(numpy):
    f = 't/simple.txt'
    table = open(f).readlines()
    testfile = get_testfiles(f)
    data = asciitable.read(table, numpy=numpy, **testfile['opts'])
    assert_equal(data.dtype.names, testfile['cols'])
    assert_equal(len(data), testfile['nrows'])
    
@has_numpy_and_not_has_numpy
def test_comment_lines(numpy):
    table = asciitable.get_reader(Reader=asciitable.RdbReader, numpy=numpy)
    data = table.read('t/apostrophe.rdb')
    assert_equal(table.comment_lines, ['# first comment', '  # second comment'])

@has_numpy_and_not_has_numpy
def test_fill_values(numpy):
    f = 't/fill_values.txt'
    testfile = get_testfiles(f)
    data = asciitable.read(f, numpy=numpy, fill_values=('a','1'), **testfile['opts'])
    if numpy:
        assert_true((data.mask['a']==[False,True]).all())
        assert_true((data.data['a']==[1,1]).all())
        assert_true((data.mask['b']==[False,True]).all())
        assert_true((data.data['b']==[2,1]).all())
        
    else:
        assert_equal(data['a'],[1,1])
        assert_equal(data['b'],[2,1])
        
@has_numpy_and_not_has_numpy
def test_fill_values_col(numpy):
    f = 't/fill_values.txt'
    testfile = get_testfiles(f)
    data = asciitable.read(f, numpy=numpy, fill_values=('a','1', 'b'), **testfile['opts'])
    check_fill_values(numpy, data)

@has_numpy_and_not_has_numpy
def test_fill_values_include_names(numpy):
    f = 't/fill_values.txt'
    testfile = get_testfiles(f)
    data = asciitable.read(f, numpy=numpy, fill_values=('a','1'),
                           fill_include_names = ['b'], **testfile['opts'])
    check_fill_values(numpy, data)
        
@has_numpy_and_not_has_numpy
def test_fill_values_exclude_names(numpy):
    f = 't/fill_values.txt'
    testfile = get_testfiles(f)
    data = asciitable.read(f, numpy=numpy, fill_values=('a','1'),
                           fill_exclude_names = ['a'], **testfile['opts'])
    check_fill_values(numpy, data)

def check_fill_values(numpy, data):
    """compare array column by column with expectation """
    if numpy:
        assert_true((data.mask['a']==[False,False]).all())
        assert_true((data.data['a']==['1','a']).all())
        assert_true((data.mask['b']==[False,True]).all())
        assert_true((data.data['b']==[2,1]).all())        
    else:
        assert_equal(data['a'],['1','a'])
        assert_equal(data['b'],[2,1])

@has_numpy_and_not_has_numpy
def test_fill_values_list(numpy):
    f = 't/fill_values.txt'
    testfile = get_testfiles(f)
    data = asciitable.read(f, numpy=numpy, fill_values=[('a','42'),('1','42','a')],
                           **testfile['opts'])
    if numpy:
        assert_true((data.data['a']==[42,42]).all())
    else:
        assert_equal(data['a'],[42,42])

@has_numpy_and_not_has_numpy
def test_masking_Cds(numpy):
    f = 't/cds.dat'
    testfile = get_testfiles(f)
    data = asciitable.read(f, numpy=numpy, 
                           **testfile['opts'])
    if numpy:
        assert_true(data['AK'].mask[0])
        assert_true(not data['Fit'].mask[0])
    else:
        assert_true(isnan(data['AK'][0]))
        assert_true(not isnan(data['Fit'][0]))

@has_numpy_and_not_has_numpy
def test_set_guess_kwarg(numpy):
    """Read a file using guess with one of the typical guess_kwargs explicitly set."""
    data = asciitable.read('t/space_delim_no_header.dat', numpy=numpy, 
                           delimiter=',', guess=True)
    assert(data.dtype.names == ('1 3.4 hello', ))
    assert(len(data) == 1)

@has_numpy_and_not_has_numpy
@raises(asciitable.InconsistentTableError)
def test_read_rdb_wrong_type(numpy):
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
         'opts': {'Reader': asciitable.Cds}},
        {'cols': ('a', 'b', 'c'),
         'name': 't/commented_header.dat',
         'nrows': 2,
         'opts': {'Reader': asciitable.CommentedHeader}},
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
        {'cols': ('ra', 'dec', 'sai', 'v2', 'sptype'),
         'name': 't/ipac.dat',
         'nrows': 2,
         'opts': {'Reader': asciitable.Ipac}},
        {'cols': ('',
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
        {'cols': ('ra', 'dec', 'sai', 'v2', 'sptype'),
         'name': 't/no_data_ipac.dat',
         'nrows': 0,
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
         'nrows': 1172,
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
    
