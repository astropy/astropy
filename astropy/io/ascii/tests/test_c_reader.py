# Licensed under a 3-clause BSD style license - see LICENSE.rst

from ....table import Table, MaskedColumn
from ... import ascii
from ...ascii.core import ParameterError
from ...ascii.cparser import CParserError
from ..fastbasic import FastBasic, FastCsv
from .common import assert_equal, assert_almost_equal, assert_true
from ....tests.helper import pytest
from cStringIO import StringIO
import numpy as np
from numpy import ma
from ....extern import six

def assert_table_equal(t1, t2):
	assert_equal(len(t1), len(t2))
	assert_equal(t1.colnames, t2.colnames)
	for name in t1.colnames:
		if not isinstance(t1[name], MaskedColumn):
			for i, el in enumerate(t1[name]):
				if not isinstance(el, six.string_types) and np.isnan(el):
					assert_true(not isinstance(t2[name][i], six.string_types) and np.isnan(t2[name][i]))
				elif isinstance(el, six.string_types):
					assert_equal(el, t2[name][i])
				else:
					assert_almost_equal(el, t2[name][i])

def read_basic(table, **kwargs):
	reader = FastBasic(**kwargs)
	t1 = reader.read(table)
	t2 = ascii.read(table, format='fast_basic', guess=False, **kwargs)
	t3 = ascii.read(table, format='basic', guess=False, **kwargs)
	assert_table_equal(t1, t2)
	assert_table_equal(t2, t3)
	return t1

def read_csv(table, **kwargs):
	reader = FastCsv(**kwargs)
	t1 = reader.read(table)
	t2 = ascii.read(table, format='fast_csv', guess=False, **kwargs)
	t3 = ascii.read(table, format='csv', guess=False, **kwargs)
	assert_table_equal(t1, t2)
	assert_table_equal(t2, t3)
	return t1

def test_simple_data():
	"""
	Make sure the fast reader works with basic input data.
	"""
	table = read_basic(StringIO("A B C\n1 2 3\n4 5 6"))
	expected = Table([[1, 4], [2, 5], [3, 6]], names=('A', 'B', 'C'))
	assert_table_equal(table, expected)

def test_supplied_names():
	"""
	If passed as a parameter, names should replace any
	column names found in the header.
	"""
	table = read_basic(StringIO("A B C\n1 2 3\n4 5 6"), names=('X', 'Y', 'Z'))
	expected = Table([[1, 4], [2, 5], [3, 6]], names=('X', 'Y', 'Z'))
	assert_table_equal(table, expected)

def test_no_header():
	"""
	The header should not be read when header_start=None. Unless names is
	passed, the column names should be auto-generated.
	"""
	table = read_basic(StringIO("A B C\n1 2 3\n4 5 6"), header_start=None, data_start=0)
	expected = Table([['A', '1', '4'], ['B', '2', '5'], ['C', '3', '6']], names=('col1', 'col2', 'col3'))
	assert_table_equal(table, expected)

def test_no_header_supplied_names():
	"""
	If header_start=None and names is passed as a parameter, header
	data should not be read and names should be used instead.
	"""
	table = read_basic(StringIO("A B C\n1 2 3\n4 5 6"), header_start=None, data_start=0,
					   names=('X', 'Y', 'Z'))
	expected = Table([['A', '1', '4'], ['B', '2', '5'], ['C', '3', '6']], names=('X', 'Y', 'Z'))
	assert_table_equal(table, expected)

def test_comment():
	"""
	Make sure that line comments are ignored by the C reader.
	"""
	table = read_basic(StringIO("# comment\nA B C\n# another comment\n1 2 3\n4 5 6"))
	expected = Table([[1, 4], [2, 5], [3, 6]], names=('A', 'B', 'C'))
	assert_table_equal(table, expected)

def test_empty_lines():
	"""
	Make sure that empty lines are ignored by the C reader.
	"""
	table = read_basic(StringIO("\n\nA B C\n1 2 3\n\n\n4 5 6\n\n\n\n"))
	expected = Table([[1, 4], [2, 5], [3, 6]], names=('A', 'B', 'C'))
	assert_table_equal(table, expected)

def test_lstrip_whitespace():
	"""
	Test to make sure the whitespace ignores whitespace at the beginning of fields.
	"""
	text = """
     1,  2,   \t3
 A,\t\t B,  C
  a, b,   c
  
   """
	table = read_basic(StringIO(text), delimiter=',')
	expected = Table([['A', 'a'], ['B', 'b'], ['C', 'c']], names=('1', '2', '3'))
	assert_table_equal(table, expected)

def test_conversion():
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
	table = read_basic(StringIO(text))
	assert_true(np.issubdtype(table['A'].dtype, np.float_))
	assert_true(np.issubdtype(table['B'].dtype, np.object_))
	print(table['C'].dtype) # hm?
	assert_true(np.issubdtype(table['C'].dtype, np.int_))
	assert_true(np.issubdtype(table['D'].dtype, np.float_))
	assert_true(np.issubdtype(table['E'].dtype, np.object_))

def test_delimiter():
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
		table = read_basic(StringIO(text.replace(' ', sep)), delimiter=sep)
		assert_table_equal(table, expected)

def test_include_names():
	"""
	If include_names is not None, the parser should read only those columns in include_names.
	"""
	table = read_basic(StringIO("A B C D\n1 2 3 4\n5 6 7 8"), include_names=['A', 'D'])
	expected = Table([[1, 5], [4, 8]], names=('A', 'D'))
	assert_table_equal(table, expected)

def test_exclude_names():
	"""
	If exclude_names is not None, the parser should exclude the columns in exclude_names.
	"""
	table = read_basic(StringIO("A B C D\n1 2 3 4\n5 6 7 8"), exclude_names=['A', 'D'])
	expected = Table([[2, 6], [3, 7]], names=('B', 'C'))
	assert_table_equal(table, expected)

def test_include_exclude_names():
	"""
	Make sure that include_names is applied before exclude_names if both are specified.
	"""
	text = """
A B C D E F G H
1 2 3 4 5 6 7 8
9 10 11 12 13 14 15 16
"""
	table = read_basic(StringIO(text), include_names=['A', 'B', 'D', 'F', 'H'],
					exclude_names=['B', 'F'])
	expected = Table([[1, 9], [4, 12], [8, 16]], names=('A', 'D', 'H'))
	assert_table_equal(table, expected)

def test_quoted_fields():
	"""
	The character quotechar (default '"') should denote the start of a field which can
	contain the field delimiter and newlines.
	"""
	text = """
"A B" C D
1.5 2.1 -37.1
a b "   c
 d"
"""
	table = read_basic(StringIO(text))
	expected = Table([['1.5', 'a'], ['2.1', 'b'], ['-37.1', 'cd']], names=('A B', 'C', 'D'))
	assert_table_equal(table, expected)
	table = read_basic(StringIO(text.replace('"', "'")), quotechar="'")
	assert_table_equal(table, expected)

def test_invalid_parameters():
	"""
	Make sure the C reader raises a ParameterError if passed parameters it can't handle.
	"""
	with pytest.raises(ParameterError):
		table = FastBasic(delimiter=',,').read(StringIO('1 2 3\n4 5 6'))
	with pytest.raises(ParameterError):
		table = FastBasic(comment='##').read(StringIO('1 2 3\n4 5 6'))
	with pytest.raises(ParameterError):
		table = FastBasic(data_start=None).read(StringIO('1 2 3\n4 5 6'))
	with pytest.raises(ParameterError):
		table = FastBasic(quotechar='""').read(StringIO('1 2 3\n4 5 6'))
	with pytest.raises(ParameterError):
		int_converter = ascii.convert_numpy(np.uint)
		converters = dict((i + 1, ascii.convert_numpy(np.uint)) for i in range(3))
		table = FastBasic(converters=converters).read(StringIO('1 2 3\n4 5 6'))
	with pytest.raises(ParameterError):
		table = FastBasic(Outputter=ascii.TableOutputter).read(StringIO('1 2 3\n4 5 6'))
	with pytest.raises(ParameterError):
		table = FastBasic(Inputter=ascii.ContinuationLinesInputter).read(StringIO('1 2 3\n4 5 6'))
	for arg in ('header_Splitter', 'data_Splitter'):
		with pytest.raises(ParameterError):
			table = FastBasic(**{arg: ascii.DefaultSplitter}).read(StringIO('1 2 3\n4 5 6'))
	# unexpected argument
	with pytest.raises(TypeError):
		table = FastBasic(foo=7).read(StringIO('1 2 3\n4 5 6'))

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
		table = FastBasic().read(StringIO(text))
	assert 'CParserError: an error occurred while tokenizing data: too many ' \
		'columns found in line 3 of data' in str(e)

def test_not_enough_cols():
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
	table = read_csv(StringIO(text))
	assert table['B'][1] is not ma.masked
	assert table['C'][1] is ma.masked

	with pytest.raises(CParserError) as e:
		table = FastBasic(delimiter=',').read(StringIO(text)) #TODO: maybe make error types the same as the old ones

def test_data_end():
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
	table = read_basic(StringIO(text), data_end=3)
	expected = Table([[1, 4], [2, 5], [3, 6]], names=('A', 'B', 'C'))
	assert_table_equal(table, expected)

	# data_end supports negative indexing
	table = read_basic(StringIO(text), data_end=-2)
	assert_table_equal(table, expected)

def test_fill_values():
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
	table = read_basic(StringIO(text), delimiter=',')
	# The empty value in row A should become a masked '0'
	assert isinstance(table['A'], MaskedColumn)
	assert table['A'][0] is ma.masked
	# '0' rather than 0 because there is a string in the column
	assert_equal(table['A'].data.data[0], '0')
	assert table['A'][1] is not ma.masked

	table = read_basic(StringIO(text), delimiter=',', fill_values=('-999', '0'))
	assert isinstance(table['B'], MaskedColumn)
	assert table['A'][0] is not ma.masked # empty value unaffected
	assert table['C'][2] is not ma.masked # -9999 is not an exact match
	assert table['B'][1] is ma.masked
	# Integer because the rest of the column contains integer data
	assert_equal(table['B'].data.data[1], 0)
	assert table['B'][0] is not ma.masked

	table = read_basic(StringIO(text), delimiter=',', fill_values=[])
	# None of the columns should be masked
	for name in 'ABC':
		assert not isinstance(table[name], MaskedColumn)
	
	table = read_basic(StringIO(text), delimiter=',', fill_values=[('', '0', 'A'),
								('nan', '999', 'A', 'C')])
	assert table['B'][3] is not ma.masked # nan filling skips column B
	assert table['A'][0] is ma.masked
	assert table['A'][2] is ma.masked
	assert_equal(table['A'].data.data[0], '0')
	assert_equal(table['A'].data.data[2], '999')
	assert table['C'][0] is ma.masked
	assert_almost_equal(table['C'].data.data[0], 999.0)
	assert_almost_equal(table['C'][1], -3.4) # column is still of type float

def test_fill_include_exclude_names():
	"""
	fill_include_names and fill_exclude_names should filter missing/empty value handling
	in the same way that include_names and exclude_names filter output columns.
	"""
	text = """
A, B, C
, 1, 2
3, , 4
5, 5,
""" #TODO: fix issue with empty final field (read_basic doesn't work here)
	table = read_csv(StringIO(text), fill_include_names=['A', 'B'])
	assert table['A'][0] is ma.masked
	assert table['B'][1] is ma.masked
	assert table['C'][2] is not ma.masked # C not in fill_include_names

	table = read_csv(StringIO(text), fill_exclude_names=['A', 'B'])
	assert table['C'][2] is ma.masked
	assert table['A'][0] is not ma.masked
	assert table['B'][1] is not ma.masked # A and B excluded from fill handling

	table = read_csv(StringIO(text), fill_include_names=['A', 'B'],
					   fill_exclude_names=['B'])
	assert table['A'][0] is ma.masked
	assert table['B'][1] is not ma.masked # fill_exclude_names applies after fill_include_names
	assert table['C'][2] is not ma.masked

def test_many_rows():
	"""
	Make sure memory reallocation works okay when the number of rows
	is large (so that each column string is longer than INITIAL_COL_SIZE).
	"""
	text = 'A B C\n'
	for i in range(500): # create 500 rows
		text += ' '.join([str(i) for i in range(3)])
		text += '\n'

	table = read_basic(StringIO(text))
	expected = Table([[0] * 500, [1] * 500, [2] * 500], names=('A', 'B', 'C'))
	assert_table_equal(table, expected)

def test_many_columns():
	"""
	Make sure memory reallocation works okay when the number of columns
	is large (so that each hedaer string is longer than INITIAL_HEADER_SIZE).
	"""
	# create a string with 500 columns and two data rows
	text = ' '.join([str(i) for i in range(500)])
	text += ('\n' + text + '\n' + text)
	table = read_basic(StringIO(text))
	expected = Table([[i, i] for i in range(500)], names=[str(i) for i in range(500)])
	assert_table_equal(table, expected)
