# Licensed under a 3-clause BSD style license - see LICENSE.rst

from ....table import Table
from ... import ascii
from ..fastbasic import FastBasic
from .common import assert_equal, assert_true
from cStringIO import StringIO
import numpy as np

def assert_table_equal(t1, t2):
	assert_equal(len(t1), len(t2))
	assert_equal(t1.colnames, t2.colnames)
	for name in t1.colnames:
		for i, el in enumerate(t1[name]):
			assert_equal(el, t2[name][i])

def read_basic(table, **kwargs):
	reader = FastBasic(**kwargs)
	t1 = reader.read(table)
	t2 = ascii.read(table, format='fast_basic', guess=False, **kwargs)
	t3 = ascii.read(table, format='basic', guess=False, **kwargs)
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
	assert_true(np.issubdtype(table['B'].dtype, np.str_))
	assert_true(np.issubdtype(table['C'].dtype, np.int_))
	assert_true(np.issubdtype(table['D'].dtype, np.float_))
	assert_true(np.issubdtype(table['E'].dtype, np.str_))

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
a b "c
 d"
"""
	table = read_basic(StringIO(text))
	expected = Table([['1.5', 'a'], ['2.1', 'b'], ['-37.1', 'cd']], names=('A B', 'C', 'D'))
	assert_table_equal(table, expected)
	table = read_basic(StringIO(text.replace('"', "'")), quotechar="'")
	assert_table_equal(table, expected)
