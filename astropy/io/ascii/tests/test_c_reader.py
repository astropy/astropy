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

