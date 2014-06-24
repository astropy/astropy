# Licensed under a 3-clause BSD style license - see LICENSE.rst

from ....table import Table
from ... import ascii
from ..fastbasic import FastBasic
from .common import assert_equal
from cStringIO import StringIO

def read_basic(table, **kwargs):
	reader = FastBasic(**kwargs)
	t1 = reader.read(table)
	t2 = ascii.read(table, format='fast_basic', guess=False, **kwargs)
	assert_equal(t1.colnames, t2.colnames)
	return t1

def test_simple_data():
	"""
	Make sure the fast reader works with basic input data.
	"""
	table = read_basic(StringIO("A B C\n1 2 3\n4 5 6"))
	assert_equal(table.colnames, ['A', 'B', 'C'])
	assert_equal(table['A'][0], 1)
	assert_equal(table['B'][1], 5)
	assert_equal(len(table), 2)

def test_supplied_names():
	"""
	If passed as a parameter, names should replace any
	column names found in the header.
	"""
	table = read_basic(StringIO("A B C\n1 2 3\n4 5 6"), names=('X', 'Y', 'Z'))
	assert_equal(table.colnames, ['X', 'Y', 'Z'])
	assert_equal(len(table), 2)

def test_no_header():
	"""
	The header should not be read when header_start=None. Unless names is
	passed, the column names should be auto-generated.
	"""
	table = read_basic(StringIO("A B C\n1 2 3\n4 5 6"), header_start=None, data_start=0)
	assert_equal(table.colnames, ['col1', 'col2', 'col3'])
	assert_equal(table['col1'][0], 'A')
	assert_equal(table['col2'][1], '2')
	assert_equal(len(table), 3)

def test_no_header_supplied_names():
	"""
	If header_start=None and names is passed as a parameter, header
	data should not be read and names should be used instead.
	"""
	table = read_basic(StringIO("A B C\n1 2 3\n4 5 6"), header_start=None, data_start=0,
					   names=('X', 'Y', 'Z'))
	assert_equal(table.colnames, ['X', 'Y', 'Z'])
	assert_equal(table['X'][0], 'A')
	assert_equal(table['Z'][2], '6')
	assert_equal(len(table), 3)

def test_comment():
	"""
	Make sure that line comments are ignored by the C reader.
	"""
	table = read_basic(StringIO("# comment\nA B C\n# another comment\n1 2 3\n4 5 6"))
	assert_equal(table.colnames, ['A', 'B', 'C'])
	assert_equal(table['A'][0], 1)
	assert_equal(table['B'][1], 5)
	assert_equal(len(table), 2)

def test_empty_lines():
	"""
	Make sure that empty lines are ignored by the C reader.
	"""
	table = read_basic(StringIO("\n\nA B C\n1 2 3\n\n\n4 5 6\n\n\n\n"))
	assert_equal(table.colnames, ['A', 'B', 'C'])
	assert_equal(table['A'][0], 1)
	assert_equal(table['B'][1], 5)
	assert_equal(len(table), 2)
