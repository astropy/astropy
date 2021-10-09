# -*- coding: utf-8 -*-

# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module tests some of the methods related to the ``HTML``
reader/writer and aims to document its functionality.

Requires `BeautifulSoup <http://www.crummy.com/software/BeautifulSoup/>`_
to be installed.
"""

from io import StringIO

from astropy.io.ascii import html
from astropy.io.ascii import core
from astropy.table import Table

import pytest
import numpy as np

from .common import setup_function, teardown_function  # noqa
from astropy.io import ascii

from astropy.utils.compat.optional_deps import HAS_BLEACH, HAS_BS4  # noqa

if HAS_BS4:
    from bs4 import BeautifulSoup, FeatureNotFound


@pytest.mark.skipif('not HAS_BS4')
def test_soupstring():
    """
    Test to make sure the class SoupString behaves properly.
    """

    soup = BeautifulSoup('<html><head></head><body><p>foo</p></body></html>',
                         'html.parser')
    soup_str = html.SoupString(soup)
    assert isinstance(soup_str, str)
    assert isinstance(soup_str, html.SoupString)
    assert soup_str == '<html><head></head><body><p>foo</p></body></html>'
    assert soup_str.soup is soup


def test_listwriter():
    """
    Test to make sure the class ListWriter behaves properly.
    """

    lst = []
    writer = html.ListWriter(lst)

    for i in range(5):
        writer.write(i)
    for ch in 'abcde':
        writer.write(ch)

    assert lst == [0, 1, 2, 3, 4, 'a', 'b', 'c', 'd', 'e']


@pytest.mark.skipif('not HAS_BS4')
def test_identify_table():
    """
    Test to make sure that identify_table() returns whether the
    given BeautifulSoup tag is the correct table to process.
    """

    # Should return False on non-<table> tags and None
    soup = BeautifulSoup('<html><body></body></html>', 'html.parser')
    assert html.identify_table(soup, {}, 0) is False
    assert html.identify_table(None, {}, 0) is False

    soup = BeautifulSoup('<table id="foo"><tr><th>A</th></tr><tr>'
                         '<td>B</td></tr></table>', 'html.parser').table
    assert html.identify_table(soup, {}, 2) is False
    assert html.identify_table(soup, {}, 1) is True  # Default index of 1

    # Same tests, but with explicit parameter
    assert html.identify_table(soup, {'table_id': 2}, 1) is False
    assert html.identify_table(soup, {'table_id': 1}, 1) is True

    # Test identification by string ID
    assert html.identify_table(soup, {'table_id': 'bar'}, 1) is False
    assert html.identify_table(soup, {'table_id': 'foo'}, 1) is True


@pytest.mark.skipif('not HAS_BS4')
def test_missing_data():
    """
    Test reading a table with missing data
    """
    # First with default where blank => '0'
    table_in = ['<table>',
                '<tr><th>A</th></tr>',
                '<tr><td></td></tr>',
                '<tr><td>1</td></tr>',
                '</table>']
    dat = Table.read(table_in, format='ascii.html')
    assert dat.masked is False
    assert np.all(dat['A'].mask == [True, False])
    assert dat['A'].dtype.kind == 'i'

    # Now with a specific value '...' => missing
    table_in = ['<table>',
                '<tr><th>A</th></tr>',
                '<tr><td>...</td></tr>',
                '<tr><td>1</td></tr>',
                '</table>']
    dat = Table.read(table_in, format='ascii.html', fill_values=[('...', '0')])
    assert dat.masked is False
    assert np.all(dat['A'].mask == [True, False])
    assert dat['A'].dtype.kind == 'i'


@pytest.mark.skipif('not HAS_BS4')
def test_rename_cols():
    """
    Test reading a table and renaming cols
    """
    table_in = ['<table>',
                '<tr><th>A</th> <th>B</th></tr>',
                '<tr><td>1</td><td>2</td></tr>',
                '</table>']

    # Swap column names
    dat = Table.read(table_in, format='ascii.html', names=['B', 'A'])
    assert dat.colnames == ['B', 'A']
    assert len(dat) == 1

    # Swap column names and only include A (the renamed version)
    dat = Table.read(table_in, format='ascii.html', names=['B', 'A'], include_names=['A'])
    assert dat.colnames == ['A']
    assert len(dat) == 1
    assert np.all(dat['A'] == 2)


@pytest.mark.skipif('not HAS_BS4')
def test_no_names():
    """
    Test reading a table with no column header
    """
    table_in = ['<table>',
                '<tr><td>1</td></tr>',
                '<tr><td>2</td></tr>',
                '</table>']
    dat = Table.read(table_in, format='ascii.html')
    assert dat.colnames == ['col1']
    assert len(dat) == 2

    dat = Table.read(table_in, format='ascii.html', names=['a'])
    assert dat.colnames == ['a']
    assert len(dat) == 2


@pytest.mark.skipif('not HAS_BS4')
def test_identify_table_fail():
    """
    Raise an exception with an informative error message if table_id
    is not found.
    """
    table_in = ['<table id="foo"><tr><th>A</th></tr>',
                '<tr><td>B</td></tr></table>']

    with pytest.raises(core.InconsistentTableError) as err:
        Table.read(table_in, format='ascii.html', htmldict={'table_id': 'bad_id'},
                   guess=False)
    assert err.match("ERROR: HTML table id 'bad_id' not found$")

    with pytest.raises(core.InconsistentTableError) as err:
        Table.read(table_in, format='ascii.html', htmldict={'table_id': 3},
                   guess=False)
    assert err.match("ERROR: HTML table number 3 not found$")


@pytest.mark.skipif('not HAS_BS4')
def test_backend_parsers():
    """
    Make sure the user can specify which back-end parser to use
    and that an error is raised if the parser is invalid.
    """
    for parser in ('lxml', 'xml', 'html.parser', 'html5lib'):
        try:
            Table.read('data/html2.html', format='ascii.html',
                       htmldict={'parser': parser}, guess=False)
        except FeatureNotFound:
            if parser == 'html.parser':
                raise
            # otherwise ignore if the dependency isn't present

    # reading should fail if the parser is invalid
    with pytest.raises(FeatureNotFound):
        Table.read('data/html2.html', format='ascii.html',
                   htmldict={'parser': 'foo'}, guess=False)


@pytest.mark.skipif('HAS_BS4')
def test_htmlinputter_no_bs4():
    """
    This should return an OptionalTableImportError if BeautifulSoup
    is not installed.
    """

    inputter = html.HTMLInputter()
    with pytest.raises(core.OptionalTableImportError):
        inputter.process_lines([])


@pytest.mark.skipif('not HAS_BS4')
def test_htmlinputter():
    """
    Test to ensure that HTMLInputter correctly converts input
    into a list of SoupStrings representing table elements.
    """

    f = 'data/html.html'
    with open(f) as fd:
        table = fd.read()

    inputter = html.HTMLInputter()
    inputter.html = {}

    # In absence of table_id, defaults to the first table
    expected = ['<tr><th>Column 1</th><th>Column 2</th><th>Column 3</th></tr>',
                '<tr><td>1</td><td>a</td><td>1.05</td></tr>',
                '<tr><td>2</td><td>b</td><td>2.75</td></tr>',
                '<tr><td>3</td><td>c</td><td>-1.25</td></tr>']
    assert [str(x) for x in inputter.get_lines(table)] == expected

    # Should raise an InconsistentTableError if the table is not found
    inputter.html = {'table_id': 4}
    with pytest.raises(core.InconsistentTableError):
        inputter.get_lines(table)

    # Identification by string ID
    inputter.html['table_id'] = 'second'
    expected = ['<tr><th>Column A</th><th>Column B</th><th>Column C</th></tr>',
                '<tr><td>4</td><td>d</td><td>10.5</td></tr>',
                '<tr><td>5</td><td>e</td><td>27.5</td></tr>',
                '<tr><td>6</td><td>f</td><td>-12.5</td></tr>']
    assert [str(x) for x in inputter.get_lines(table)] == expected

    # Identification by integer index
    inputter.html['table_id'] = 3
    expected = ['<tr><th>C1</th><th>C2</th><th>C3</th></tr>',
                '<tr><td>7</td><td>g</td><td>105.0</td></tr>',
                '<tr><td>8</td><td>h</td><td>275.0</td></tr>',
                '<tr><td>9</td><td>i</td><td>-125.0</td></tr>']
    assert [str(x) for x in inputter.get_lines(table)] == expected


@pytest.mark.skipif('not HAS_BS4')
def test_htmlsplitter():
    """
    Test to make sure that HTMLSplitter correctly inputs lines
    of type SoupString to return a generator that gives all
    header and data elements.
    """

    splitter = html.HTMLSplitter()

    lines = [html.SoupString(BeautifulSoup('<table><tr><th>Col 1</th><th>Col 2</th></tr></table>',
                                           'html.parser').tr),
             html.SoupString(BeautifulSoup('<table><tr><td>Data 1</td><td>Data 2</td></tr></table>',
                                           'html.parser').tr)]
    expected_data = [['Col 1', 'Col 2'], ['Data 1', 'Data 2']]
    assert list(splitter(lines)) == expected_data

    # Make sure the presence of a non-SoupString triggers a TypeError
    lines.append('<tr><td>Data 3</td><td>Data 4</td></tr>')
    with pytest.raises(TypeError):
        list(splitter(lines))

    # Make sure that passing an empty list triggers an error
    with pytest.raises(core.InconsistentTableError):
        list(splitter([]))


@pytest.mark.skipif('not HAS_BS4')
def test_htmlheader_start():
    """
    Test to ensure that the start_line method of HTMLHeader
    returns the first line of header data. Uses t/html.html
    for sample input.
    """

    f = 'data/html.html'
    with open(f) as fd:
        table = fd.read()

    inputter = html.HTMLInputter()
    inputter.html = {}
    header = html.HTMLHeader()

    lines = inputter.get_lines(table)
    assert str(lines[header.start_line(lines)]) == \
        '<tr><th>Column 1</th><th>Column 2</th><th>Column 3</th></tr>'
    inputter.html['table_id'] = 'second'
    lines = inputter.get_lines(table)
    assert str(lines[header.start_line(lines)]) == \
        '<tr><th>Column A</th><th>Column B</th><th>Column C</th></tr>'
    inputter.html['table_id'] = 3
    lines = inputter.get_lines(table)
    assert str(lines[header.start_line(lines)]) == \
        '<tr><th>C1</th><th>C2</th><th>C3</th></tr>'

    # start_line should return None if no valid header is found
    lines = [html.SoupString(BeautifulSoup('<table><tr><td>Data</td></tr></table>',
                                           'html.parser').tr),
             html.SoupString(BeautifulSoup('<p>Text</p>', 'html.parser').p)]
    assert header.start_line(lines) is None

    # Should raise an error if a non-SoupString is present
    lines.append('<tr><th>Header</th></tr>')
    with pytest.raises(TypeError):
        header.start_line(lines)


@pytest.mark.skipif('not HAS_BS4')
def test_htmldata():
    """
    Test to ensure that the start_line and end_lines methods
    of HTMLData returns the first line of table data. Uses
    t/html.html for sample input.
    """

    f = 'data/html.html'
    with open(f) as fd:
        table = fd.read()

    inputter = html.HTMLInputter()
    inputter.html = {}
    data = html.HTMLData()

    lines = inputter.get_lines(table)
    assert str(lines[data.start_line(lines)]) == \
        '<tr><td>1</td><td>a</td><td>1.05</td></tr>'
    # end_line returns the index of the last data element + 1
    assert str(lines[data.end_line(lines) - 1]) == \
        '<tr><td>3</td><td>c</td><td>-1.25</td></tr>'

    inputter.html['table_id'] = 'second'
    lines = inputter.get_lines(table)
    assert str(lines[data.start_line(lines)]) == \
        '<tr><td>4</td><td>d</td><td>10.5</td></tr>'
    assert str(lines[data.end_line(lines) - 1]) == \
        '<tr><td>6</td><td>f</td><td>-12.5</td></tr>'

    inputter.html['table_id'] = 3
    lines = inputter.get_lines(table)
    assert str(lines[data.start_line(lines)]) == \
        '<tr><td>7</td><td>g</td><td>105.0</td></tr>'
    assert str(lines[data.end_line(lines) - 1]) == \
        '<tr><td>9</td><td>i</td><td>-125.0</td></tr>'

    # start_line should raise an error if no table data exists
    lines = [html.SoupString(BeautifulSoup('<div></div>', 'html.parser').div),
             html.SoupString(BeautifulSoup('<p>Text</p>', 'html.parser').p)]
    with pytest.raises(core.InconsistentTableError):
        data.start_line(lines)

    # end_line should return None if no table data exists
    assert data.end_line(lines) is None

    # Should raise an error if a non-SoupString is present
    lines.append('<tr><td>Data</td></tr>')
    with pytest.raises(TypeError):
        data.start_line(lines)
    with pytest.raises(TypeError):
        data.end_line(lines)


def test_multicolumn_write():
    """
    Test to make sure that the HTML writer writes multidimensional
    columns (those with iterable elements) using the colspan
    attribute of <th>.
    """

    col1 = [1, 2, 3]
    col2 = [(1.0, 1.0), (2.0, 2.0), (3.0, 3.0)]
    col3 = [('a', 'a', 'a'), ('b', 'b', 'b'), ('c', 'c', 'c')]
    table = Table([col1, col2, col3], names=('C1', 'C2', 'C3'))
    expected = """\
<html>
 <head>
  <meta charset="utf-8"/>
  <meta content="text/html;charset=UTF-8" http-equiv="Content-type"/>
 </head>
 <body>
  <table>
   <thead>
    <tr>
     <th>C1</th>
     <th colspan="2">C2</th>
     <th colspan="3">C3</th>
    </tr>
   </thead>
   <tr>
    <td>1</td>
    <td>1.0</td>
    <td>1.0</td>
    <td>a</td>
    <td>a</td>
    <td>a</td>
   </tr>
   <tr>
    <td>2</td>
    <td>2.0</td>
    <td>2.0</td>
    <td>b</td>
    <td>b</td>
    <td>b</td>
   </tr>
   <tr>
    <td>3</td>
    <td>3.0</td>
    <td>3.0</td>
    <td>c</td>
    <td>c</td>
    <td>c</td>
   </tr>
  </table>
 </body>
</html>
    """
    out = html.HTML().write(table)[0].strip()
    assert out == expected.strip()


@pytest.mark.skipif('not HAS_BLEACH')
def test_multicolumn_write_escape():
    """
    Test to make sure that the HTML writer writes multidimensional
    columns (those with iterable elements) using the colspan
    attribute of <th>.
    """

    col1 = [1, 2, 3]
    col2 = [(1.0, 1.0), (2.0, 2.0), (3.0, 3.0)]
    col3 = [('<a></a>', '<a></a>', 'a'), ('<b></b>', 'b', 'b'), ('c', 'c', 'c')]
    table = Table([col1, col2, col3], names=('C1', 'C2', 'C3'))
    expected = """\
<html>
 <head>
  <meta charset="utf-8"/>
  <meta content="text/html;charset=UTF-8" http-equiv="Content-type"/>
 </head>
 <body>
  <table>
   <thead>
    <tr>
     <th>C1</th>
     <th colspan="2">C2</th>
     <th colspan="3">C3</th>
    </tr>
   </thead>
   <tr>
    <td>1</td>
    <td>1.0</td>
    <td>1.0</td>
    <td><a></a></td>
    <td><a></a></td>
    <td>a</td>
   </tr>
   <tr>
    <td>2</td>
    <td>2.0</td>
    <td>2.0</td>
    <td><b></b></td>
    <td>b</td>
    <td>b</td>
   </tr>
   <tr>
    <td>3</td>
    <td>3.0</td>
    <td>3.0</td>
    <td>c</td>
    <td>c</td>
    <td>c</td>
   </tr>
  </table>
 </body>
</html>
    """
    out = html.HTML(htmldict={'raw_html_cols': 'C3'}).write(table)[0].strip()
    assert out == expected.strip()


def test_write_no_multicols():
    """
    Test to make sure that the HTML writer will not use
    multi-dimensional columns if the multicol parameter
    is False.
    """

    col1 = [1, 2, 3]
    col2 = [(1.0, 1.0), (2.0, 2.0), (3.0, 3.0)]
    col3 = [('a', 'a', 'a'), ('b', 'b', 'b'), ('c', 'c', 'c')]
    table = Table([col1, col2, col3], names=('C1', 'C2', 'C3'))
    expected = """\
<html>
 <head>
  <meta charset="utf-8"/>
  <meta content="text/html;charset=UTF-8" http-equiv="Content-type"/>
 </head>
 <body>
  <table>
   <thead>
    <tr>
     <th>C1</th>
     <th>C2</th>
     <th>C3</th>
    </tr>
   </thead>
   <tr>
    <td>1</td>
    <td>1.0 .. 1.0</td>
    <td>a .. a</td>
   </tr>
   <tr>
    <td>2</td>
    <td>2.0 .. 2.0</td>
    <td>b .. b</td>
   </tr>
   <tr>
    <td>3</td>
    <td>3.0 .. 3.0</td>
    <td>c .. c</td>
   </tr>
  </table>
 </body>
</html>
    """
    assert html.HTML({'multicol': False}).write(table)[0].strip() == \
        expected.strip()


@pytest.mark.skipif('not HAS_BS4')
def test_multicolumn_read():
    """
    Test to make sure that the HTML reader inputs multidimensional
    columns (those with iterable elements) using the colspan
    attribute of <th>.

    Ensure that any string element within a multidimensional column
    casts all elements to string prior to type conversion operations.
    """

    table = Table.read('data/html2.html', format='ascii.html')
    str_type = np.dtype((str, 21))
    expected = Table(np.array([(['1', '2.5000000000000000001'], 3),
                               (['1a', '1'], 3.5)],
                              dtype=[('A', str_type, (2,)), ('B', '<f8')]))
    assert np.all(table == expected)


@pytest.mark.skipif('not HAS_BLEACH')
def test_raw_html_write():
    """
    Test that columns can contain raw HTML which is not escaped.
    """
    t = Table([['<em>x</em>'], ['<em>y</em>']], names=['a', 'b'])

    # One column contains raw HTML (string input)
    out = StringIO()
    t.write(out, format='ascii.html', htmldict={'raw_html_cols': 'a'})
    expected = """\
   <tr>
    <td><em>x</em></td>
    <td>&lt;em&gt;y&lt;/em&gt;</td>
   </tr>"""
    assert expected in out.getvalue()

    # One column contains raw HTML (list input)
    out = StringIO()
    t.write(out, format='ascii.html', htmldict={'raw_html_cols': ['a']})
    assert expected in out.getvalue()

    # Two columns contains raw HTML (list input)
    out = StringIO()
    t.write(out, format='ascii.html', htmldict={'raw_html_cols': ['a', 'b']})
    expected = """\
   <tr>
    <td><em>x</em></td>
    <td><em>y</em></td>
   </tr>"""
    assert expected in out.getvalue()


@pytest.mark.skipif('not HAS_BLEACH')
def test_raw_html_write_clean():
    """
    Test that columns can contain raw HTML which is not escaped.
    """
    import bleach  # noqa

    t = Table([['<script>x</script>'], ['<p>y</p>'], ['<em>y</em>']], names=['a', 'b', 'c'])

    # Confirm that <script> and <p> get escaped but not <em>
    out = StringIO()
    t.write(out, format='ascii.html', htmldict={'raw_html_cols': t.colnames})
    expected = """\
   <tr>
    <td>&lt;script&gt;x&lt;/script&gt;</td>
    <td>&lt;p&gt;y&lt;/p&gt;</td>
    <td><em>y</em></td>
   </tr>"""
    assert expected in out.getvalue()

    # Confirm that we can whitelist <p>
    out = StringIO()
    t.write(out, format='ascii.html',
            htmldict={'raw_html_cols': t.colnames,
                      'raw_html_clean_kwargs': {'tags': bleach.ALLOWED_TAGS + ['p']}})
    expected = """\
   <tr>
    <td>&lt;script&gt;x&lt;/script&gt;</td>
    <td><p>y</p></td>
    <td><em>y</em></td>
   </tr>"""
    assert expected in out.getvalue()


def test_write_table_html_fill_values():
    """
    Test that passing fill_values should replace any matching row
    """
    buffer_output = StringIO()
    t = Table([[1], [2]], names=('a', 'b'))
    ascii.write(t, buffer_output, fill_values=('1', 'Hello world'),
                format='html')

    t_expected = Table([['Hello world'], [2]], names=('a', 'b'))
    buffer_expected = StringIO()
    ascii.write(t_expected, buffer_expected, format='html')

    assert buffer_output.getvalue() == buffer_expected.getvalue()


def test_write_table_html_fill_values_optional_columns():
    """
    Test that passing optional column in fill_values should only replace
    matching columns
    """
    buffer_output = StringIO()
    t = Table([[1], [1]], names=('a', 'b'))
    ascii.write(t, buffer_output, fill_values=('1', 'Hello world', 'b'),
                format='html')

    t_expected = Table([[1], ['Hello world']], names=('a', 'b'))
    buffer_expected = StringIO()
    ascii.write(t_expected, buffer_expected, format='html')

    assert buffer_output.getvalue() == buffer_expected.getvalue()


def test_write_table_html_fill_values_masked():
    """
    Test that passing masked values in fill_values should only replace
    masked columns or values
    """
    buffer_output = StringIO()
    t = Table([[1], [1]], names=('a', 'b'), masked=True, dtype=('i4', 'i8'))
    t['a'] = np.ma.masked
    ascii.write(t, buffer_output, fill_values=(ascii.masked, 'TEST'),
                format='html')

    t_expected = Table([['TEST'], [1]], names=('a', 'b'))
    buffer_expected = StringIO()
    ascii.write(t_expected, buffer_expected, format='html')

    assert buffer_output.getvalue() == buffer_expected.getvalue()


def test_multicolumn_table_html_fill_values():
    """
    Test to make sure that the HTML writer writes multidimensional
    columns with correctly replaced fill_values.
    """
    col1 = [1, 2, 3]
    col2 = [(1.0, 1.0), (2.0, 2.0), (3.0, 3.0)]
    col3 = [('a', 'a', 'a'), ('b', 'b', 'b'), ('c', 'c', 'c')]

    buffer_output = StringIO()
    t = Table([col1, col2, col3], names=('C1', 'C2', 'C3'))
    ascii.write(t, buffer_output, fill_values=('a', 'z'),
                format='html')

    col1 = [1, 2, 3]
    col2 = [(1.0, 1.0), (2.0, 2.0), (3.0, 3.0)]
    col3 = [('z', 'z', 'z'), ('b', 'b', 'b'), ('c', 'c', 'c')]

    buffer_expected = StringIO()
    t_expected = Table([col1, col2, col3], names=('C1', 'C2', 'C3'))
    ascii.write(t_expected, buffer_expected, format='html')

    assert buffer_output.getvalue() == buffer_expected.getvalue()


def test_multi_column_write_table_html_fill_values_masked():
    """
    Test that passing masked values in fill_values should only replace
    masked columns or values for multidimensional tables
    """
    buffer_output = StringIO()
    t = Table([[1, 2, 3, 4], ['--', 'a', '--', 'b']], names=('a', 'b'), masked=True)
    t['a'][0:2] = np.ma.masked
    t['b'][0:2] = np.ma.masked
    ascii.write(t, buffer_output, fill_values=[(ascii.masked, 'MASKED')],
                format='html')

    t_expected = Table([['MASKED', 'MASKED', 3, 4], [
                       'MASKED', 'MASKED', '--', 'b']], names=('a', 'b'))
    buffer_expected = StringIO()
    ascii.write(t_expected, buffer_expected, format='html')
    print(buffer_expected.getvalue())

    assert buffer_output.getvalue() == buffer_expected.getvalue()


@pytest.mark.skipif('not HAS_BS4')
def test_read_html_unicode():
    """
    Test reading an HTML table with unicode values
    """
    table_in = ['<table>',
                '<tr><td>&#x0394;</td></tr>',
                '<tr><td>Δ</td></tr>',
                '</table>']
    dat = Table.read(table_in, format='ascii.html')
    assert np.all(dat['col1'] == ['Δ', 'Δ'])
