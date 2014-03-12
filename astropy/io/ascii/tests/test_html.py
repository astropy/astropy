# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module tests some of the methods related to the ``HTML``
reader/writer and aims to document its functionality.

Requires `BeautifulSoup <http://www.crummy.com/software/BeautifulSoup/>`_
to be installed.
"""

from .. import html
from .. import core
from ....table import Table

from ....tests.helper import pytest
from ....extern.six.moves import zip as izip
from .common import (raises, assert_equal, assert_almost_equal,
                     assert_true, setup_function, teardown_function)

try:
    from itertools import izip
except ImportError:
    izip = zip

# Check to see if the BeautifulSoup dependency is present.

try:
    from bs4 import BeautifulSoup
    HAS_BEAUTIFUL_SOUP = True
except ImportError:
    HAS_BEAUTIFUL_SOUP = False

@pytest.mark.skipif('not HAS_BEAUTIFUL_SOUP')
def test_soupstring():
    """
    Test to make sure the class SoupString behaves properly.
    """
    
    soup = BeautifulSoup('<html><body><p>foo</p></body></html>')
    soup_str = html.SoupString(soup)
    assert isinstance(soup_str, str)
    assert isinstance(soup_str, html.SoupString)
    assert soup_str == '<html><body><p>foo</p></body></html>'
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

@pytest.mark.skipif('not HAS_BEAUTIFUL_SOUP')
def test_correct_table():
    """
    Test to make sure that correct_table() returns whether the
    given BeautifulSoup tag is the correct table to process.
    """

    # Should return False on non-tables
    soup = BeautifulSoup('<html><body></body></html>')
    assert html.correct_table(soup, {}, 0) is False

    soup = BeautifulSoup('<table id="foo"><tr><th>A</th></tr><tr>' \
                         '<td>B</td></tr></table>').table
    assert html.correct_table(soup, {}, 2) is False
    assert html.correct_table(soup, {}, 1) is True # Default index of 1

    # Same tests, but with explicit parameter
    assert html.correct_table(soup, {'table_id': 2}, 1) is False
    assert html.correct_table(soup, {'table_id': 1}, 1) is True

    # Test identification by string ID
    assert html.correct_table(soup, {'table_id': 'bar'}, 1) is False
    assert html.correct_table(soup, {'table_id': 'foo'}, 1) is True

@pytest.mark.skipif('HAS_BEAUTIFUL_SOUP')
def test_htmlinputter_no_bs4():
    """
    This should return an InconsistentTableError if BeautifulSoup
    is not installed.
    """

    inputter = html.HTMLInputter()
    with pytest.raises(core.InconsistentTableError):
        inputter.process_lines([])
    
@pytest.mark.skipif('not HAS_BEAUTIFUL_SOUP')
def test_htmlinputter():
    """
    Test to ensure that HTMLInputter correctly converts input
    into a list of valid SoupString objects.
    """

    inputter = html.HTMLInputter()
    
    lines = ['<html>', '<body>', '<a href="http://www.astropy.org/">Link',
             '</a><br></br>', '</body>', '<!--Comment-->', '</html>']
    soup_list = inputter.process_lines(lines)
    expected_soups = ['<html><body><a href="http://www.astropy.org/">Link</a><br/>' \
        '</body><!--Comment--></html>', '<body><a href="http://www.astropy.org/">Link' \
        '</a><br/></body>', '<a href="http://www.astropy.org/">Link</a>',
        'Link', '<br/>']
    for soup, expected_soup in izip(soup_list, expected_soups):
        assert str(soup).replace('\n', '') == expected_soup

@pytest.mark.skipif('not HAS_BEAUTIFUL_SOUP')
def test_htmlsplitter():
    """
    Test to make sure that HTMLSplitter correctly inputs lines
    of type SoupString to return a generator that gives all
    header and data elements.
    """

    splitter = html.HTMLSplitter()

    lines = [html.SoupString(BeautifulSoup('<a href="http://www.astropy.org/">Link</a>').a),
                html.SoupString(BeautifulSoup('<tr><th>Col 1</th><th>Col 2</th></tr>').tr),
                html.SoupString(BeautifulSoup('<tr><td>Data 1</td><td>Data 2</td></tr>').tr)]
    expected_data = [['Col 1', 'Col 2'], ['Data 1', 'Data 2']]
    assert list(splitter(lines)) == expected_data

    # Make sure the presence of a non-SoupString triggers a TypeError
    lines.append('<tr><td>Data 3</td><td>Data 4</td></tr>')
    with pytest.raises(TypeError):
        list(splitter(lines))

    # Make sure that not having a <tr> tag in lines triggers an error
    lines = [html.SoupString(BeautifulSoup('<a href="http://www.astropy.org/">Link</a>').a),
                html.SoupString(BeautifulSoup('<p>Text</p>').p)]
    with pytest.raises(core.InconsistentTableError):
        list(splitter(lines))

@pytest.mark.skipif('not HAS_BEAUTIFUL_SOUP')
def test_htmlheader_start():
    """
    Test to ensure that the start_line method of HTMLHeader
    returns the first line of header data. Uses t/html.html
    for sample input.
    """

    f = 't/html.html'
    with open(f) as fd:
        table = fd.read()

    inputter = html.HTMLInputter()
    
    # Create a list of SoupStrings from raw text
    lines = inputter.get_lines(table)
    header = html.HTMLHeader()
    header.html = {}
    
    # In absence of table_id, defaults to the first line of the first table
    assert str(lines[header.start_line(lines)]) == \
           '<tr><th>Column 1</th><th>Column 2</th><th>Column 3</th></tr>'
    header.html = {'table_id': 'second'}
    assert str(lines[header.start_line(lines)]) == \
           '<tr><th>Column A</th><th>Column B</th><th>Column C</th></tr>'
    header.html = {'table_id': 3}
    assert str(lines[header.start_line(lines)]) == \
           '<tr><th>C1</th><th>C2</th><th>C3</th></tr>'

    # Should raise an error if the desired table is not found
    header.html = {'table_id': 4}
    with pytest.raises(core.InconsistentTableError):
        header.start_line(lines)

    # Should raise an error if a non-SoupString is present
    lines.append('<tr><th>Header</th></tr>')
    with pytest.raises(TypeError):
        header.start_line(lines)

@pytest.mark.skipif('not HAS_BEAUTIFUL_SOUP')
def test_htmldata():
    """
    Test to ensure that the start_line and end_lines methods
    of HTMLData returns the first line of table data. Uses
    t/html.html for sample input.
    """

    f = 't/html.html'
    with open(f) as fd:
        table = fd.read()

    inputter = html.HTMLInputter()

    # Create a list of SoupStrings from raw text
    lines = inputter.get_lines(table)
    data = html.HTMLData()
    data.html = {}

    # In absence of table_id, defaults to the first table
    assert str(lines[data.start_line(lines)]) == \
           '<tr><td>1</td><td>a</td><td>1.05</td></tr>'
    # end_line returns the index of the last data element
    assert str(lines[data.end_line(lines)]) == '<td>3</td>'
    
    data.html = {'table_id': 'second'}
    assert str(lines[data.start_line(lines)]) == \
           '<tr><td>4</td><td>d</td><td>10.5</td></tr>'
    assert str(lines[data.end_line(lines)]) == '<td>6</td>'
    
    data.html = {'table_id': 3}
    assert str(lines[data.start_line(lines)]) == \
           '<tr><td>7</td><td>g</td><td>105.0</td></tr>'
    assert str(lines[data.end_line(lines)]) == '<td>9</td>'
    
    # start_line should raise an error if the desired table is not found
    data.html = {'table_id': 'foo'}
    with pytest.raises(core.InconsistentTableError):
        data.start_line(lines)

    # end_line should return None if the desired table is not found
    assert data.end_line(lines) is None

    # Should raise an error if a non-SoupString is present
    lines.append('<tr><td>Data</td></tr>')
    with pytest.raises(TypeError):
        data.start_line(lines)
    with pytest.raises(TypeError):
        data.end_line(lines)

def test_multicolumn_write():
    """
    Test to make sure that the HTML writer writes multimensional
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
   <tr>
    <th>C1</th>
    <th colspan="2">C2</th>
    <th colspan="3">C3</th>
   </tr>
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
    assert html.HTML().write(table)[0].strip() == expected.strip()

@pytest.mark.skipif('not HAS_BEAUTIFUL_SOUP')
def test_multicolumn_read():
    """
    Test to make sure that the HTML reader inputs multimensional
    columns (those with iterable elements) using the colspan
    attribute of <th>.
    """

    table = Table.read('t/html2.html', format='ascii.html')
    col1 = [(1, 2)]
    col2 = [3]
    expected = Table([col1, col2], names=('A', 'B'))
    assert table == expected
