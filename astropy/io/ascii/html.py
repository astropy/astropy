# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""An extensible HTML table reader and writer.

html.py:
  Classes to read and write HTML tables

This requires `BeautifulSoup
        <http://www.crummy.com/software/BeautifulSoup/>`_ to be installed.
"""

from __future__ import absolute_import, division, print_function
from ...extern import six

from . import core
from ...table import Table
from ...utils.xml import writer

class SoupString(str):
    """
    Allows for strings to hold BeautifulSoup data.
    """

    def __new__(cls, *args, **kwargs):
        return str.__new__(cls, *args, **kwargs)

    def __init__(self, val):
        self.soup = val

class ListWriter:
    """
    Allows for XMLWriter to write to a list instead of a file.
    """

    def __init__(self, out):
        self.out = out

    def write(self, data):
        self.out.append(data)

def correct_table(soup, htmldict, numtable):
    """
    Checks whether the given BeautifulSoup tag is the table
    the user intends to process.
    """

    if soup.name != 'table':
        return False
    if 'table_id' not in htmldict:
        return numtable == 1

    table_id = htmldict['table_id']

    if isinstance(table_id, six.string_types):
        return 'id' in soup.attrs and soup['id'] == table_id
    elif isinstance(table_id, int):
        return table_id == numtable
    
    # Return False if an invalid parameter is given
    return False
    

class HTMLInputter(core.BaseInputter):
    """
    Input lines of HTML in a valid form.

    This requires `BeautifulSoup
        <http://www.crummy.com/software/BeautifulSoup/>`_ to be installed.
    """

    def process_lines(self, lines):
        """
        Convert the given input into a list of SoupString objects
        for further processing.
        """
        
        try:
            from bs4 import BeautifulSoup
            from bs4.element import Comment
        except ImportError:
            raise core.InconsistentTableError('BeautifulSoup must be '
                                        'installed to read HTML tables')
        
        soup = BeautifulSoup('\n'.join(lines))
        soup_list = []
        for x in soup.descendants: # Navigate down HTML hierarchy
            # Remove all blank elements and comments
            if str(x).strip() and not isinstance(x, Comment):
                soup_obj = SoupString(x)
                soup_list.append(soup_obj)
        return soup_list
        
class HTMLSplitter(core.BaseSplitter):
    """
    Split HTML table data.
    """

    def __call__(self, lines):
        """
        Return HTML data from lines as a generator.
        """
        data_found = False
        for line in lines:
            if not isinstance(line, SoupString):
                raise TypeError('HTML lines should be of type SoupString')
            soup = line.soup
            if soup.name == 'tr':
                data_found = True
                header_elements = soup.find_all('th')
                if header_elements:
                    yield [el.string.strip() for el in header_elements]
                data_elements = soup.find_all('td')
                if data_elements:
                    yield [el.string.strip() for el in data_elements]
        if not data_found:
            raise core.InconsistentTableError('HTML tables must contain data '
                                              'in a <table> tag')

class HTMLHeader(core.BaseHeader):
    def start_line(self, lines):
        """
        Return the line number at which header data begins.
        """
        tables = 0
        
        for i, line in enumerate(lines):
            if not isinstance(line, SoupString):
                raise TypeError('HTML lines should be of type SoupString')
            soup = line.soup
            if soup.name == 'table':
                tables += 1
            elif soup.name == 'tr' and correct_table(soup.parent, self.html, tables):
                return i
        
        raise core.InconsistentTableError('No start line found for HTML header')

class HTMLData(core.BaseData):
    def start_line(self, lines):
        """
        Return the line number at which table data begins.
        """
        tables = 0
        
        for i, line in enumerate(lines):
            if not isinstance(line, SoupString):
                raise TypeError('HTML lines should be of type SoupString')
            soup = line.soup
            
            if soup.name == 'table':
                tables += 1
            elif soup.name == 'tr' and correct_table(soup.parent, self.html, tables) \
                                                and soup.td is not None:
                if soup.th is not None:
                    raise core.InconsistentTableError('HTML tables cannot '
                                'have headings and data in the same row')
                return i
        
        raise core.InconsistentTableError('No start line found for HTML data')
    def end_line(self, lines):
        """
        Return the line number at which table data ends.
        """
        last_index = -1
        tables = 0
        
        for i, line in enumerate(lines):
            if not isinstance(line, SoupString):
                raise TypeError('HTML lines should be of type SoupString')
            soup = line.soup
            if soup.name == 'table':
                tables += 1
            elif soup.name == 'tr' and correct_table(soup.parent, self.html, tables):
                last_index = i

        if last_index == -1:
            return None
        return last_index + 1

class HTML(core.BaseReader):
    """
    Read and write HTML tables.

    In order to customize input and output, a dict of parameters may
    be passed to this class holding specific customizations.

    **htmldict** : Dictionary of parameters for HTML input/output.

        * css : Customized styling
            If present, this parameter will be included in a <style>
            tag and will define stylistic attributes of the output.

        * table_id : ID for the input table
            If a string, this defines the HTML id of the table to be processed.
            If an integer, this specificies the index of the input table in the
            available tables. Unless this parameter is given, the reader will
            use the first table found in the input file.
    
    """
    
    _format_name = 'html'
    _io_registry_format_aliases = ['html']
    _io_registry_suffix = '.html'
    _description = 'HTML table'

    def __init__(self, htmldict={}):
        """
        Initialize classes for HTML reading and writing.
        """
        core.BaseReader.__init__(self)
        self.inputter = HTMLInputter()
        self.header = HTMLHeader()
        self.data = HTMLData()
        self.header.splitter = HTMLSplitter()
        self.header.inputter = HTMLInputter()
        self.data.splitter = HTMLSplitter()
        self.data.inputter = HTMLInputter()
        self.data.header = self.header
        self.header.data = self.data
        self.html = htmldict
        self.header.html = self.html
        self.data.html = self.html

    def write(self, table):
        """
        Return data in ``table`` converted to HTML as a list of strings.
        """

        cols = list(six.itervalues(table.columns))
        lines = []

        # Use XMLWriter to output HTML to lines
        w = writer.XMLWriter(ListWriter(lines))
        
        with w.tag('html'):
            with w.tag('head'):
                # Declare encoding and set CSS style for table
                with w.tag('meta', attrib={'charset':'utf-8'}):
                    pass
                with w.tag('meta', attrib={'http-equiv':'Content-type',
                                    'content':'text/html;charset=UTF-8'}):
                    pass
                if 'css' in self.html:
                    with w.tag('style'):
                        w.data(self.html['css'])
            with w.tag('body'):
                with w.tag('table'):
                    with w.tag('tr'):
                        for col in cols:
                            with w.tag('th'):
                                w.data(col.name)
                    col_str_iters = [col.iter_str_vals() for col in cols]
                    for row in zip(*col_str_iters):
                        with w.tag('tr'):
                            for el in row:
                               with w.tag('td'):
                                   w.data(el)

        # Fixes XMLWriter's insertion of unwanted line breaks
        return [''.join(lines)]
