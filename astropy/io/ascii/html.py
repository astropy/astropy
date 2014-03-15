# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""An extensible HTML table reader and writer.

html.py:
  Classes to read and write HTML tables

`BeautifulSoup <http://www.crummy.com/software/BeautifulSoup/>`_
must be installed to read HTML tables.
"""

from __future__ import absolute_import, division, print_function
from ...extern import six
from ...extern.six.moves import zip as izip

from . import core
from ...table import Table, Column
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

def identify_table(soup, htmldict, numtable):
    """
    Checks whether the given BeautifulSoup tag is inside the table
    the user intends to process.
    """

    for parent in soup.parents:
        if parent.name == 'table':
            table_soup = parent
            break
    else: # Tag is not inside any table
        return False
    if 'table_id' not in htmldict:
        return numtable == 1

    table_id = htmldict['table_id']

    if isinstance(table_id, six.string_types):
        return 'id' in table_soup.attrs and table_soup['id'] == table_id
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
            raise core.OptionalTableImportError('BeautifulSoup must be '
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
                    # Return multicolumns as tuples for HTMLHeader handling
                    yield [(el.string.strip(), el['colspan']) if el.has_attr('colspan')
                           else el.string.strip() for el in header_elements]
                data_elements = soup.find_all('td')
                if data_elements:
                    yield [el.string.strip() for el in data_elements]
        if not data_found:
            raise core.InconsistentTableError('HTML tables must contain data '
                                              'in a <table> tag')

class HTMLOutputter(core.TableOutputter):
    """
    Output the HTML data as an ``astropy.table.Table`` object.

    This subclass allows for the final table to contain
    multidimensional columns (defined using the colspan attribute
    of <th>).
    """

    def __call__(self, cols, meta):
        """
        Process the data in multidimensional columns.
        """
        self._convert_vals(cols)
        new_cols = []
        col_num = 0

        while col_num < len(cols):
            col = cols[col_num]
            if hasattr(col, 'colspan'):
                # Join elements of spanned columns together into tuples
                data = [tuple([cols[i].data[row] for i in range(col_num,
                    col_num + col.colspan)]) for row in range(len(col.data))]
                new_col = core.Column(col.name)
                new_col.data = data
                new_cols.append(new_col)
                col_num += col.colspan
            else:
                new_cols.append(col)
                col_num += 1

        return core.TableOutputter.__call__(self, new_cols, meta)
        

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
            elif soup.name == 'tr' and identify_table(soup, self.html, tables) \
                                             and soup.th is not None:
                return i

        return None

    def _set_cols_from_names(self):
        """
        Set columns from header names, handling multicolumns appropriately.
        """
        self.cols = []
        new_names = []

        for name in self.names:
            if isinstance(name, tuple):
                col = core.Column(name=name[0])
                col.colspan = int(name[1])
                self.cols.append(col)
                new_names.append(name[0])
                for i in range(1, int(name[1])):
                    # Add dummy columns
                    self.cols.append(core.Column(''))
                    new_names.append('')
            else:
                self.cols.append(core.Column(name=name))
                new_names.append(name)

        self.names = new_names
    

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
            elif soup.name == 'tr' and identify_table(soup, self.html, tables) \
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
            elif soup.name == 'tr' and identify_table(soup, self.html, tables):
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

        * multicol : Use multi-dimensional columns for output
            The writer will output tuples as elements of multi-dimensional
            columns if this parameter is true, and if not then it will
            use the syntax 1.36583e-13 .. 1.36583e-13 for output. If not
            present, this parameter will be true by default.
    
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
        if 'multicol' not in htmldict:
            self.html['multicol'] = True
        self.header.html = self.html
        self.data.html = self.html

    def read(self, table):
        """
        Read the ``table`` in HTML format and return a resulting ``Table``.
        """

        self.outputter = HTMLOutputter()
        return core.BaseReader.read(self, table)
    
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
                            if len(col.shape) > 1 and self.html['multicol']:
                                # Set colspan attribute for multicolumns
                                w.start('th', colspan=col.shape[1])
                            else:
                                w.start('th')
                            w.data(col.name.strip())
                            w.end(indent=False)
                    col_str_iters = []
                    for col in cols:
                        if len(col.shape) > 1 and self.html['multicol']:
                            span = col.shape[1]
                            for i in range(span):
                                # Split up multicolumns into separate columns
                                new_col = Column([el[i] for el in col])
                                col_str_iters.append(new_col.iter_str_vals())
                        else:
                            col_str_iters.append(col.iter_str_vals())

                    for row in izip(*col_str_iters):
                        with w.tag('tr'):
                            for el in row:
                                w.start('td')
                                w.data(el.strip())
                                w.end(indent=False)

        # Fixes XMLWriter's insertion of unwanted line breaks
        return [''.join(lines)]
