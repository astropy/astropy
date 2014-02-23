# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package contains functions for reading and writing HTML files
containing table data, meant to be used as readers/writers in
`astropy.table`. See :ref:`table_io` for more details.
"""

from __future__ import absolute_import

from ...table import Table
from ...utils.xml import writer

__all__ = ['html_write', 'html_read', 'html_identify']

def html_write(table, filename, clobber=False):
    """
    Writes the table data to filename in HTML.
    """
    if not clobber and os.path.exists(filename):
        return
    html_file = open(filename, 'w')
    xml_writer = writer.XMLWriter(html_file)
    with xml_writer.tag('html'):
        with xml_writer.tag('head'):
            with xml_writer.tag('style'):
                xml_writer.data(
                    'table,th,td{border:1px solid black;'
                    'border-collapse:collapse;}'
                    'th,td{padding:5px;}')
        with xml_writer.tag('table'):
            with xml_writer.tag('tr'):
                for colname in table.colnames:
                    with xml_writer.tag('th'):
                        xml_writer.data(colname)
            for row_num in range(len(table)):
                with xml_writer.tag('tr'):
                    for colname in table.colnames:
                        with xml_writer.tag('td'):
                            xml_writer.data(str(table[colname][row_num]))
    html_file.close()

def html_read(fileorname):
    """
    Reads HTML input from the given file into a table.

    Parameters
    ----------
    fileorname : str or `file`-like object
        The file name or file from which to extract table data.

    This requires `BeautifulSoup
    <http://www.crummy.com/software/BeautifulSoup/>`_ to be installed.
    """

    from bs4 import BeautifulSoup

    if isinstance(fileorname, basestring):
        html_file = open(fileorname)
        html = html_file.read()
        html_file.close()
    else:
        html = fileorname.read()
    soup = BeautifulSoup(html)
    table = soup.find('table')
    colnames = []
    grid = []
    for row in table.findAll('tr'):
        row_list = []
        for header_el in row.findAll('th'):
            colnames.append(header_el.string.replace('\n', ''))
        for el in row.findAll('td'):
            row_list.append(el.string.replace('\n', ''))
        if len(row_list) > 0:
            grid.append(row_list)

    # Reverse rows and columns as required by Table.__init__
    transpose = [[grid[j][i] for j in range(len(grid))]
                 for i in range(len(grid[0]))]
    return Table(transpose, names=colnames)

def html_identify(origin, filepath, fileobj, *args, **kwargs):
    """
    Determines whether the given filename is an HTML file.
    """
    if fileobj is not None:
        filepath = fileobj.name
    if filepath is not None:
        return filepath.endswith('.html') or filepath.endswith('.htm')
    return False
