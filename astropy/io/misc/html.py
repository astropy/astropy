from __future__ import absolute_import

from ...table import Table
from ...utils.xml import writer
from .. import registry

def html_write(table, filename, clobber=False):
    """
    Writes the table data to filename in HTML.
    """
    if not clobber and os.path.exists(filename):
        return
    html_file = open(filename, 'w')
    xml_writer = writer.XMLWriter(html_file)
    with xml_writer.tag('html'):
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
                            
def html_identify(origin, path, fileobj, *args, **kwargs):
    """
    Determines whether the given filename is an HTML file.
    """
    return path.lower().split('.')[-1] in ['htm', 'html']

registry.register_identifier('html', Table, html_identify)
registry.register_writer('html', Table, html_write)
