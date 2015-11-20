# Licensed under a 3-clause BSD style license - see LICENSE.rst

from os.path import abspath, dirname, join

from .table import Table

from ..io import registry as io_registry
from .. import config as _config
from .. import extern


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.table.jsviewer`.
    """

    jquery_url = _config.ConfigItem(
        'https://code.jquery.com/jquery-1.11.3.min.js',
        'The URL to the jquery library.')

    datatables_url = _config.ConfigItem(
        'https://cdn.datatables.net/1.10.9/js/jquery.dataTables.min.js',
        'The URL to the jquery datatables library.')

    css_urls = _config.ConfigItem(
        ['https://cdn.datatables.net/1.10.9/css/jquery.dataTables.css'],
        'The URLs to the css file(s) to include.', cfgtype='list')


conf = Conf()


JQUERY_URL = _config.ConfigAlias(
    '0.4', 'JQUERY_URL', 'jquery_url',
    'astropy.table.jsviewer', 'astropy.table.jsviewer')

DATATABLES_URL = _config.ConfigAlias(
    '0.4', 'DATATABLES_URL', 'datatables_url',
    'astropy.table.jsviewer', 'astropy.table.jsviewer')

EXTERN_JS_DIR = abspath(join(dirname(extern.__file__), 'js'))
EXTERN_CSS_DIR = abspath(join(dirname(extern.__file__), 'css'))

IPYNB_JS_SCRIPT = """
<script>
require.config({{paths: {{
    datatables: '{datatables_url}'
}}}});
require(["datatables"], function(){{
    console.log("$('#{tid}').dataTable()");
    $('#{tid}').dataTable({{
        "iDisplayLength": {display_length},
        "aLengthMenu": {display_length_menu},
        "pagingType": "full_numbers"
    }});
}});
</script>
"""

HTML_JS_SCRIPT = """
$(document).ready(function() {{
    $('#{tid}').dataTable({{
     "iDisplayLength": {display_length},
     "aLengthMenu": {display_length_menu},
     "pagingType": "full_numbers"
    }});
}} );
"""


# Default CSS for the JSViewer writer
DEFAULT_CSS = """\
body {font-family: sans-serif;}
table.dataTable {width: auto !important; margin: 0 !important;}
.dataTables_filter, .dataTables_paginate {float: left !important; margin-left:1em}
"""


# Default CSS used when rendering a table in the IPython notebook
DEFAULT_CSS_NB = """\
table.dataTable {clear: both; width: auto !important; margin: 0 !important;}
.dataTables_info, .dataTables_length, .dataTables_filter, .dataTables_paginate{
display: inline-block; margin-right: 1em; }
.paginate_button { margin-right: 5px; }
"""


class JSViewer(object):
    """Provides an interactive HTML export of a Table.

    This class provides an interface to the `DataTables
    <http://datatables.net/>`_ library, which allow to visualize interactively
    an HTML table. It is used by the `~astropy.table.Table.show_in_browser`
    method.

    Parameters
    ----------
    use_local_files : bool, optional
        Use local files or a CDN for JavaScript libraries. Default False.
    display_length : int, optional
        Number or rows to show. Default to 50.

    """

    def __init__(self, use_local_files=False, display_length=50):
        self._use_local_files = use_local_files
        self.display_length_menu = [[10, 25, 50, 100, 500, 1000, -1],
                                    [10, 25, 50, 100, 500, 1000, "All"]]
        self.display_length = display_length
        for L in self.display_length_menu:
            if display_length not in L:
                L.insert(0, display_length)

    @property
    def jquery_urls(self):
        if self._use_local_files:
            return ['file://' + join(EXTERN_JS_DIR, 'jquery-1.11.3.min.js'),
                    'file://' + join(EXTERN_JS_DIR, 'jquery.dataTables.min.js')]
        else:
            return [conf.jquery_url, conf.datatables_url]

    @property
    def css_urls(self):
        if self._use_local_files:
            return ['file://' + join(EXTERN_CSS_DIR,
                                     'jquery.dataTables.css')]
        else:
            return conf.css_urls

    def _jstable_file(self):
        if self._use_local_files:
            return 'file://' + join(EXTERN_JS_DIR, 'jquery.dataTables.min')
        else:
            return conf.datatables_url[:-3]

    def ipynb(self, table_id, css=None):
        html = '<style>{0}</style>'.format(css if css is not None
                                           else DEFAULT_CSS_NB)
        html += IPYNB_JS_SCRIPT.format(
            display_length=self.display_length,
            display_length_menu=self.display_length_menu,
            datatables_url=self._jstable_file(),
            tid=table_id)
        return html

    def html_js(self, table_id='table0'):
        return HTML_JS_SCRIPT.format(
            display_length=self.display_length,
            display_length_menu=self.display_length_menu,
            tid=table_id).strip()


def write_table_jsviewer(table, filename, table_id=None, max_lines=5000,
                         table_class="display compact", jskwargs=None,
                         css=DEFAULT_CSS):
    if table_id is None:
        table_id = 'table{id}'.format(id=id(table))

    jskwargs = jskwargs or {}
    jsv = JSViewer(**jskwargs)

    htmldict = {
        'table_id': table_id,
        'table_class': table_class,
        'css': css,
        'cssfiles': jsv.css_urls,
        'jsfiles': jsv.jquery_urls,
        'js':  jsv.html_js(table_id=table_id)
    }

    if max_lines < len(table):
        table = table[:max_lines]
    table.write(filename, format='html', htmldict=htmldict)

io_registry.register_writer('jsviewer', Table, write_table_jsviewer)
