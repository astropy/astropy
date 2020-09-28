# Licensed under a 3-clause BSD style license - see LICENSE.rst

from os.path import abspath, dirname, join

from .table import Table

import astropy.io.registry as io_registry
import astropy.config as _config
from astropy import extern


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.table.jsviewer`.
    """

    jquery_url = _config.ConfigItem(
        'https://code.jquery.com/jquery-3.1.1.min.js',
        'The URL to the jquery library.')

    datatables_url = _config.ConfigItem(
        'https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min.js',
        'The URL to the jquery datatables library.')

    css_urls = _config.ConfigItem(
        ['https://cdn.datatables.net/1.10.12/css/jquery.dataTables.css'],
        'The URLs to the css file(s) to include.', cfgtype='string_list')


conf = Conf()


EXTERN_JS_DIR = abspath(join(dirname(extern.__file__), 'jquery', 'data', 'js'))
EXTERN_CSS_DIR = abspath(join(dirname(extern.__file__), 'jquery', 'data', 'css'))

_SORTING_SCRIPT_PART_1 = """
var astropy_sort_num = function(a, b) {{
    var a_num = parseFloat(a);
    var b_num = parseFloat(b);

    if (isNaN(a_num) && isNaN(b_num))
        return ((a < b) ? -1 : ((a > b) ? 1 : 0));
    else if (!isNaN(a_num) && !isNaN(b_num))
        return ((a_num < b_num) ? -1 : ((a_num > b_num) ? 1 : 0));
    else
        return isNaN(a_num) ? -1 : 1;
}}
"""

_SORTING_SCRIPT_PART_2 = """
jQuery.extend( jQuery.fn.dataTableExt.oSort, {{
    "optionalnum-asc": astropy_sort_num,
    "optionalnum-desc": function (a,b) {{ return -astropy_sort_num(a, b); }}
}});
"""

IPYNB_JS_SCRIPT = """
<script>
%(sorting_script1)s
require.config({{paths: {{
    datatables: '{datatables_url}'
}}}});
require(["datatables"], function(){{
    console.log("$('#{tid}').dataTable()");
    %(sorting_script2)s
    $('#{tid}').dataTable({{
        order: [],
        pageLength: {display_length},
        lengthMenu: {display_length_menu},
        pagingType: "full_numbers",
        columnDefs: [{{targets: {sort_columns}, type: "optionalnum"}}]
    }});
}});
</script>
""" % dict(sorting_script1=_SORTING_SCRIPT_PART_1,
           sorting_script2=_SORTING_SCRIPT_PART_2)

HTML_JS_SCRIPT = _SORTING_SCRIPT_PART_1 + _SORTING_SCRIPT_PART_2 + """
$(document).ready(function() {{
    $('#{tid}').dataTable({{
        order: [],
        pageLength: {display_length},
        lengthMenu: {display_length_menu},
        pagingType: "full_numbers",
        columnDefs: [{{targets: {sort_columns}, type: "optionalnum"}}]
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


class JSViewer:
    """Provides an interactive HTML export of a Table.

    This class provides an interface to the `DataTables
    <https://datatables.net/>`_ library, which allow to visualize interactively
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
            return ['file://' + join(EXTERN_JS_DIR, 'jquery-3.1.1.min.js'),
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

    def ipynb(self, table_id, css=None, sort_columns='[]'):
        html = f'<style>{css if css is not None else DEFAULT_CSS_NB}</style>'
        html += IPYNB_JS_SCRIPT.format(
            display_length=self.display_length,
            display_length_menu=self.display_length_menu,
            datatables_url=self._jstable_file(),
            tid=table_id, sort_columns=sort_columns)
        return html

    def html_js(self, table_id='table0', sort_columns='[]'):
        return HTML_JS_SCRIPT.format(
            display_length=self.display_length,
            display_length_menu=self.display_length_menu,
            tid=table_id, sort_columns=sort_columns).strip()


def write_table_jsviewer(table, filename, table_id=None, max_lines=5000,
                         table_class="display compact", jskwargs=None,
                         css=DEFAULT_CSS, htmldict=None, overwrite=False):
    if table_id is None:
        table_id = f'table{id(table)}'

    jskwargs = jskwargs or {}
    jsv = JSViewer(**jskwargs)

    sortable_columns = [i for i, col in enumerate(table.columns.values())
                        if col.info.dtype.kind in 'iufc']
    html_options = {
        'table_id': table_id,
        'table_class': table_class,
        'css': css,
        'cssfiles': jsv.css_urls,
        'jsfiles': jsv.jquery_urls,
        'js': jsv.html_js(table_id=table_id, sort_columns=sortable_columns)
    }
    if htmldict:
        html_options.update(htmldict)

    if max_lines < len(table):
        table = table[:max_lines]
    table.write(filename, format='html', htmldict=html_options,
                overwrite=overwrite)


io_registry.register_writer('jsviewer', Table, write_table_jsviewer)
