# Licensed under a 3-clause BSD style license - see LICENSE.rst
# TODO: need to download some images used by the jquery-ui css file:
# images/ui-icons_888888_256x240.png

import os
import glob

from .table import Table

from ..io.registry import BaseIO
from .. import config as _config
from .. import extern


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.table.jsviewer`.
    """

    jquery_url = _config.ConfigItem(
        'http://code.jquery.com/jquery-1.11.1.min.js',
        'The URL to the jquery library.')

    datatables_url = _config.ConfigItem(
        'http://cdn.datatables.net/1.10.2/js/jquery.dataTables.min.js',
        'The URL to the jquery datatables library.')

    css_urls = _config.ConfigItem(
        ['https://code.jquery.com/ui/1.11.1/themes/overcast/jquery-ui.css'],
        'The URLs to the css file(s) to include.', cfgtype='list')


conf = Conf()


JQUERY_URL = _config.ConfigAlias(
    '0.4', 'JQUERY_URL', 'jquery_url',
    'astropy.table.jsviewer', 'astropy.table.jsviewer')

DATATABLES_URL = _config.ConfigAlias(
    '0.4', 'DATATABLES_URL', 'datatables_url',
    'astropy.table.jsviewer', 'astropy.table.jsviewer')


DATA_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
EXTERN_JS_DIR = os.path.abspath(os.path.join(os.path.dirname(extern.__file__), 'js'))


IPYNB_JS_SCRIPT = """
<script>
    function html_repr_full() {{
        var kernel = IPython.notebook.kernel;
        var button = $("#MakeTableBrowseable{tid}");
        var tablename = button.parents()[4].getElementsByClassName("input_area")[0].innerText;
        tablename = tablename.replace(/\s+/g, '');
        var command = "print ''.join(" + tablename + ".pformat(html=True, max_lines=1000, max_width=1000, table_id={tid}))";
        console.log(command);
        var result = kernel.execute(command, {{'output': callback}}, {{silent:false}});
        console.log(result);
    }}
    function callback(output_type, out) {{
        console.log(output_type);
        console.log(out);
        var button = $("#MakeTableBrowseable{tid}");
        button[0].parentNode.innerHTML = out.data;
        return out.data;
    }}
    function make_table_browseable() {{
        console.log("$('#{tid}').dataTable()");
        $('#{tid}').dataTable({{
             "iDisplayLength": {display_length},
             "aLengthMenu": {display_length_menu},
             "bJQueryUI": true,
             "sPaginationType": "full_numbers"
        }});
    }}
    function replace_table() {{
        html_repr_full();
        make_table_browseable();
    }}
</script>
<button id='MakeTableBrowseable{tid}' onclick="make_table_browseable()">Make Table Browseable</button>
"""


HTML_JS_SCRIPT = """
$(document).ready(function() {{
    $('#{tid}').dataTable({{
     "iDisplayLength": {display_length},
     "aLengthMenu": {display_length_menu},
     "bJQueryUI": true,
     "sPaginationType": "full_numbers"
    }});
}} );
"""


class JSViewer(object):
    def __init__(self,
                 use_local_files=False,
                 display_length=50):
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
            return ['file://' + os.path.join(EXTERN_JS_DIR, 'jquery-1.11.0.js'),
                    'file://' + os.path.join(EXTERN_JS_DIR, 'jquery.dataTables.js')]
        else:
            return [conf.jquery_url, conf.datatables_url]

    @property
    def css_urls(self):
        if self._use_local_files:
            return ["file://" + filename for filename in glob.glob(os.path.join(DATA_PATH, '*.css'))]
        else:
            return conf.css_urls

    def _jstable_file(self):
        # downloaded from http://datatables.net/download/build/
        datatables_url = conf.datatables_url
        if not datatables_url:
            datatables_url = 'file://' + os.path.abspath(
                os.path.join(
                    os.path.dirname(extern.__file__), 'js',
                    'jquery.dataTables.js'))
        return '<script class="jsbin" src="{0}"></script>'.format(
            datatables_url)

    def _css_files(self):
        return [
            '<link rel="stylesheet" href="{css}" type="text/css">'.format(css=css)
            for css in self.css_urls]

    def ipynb(self, table_id):
        js = self._css_files()
        js.append(self._jstable_file())
        js.append(IPYNB_JS_SCRIPT.format(
            display_length=self.display_length,
            display_length_menu=self.display_length_menu,
            tid=table_id,
            data_path="file://"+DATA_PATH))
        return js

    def html_js(self, table_id='table0'):
        return HTML_JS_SCRIPT.format(display_length=self.display_length,
                                     display_length_menu=self.display_length_menu,
                                     tid=table_id).strip()


class JSViewerTableIO(BaseIO):

    _format_name = 'jsviewer'
    _supported_class = Table

    @staticmethod
    def write(table, filename, table_id=None,
              css="table,th,td,tr,tbody {border: 1px solid black; border-collapse: collapse;}",
              max_lines=5000,
              jskwargs={}):

        if table_id is None:
            table_id = 'table{id}'.format(id=id(table))

        jsv = JSViewer(**jskwargs)

        htmldict = {}
        htmldict['table_id'] = table_id
        htmldict['css'] = css
        htmldict['cssfiles'] = jsv.css_urls
        htmldict['jsfiles'] = jsv.jquery_urls
        htmldict['js'] =  jsv.html_js(table_id=table_id)

        table.write(filename, format='ascii.html', htmldict=htmldict)
