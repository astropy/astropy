# Licensed under a 3-clause BSD style license - see LICENSE.rst
# TODO: need to download some images used by the jquery-ui css file:
# images/ui-icons_888888_256x240.png
import os

from .. import config as _config
from .. import extern


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.table.jsviewer`.
    """

    jquery_url = _config.ConfigItem(
        '',
        'The URL to the jquery library.')

    datatables_url = _config.ConfigItem(
        '',
        'The URL to the jquery datatables library.')

conf = Conf()


JQUERY_URL = _config.ConfigAlias(
    '0.4', 'JQUERY_URL', 'jquery_url',
    'astropy.table.jsviewer', 'astropy.table.jsviewer')

DATATABLES_URL = _config.ConfigAlias(
    '0.4', 'DATATABLES_URL', 'datatables_url',
    'astropy.table.jsviewer', 'astropy.table.jsviewer')


data_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')


ipynb_js_script = """
<script>
    function html_repr_full() {{
        var kernel = IPython.notebook.kernel;
        var button = $("#MakeTableBrowseable{tid}");
        var tablename = button.parents()[4].getElementsByClassName("input_area")[0].innerText;
        tablename = tablename.replace(/\s+/g, '');
        var command = "print ''.join(" + tablename + ".pformat(html=True, max_lines=1000, max_width=1000, tableid={tid}))";
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


commandline_js_script = """
<script>
    $(document).ready(function() {{
        $('#{tid}').dataTable({{
         "iDisplayLength": {display_length},
         "aLengthMenu": {display_length_menu},
         "bJQueryUI": true,
         "sPaginationType": "full_numbers"
        }});
    }} );
</script>
"""


class JSViewer(object):
    def __init__(self,
                 css_files=['jquery-ui.css', 'demo_page.css', 'demo_table.css'],
                 display_length=50):
        self.css_urls = ["file://" + os.path.join(data_path, c) for c in css_files]
        self.display_length_menu = [[10, 25, 50, 100, 500, 1000, -1],
                                    [10, 25, 50, 100, 500, 1000, "All"]]
        self.display_length = display_length
        for L in self.display_length_menu:
            if display_length not in L:
                L.insert(0, display_length)

    def _jquery_file(self):
        jquery_url = conf.jquery_url
        if not jquery_url:
            jquery_url = 'file://' + os.path.abspath(
                os.path.join(
                    os.path.dirname(extern.__file__), 'js', 'jquery-1.11.0.js'))
        return '<script src="{0}"></script>'.format(jquery_url)

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

    def ipynb(self, tableid):
        js = self._css_files()
        js.append(self._jstable_file())
        js.append(ipynb_js_script.format(
            display_length=self.display_length,
            display_length_menu=self.display_length_menu,
            tid=tableid,
            data_path="file://"+data_path))
        return js

    def command_line(self, tableid='table0'):
        js = self._css_files()
        js.append(self._jquery_file())
        js.append(self._jstable_file())

        js.append(commandline_js_script.format(
            display_length=self.display_length,
            display_length_menu=self.display_length_menu,
            tid=tableid))
        return js
