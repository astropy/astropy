# Licensed under a 3-clause BSD style license - see LICENSE.rst

from pathlib import Path
from warnings import warn

import astropy.config as _config
import astropy.io.registry as io_registry
from astropy import extern

from .table import Table


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.table.jsviewer`.
    """

    jquery_url = _config.ConfigItem(
        "https://code.jquery.com/jquery-3.6.0.min.js", "The URL to the jquery library."
    )

    datatables_url = _config.ConfigItem(
        "https://cdn.datatables.net/2.1.8/js/dataTables.min.js",
        "The URL to the jquery datatables library.",
    )

    css_urls = _config.ConfigItem(
        ["https://cdn.datatables.net/2.1.8/css/dataTables.dataTables.min.css"],
        "The URLs to the css file(s) to include.",
        cfgtype="string_list",
    )


conf = Conf()

_TABLE_DIR = Path(extern.__file__).parent
EXTERN_JS_DIR = _TABLE_DIR.joinpath("jquery", "data", "js").resolve()
EXTERN_CSS_DIR = _TABLE_DIR.joinpath("jquery", "data", "css").resolve()

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
""" % dict(  # noqa: UP031
    sorting_script1=_SORTING_SCRIPT_PART_1, sorting_script2=_SORTING_SCRIPT_PART_2
)

HTML_JS_SCRIPT = (
    _SORTING_SCRIPT_PART_1
    + _SORTING_SCRIPT_PART_2
    + """
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
)


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
        if use_local_files:
            warn(
                "`use_local_files` is deprecated and has no effect; for security reasons no static versions of the required js libraries are included in astropy.",
                DeprecationWarning,
            )
        self.display_length_menu = [
            [10, 25, 50, 100, 500, 1000, -1],
            [10, 25, 50, 100, 500, 1000, "All"],
        ]
        self.display_length = display_length
        for L in self.display_length_menu:
            if display_length not in L:
                L.insert(0, display_length)

    @property
    def jquery_urls(self):
        return [conf.jquery_url, conf.datatables_url]

    @property
    def css_urls(self):
        return conf.css_urls

    def _jstable_file(self):
        return conf.datatables_url[:-3]

    def ipynb(self, table_id, css=None, sort_columns="[]"):
        html = f"<style>{css if css is not None else DEFAULT_CSS_NB}</style>"
        html += IPYNB_JS_SCRIPT.format(
            display_length=self.display_length,
            display_length_menu=self.display_length_menu,
            datatables_url=self._jstable_file(),
            tid=table_id,
            sort_columns=sort_columns,
        )
        return html

    def html_js(self, table_id="table0", sort_columns="[]"):
        return HTML_JS_SCRIPT.format(
            display_length=self.display_length,
            display_length_menu=self.display_length_menu,
            tid=table_id,
            sort_columns=sort_columns,
        ).strip()


def write_table_jsviewer(
    table,
    filename,
    table_id=None,
    max_lines=5000,
    table_class="display compact",
    jskwargs=None,
    css=DEFAULT_CSS,
    htmldict=None,
    overwrite=False,
):
    """
    Write an Astropy Table to an HTML file with JavaScript viewer.

    This function uses the JSViewer class to generate the necessary JavaScript
    and CSS for displaying the table interactively in a web browser.

    Parameters
    ----------
    table : Table
        The Astropy Table to be written to an HTML file.
    filename : str, Path
        The name of the output HTML file.
    table_id : str, optional
        The HTML id attribute for the table. Defaults to ``f"table({id(table)}"``.
    max_lines : int, optional
        The maximum number of lines to include in the output table. Default is 5000.
    table_class : str, optional
        The CSS class for the table. Default is "display compact".
    jskwargs : dict, optional
        Additional keyword arguments to pass to the JSViewer.
    css : str, optional
        CSS styles to include in the HTML file. Default is `DEFAULT_CSS`.
    htmldict : dict, optional
        Additional HTML options passed to :class:`~astropy.io.ascii.HTML`.
    overwrite : bool, optional
        If True, overwrite the output file if it exists. Default is False.

    Returns
    -------
    None
    """
    if table_id is None:
        table_id = f"table{id(table)}"

    jskwargs = jskwargs or {}
    jsv = JSViewer(**jskwargs)

    sortable_columns = [
        i
        for i, col in enumerate(table.columns.values())
        if col.info.dtype.kind in "iufc"
    ]
    html_options = {
        "table_id": table_id,
        "table_class": table_class,
        "css": css,
        "cssfiles": jsv.css_urls,
        "jsfiles": jsv.jquery_urls,
        "js": jsv.html_js(table_id=table_id, sort_columns=sortable_columns),
    }
    if htmldict:
        html_options.update(htmldict)

    if max_lines < len(table):
        table = table[:max_lines]
    table.write(filename, format="html", htmldict=html_options, overwrite=overwrite)


io_registry.register_writer("jsviewer", Table, write_table_jsviewer)
