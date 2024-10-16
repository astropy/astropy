import textwrap
from os.path import abspath, dirname, join

import pytest

from astropy import extern
from astropy.coordinates import SkyCoord
from astropy.table.table import Table
from astropy.time import Time
from astropy.utils.compat.optional_deps import HAS_BLEACH, HAS_IPYDATAGRID, HAS_IPYTHON
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.utils.misc import _NOT_OVERWRITING_MSG_MATCH

EXTERN_DIR = abspath(join(dirname(extern.__file__), "jquery", "data"))
JQUERY_MIN_JS = "jquery-3.6.0.min.js"


REFERENCE = """
<html>
 <head>
  <meta charset="utf-8"/>
  <meta content="text/html;charset=UTF-8" http-equiv="Content-type"/>
  <style>
body {font-family: sans-serif;}
table.dataTable {width: auto !important; margin: 0 !important;}
.dataTables_filter, .dataTables_paginate {float: left !important; margin-left:1em}
  </style>
  <link href="%(datatables_css_url)s" rel="stylesheet" type="text/css"/>
  <script src="%(jquery_url)s">
  </script>
  <script src="%(datatables_js_url)s">
  </script>
 </head>
 <body>
  <script>
var astropy_sort_num = function(a, b) {
    var a_num = parseFloat(a);
    var b_num = parseFloat(b);

    if (isNaN(a_num) && isNaN(b_num))
        return ((a < b) ? -1 : ((a > b) ? 1 : 0));
    else if (!isNaN(a_num) && !isNaN(b_num))
        return ((a_num < b_num) ? -1 : ((a_num > b_num) ? 1 : 0));
    else
        return isNaN(a_num) ? -1 : 1;
}

jQuery.extend( jQuery.fn.dataTableExt.oSort, {
    "optionalnum-asc": astropy_sort_num,
    "optionalnum-desc": function (a,b) { return -astropy_sort_num(a, b); }
});

$(document).ready(function() {
    $('#%(table_id)s').dataTable({
        order: [],
        pageLength: %(length)s,
        lengthMenu: [[%(display_length)s, -1], [%(display_length)s, 'All']],
        pagingType: "full_numbers",
        columnDefs: [{targets: [0], type: "optionalnum"}]
    });
} );  </script>
  <table class="%(table_class)s" id="%(table_id)s">
   <thead>
    <tr>
     <th>a</th>
     <th>b</th>
    </tr>
   </thead>
%(lines)s
  </table>
 </body>
</html>
"""

TPL = "   <tr>\n    <td>{0}</td>\n    <td>{1}</td>\n   </tr>"


def format_lines(col1, col2):
    col1_format = getattr(col1.info, "default_format", lambda x: x)
    col2_format = getattr(col2.info, "default_format", lambda x: x)
    return "\n".join(
        TPL.format(col1_format(v1), col2_format(v2)) for v1, v2 in zip(col1, col2)
    )


def test_write_jsviewer_default(tmp_path):
    t = Table()
    t["a"] = [1, 2, 3, 4, 5]
    t["b"] = ["a", "b", "c", "d", "e"]
    t["a"].unit = "m"

    tmpfile = tmp_path / "test.html"

    t.write(tmpfile, format="jsviewer")
    ref = REFERENCE % dict(
        lines=format_lines(t["a"], t["b"]),
        table_class="display compact",
        table_id=f"table{id(t)}",
        length="50",
        display_length="10, 25, 50, 100, 500, 1000",
        datatables_css_url=(
            "https://cdn.datatables.net/1.10.12/css/jquery.dataTables.css"
        ),
        datatables_js_url=(
            "https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min.js"
        ),
        jquery_url="https://code.jquery.com/" + JQUERY_MIN_JS,
    )
    with open(tmpfile) as f:
        assert f.read().strip() == ref.strip()


def test_write_jsviewer_overwrite(tmp_path):
    t = Table()
    t["a"] = [1, 2, 3, 4, 5]
    t["b"] = ["a", "b", "c", "d", "e"]
    t["a"].unit = "m"
    tmpfile = tmp_path / "test.html"

    # normal write
    t.write(tmpfile, format="jsviewer")
    # errors on overwrite
    with pytest.raises(OSError, match=_NOT_OVERWRITING_MSG_MATCH):
        t.write(tmpfile, format="jsviewer")
    # unless specified
    t.write(tmpfile, format="jsviewer", overwrite=True)


@pytest.mark.parametrize(
    "mixin",
    [
        Time(["J2000", "J2001"]),
        Time([50000.0, 50001.0001], format="mjd"),
        SkyCoord(ra=[100.0, 110.0], dec=[-10.0, 10.0], unit="deg"),
    ],
)
def test_write_jsviewer_mixin(tmp_path, mixin):
    t = Table()
    t["a"] = [1, 2]
    t["b"] = mixin
    t["a"].unit = "m"

    tmpfile = tmp_path / "test.html"

    t.write(tmpfile, format="jsviewer")
    ref = REFERENCE % dict(
        lines=format_lines(t["a"], t["b"]),
        table_class="display compact",
        table_id=f"table{id(t)}",
        length="50",
        display_length="10, 25, 50, 100, 500, 1000",
        datatables_css_url=(
            "https://cdn.datatables.net/1.10.12/css/jquery.dataTables.css"
        ),
        datatables_js_url=(
            "https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min.js"
        ),
        jquery_url="https://code.jquery.com/" + JQUERY_MIN_JS,
    )
    with open(tmpfile) as f:
        assert f.read().strip() == ref.strip()


@pytest.mark.skipif(not HAS_BLEACH, reason="requires bleach")
def test_write_jsviewer_options(tmp_path):
    t = Table()
    t["a"] = [1, 2, 3, 4, 5]
    t["b"] = ["<b>a</b>", "b", "c", "d", "e"]
    t["a"].unit = "m"

    tmpfile = tmp_path / "test.html"
    t.write(
        tmpfile,
        format="jsviewer",
        table_id="test",
        max_lines=3,
        jskwargs={"display_length": 5},
        table_class="display hover",
        htmldict=dict(raw_html_cols="b"),
    )

    ref = REFERENCE % dict(
        lines=format_lines(t["a"][:3], t["b"][:3]),
        table_class="display hover",
        table_id="test",
        length="5",
        display_length="5, 10, 25, 50, 100, 500, 1000",
        datatables_css_url=(
            "https://cdn.datatables.net/1.10.12/css/jquery.dataTables.css"
        ),
        datatables_js_url=(
            "https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min.js"
        ),
        jquery_url="https://code.jquery.com/" + JQUERY_MIN_JS,
    )
    with open(tmpfile) as f:
        assert f.read().strip() == ref.strip()


def test_write_jsviewer_local(tmp_path):
    t = Table()
    t["a"] = [1, 2, 3, 4, 5]
    t["b"] = ["a", "b", "c", "d", "e"]
    t["a"].unit = "m"

    tmpfile = tmp_path / "test.html"

    t.write(
        tmpfile, format="jsviewer", table_id="test", jskwargs={"use_local_files": True}
    )
    ref = REFERENCE % dict(
        lines=format_lines(t["a"], t["b"]),
        table_class="display compact",
        table_id="test",
        length="50",
        display_length="10, 25, 50, 100, 500, 1000",
        datatables_css_url="file://" + join(EXTERN_DIR, "css", "jquery.dataTables.css"),
        datatables_js_url="file://"
        + join(EXTERN_DIR, "js", "jquery.dataTables.min.js"),
        jquery_url="file://" + join(EXTERN_DIR, "js", JQUERY_MIN_JS),
    )
    with open(tmpfile) as f:
        assert f.read().strip() == ref.strip()


@pytest.mark.skipif(not HAS_IPYTHON, reason="requires IPython")
def test_show_in_notebook_classic():
    t = Table()
    t["a"] = [1, 2, 3, 4, 5]
    t["b"] = ["b", "c", "a", "d", "e"]

    with pytest.warns(AstropyDeprecationWarning):
        htmlstr_windx = t.show_in_notebook(
            backend="classic"
        ).data  # should default to 'idx'
        htmlstr_windx_named = t.show_in_notebook(
            backend="classic", show_row_index="realidx"
        ).data
        htmlstr_woindx = t.show_in_notebook(
            backend="classic", show_row_index=False
        ).data

    assert (
        textwrap.dedent(
            """
    <thead><tr><th>idx</th><th>a</th><th>b</th></tr></thead>
    <tr><td>0</td><td>1</td><td>b</td></tr>
    <tr><td>1</td><td>2</td><td>c</td></tr>
    <tr><td>2</td><td>3</td><td>a</td></tr>
    <tr><td>3</td><td>4</td><td>d</td></tr>
    <tr><td>4</td><td>5</td><td>e</td></tr>
    """
        ).strip()
        in htmlstr_windx
    )

    assert (
        "<thead><tr><th>realidx</th><th>a</th><th>b</th></tr></thead>"
        in htmlstr_windx_named
    )

    assert "<thead><tr><th>a</th><th>b</th></tr></thead>" in htmlstr_woindx


@pytest.mark.skipif(not HAS_IPYDATAGRID, reason="requires ipydatagrid")
# https://github.com/bqplot/bqplot/issues/1624 and such
@pytest.mark.filterwarnings(r"ignore:((.|\n)*)traitlets((.|\n)*):DeprecationWarning")
def test_show_in_notebook_ipydatagrid():
    from ipydatagrid import DataGrid

    t = Table()
    dg = t.show_in_notebook()
    assert isinstance(dg, DataGrid)


def test_show_in_notebook_invalid_backend():
    t = Table()
    with pytest.raises(NotImplementedError, match=".* backend is not supported"):
        t.show_in_notebook(backend="foo")
