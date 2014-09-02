from ..table import Table

REFERENCE = """
<html>
 <head>
  <meta charset="utf-8"/>
  <meta content="text/html;charset=UTF-8" http-equiv="Content-type"/>
  <style>
table,th,td,tr,tbody {border: 1px solid black; border-collapse: collapse;}  </style>
  <link href="https://code.jquery.com/ui/1.11.1/themes/overcast/jquery-ui.css" rel="stylesheet" type="text/css"/>
  <script src="http://code.jquery.com/jquery-1.11.1.min.js">
  </script>
  <script src="http://cdn.datatables.net/1.10.2/js/jquery.dataTables.min.js">
  </script>
 </head>
 <body>
  <script>
$(document).ready(function() {
    $('#test').dataTable({
     "iDisplayLength": 50,
     "aLengthMenu": [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
     "bJQueryUI": true,
     "sPaginationType": "full_numbers"
    });
} );  </script>
  <table id="test">
   <thead>
    <tr>
     <th>a</th>
     <th>b</th>
    </tr>
   </thead>
   <tr>
    <td>1</td>
    <td>a</td>
   </tr>
   <tr>
    <td>2</td>
    <td>b</td>
   </tr>
   <tr>
    <td>3</td>
    <td>c</td>
   </tr>
  </table>
 </body>
</html>
"""


def test_write_jsviewer(tmpdir):

    t = Table()
    t['a'] = [1,2,3]
    t['b'] = ['a','b','c']
    t['a'].unit = 'm'

    tmpfile = tmpdir.join('test.html').strpath

    t.write(tmpfile, format='jsviewer', table_id='test')

    with open(tmpfile) as f:
        content = f.read()

    assert content.strip() == REFERENCE.strip()
