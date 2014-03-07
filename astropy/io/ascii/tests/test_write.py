# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

import copy

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

from ... import ascii as asciitable

from .common import setup_function, teardown_function

test_defs = [
    dict(kwargs=dict(),
         out="""\
ID XCENTER YCENTER MAG MERR MSKY NITER SHARPNESS CHI PIER PERROR
14 138.538 256.405 15.461 0.003 34.85955 4 -0.032 0.802 0 No_error
18 18.114 280.170 22.329 0.206 30.12784 4 -2.544 1.104 0 No_error
"""
         ),
    dict(kwargs=dict(delimiter=None),
         out="""\
ID XCENTER YCENTER MAG MERR MSKY NITER SHARPNESS CHI PIER PERROR
14 138.538 256.405 15.461 0.003 34.85955 4 -0.032 0.802 0 No_error
18 18.114 280.170 22.329 0.206 30.12784 4 -2.544 1.104 0 No_error
"""
         ),
    dict(kwargs=dict(formats={'XCENTER': '%12.1f',
                              'YCENTER': '{0:.1f}'},
                     include_names=['XCENTER', 'YCENTER'],
                     strip_whitespace=False),
         out="""\
XCENTER YCENTER
"       138.5" 256.4
"        18.1" 280.2
"""
         ),
    dict(kwargs=dict(Writer=asciitable.Rdb, exclude_names=['CHI']),
         out="""\
ID	XCENTER	YCENTER	MAG	MERR	MSKY	NITER	SHARPNESS	PIER	PERROR
N	N	N	N	N	N	N	N	N	S
14	138.538	256.405	15.461	0.003	34.85955	4	-0.032	0	No_error
18	18.114	280.170	22.329	0.206	30.12784	4	-2.544	0	No_error
"""
         ),
    dict(kwargs=dict(Writer=asciitable.Tab),
         out="""\
ID	XCENTER	YCENTER	MAG	MERR	MSKY	NITER	SHARPNESS	CHI	PIER	PERROR
14	138.538	256.405	15.461	0.003	34.85955	4	-0.032	0.802	0	No_error
18	18.114	280.170	22.329	0.206	30.12784	4	-2.544	1.104	0	No_error
"""
         ),
    dict(kwargs=dict(Writer=asciitable.Csv),
         out="""\
ID,XCENTER,YCENTER,MAG,MERR,MSKY,NITER,SHARPNESS,CHI,PIER,PERROR
14,138.538,256.405,15.461,0.003,34.85955,4,-0.032,0.802,0,No_error
18,18.114,280.170,22.329,0.206,30.12784,4,-2.544,1.104,0,No_error
"""
         ),
    dict(kwargs=dict(Writer=asciitable.NoHeader),
         out="""\
14 138.538 256.405 15.461 0.003 34.85955 4 -0.032 0.802 0 No_error
18 18.114 280.170 22.329 0.206 30.12784 4 -2.544 1.104 0 No_error
"""
         ),
    dict(kwargs=dict(Writer=asciitable.CommentedHeader),
         out="""\
# ID XCENTER YCENTER MAG MERR MSKY NITER SHARPNESS CHI PIER PERROR
14 138.538 256.405 15.461 0.003 34.85955 4 -0.032 0.802 0 No_error
18 18.114 280.170 22.329 0.206 30.12784 4 -2.544 1.104 0 No_error
"""
         ),
    dict(kwargs=dict(Writer=asciitable.CommentedHeader, comment='&'),
         out="""\
&ID XCENTER YCENTER MAG MERR MSKY NITER SHARPNESS CHI PIER PERROR
14 138.538 256.405 15.461 0.003 34.85955 4 -0.032 0.802 0 No_error
18 18.114 280.170 22.329 0.206 30.12784 4 -2.544 1.104 0 No_error
"""
         ),
    dict(kwargs=dict(Writer=asciitable.Latex),
         out="""\
\\begin{table}
\\begin{tabular}{ccccccccccc}
ID & XCENTER & YCENTER & MAG & MERR & MSKY & NITER & SHARPNESS & CHI & PIER & PERROR \\\\
14 & 138.538 & 256.405 & 15.461 & 0.003 & 34.85955 & 4 & -0.032 & 0.802 & 0 & No_error \\\\
18 & 18.114 & 280.170 & 22.329 & 0.206 & 30.12784 & 4 & -2.544 & 1.104 & 0 & No_error \\\\
\\end{tabular}
\\end{table}
"""
         ),
    dict(kwargs=dict(Writer=asciitable.AASTex),
         out="""\
\\begin{deluxetable}{ccccccccccc}
\\tablehead{\\colhead{ID} & \\colhead{XCENTER} & \\colhead{YCENTER} & \\colhead{MAG} & \\colhead{MERR} & \\colhead{MSKY} & \\colhead{NITER} & \\colhead{SHARPNESS} & \\colhead{CHI} & \\colhead{PIER} & \\colhead{PERROR}}
\\startdata
14 & 138.538 & 256.405 & 15.461 & 0.003 & 34.85955 & 4 & -0.032 & 0.802 & 0 & No_error \\\\
18 & 18.114 & 280.170 & 22.329 & 0.206 & 30.12784 & 4 & -2.544 & 1.104 & 0 & No_error \\\\
\\enddata
\\end{deluxetable}
"""
         ),
    dict(
        kwargs=dict(Writer=asciitable.AASTex, caption='Mag values \\label{tab1}', latexdict={
                    'units': {'MAG': '[mag]', 'XCENTER': '[pixel]'}, 'tabletype': 'deluxetable*',
                    'tablealign':'htpb'}),
        out="""\
\\begin{deluxetable*}{ccccccccccc}[htpb]
\\tablecaption{Mag values \\label{tab1}}
\\tablehead{\\colhead{ID} & \\colhead{XCENTER} & \\colhead{YCENTER} & \\colhead{MAG} & \\colhead{MERR} & \\colhead{MSKY} & \\colhead{NITER} & \\colhead{SHARPNESS} & \\colhead{CHI} & \\colhead{PIER} & \\colhead{PERROR}\\\\ \\colhead{ } & \\colhead{[pixel]} & \\colhead{ } & \\colhead{[mag]} & \\colhead{ } & \\colhead{ } & \\colhead{ } & \\colhead{ } & \\colhead{ } & \\colhead{ } & \\colhead{ }}
\\startdata
14 & 138.538 & 256.405 & 15.461 & 0.003 & 34.85955 & 4 & -0.032 & 0.802 & 0 & No_error \\\\
18 & 18.114 & 280.170 & 22.329 & 0.206 & 30.12784 & 4 & -2.544 & 1.104 & 0 & No_error \\\\
\\enddata
\\end{deluxetable*}
"""
    ),
    dict(
        kwargs=dict(Writer=asciitable.Latex, caption='Mag values \\label{tab1}', latexdict={'preamble': '\\begin{center}', 'tablefoot': '\\end{center}', 'data_end': [
                    '\\hline', '\\hline'], 'units':{'MAG': '[mag]', 'XCENTER': '[pixel]'},
                    'tabletype': 'table*',
                    'tablealign': 'h'},
                    col_align='|lcccccccccc|'),
        out="""\
\\begin{table*}[h]
\\begin{center}
\\caption{Mag values \\label{tab1}}
\\begin{tabular}{|lcccccccccc|}
ID & XCENTER & YCENTER & MAG & MERR & MSKY & NITER & SHARPNESS & CHI & PIER & PERROR \\\\
 & [pixel] &  & [mag] &  &  &  &  &  &  &  \\\\
14 & 138.538 & 256.405 & 15.461 & 0.003 & 34.85955 & 4 & -0.032 & 0.802 & 0 & No_error \\\\
18 & 18.114 & 280.170 & 22.329 & 0.206 & 30.12784 & 4 & -2.544 & 1.104 & 0 & No_error \\\\
\\hline
\\hline
\\end{tabular}
\\end{center}
\\end{table*}
"""
    ),
    dict(kwargs=dict(Writer=asciitable.Latex, latexdict=asciitable.latexdicts['template']),
         out="""\
\\begin{tabletype}[tablealign]
preamble
\\caption{caption}
\\begin{tabular}{col_align}
header_start
ID & XCENTER & YCENTER & MAG & MERR & MSKY & NITER & SHARPNESS & CHI & PIER & PERROR \\\\
 &  &  &  &  &  &  &  &  &  &  \\\\
header_end
data_start
14 & 138.538 & 256.405 & 15.461 & 0.003 & 34.85955 & 4 & -0.032 & 0.802 & 0 & No_error \\\\
18 & 18.114 & 280.170 & 22.329 & 0.206 & 30.12784 & 4 & -2.544 & 1.104 & 0 & No_error \\\\
data_end
\\end{tabular}
tablefoot
\\end{tabletype}
"""
         ),
    dict(kwargs=dict(Writer=asciitable.HTML, htmldict={'css':'table,th,td{border:1px solid black;'}),
         out="""\
<html>
 <head>
  <meta charset="utf-8"/>
  <meta content="text/html;charset=UTF-8" http-equiv="Content-type"/>
  <style>
table,th,td{border:1px solid black;  </style>
 </head>
 <body>
  <table>
   <tr>
    <th>
ID    </th>
    <th>
XCENTER    </th>
    <th>
YCENTER    </th>
    <th>
MAG    </th>
    <th>
MERR    </th>
    <th>
MSKY    </th>
    <th>
NITER    </th>
    <th>
SHARPNESS    </th>
    <th>
CHI    </th>
    <th>
PIER    </th>
    <th>
PERROR    </th>
   </tr>
   <tr>
    <td>
14           </td>
    <td>
138.538       </td>
    <td>
256.405       </td>
    <td>
15.461          </td>
    <td>
0.003             </td>
    <td>
34.85955           </td>
    <td>
4         </td>
    <td>
-0.032                     </td>
    <td>
0.802           </td>
    <td>
0         </td>
    <td>
No_error         </td>
   </tr>
   <tr>
    <td>
18           </td>
    <td>
18.114        </td>
    <td>
280.170       </td>
    <td>
22.329          </td>
    <td>
0.206             </td>
    <td>
30.12784           </td>
    <td>
4         </td>
    <td>
-2.544                     </td>
    <td>
1.104           </td>
    <td>
0         </td>
    <td>
No_error         </td>
   </tr>
  </table>
 </body>
</html>
"""
         ),
    dict(kwargs=dict(Writer=asciitable.Ipac),
         out="""\
\\MERGERAD='INDEF'
\\IRAF='NOAO/IRAFV2.10EXPORT'
\\USER=''
\\HOST='tucana'
\\DATE='05-28-93'
\\TIME='14:46:13'
\\PACKAGE='daophot'
\\TASK='nstar'
\\IMAGE='test'
\\GRPFILE='test.psg.1'
\\PSFIMAGE='test.psf.1'
\\NSTARFILE='test.nst.1'
\\REJFILE='"hello world"'
\\SCALE='1.'
\\DATAMIN='50.'
\\DATAMAX='24500.'
\\GAIN='1.'
\\READNOISE='0.'
\\OTIME='00:07:59.0'
\\XAIRMASS='1.238106'
\\IFILTER='V'
\\RECENTER='yes'
\\FITSKY='no'
\\PSFMAG='16.594'
\\PSFRAD='5.'
\\FITRAD='3.'
\\MAXITER='50'
\\MAXGROUP='60'
\\FLATERROR='0.75'
\\PROFERROR='5.'
\\CLIPEXP='6'
\\CLIPRANGE='2.5'
|       ID|   XCENTER|   YCENTER|         MAG|          MERR|           MSKY| NITER|              SHARPNESS|         CHI|  PIER|       PERROR|
|     long|    double|    double|      double|        double|         double|  long|                 double|      double|  long|         char|
|         |    pixels|    pixels|  magnitudes|    magnitudes|         counts|      |                       |            |      |      perrors|
|     null|      null|      null|        null|          null|           null|  null|                   null|        null|  null|         null|
 14        138.538    256.405    15.461       0.003          34.85955        4      -0.032                  0.802        0      No_error
 18        18.114     280.170    22.329       0.206          30.12784        4      -2.544                  1.104        0      No_error
"""
         ),
]

tab_to_fill = ['a b c', '1 2 3', '1 1 3']

test_defs_fill_value = [
    dict(kwargs=dict(),
         out="""\
a b c
1 2 3
1 1 3
"""
         ),
    dict(kwargs=dict(fill_values=('1', 'w')),
         out="""\
a b c
w 2 3
w w 3
"""
         ),
    dict(kwargs=dict(fill_values=('1', 'w', 'b')),
         out="""\
a b c
1 2 3
1 w 3
"""
         ),
    dict(kwargs=dict(fill_values = ('1', 'w'),
                     fill_include_names = ['b']),
         out="""\
a b c
1 2 3
1 w 3
"""
         ),
    dict(kwargs=dict(fill_values = ('1', 'w'),
                     fill_exclude_names = ['a']),
         out="""\
a b c
1 2 3
1 w 3
"""
         ),
    dict(kwargs=dict(fill_values = ('1', 'w'),
                     fill_include_names = ['a'],
                     fill_exclude_names = ['a', 'b']),
         out="""\
a b c
1 2 3
1 1 3
"""
         ),
    dict(kwargs=dict(fill_values = [('1', 'w')],
                     formats={'a': '%4.2f'}),
         out="""\
a b c
1.00 2 3
1.00 w 3
"""
         ),
]


def check_write_table(test_def, table):
    out = StringIO()
    asciitable.write(table, out, **test_def['kwargs'])
    print('Expected:\n%s' % test_def['out'])
    print('Actual:\n%s' % out.getvalue())
    assert [x.strip() for x in out.getvalue().strip().splitlines()] == [
        x.strip() for x in test_def['out'].strip().splitlines()]


def check_write_table_via_table(test_def, table):
    out = StringIO()

    test_def = copy.deepcopy(test_def)
    if 'Writer' in test_def['kwargs']:
        format = 'ascii.{0}'.format(test_def['kwargs']['Writer']._format_name)
        del test_def['kwargs']['Writer']
    else:
        format = 'ascii'

    table.write(out, format=format, **test_def['kwargs'])
    print('Expected:\n%s' % test_def['out'])
    print('Actual:\n%s' % out.getvalue())
    assert [x.strip() for x in out.getvalue().strip().splitlines()] == [
        x.strip() for x in test_def['out'].strip().splitlines()]


def test_write_table():
    table = asciitable.get_reader(Reader=asciitable.Daophot)
    data = table.read('t/daophot.dat')

    for test_def in test_defs:
        check_write_table(test_def, data)
        check_write_table_via_table(test_def, data)


def test_write_fill_values():
    data = asciitable.read(tab_to_fill)

    for test_def in test_defs_fill_value:
        check_write_table(test_def, data)
