# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

import os
import copy

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

from ... import ascii
from .... import table
from ....tests.helper import pytest
from .... import units

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
    dict(kwargs=dict(Writer=ascii.Rdb, exclude_names=['CHI']),
         out="""\
ID	XCENTER	YCENTER	MAG	MERR	MSKY	NITER	SHARPNESS	PIER	PERROR
N	N	N	N	N	N	N	N	N	S
14	138.538	256.405	15.461	0.003	34.85955	4	-0.032	0	No_error
18	18.114	280.170	22.329	0.206	30.12784	4	-2.544	0	No_error
"""
         ),
    dict(kwargs=dict(Writer=ascii.Tab),
         out="""\
ID	XCENTER	YCENTER	MAG	MERR	MSKY	NITER	SHARPNESS	CHI	PIER	PERROR
14	138.538	256.405	15.461	0.003	34.85955	4	-0.032	0.802	0	No_error
18	18.114	280.170	22.329	0.206	30.12784	4	-2.544	1.104	0	No_error
"""
         ),
    dict(kwargs=dict(Writer=ascii.Csv),
         out="""\
ID,XCENTER,YCENTER,MAG,MERR,MSKY,NITER,SHARPNESS,CHI,PIER,PERROR
14,138.538,256.405,15.461,0.003,34.85955,4,-0.032,0.802,0,No_error
18,18.114,280.170,22.329,0.206,30.12784,4,-2.544,1.104,0,No_error
"""
         ),
    dict(kwargs=dict(Writer=ascii.NoHeader),
         out="""\
14 138.538 256.405 15.461 0.003 34.85955 4 -0.032 0.802 0 No_error
18 18.114 280.170 22.329 0.206 30.12784 4 -2.544 1.104 0 No_error
"""
         ),
    dict(kwargs=dict(Writer=ascii.CommentedHeader),
         out="""\
# ID XCENTER YCENTER MAG MERR MSKY NITER SHARPNESS CHI PIER PERROR
14 138.538 256.405 15.461 0.003 34.85955 4 -0.032 0.802 0 No_error
18 18.114 280.170 22.329 0.206 30.12784 4 -2.544 1.104 0 No_error
"""
         ),
    dict(kwargs=dict(Writer=ascii.CommentedHeader, comment='&'),
         out="""\
&ID XCENTER YCENTER MAG MERR MSKY NITER SHARPNESS CHI PIER PERROR
14 138.538 256.405 15.461 0.003 34.85955 4 -0.032 0.802 0 No_error
18 18.114 280.170 22.329 0.206 30.12784 4 -2.544 1.104 0 No_error
"""
         ),
    dict(kwargs=dict(Writer=ascii.Latex),
         out="""\
\\begin{table}
\\begin{tabular}{ccccccccccc}
ID & XCENTER & YCENTER & MAG & MERR & MSKY & NITER & SHARPNESS & CHI & PIER & PERROR \\\\
 & pixels & pixels & magnitudes & magnitudes & counts &  &  &  &  & perrors \\\\
14 & 138.538 & 256.405 & 15.461 & 0.003 & 34.85955 & 4 & -0.032 & 0.802 & 0 & No_error \\\\
18 & 18.114 & 280.170 & 22.329 & 0.206 & 30.12784 & 4 & -2.544 & 1.104 & 0 & No_error \\\\
\\end{tabular}
\\end{table}
"""
         ),
    dict(kwargs=dict(Writer=ascii.AASTex),
         out="""\
\\begin{deluxetable}{ccccccccccc}
\\tablehead{\\colhead{ID} & \\colhead{XCENTER} & \\colhead{YCENTER} & \\colhead{MAG} & \\colhead{MERR} & \\colhead{MSKY} & \\colhead{NITER} & \\colhead{SHARPNESS} & \\colhead{CHI} & \\colhead{PIER} & \\colhead{PERROR}\\\\ \\colhead{ } & \\colhead{pixels} & \\colhead{pixels} & \\colhead{magnitudes} & \\colhead{magnitudes} & \\colhead{counts} & \\colhead{ } & \\colhead{ } & \\colhead{ } & \\colhead{ } & \\colhead{perrors}}
\\startdata
14 & 138.538 & 256.405 & 15.461 & 0.003 & 34.85955 & 4 & -0.032 & 0.802 & 0 & No_error \\\\
18 & 18.114 & 280.170 & 22.329 & 0.206 & 30.12784 & 4 & -2.544 & 1.104 & 0 & No_error \\\\
\\enddata
\\end{deluxetable}
"""
         ),
    dict(
        kwargs=dict(Writer=ascii.AASTex, caption='Mag values \\label{tab1}', latexdict={
                    'units': {'MAG': '[mag]', 'XCENTER': '[pixel]'}, 'tabletype': 'deluxetable*',
                    'tablealign': 'htpb'}),
        out="""\
\\begin{deluxetable*}{ccccccccccc}[htpb]
\\tablecaption{Mag values \\label{tab1}}
\\tablehead{\\colhead{ID} & \\colhead{XCENTER} & \\colhead{YCENTER} & \\colhead{MAG} & \\colhead{MERR} & \\colhead{MSKY} & \\colhead{NITER} & \\colhead{SHARPNESS} & \\colhead{CHI} & \\colhead{PIER} & \\colhead{PERROR}\\\\ \\colhead{ } & \\colhead{[pixel]} & \\colhead{pixels} & \\colhead{[mag]} & \\colhead{magnitudes} & \\colhead{counts} & \\colhead{ } & \\colhead{ } & \\colhead{ } & \\colhead{ } & \\colhead{perrors}}
\\startdata
14 & 138.538 & 256.405 & 15.461 & 0.003 & 34.85955 & 4 & -0.032 & 0.802 & 0 & No_error \\\\
18 & 18.114 & 280.170 & 22.329 & 0.206 & 30.12784 & 4 & -2.544 & 1.104 & 0 & No_error \\\\
\\enddata
\\end{deluxetable*}
"""
    ),
    dict(
        kwargs=dict(Writer=ascii.Latex, caption='Mag values \\label{tab1}',
                    latexdict={'preamble': '\\begin{center}', 'tablefoot': '\\end{center}',
                               'data_end': ['\\hline', '\\hline'],
                               'units':{'MAG': '[mag]', 'XCENTER': '[pixel]'},
                    'tabletype': 'table*',
                    'tablealign': 'h'},
                    col_align='|lcccccccccc|'),
        out="""\
\\begin{table*}[h]
\\begin{center}
\\caption{Mag values \\label{tab1}}
\\begin{tabular}{|lcccccccccc|}
ID & XCENTER & YCENTER & MAG & MERR & MSKY & NITER & SHARPNESS & CHI & PIER & PERROR \\\\
 & [pixel] & pixels & [mag] & magnitudes & counts &  &  &  &  & perrors \\\\
14 & 138.538 & 256.405 & 15.461 & 0.003 & 34.85955 & 4 & -0.032 & 0.802 & 0 & No_error \\\\
18 & 18.114 & 280.170 & 22.329 & 0.206 & 30.12784 & 4 & -2.544 & 1.104 & 0 & No_error \\\\
\\hline
\\hline
\\end{tabular}
\\end{center}
\\end{table*}
"""
    ),
    dict(kwargs=dict(Writer=ascii.Latex, latexdict=ascii.latexdicts['template']),
         out="""\
\\begin{tabletype}[tablealign]
preamble
\\caption{caption}
\\begin{tabular}{col_align}
header_start
ID & XCENTER & YCENTER & MAG & MERR & MSKY & NITER & SHARPNESS & CHI & PIER & PERROR \\\\
 & pixels & pixels & magnitudes & magnitudes & counts &  &  &  &  & perrors \\\\
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
    dict(kwargs=dict(Writer=ascii.HTML, htmldict={'css': 'table,th,td{border:1px solid black;'}),
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
   <thead>
    <tr>
     <th>ID</th>
     <th>XCENTER</th>
     <th>YCENTER</th>
     <th>MAG</th>
     <th>MERR</th>
     <th>MSKY</th>
     <th>NITER</th>
     <th>SHARPNESS</th>
     <th>CHI</th>
     <th>PIER</th>
     <th>PERROR</th>
    </tr>
   </thead>
   <tr>
    <td>14</td>
    <td>138.538</td>
    <td>256.405</td>
    <td>15.461</td>
    <td>0.003</td>
    <td>34.85955</td>
    <td>4</td>
    <td>-0.032</td>
    <td>0.802</td>
    <td>0</td>
    <td>No_error</td>
   </tr>
   <tr>
    <td>18</td>
    <td>18.114</td>
    <td>280.170</td>
    <td>22.329</td>
    <td>0.206</td>
    <td>30.12784</td>
    <td>4</td>
    <td>-2.544</td>
    <td>1.104</td>
    <td>0</td>
    <td>No_error</td>
   </tr>
  </table>
 </body>
</html>
"""
         ),
    dict(kwargs=dict(Writer=ascii.Ipac),
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

test_defs_no_data = [
    dict(kwargs=dict(Writer=ascii.Ipac),
         out="""\
\\ This is an example of a valid comment.
\\ The 2nd data line is used to verify the exact column parsing
\\ (unclear if this is a valid for the IPAC format)
\\catalog='sao'
\\date='Wed Sp 20 09:48:36 1995'
\\mykeyword='Another way for defining keyvalue string'
|    ra|   dec| sai|    v2|sptype|
|double|double|long|double|  char|
|  unit|  unit|unit|  unit|  ergs|
|  null|  null|null|  null|  null|
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

test_def_masked_fill_value = [
    dict(kwargs=dict(),
         out="""\
a b c
-- 2 3
1 1 --
"""
         ),
    dict(kwargs=dict(fill_values=[('1', 'w'), (ascii.masked, 'X')]),
         out="""\
a b c
X 2 3
w w X
"""
         ),
    dict(kwargs=dict(fill_values=[('1', 'w'), (ascii.masked, 'XXX')],
                     formats={'a': '%4.1f'}),
         out="""\
a b c
XXX 2 3
1.0 w XXX
"""
         ),
    dict(kwargs=dict(Writer=ascii.Csv),
         out="""\
a,b,c
,2,3
1,1,
"""
         ),
]


def check_write_table(test_def, table, fast_writer):
    out = StringIO()
    try:
        ascii.write(table, out, fast_writer=fast_writer, **test_def['kwargs'])
    except ValueError as e: # if format doesn't have a fast writer, ignore
        if not 'not in the list of formats with fast writers' in str(e):
            raise e
        return
    print('Expected:\n%s' % test_def['out'])
    print('Actual:\n%s' % out.getvalue())
    assert [x.strip() for x in out.getvalue().strip().splitlines()] == [
        x.strip() for x in test_def['out'].strip().splitlines()]


def check_write_table_via_table(test_def, table, fast_writer):
    out = StringIO()

    test_def = copy.deepcopy(test_def)
    if 'Writer' in test_def['kwargs']:
        format = 'ascii.{0}'.format(test_def['kwargs']['Writer']._format_name)
        del test_def['kwargs']['Writer']
    else:
        format = 'ascii'

    try:
        table.write(out, format=format,  fast_writer=fast_writer, **test_def['kwargs'])
    except ValueError as e: # if format doesn't have a fast writer, ignore
        if not 'not in the list of formats with fast writers' in str(e):
            raise e
        return
    print('Expected:\n%s' % test_def['out'])
    print('Actual:\n%s' % out.getvalue())
    assert [x.strip() for x in out.getvalue().strip().splitlines()] == [
        x.strip() for x in test_def['out'].strip().splitlines()]


@pytest.mark.parametrize("fast_writer", [True, False])
def test_write_table(fast_writer):
    table = ascii.get_reader(Reader=ascii.Daophot)
    data = table.read('t/daophot.dat')

    for test_def in test_defs:
        check_write_table(test_def, data, fast_writer)
        check_write_table_via_table(test_def, data, fast_writer)


@pytest.mark.parametrize("fast_writer", [True, False])
def test_write_fill_values(fast_writer):
    data = ascii.read(tab_to_fill)

    for test_def in test_defs_fill_value:
        check_write_table(test_def, data, fast_writer)


@pytest.mark.parametrize("fast_writer", [True, False])
def test_write_fill_masked_different(fast_writer):
    '''see discussion in #2255'''
    data = ascii.read(tab_to_fill)
    data = table.Table(data, masked=True)
    data['a'].mask = [True, False]
    data['c'].mask = [False, True]

    for test_def in test_def_masked_fill_value:
        check_write_table(test_def, data, fast_writer)


@pytest.mark.parametrize("fast_writer", [True, False])
def test_write_no_data_ipac(fast_writer):
    """Write an IPAC table that contains no data."""
    table = ascii.get_reader(Reader=ascii.Ipac)
    data = table.read('t/no_data_ipac.dat')

    for test_def in test_defs_no_data:
        check_write_table(test_def, data, fast_writer)
        check_write_table_via_table(test_def, data, fast_writer)

@pytest.mark.parametrize("fast_writer", [True, False])
def test_write_comments(fast_writer):
    """Write comments in output originally read by io.ascii."""
    data = ascii.read('#c1\n  # c2\t\na,b,c\n#  c3\n1,2,3')
    out = StringIO()
    ascii.write(data, out, format='basic', fast_writer=fast_writer)
    expected = ['# c1', '# c2', '# c3', 'a b c', '1 2 3']
    assert out.getvalue().splitlines() == expected
    # header comes before comments for commented-header
    out = StringIO()
    ascii.write(data, out, format='commented_header', fast_writer=fast_writer)
    expected = ['# a b c', '# c1', '# c2', '# c3', '1 2 3']
    assert out.getvalue().splitlines() == expected

@pytest.mark.parametrize("fast_writer", [True, False])
def test_strip_names(fast_writer):
    """Names should be stripped of whitespace by default."""
    data = table.Table([[1], [2], [3]], names=(' A', 'B ', ' C '))
    out = StringIO()
    ascii.write(data, out, format='csv', fast_writer=fast_writer)
    assert out.getvalue().splitlines()[0] == 'A,B,C'


@pytest.mark.skipif(os.environ.get('APPVEYOR'), reason="fails on AppVeyor")
def test_latex_units():
    """
    Check to make sure that Latex and AASTex writers attempt to fall
    back on the **unit** attribute of **Column** if the supplied
    **latexdict** does not specify units.
    """
    t = table.Table([table.Column(name='date', data=['a','b']),
               table.Column(name='NUV exp.time', data=[1,2])])
    latexdict = copy.deepcopy(ascii.latexdicts['AA'])
    latexdict['units'] = {'NUV exp.time':'s'}
    out = StringIO()
    expected = '''\
\\begin{table}{cc}
\\tablehead{\\colhead{date} & \\colhead{NUV exp.time}\\\\ \\colhead{ } & \\colhead{s}}
\\startdata
a & 1 \\\\
b & 2 \\\\
\\enddata
\\end{table}
'''
    ascii.write(t, out, format='aastex', latexdict=latexdict)
    assert out.getvalue() == expected
    # use unit attribute instead
    t['NUV exp.time'].unit = units.s
    t['date'].unit = units.yr
    out = StringIO()
    ascii.write(t, out, format='aastex', latexdict=ascii.latexdicts['AA'])
    assert out.getvalue() == expected.replace(
        'colhead{s}', 'colhead{$\mathrm{s}$}').replace(
        'colhead{ }', 'colhead{$\mathrm{yr}$}')

