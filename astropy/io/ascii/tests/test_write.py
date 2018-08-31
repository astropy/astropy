# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

import os
import copy
from itertools import chain

import pytest
import numpy as np

from ....extern.six.moves import cStringIO as StringIO
from ... import ascii
from .... import table
from ....table.table_helpers import simple_table
from ....tests.helper import catch_warnings
from ....utils.exceptions import AstropyWarning, AstropyDeprecationWarning
from .... import units

from .common import setup_function, teardown_function

# Check to see if the BeautifulSoup dependency is present.
try:
    from bs4 import BeautifulSoup, FeatureNotFound
    HAS_BEAUTIFUL_SOUP = True
except ImportError:
    HAS_BEAUTIFUL_SOUP = False

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
ID\tXCENTER\tYCENTER\tMAG\tMERR\tMSKY\tNITER\tSHARPNESS\tPIER\tPERROR
N\tN\tN\tN\tN\tN\tN\tN\tN\tS
14\t138.538\t256.405\t15.461\t0.003\t34.85955\t4\t-0.032\t0\tNo_error
18\t18.114\t280.170\t22.329\t0.206\t30.12784\t4\t-2.544\t0\tNo_error
"""
         ),
    dict(kwargs=dict(Writer=ascii.Tab),
         out="""\
ID\tXCENTER\tYCENTER\tMAG\tMERR\tMSKY\tNITER\tSHARPNESS\tCHI\tPIER\tPERROR
14\t138.538\t256.405\t15.461\t0.003\t34.85955\t4\t-0.032\t0.802\t0\tNo_error
18\t18.114\t280.170\t22.329\t0.206\t30.12784\t4\t-2.544\t1.104\t0\tNo_error
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
18 & 18.114 & 280.170 & 22.329 & 0.206 & 30.12784 & 4 & -2.544 & 1.104 & 0 & No_error
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
18 & 18.114 & 280.170 & 22.329 & 0.206 & 30.12784 & 4 & -2.544 & 1.104 & 0 & No_error
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
    dict(kwargs=dict(Writer=ascii.Latex, latexdict={'tabletype': None}),
         out="""\
\\begin{tabular}{ccccccccccc}
ID & XCENTER & YCENTER & MAG & MERR & MSKY & NITER & SHARPNESS & CHI & PIER & PERROR \\\\
 & pixels & pixels & magnitudes & magnitudes & counts &  &  &  &  & perrors \\\\
14 & 138.538 & 256.405 & 15.461 & 0.003 & 34.85955 & 4 & -0.032 & 0.802 & 0 & No_error \\\\
18 & 18.114 & 280.170 & 22.329 & 0.206 & 30.12784 & 4 & -2.544 & 1.104 & 0 & No_error \\\\
\\end{tabular}
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
    dict(kwargs=dict(fill_values=('1', 'w'),
                     fill_include_names=['b']),
         out="""\
a b c
1 2 3
1 w 3
"""
         ),
    dict(kwargs=dict(fill_values=('1', 'w'),
                     fill_exclude_names=['a']),
         out="""\
a b c
1 2 3
1 w 3
"""
         ),
    dict(kwargs=dict(fill_values=('1', 'w'),
                     fill_include_names=['a'],
                     fill_exclude_names=['a', 'b']),
         out="""\
a b c
1 2 3
1 1 3
"""
         ),
    dict(kwargs=dict(fill_values=[('1', 'w')],
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
"" 2 3
1 1 ""
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
    except ValueError as e:  # if format doesn't have a fast writer, ignore
        if 'not in the list of formats with fast writers' not in str(e):
            raise e
        return
    print('Expected:\n{}'.format(test_def['out']))
    print('Actual:\n{}'.format(out.getvalue()))
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
        table.write(out, format=format, fast_writer=fast_writer, **test_def['kwargs'])
    except ValueError as e:  # if format doesn't have a fast writer, ignore
        if 'not in the list of formats with fast writers' not in str(e):
            raise e
        return
    print('Expected:\n{}'.format(test_def['out']))
    print('Actual:\n{}'.format(out.getvalue()))
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


def test_write_invalid_toplevel_meta_ipac():
    """Write an IPAC table that contains no data but has invalid (incorrectly
    specified) metadata stored in the top-level metadata and therefore should
    raise a warning, and check that the warning has been raised"""
    table = ascii.get_reader(Reader=ascii.Ipac)
    data = table.read('t/no_data_ipac.dat')
    data.meta['blah'] = 'extra'

    with catch_warnings(AstropyWarning) as ASwarn:
        out = StringIO()
        data.write(out, format='ascii.ipac')
    assert len(ASwarn) == 1
    assert "were not written" in str(ASwarn[0].message)


def test_write_invalid_keyword_meta_ipac():
    """Write an IPAC table that contains no data but has invalid (incorrectly
    specified) metadata stored appropriately in the ``keywords`` section
    of the metadata but with invalid format and therefore should raise a
    warning, and check that the warning has been raised"""
    table = ascii.get_reader(Reader=ascii.Ipac)
    data = table.read('t/no_data_ipac.dat')
    data.meta['keywords']['blah'] = 'invalid'

    with catch_warnings(AstropyWarning) as ASwarn:
        out = StringIO()
        data.write(out, format='ascii.ipac')
    assert len(ASwarn) == 1
    assert "has been skipped" in str(ASwarn[0].message)


def test_write_valid_meta_ipac():
    """Write an IPAC table that contains no data and has *correctly* specified
    metadata.  No warnings should be issued"""
    table = ascii.get_reader(Reader=ascii.Ipac)
    data = table.read('t/no_data_ipac.dat')
    data.meta['keywords']['blah'] = {'value': 'invalid'}

    with catch_warnings(AstropyWarning) as ASwarn:
        out = StringIO()
        data.write(out, format='ascii.ipac')
    assert len(ASwarn) == 0


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

    # setting comment=False should disable comment writing
    out = StringIO()
    ascii.write(data, out, format='basic', comment=False, fast_writer=fast_writer)
    expected = ['a b c', '1 2 3']
    assert out.getvalue().splitlines() == expected


@pytest.mark.parametrize("fast_writer", [True, False])
@pytest.mark.parametrize("fmt", ['%0.1f', '.1f', '0.1f', '{0:0.1f}'])
def test_write_format(fast_writer, fmt):
    """Check different formats for a column."""
    data = ascii.read('#c1\n  # c2\t\na,b,c\n#  c3\n1.11,2.22,3.33')
    out = StringIO()
    expected = ['# c1', '# c2', '# c3', 'a b c', '1.1 2.22 3.33']
    data['a'].format = fmt
    ascii.write(data, out, format='basic', fast_writer=fast_writer)
    assert out.getvalue().splitlines() == expected


@pytest.mark.parametrize("fast_writer", [True, False])
def test_strip_names(fast_writer):
    """Names should be stripped of whitespace by default."""
    data = table.Table([[1], [2], [3]], names=(' A', 'B ', ' C '))
    out = StringIO()
    ascii.write(data, out, format='csv', fast_writer=fast_writer)
    assert out.getvalue().splitlines()[0] == 'A,B,C'


def test_latex_units():
    """
    Check to make sure that Latex and AASTex writers attempt to fall
    back on the **unit** attribute of **Column** if the supplied
    **latexdict** does not specify units.
    """
    t = table.Table([table.Column(name='date', data=['a', 'b']),
               table.Column(name='NUV exp.time', data=[1, 2])])
    latexdict = copy.deepcopy(ascii.latexdicts['AA'])
    latexdict['units'] = {'NUV exp.time': 's'}
    out = StringIO()
    expected = '''\
\\begin{table}{cc}
\\tablehead{\\colhead{date} & \\colhead{NUV exp.time}\\\\ \\colhead{ } & \\colhead{s}}
\\startdata
a & 1 \\\\
b & 2
\\enddata
\\end{table}
'''.replace('\n', os.linesep)

    ascii.write(t, out, format='aastex', latexdict=latexdict)
    assert out.getvalue() == expected
    # use unit attribute instead
    t['NUV exp.time'].unit = units.s
    t['date'].unit = units.yr
    out = StringIO()
    ascii.write(t, out, format='aastex', latexdict=ascii.latexdicts['AA'])
    assert out.getvalue() == expected.replace(
        'colhead{s}', r'colhead{$\mathrm{s}$}').replace(
        'colhead{ }', r'colhead{$\mathrm{yr}$}')


@pytest.mark.parametrize("fast_writer", [True, False])
def test_commented_header_comments(fast_writer):
    """
    Test the fix for #3562 with confusing exception using comment=False
    for the commented_header writer.
    """
    t = table.Table([[1, 2]])
    with pytest.raises(ValueError) as err:
        out = StringIO()
        ascii.write(t, out, format='commented_header', comment=False,
                    fast_writer=fast_writer)
    assert "for the commented_header writer you must supply a string" in str(err.value)


@pytest.mark.parametrize("fast_writer", [True, False])
def test_byte_string_output(fast_writer):
    """
    Test the fix for #4350 where byte strings were output with a
    leading `b` on Py3.
    """
    t = table.Table([['Hello', 'World']], dtype=['S10'])
    out = StringIO()
    ascii.write(t, out, fast_writer=fast_writer)
    assert out.getvalue().splitlines() == ['col0', 'Hello', 'World']


@pytest.mark.parametrize('names, include_names, exclude_names, formats, issues_warning', [
    (['x', 'y'], ['x', 'y'], ['x'], {'x': '%d', 'y': '%f'}, True),
    (['x', 'y'], ['x', 'y'], ['y'], {'x': '%d'}, False),
    (['x', 'y'], ['x', 'y'], [], {'p': '%d', 'q': '%f'}, True),
    (['x', 'y'], ['x', 'y'], [], {'z': '%f'}, True),
    (['x', 'y'], ['x', 'y'], [], {'x': '%d'}, False),
    (['x', 'y'], ['x', 'y'], [], {'p': '%d', 'y': '%f'}, True),
    (['x', 'y'], ['x', 'y'], [], {}, False)
])
def test_names_with_formats(names, include_names, exclude_names, formats, issues_warning):
    """Test for #4508."""
    t = table.Table([[1, 2, 3], [4.1, 5.2, 6.3]])
    with catch_warnings(AstropyWarning) as ASwarn:
        out = StringIO()
        ascii.write(t, out, names=names, include_names=include_names,
        exclude_names=exclude_names, formats=formats)
    assert (issues_warning == (len(ASwarn) == 1))


@pytest.mark.parametrize('formats, issues_warning', [
    ({'p': '%d', 'y': '%f'}, True),
    ({'x': '%d', 'y': '%f'}, True),
    ({'z': '%f'}, True),
    ({}, False)
])
def test_columns_names_with_formats(formats, issues_warning):
    """Test the fix for #4508."""
    t = table.Table([[1, 2, 3], [4.1, 5.2, 6.3]])
    with catch_warnings(AstropyWarning) as ASwarn:
        out = StringIO()
        ascii.write(t, out, formats=formats)
    assert (issues_warning == (len(ASwarn) == 1))


@pytest.mark.parametrize("fast_writer", [True, False])
def test_write_quoted_empty_field(fast_writer):
    """
    Test the fix for #4350 where byte strings were output with a
    leading `b` on Py3.
    """
    t = table.Table([['Hello', ''], ['', '']], dtype=['S10', 'S10'])
    out = StringIO()
    ascii.write(t, out, fast_writer=fast_writer)
    assert out.getvalue().splitlines() == ['col0 col1', 'Hello ""', '"" ""']

    out = StringIO()
    ascii.write(t, out, fast_writer=fast_writer, delimiter=',')
    assert out.getvalue().splitlines() == ['col0,col1', 'Hello,', ',']


@pytest.mark.parametrize("format", ['ascii', 'csv', 'html', 'latex',
                                    'ascii.fixed_width', 'html'])
@pytest.mark.parametrize("fast_writer", [True, False])
def test_write_overwrite_ascii(format, fast_writer, tmpdir):
    """Test overwrite argument for various ASCII writers"""
    filename = tmpdir.join("table-tmp.dat").strpath
    with open(filename, 'w'):
        # create empty file
        pass
    t = table.Table([['Hello', ''], ['', '']], dtype=['S10', 'S10'])

    with pytest.raises(IOError) as err:
        t.write(filename, overwrite=False, format=format,
                fast_writer=fast_writer)
    assert str(err.value).endswith('already exists')

    with catch_warnings(AstropyDeprecationWarning) as warning:
        t.write(filename, format=format, fast_writer=fast_writer)
    assert len(warning) == 1
    assert str(warning[0].message).endswith(
        "Automatically overwriting ASCII files is deprecated. "
        "Use the argument 'overwrite=True' in the future.")

    t.write(filename, overwrite=True, format=format,
            fast_writer=fast_writer)

    # If the output is a file object, overwrite is ignored
    with open(filename, 'w') as fp:
        t.write(fp, format=format,
                fast_writer=fast_writer)
        t.write(fp, overwrite=False, format=format,
                fast_writer=fast_writer)
        t.write(fp, overwrite=True, format=format,
                fast_writer=fast_writer)


fmt_name_classes = list(chain(ascii.core.FAST_CLASSES.items(),
                              ascii.core.FORMAT_CLASSES.items()))


@pytest.mark.parametrize("fmt_name_class", fmt_name_classes)
def test_roundtrip_masked(fmt_name_class):
    """
    Round trip a simple masked table through every writable format and confirm
    that reading back gives the same result.
    """
    fmt_name, fmt_cls = fmt_name_class

    if not getattr(fmt_cls, '_io_registry_can_write', True):
        return

    # Skip tests for fixed_width or HTML without bs4
    if ((fmt_name == 'html' and not HAS_BEAUTIFUL_SOUP)
            or fmt_name == 'fixed_width'):
        return

    t = simple_table(masked=True)

    out = StringIO()
    fast = fmt_name in ascii.core.FAST_CLASSES
    try:
        ascii.write(t, out, format=fmt_name, fast_writer=fast)
    except ImportError:  # Some failed dependency, e.g. PyYAML, skip test
        return

    # No-header formats need to be told the column names
    kwargs = {'names': t.colnames} if 'no_header' in fmt_name else {}

    t2 = ascii.read(out.getvalue(), format=fmt_name, fast_reader=fast, guess=False, **kwargs)

    assert t.colnames == t2.colnames
    for col, col2 in zip(t.itercols(), t2.itercols()):
        assert col.dtype.kind == col2.dtype.kind
        assert np.all(col == col2)
