# Licensed under a 3-clause BSD style license - see LICENSE.rst
from distutils import version
import numpy as np

from ...tests.helper import pytest
from ... import table
from ...table import pprint

BIG_WIDE_ARR = np.arange(2000, dtype=np.float).reshape(100, 20)
SMALL_ARR = np.arange(12, dtype=np.int).reshape(4, 3)


class MaskedTable(table.Table):
    def __init__(self, *args, **kwargs):
        kwargs['masked'] = True
        table.Table.__init__(self, *args, **kwargs)


# Fixture to run all tests for both an unmasked (ndarray) and masked
# (MaskedArray) column.
@pytest.fixture(params=[False, True])
def table_type(request):
    # return MaskedTable if request.param else table.Table
    try:
        request.param
        return MaskedTable
    except AttributeError:
        return table.Table


@pytest.mark.usefixtures('table_type')
class TestMultiD():

    def test_multidim(self, table_type):
        """Test printing with multidimensional column"""
        arr = [np.array([[1, 2],
                         [10, 20]]),
               np.array([[3, 4],
                         [30, 40]]),
               np.array([[5, 6],
                         [50, 60]])]
        t = table_type(arr)
        lines = t.pformat()
        print lines
        assert lines == ['col0 [2] col1 [2] col2 [2]',
                         '-------- -------- --------',
                         '  1 .. 2   3 .. 4   5 .. 6',
                         '10 .. 20 30 .. 40 50 .. 60']

        lines = t.pformat(html=True)
        assert lines == ['<table>',
                         '<tr><th>col0 [2]</th><th>col1 [2]</th><th>col2 [2]</th></tr>',
                         '<tr><td>1 .. 2</td><td>3 .. 4</td><td>5 .. 6</td></tr>',
                         '<tr><td>10 .. 20</td><td>30 .. 40</td><td>50 .. 60</td></tr>',
                         '</table>']
        assert t._repr_html_() == ('<table><tr><th>col0 [2]</th><th>col1 [2]</th><th>col2 [2]'
                                   '</th></tr><tr><td>1 .. 2</td><td>3 .. 4</td><td>5 .. 6</td>'
                                   '</tr><tr><td>10 .. 20</td><td>30 .. 40</td><td>50 .. 60</td>'
                                   '</tr></table>')

        t = table_type([arr])
        lines = t.pformat()
        print lines
        assert lines == ['col0 [2,2]',
                         '----------',
                         '   1 .. 20',
                         '   3 .. 40',
                         '   5 .. 60']


def test_html_escaping():
    t = table.Table([('<script>alert("gotcha");</script>', 2, 3)])
    assert t._repr_html_() == (
        '<table><tr><th>col0</th></tr><tr>'
        '<td>&lt;script&gt;alert(&quot;gotcha&quot;);&lt;/script&gt;</td>'
        '</tr><tr><td>2</td></tr><tr><td>3</td></tr></table>')


@pytest.mark.usefixtures('table_type')
class TestPprint():

    def _setup(self, table_type):
        self.tb = table_type(BIG_WIDE_ARR)
        self.tb['col0'].format = '%e'
        self.tb['col1'].format = '%.6f'

        self.tb['col0'].unit = 'km**2'
        self.tb['col19'].unit = 'kg s m**-2'
        self.ts = table_type(SMALL_ARR)

    def test_format0(self, table_type):
        """Try getting screen size but fail to defaults because testing doesn't
        have access to screen (fcntl.ioctl fails).
        """
        self._setup(table_type)
        arr = np.arange(4000, dtype=np.float).reshape(100, 40)
        lines = table_type(arr).pformat()
        assert len(lines) == pprint.MAX_LINES()
        for line in lines:
            assert (len(line) > pprint.MAX_WIDTH() - 10 and
                    len(line) <= pprint.MAX_WIDTH())

    def test_format1(self, table_type):
        """Basic test of formatting"""
        self._setup(table_type)
        lines = self.tb.pformat(max_lines=8, max_width=40)
        assert lines == ['    col0         col1    ... col19 ',
                         '------------ ----------- ... ------',
                         '0.000000e+00    1.000000 ...   19.0',
                         '2.000000e+01   21.000000 ...   39.0',
                         '4.000000e+01   41.000000 ...   59.0',
                         '         ...         ... ...    ...',
                         '1.960000e+03 1961.000000 ... 1979.0',
                         '1.980000e+03 1981.000000 ... 1999.0']

    def test_format2(self, table_type):
        """Include the unit header row"""
        self._setup(table_type)
        lines = self.tb.pformat(max_lines=8, max_width=40, show_unit=True)

        print(lines)
        assert lines == ['    col0         col1    ...   col19  ',
                         '    km2                  ... kg s / m2',
                         '------------ ----------- ... ---------',
                         '0.000000e+00    1.000000 ...      19.0',
                         '2.000000e+01   21.000000 ...      39.0',
                         '         ...         ... ...       ...',
                         '1.960000e+03 1961.000000 ...    1979.0',
                         '1.980000e+03 1981.000000 ...    1999.0']

    def test_format3(self, table_type):
        """Do not include the name header row"""
        self._setup(table_type)
        lines = self.tb.pformat(max_lines=8, max_width=40, show_name=False)
        assert lines == ['0.000000e+00    1.000000 ...   19.0',
                         '2.000000e+01   21.000000 ...   39.0',
                         '4.000000e+01   41.000000 ...   59.0',
                         '6.000000e+01   61.000000 ...   79.0',
                         '         ...         ... ...    ...',
                         '1.940000e+03 1941.000000 ... 1959.0',
                         '1.960000e+03 1961.000000 ... 1979.0',
                         '1.980000e+03 1981.000000 ... 1999.0']

    def test_noclip(self, table_type):
        """Basic table print"""
        self._setup(table_type)
        lines = self.ts.pformat(max_lines=-1, max_width=-1)
        assert lines == ['col0 col1 col2',
                         '---- ---- ----',
                         '   0    1    2',
                         '   3    4    5',
                         '   6    7    8',
                         '   9   10   11']

    def test_clip1(self, table_type):
        """max lines below hard limit of 6
        """
        self._setup(table_type)
        lines = self.ts.pformat(max_lines=3, max_width=-1)
        assert lines == ['col0 col1 col2',
                         '---- ---- ----',
                         '   0    1    2',
                         '   3    4    5',
                         '   6    7    8',
                         '   9   10   11']

    def test_clip2(self, table_type):
        """max lines below hard limit of 6 and output longer than 6
        """
        self._setup(table_type)
        lines = self.ts.pformat(max_lines=3, max_width=-1, show_unit=True)
        assert lines == ['col0 col1 col2',
                         '              ',
                         '---- ---- ----',
                         '   0    1    2',
                         ' ...  ...  ...',
                         '   9   10   11']

    def test_clip3(self, table_type):
        """Max lines below hard limit of 6 and max width below hard limit
        of 10
        """
        self._setup(table_type)
        lines = self.ts.pformat(max_lines=3, max_width=1, show_unit=True)
        assert lines == ['col0 ...',
                         '     ...',
                         '---- ...',
                         '   0 ...',
                         ' ... ...',
                         '   9 ...']

    def test_clip4(self, table_type):
        """Test a range of max_lines"""
        self._setup(table_type)
        for max_lines in range(130):
            lines = self.tb.pformat(max_lines=max_lines)
            assert len(lines) == max(6, min(102, max_lines))


@pytest.mark.usefixtures('table_type')
class TestFormat():

    def test_column_format(self, table_type):
        t = table_type([[1, 2], [3, 4]], names=('a', 'b'))
        # default (format=None)
        assert str(t['a']) == ' a \n---\n  1\n  2'

        #  Old-style that is almost new-style
        t['a'].format = '{ %4.2f }'
        assert str(t['a']) == '   a    \n--------\n{ 1.00 }\n{ 2.00 }'

        #  New-style that is almost old-style
        t['a'].format = '%{0:}'
        assert str(t['a']) == ' a \n---\n %1\n %2'

        #  New-style with extra spaces
        t['a'].format = ' {0:05d} '
        assert str(t['a']) == '   a   \n-------\n 00001 \n 00002 '

        #  New-style has precedence
        t['a'].format = '%4.2f {0:}'
        assert str(t['a']) == '   a   \n-------\n%4.2f 1\n%4.2f 2'

        #  Invalid format spec
        t['a'].format = 'fail'
        with pytest.raises(ValueError):
            str(t['a'])

    def test_column_format_with_threshold(self, table_type):
        MAX_LINES_val = pprint.MAX_LINES()
        pprint.MAX_LINES.set(6)
        t = table_type([np.arange(20)], names=['a'])
        t['a'].format = '%{0:}'
        assert str(t['a']) == ' a \n---\n %0\n %1\n...\n%19'
        t['a'].format = '{ %4.2f }'
        assert str(t['a']) == '    a    \n---------\n { 0.00 }\n' \
                              ' { 1.00 }\n      ...\n{ 19.00 }'
        pprint.MAX_LINES.set(MAX_LINES_val)

    def test_column_format_func(self, table_type):
        # run most of functions twice
        # 1) astropy.table.pprint._format_funcs gets populated
        # 2) astropy.table.pprint._format_funcs gets used

        t = table_type([[1., 2.], [3, 4]], names=('a', 'b'))

        # mathematical function
        t['a'].format = lambda x: str(x * 3.)
        assert str(t['a']) == ' a \n---\n3.0\n6.0'
        assert str(t['a']) == ' a \n---\n3.0\n6.0'

    def test_column_format_func_wrong_number_args(self, table_type):
        t = table_type([[1., 2.], [3, 4]], names=('a', 'b'))

        # function that expects wrong number of arguments
        def func(a, b):
            pass

        t['a'].format = func
        with pytest.raises(ValueError):
            str(t['a'])

    def test_column_format_func_multiD(self, table_type):
        arr = [np.array([[1, 2],
                         [10, 20]])]
        t = table_type(arr, names=['a'])

        # mathematical function
        t['a'].format = lambda x: str(x * 3.)
        outstr = '   a [2]    \n------------\n  3.0 .. 6.0\n30.0 .. 60.0'
        assert str(t['a']) == outstr
        assert str(t['a']) == outstr

    def test_column_format_func_not_str(self, table_type):
        t = table_type([[1., 2.], [3, 4]], names=('a', 'b'))

        # mathematical function
        t['a'].format = lambda x: x * 3
        with pytest.raises(ValueError):
            str(t['a'])
