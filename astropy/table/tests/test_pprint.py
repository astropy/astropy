# Licensed under a 3-clause BSD style license - see LICENSE.rst
from distutils import version
import numpy as np

from ...tests.helper import pytest
from ... import table
from ...table import pprint

BIG_WIDE_ARR = np.arange(2000, dtype=np.float).reshape(100, 20)
SMALL_ARR = np.arange(12, dtype=np.int).reshape(4, 3)

numpy_lt_1p5 = version.LooseVersion(np.__version__) < version.LooseVersion('1.5')

# Dummy init of Table for pyflakes and to be sure test fixture is working
Table = None


class MaskedTable(table.Table):
    def __init__(self, *args, **kwargs):
        kwargs['masked'] = True
        table.Table.__init__(self, *args, **kwargs)


# Fixture to run all tests for both an unmasked (ndarray) and masked (MaskedArray) column.
@pytest.fixture(params=[False] if numpy_lt_1p5 else [False, True])
def set_global_Table(request):
    global Table
    Table = MaskedTable if request.param else table.Table


@pytest.mark.usefixtures('set_global_Table')
class TestMultiD():

    def test_multidim(self):
        """Test printing with multidimensional column"""
        arr = [np.array([[1, 2],
                         [10, 20]]),
               np.array([[3, 4],
                         [30, 40]]),
               np.array([[5, 6],
                         [50, 60]])]
        t = Table(arr)
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

        t = Table([arr])
        lines = t.pformat()
        print lines
        assert lines == ['col0 [2,2]',
                         '----------',
                         '   1 .. 20',
                         '   3 .. 40',
                         '   5 .. 60']


@pytest.mark.usefixtures('set_global_Table')
class TestPprint():

    def setup_method(self, method):
        self.tb = Table(BIG_WIDE_ARR)
        self.tb['col0'].format = '%e'
        self.tb['col1'].format = '%.6f'
        self.tb['col0'].units = 'km**2'
        self.tb['col19'].units = 'kg s m**-2'
        self.ts = Table(SMALL_ARR)

    def test_format0(self):
        """Try getting screen size but fail to defaults because testing doesn't
        have access to screen (fcntl.ioctl fails).
        """
        arr = np.arange(4000, dtype=np.float).reshape(100, 40)
        lines = Table(arr).pformat()
        assert len(lines) == pprint.MAX_LINES()
        for line in lines:
            assert (len(line) > pprint.MAX_WIDTH() - 10 and
                    len(line) <= pprint.MAX_WIDTH())

    def test_format1(self):
        """Basic test of formatting"""
        lines = self.tb.pformat(max_lines=8, max_width=40)
        assert lines == ['    col0         col1    ... col19 ',
                         '------------ ----------- ... ------',
                         '0.000000e+00    1.000000 ...   19.0',
                         '2.000000e+01   21.000000 ...   39.0',
                         '4.000000e+01   41.000000 ...   59.0',
                         '         ...         ... ...    ...',
                         '1.960000e+03 1961.000000 ... 1979.0',
                         '1.980000e+03 1981.000000 ... 1999.0']

    def test_format2(self):
        """Include the units header row"""
        lines = self.tb.pformat(max_lines=8, max_width=40, show_units=True)
        print(lines)
        assert lines == ['    col0         col1    ...   col19  ',
                         '    km2                  ... kg s / m2',
                         '------------ ----------- ... ---------',
                         '0.000000e+00    1.000000 ...      19.0',
                         '2.000000e+01   21.000000 ...      39.0',
                         '         ...         ... ...       ...',
                         '1.960000e+03 1961.000000 ...    1979.0',
                         '1.980000e+03 1981.000000 ...    1999.0']

    def test_format3(self):
        """Do not include the name header row"""
        lines = self.tb.pformat(max_lines=8, max_width=40, show_name=False)
        assert lines == ['0.000000e+00    1.000000 ...   19.0',
                         '2.000000e+01   21.000000 ...   39.0',
                         '4.000000e+01   41.000000 ...   59.0',
                         '6.000000e+01   61.000000 ...   79.0',
                         '         ...         ... ...    ...',
                         '1.940000e+03 1941.000000 ... 1959.0',
                         '1.960000e+03 1961.000000 ... 1979.0',
                         '1.980000e+03 1981.000000 ... 1999.0']

    def test_noclip(self):
        """Basic table print"""
        lines = self.ts.pformat(max_lines=-1, max_width=-1)
        assert lines == ['col0 col1 col2',
                         '---- ---- ----',
                         '   0    1    2',
                         '   3    4    5',
                         '   6    7    8',
                         '   9   10   11']

    def test_clip1(self):
        """max lines below hard limit of 6
        """
        lines = self.ts.pformat(max_lines=3, max_width=-1)
        assert lines == ['col0 col1 col2',
                         '---- ---- ----',
                         '   0    1    2',
                         '   3    4    5',
                         '   6    7    8',
                         '   9   10   11']

    def test_clip2(self):
        """max lines below hard limit of 6 and output longer than 6
        """
        lines = self.ts.pformat(max_lines=3, max_width=-1, show_units=True)
        assert lines == ['col0 col1 col2',
                         '              ',
                         '---- ---- ----',
                         '   0    1    2',
                         ' ...  ...  ...',
                         '   9   10   11']

    def test_clip3(self):
        """Max lines below hard limit of 6 and max width below hard limit
        of 10
        """
        lines = self.ts.pformat(max_lines=3, max_width=1, show_units=True)
        assert lines == ['col0 ...',
                         '     ...',
                         '---- ...',
                         '   0 ...',
                         ' ... ...',
                         '   9 ...']

    def test_clip4(self):
        """Test a range of max_lines"""
        for max_lines in range(130):
            lines = self.tb.pformat(max_lines=max_lines)
            assert len(lines) == max(6, min(102, max_lines))


@pytest.mark.usefixtures('set_global_Table')
class TestFormat():

    def test_column_format(self):
        t = Table([[1, 2], [3, 4]], names=('a', 'b'))
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

    def test_column_format_with_threshold(self):
        MAX_LINES_val = pprint.MAX_LINES()
        pprint.MAX_LINES.set(6)
        t = Table([np.arange(20)], names=['a'])
        t['a'].format = '%{0:}'
        assert str(t['a']) == ' a \n---\n %0\n %1\n...\n%19'
        t['a'].format = '{ %4.2f }'
        assert str(t['a']) == '    a    \n---------\n { 0.00 }\n' \
                              ' { 1.00 }\n      ...\n{ 19.00 }'
        pprint.MAX_LINES.set(MAX_LINES_val)

    def test_column_format_func(self):
        # run most of functions twice
        # 1) astropy.table.pprint._format_funcs gets populated
        # 2) astropy.table.pprint._format_funcs gets used
        
        t = Table([[1., 2.], [3, 4]], names=('a', 'b'))
        
        # mathematical function
        t['a'].format = lambda x: x*3.
        assert str(t['a']) == ' a \n---\n3.0\n6.0'
        assert str(t['a']) == ' a \n---\n3.0\n6.0'

        #function that expects wrong number of arguments
        def func(a,b):
            pass

        t['a'].format = func
        with pytest.raises(ValueError):
            str(t['a'])

    def test_column_format_func_multiD(self):
        arr = [np.array([[1, 2],
                         [10, 20]])]
        t = Table(arr, names = ['a'])

        # mathematical function
        t['a'].format = lambda x: x*3.
        outstr = '   a [2]    \n------------\n  3.0 .. 6.0\n30.0 .. 60.0'
        assert str(t['a']) == outstr
        assert str(t['a']) == outstr
        
